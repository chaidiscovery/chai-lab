# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
from collections import defaultdict
from pathlib import Path

import antipickle
import torch
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdmolops

# for some reason calling Chem.rdDetermineBonds doesnt work
from rdkit.Chem.rdDetermineBonds import DetermineBonds
from rdkit.Geometry import Point3D
from rdkit.rdBase import BlockLogs
from tqdm import tqdm

from chai_lab.data.parsing.structure.residue import ConformerData
from chai_lab.data.residue_constants import (
    new_ligand_residue_name,
    standard_residue_pdb_codes,
)
from chai_lab.utils import paths
from chai_lab.utils.pickle import TorchAntipickleAdapter
from chai_lab.utils.timeout import timeout

# important to set this flag otherwise atom properties such as
# "name" will be lost when pickling
# See https://github.com/rdkit/rdkit/issues/1320
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

logger = logging.getLogger(__name__)


class RefConformerGenerator:
    def __init__(
        self,
        leaving_atoms_cache_file: str | None = None,
    ):
        """
        N.B. in almost all cases, you want to use RefConformerGenerator.make() rather
        than initializing the object directly, since constructor the conformer generator
        is expensive, and we want to cache the result.

        Caches idealized 3D coordinates and list of atoms for residues that exist in the PDB
        This is needed to create empty atom coordinates and mask for missing residues
        and ensure the number of tokens and atoms is the same for chains with the same entity_id
        """
        # Mapping of molecule names to (atom_names, leaving_atoms); leaving atoms
        # correspond to True. See the following file for how this was constructed:
        # src/scripts/small_molecule_preprocess/leaving_atoms.py
        self.leaving_atoms: dict[str, tuple[list[str], list[bool]]] = dict()
        if leaving_atoms_cache_file is not None:
            self.leaving_atoms = antipickle.load(leaving_atoms_cache_file)

        # download conformers' cache if needed
        conformers_cache_file = paths.cached_conformers.get_path().as_posix()
        # load cached conformers after leaving atoms cache is generated in
        # case we need to re-generate the cache
        self.cached_conformers = self._load_apkl_conformers(conformers_cache_file)

        if new_ligand_residue_name in self.cached_conformers:
            self.cached_conformers.pop(new_ligand_residue_name)

        assert len(self.cached_conformers) > 0

    def _load_apkl_conformers(self, path: str) -> dict[str, ConformerData]:
        assert path.endswith(".apkl")
        assert Path(path).exists()
        return antipickle.load(path, adapters=_get_adapters())

    def _load_cached_conformers(self, path: str) -> dict[str, ConformerData]:
        block = BlockLogs()
        with Chem.SDMolSupplier(path) as suppl:
            mols = [m for m in suppl if m is not None]
        del block
        logger.info(f"Loaded {len(mols)} cached conformers")

        residues_dict = {
            m.GetProp("_Name"): self._load_ref_conformer_from_rdkit(m)
            for m in tqdm(mols)
        }

        # check at least standard residues were loaded
        # otherwise missing protein residues cannot be handled
        for res_name in standard_residue_pdb_codes:
            assert (
                res_name in residues_dict
            ), f"Standard residue {res_name} should have a reference conformer loaded"

        return residues_dict

    @classmethod
    def _load_ref_conformer_from_rdkit(self, mol: Chem.Mol) -> ConformerData:
        mol = Chem.RemoveAllHs(mol)

        ref_pos = torch.tensor(mol.GetConformer().GetPositions(), dtype=torch.float)

        ref_atom_names = [atom.GetProp("name") for atom in mol.GetAtoms()]

        ref_atom_charge = torch.tensor(
            [atom.GetFormalCharge() for atom in mol.GetAtoms()], dtype=torch.int
        )
        ref_atom_element = torch.tensor(
            [atom.GetAtomicNum() for atom in mol.GetAtoms()], dtype=torch.int
        )

        bonds = [
            (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()
        ]

        symms = get_intra_res_atom_symmetries(mol)

        symmetries = (
            torch.stack([torch.tensor(x) for x in symms], dim=-1)
            if len(symms) > 0
            else torch.arange(len(ref_atom_names)).unsqueeze(-1)
        )

        return ConformerData(
            position=ref_pos,
            element=ref_atom_element,
            charge=ref_atom_charge,
            atom_names=ref_atom_names,
            bonds=bonds,
            symmetries=symmetries,
        )

    def get(self, residue_name: str) -> ConformerData | None:
        """
        Returns an rdkit reference conformer if residue is in CCD and conformer
        generation succeeded. Otherwise, returns None.

        N.B. we should _not_ add more post-processing logic to this method, since we
        call this for every residue and want cache lookups to be fast for large
        chains. If you need to modify the conformer data, do that when building the
        cache instead.
        """
        return self.cached_conformers.get(residue_name)

    def generate(self, smiles: str) -> ConformerData:
        """Generates a conformer for a ligand from its SMILES string."""
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"Invalid smiles {smiles}"

        mol_with_hs = Chem.AddHs(mol)

        params = rdDistGeom.ETKDGv3()
        params.useSmallRingTorsions = True
        params.randomSeed = 123
        params.useChirality = True
        # below params were added after facing 'Value Error: Bad Conformer id'
        # https://github.com/rdkit/rdkit/issues/1433#issuecomment-305097888
        params.maxAttempts = 10_000
        params.useRandomCoords = True

        rdDistGeom.EmbedMultipleConfs(mol_with_hs, numConfs=1, params=params)
        rdmolops.RemoveHs(mol_with_hs)

        element_counter: dict = defaultdict(int)
        for atom in mol_with_hs.GetAtoms():
            elem = atom.GetSymbol()
            element_counter[elem] += 1  # Start each counter at 1
            atom.SetProp("name", elem + str(element_counter[elem]))

        retval = self._load_ref_conformer_from_rdkit(mol_with_hs)
        retval.atom_names = [a.upper() for a in retval.atom_names]
        return retval


def _get_adapters():
    ## adapters define how antipickle should serialize unknown types
    from antipickle.adapters import DataclassAdapter

    return [TorchAntipickleAdapter(), DataclassAdapter(dict(conf=ConformerData))]


def conformer_data_to_rdkit_mol(conformer: ConformerData) -> Chem.Mol:
    """Convert ConformerData to RDKit Mol
    RDKit Molecules can be used infer bonds (often better than the PDB) and compute
    intra-residue atom symmetries.
    """

    # Create an editable molecule object and add atoms
    editable_mol = Chem.RWMol()

    # Add atoms to the molecule
    for atom_type, atom_name in zip(conformer.element, conformer.atom_names):
        atom = Chem.Atom(atom_type.item())
        atom.SetProp("name", atom_name)
        editable_mol.AddAtom(atom)

    # Create a conformer to hold the 3D coordinates
    rd_conformer = Chem.Conformer(len(conformer.element))

    # Set the coordinates for each atom
    for i, pos in enumerate(conformer.position.tolist()):
        rd_conformer.SetAtomPosition(i, Point3D(*pos))

    # Add the conformer and convert back to standard molecule instance
    editable_mol.AddConformer(rd_conformer)
    # add bonds
    mol = editable_mol.GetMol()
    mol = maybe_add_bonds(mol)
    return mol


def maybe_add_bonds(mol: Chem.Mol, timeout_after: float = 1.0) -> Chem.Mol:
    """Attempts to add bonds to a molecule. Returns original molecule if not
    successful

    The RDKit determineBonds function is known to hang for certain molecules.
    This function wraps the call in a timeout.

    """

    @timeout(timeout_after)
    def _add_bonds(mol):
        # hard-to-find function for inferring bond information
        # https://rdkit.org/docs/source/rdkit.Chem.rdDetermineBonds.html
        # We wrap this in a timeout because this function is known to hang
        # for some molecules. See Issue
        # (https://github.com/rdkit/rdkit/discussions/7289#discussioncomment-8930333)
        DetermineBonds(mol)
        return mol

    try:
        mol = _add_bonds(mol)
    except ValueError as e:
        # ValueError is caused by rdKit, e.g.
        # - "could not find valid bond ordering"
        # - "determineBondOrdering() does not work with element Os"
        logger.warning(f"Failed to determine bonds for {Chem.MolToSmiles(mol)}, {e}")
    except TimeoutError as e:
        # TimoutError is cause by bug in rdkit
        logger.warning(f"Failed to determine bonds for {Chem.MolToSmiles(mol)}, {e}")

    return mol


def get_intra_res_atom_symmetries(
    mol: Chem.Mol, max_symmetries: int = 1000, timeout_after: float = 1.0
) -> tuple[tuple[int, ...]]:
    """Attempts to compute full set of intra-residue atom symmetries. Returns identity
    permutation of atoms if not successful"""

    @timeout(timeout_after)
    def _get_symmetries():
        return mol.GetSubstructMatches(
            mol, uniquify=False, maxMatches=max_symmetries, useChirality=False
        )

    try:
        symms = _get_symmetries()
    except TimeoutError:
        # Issues of hangup have been reported for certain ligand pairs
        # Issue(https://github.com/michellab/BioSimSpace/issues/100)
        # NOTE: this function calls MCS algorithm described in linked issue.
        symms = (tuple(range(mol.GetNumAtoms())),)

    return symms

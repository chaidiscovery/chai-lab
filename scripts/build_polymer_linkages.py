# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

"""Build peptide-linkage metadata from the wwPDB Chemical Component Dictionary."""

import argparse
import json
from pathlib import Path

import gemmi


def _select_link_atom(candidates: list[str], conventional_name: str) -> str | None:
    candidates = sorted(set(candidates))
    if conventional_name in candidates:
        return conventional_name
    return candidates[0] if len(candidates) == 1 else None


def _extract_linkage(
    block: gemmi.cif.Block,
) -> dict[str, str | list[str] | None] | None:
    component_type = gemmi.cif.as_string(block.find_value("_chem_comp.type")).upper()
    if "PEPTIDE" not in component_type or "LINKING" not in component_type:
        return None
    if component_type == "PEPTIDE-LIKE":
        return None

    terminal_tags = [
        "_chem_comp_atom.pdbx_backbone_atom_flag",
        "_chem_comp_atom.pdbx_n_terminal_atom_flag",
        "_chem_comp_atom.pdbx_c_terminal_atom_flag",
    ]
    has_terminal_flags = all(block.find_loop(tag) for tag in terminal_tags)
    atom_tags = [
        "_chem_comp_atom.atom_id",
        "_chem_comp_atom.type_symbol",
        "_chem_comp_atom.pdbx_leaving_atom_flag",
    ]
    if has_terminal_flags:
        atom_tags.extend(terminal_tags)

    atoms = {}
    for row in block.find(atom_tags):
        atoms[row.str(0)] = {
            "element": row.str(1).upper(),
            "leaving": row.str(2).upper() == "Y",
            "n_terminal": has_terminal_flags and row.str(4).upper() == "Y",
            "c_terminal": has_terminal_flags and row.str(5).upper() == "Y",
        }

    neighbors = {atom_name: [] for atom_name in atoms}
    bond_tags = ["_chem_comp_bond.atom_id_1", "_chem_comp_bond.atom_id_2"]
    for row in block.find(bond_tags):
        atom_a, atom_b = row.str(0), row.str(1)
        if atom_a in neighbors and atom_b in neighbors:
            neighbors[atom_a].append(atom_b)
            neighbors[atom_b].append(atom_a)

    if has_terminal_flags:
        n_candidates = [
            name
            for name, atom in atoms.items()
            if atom["n_terminal"] and not atom["leaving"] and atom["element"] == "N"
        ]
        c_candidates = [
            name
            for name, atom in atoms.items()
            if atom["c_terminal"] and not atom["leaving"] and atom["element"] == "C"
        ]
    else:
        n_candidates = []
        c_candidates = []
        for name, atom in atoms.items():
            if not atom["leaving"]:
                continue
            for neighbor in neighbors[name]:
                if atom["element"] == "H" and atoms[neighbor]["element"] == "N":
                    n_candidates.append(neighbor)
                if atom["element"] == "O" and atoms[neighbor]["element"] == "C":
                    c_candidates.append(neighbor)

    n_atom = _select_link_atom(n_candidates, "N")
    c_atom = _select_link_atom(c_candidates, "C")
    if n_atom is None and "N" in atoms:
        n_atom = "N"
    if c_atom is None and "C" in atoms:
        c_atom = "C"
    if "AMINO TERMINUS" in component_type:
        n_atom = None
    if "CARBOXY TERMINUS" in component_type:
        c_atom = None

    def leaving_atoms(link_atom: str | None, terminal_key: str) -> list[str]:
        if link_atom is None:
            return []
        if has_terminal_flags:
            return sorted(
                name
                for name, atom in atoms.items()
                if atom[terminal_key] and atom["leaving"] and atom["element"] != "H"
            )
        return sorted(
            name
            for name, atom in atoms.items()
            if atom["leaving"]
            and atom["element"] != "H"
            and link_atom in neighbors[name]
        )

    if n_atom is None and c_atom is None:
        return None
    return {
        "n_atom": n_atom,
        "c_atom": c_atom,
        "n_leaving_atoms": leaving_atoms(n_atom, "n_terminal"),
        "c_leaving_atoms": leaving_atoms(c_atom, "c_terminal"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("components_cif", type=Path)
    parser.add_argument("output_json", type=Path)
    args = parser.parse_args()

    linkages = {}
    for block in gemmi.cif.read(str(args.components_cif)):
        linkage = _extract_linkage(block)
        if linkage is not None:
            linkages[block.name.upper()] = linkage

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(linkages, sort_keys=True, separators=(",", ":")) + "\n"
    )
    print(f"Wrote {len(linkages)} peptide linkages to {args.output_json}")


if __name__ == "__main__":
    main()

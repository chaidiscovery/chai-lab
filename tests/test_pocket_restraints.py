"""Pocket restraints for standard and atom-tokenized protein residues."""

import pytest
import torch

from chai_lab.data.dataset.constraints.restraint_context import (
    load_manual_restraints_for_chai1,
)
from chai_lab.data.dataset.inference_dataset import Input, load_chains_from_raw
from chai_lab.data.dataset.structure.all_atom_structure_context import (
    AllAtomStructureContext,
)
from chai_lab.data.features.generators.token_pair_pocket_restraint import (
    RestraintGroup,
    TokenPairPocketRestraint,
)
from chai_lab.data.parsing.restraints import (
    PairwiseInteraction,
    PairwiseInteractionType,
)
from chai_lab.data.parsing.structure.entity_type import EntityType


@pytest.fixture(scope="module")
def modified_pocket_context():
    chains = load_chains_from_raw(
        inputs=[
            Input(
                "CKKVAVVR(TPO)(TPJ)PKSPSSAK",
                entity_type=EntityType.PROTEIN.value,
                entity_name="target",
            ),
            Input("GGG", entity_type=EntityType.PROTEIN.value, entity_name="binder"),
        ],
        entity_name_as_subchain=False,
    )
    context = AllAtomStructureContext.merge([c.structure_context for c in chains])
    return chains, context


def _generate_pocket_feature(context, restraints):
    return TokenPairPocketRestraint()._generate(
        atom_gt_coords=context.atom_gt_coords.unsqueeze(0),
        token_asym_id=context.token_asym_id.long().unsqueeze(0),
        token_residue_index=context.token_residue_index.long().unsqueeze(0),
        token_residue_names=context.token_residue_name.unsqueeze(0),
        token_subchain_id=context.subchain_id.unsqueeze(0),
        constraints=restraints,
    )


@pytest.mark.parametrize("residue_name", ["TPO", ""])
def test_atom_tokenized_pocket_retains_all_residue_tokens(
    modified_pocket_context, residue_name, caplog
):
    _, context = modified_pocket_context
    restraints = [
        RestraintGroup("B", "A", 8, residue_name, 5.5),
        RestraintGroup("B", "A", 10, "PRO", 6.0),
    ]
    feature = _generate_pocket_feature(context, restraints)[0, :, :, 0]
    target = context.token_asym_id == 1
    binder = context.token_asym_id == 2
    tpo = target & (context.token_residue_index == 8)
    pro = target & (context.token_residue_index == 10)

    assert tpo.sum().item() == 11  # The original job failed on these 11 tokens.
    assert torch.all(feature[tpo][:, binder] == 5.5)
    assert torch.all(feature[pro][:, binder] == 6.0)
    expected = (tpo | pro)[:, None] & binder[None, :]
    assert torch.all(feature[~expected] == -1)
    assert "Error" not in caplog.text


def test_numeric_ptm_pockets_load_and_preserve_other_restraints(
    modified_pocket_context, caplog
):
    chains, context = modified_pocket_context
    provided = [
        PairwiseInteraction(
            chainA="B",
            res_idxA="",
            atom_nameA="",
            chainB="A",
            res_idxB=selector,
            atom_nameB="",
            connection_type=PairwiseInteractionType.POCKET,
            max_dist_angstrom=5.5,
        )
        for selector in ("9", "10", "P11")
    ]
    restraints = load_manual_restraints_for_chai1(chains, None, provided)
    feature = _generate_pocket_feature(context, restraints.pocket_restraints)
    feature = feature[0, :, :, 0]
    target = context.token_asym_id == 1
    binder = context.token_asym_id == 2
    selected = target & torch.isin(
        context.token_residue_index, torch.tensor([8, 9, 10])
    )
    expected = selected[:, None] & binder[None, :]

    assert torch.all(feature[expected] == 5.5)
    assert torch.all(feature[~expected] == -1)
    assert "Error" not in caplog.text


@pytest.mark.parametrize(
    "residue_index,residue_name",
    [(8, "VAL"), (10, "ALA"), (99, "")],
)
def test_invalid_pocket_selector_is_not_accepted(
    modified_pocket_context, residue_index, residue_name, caplog
):
    _, context = modified_pocket_context
    feature = _generate_pocket_feature(
        context, [RestraintGroup("B", "A", residue_index, residue_name, 5.5)]
    )

    assert torch.all(feature == -1)
    assert "Error" in caplog.text


def test_no_pocket_restraints_produce_null_feature(modified_pocket_context):
    _, context = modified_pocket_context
    assert torch.all(_generate_pocket_feature(context, None) == -1)

# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from chai_lab.data.parsing.restraints import parse_pairwise_table
from chai_lab.utils.paths import repo_root


def test_loading_restraints():
    """Small test to ensure that restraints can be loaded."""
    contact_path = repo_root / "examples" / "restraints" / "contact.restraints"
    pocket_path = repo_root / "examples" / "restraints" / "pocket.restraints"

    assert len(parse_pairwise_table(contact_path)) > 0
    assert len(parse_pairwise_table(pocket_path)) > 0

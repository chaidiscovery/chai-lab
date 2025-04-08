# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from chai_lab.tools.kalign import kalign_query_to_reference


def test_all_matches():
    ref = "RKDES"
    query = "RKDDS"
    alignment = kalign_query_to_reference(ref, query)
    assert alignment is not None

    assert alignment.reference_aligned == ref
    assert alignment.query_aligned == alignment.query_a3m_line == query
    assert alignment.reference_span == (0, 4)

# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.
import inspect
import os

import torch
from typer.testing import CliRunner

from chai_lab.chai1 import run_inference, use_deterministic_kernels


def _restore(state):
    torch.use_deterministic_algorithms(state["algos"])
    torch.backends.cudnn.deterministic = state["cudnn_det"]
    torch.backends.cudnn.benchmark = state["cudnn_bench"]
    if state["env"] is None:
        os.environ.pop("CUBLAS_WORKSPACE_CONFIG", None)
    else:
        os.environ["CUBLAS_WORKSPACE_CONFIG"] = state["env"]


def _save():
    return {
        "algos": torch.are_deterministic_algorithms_enabled(),
        "cudnn_det": torch.backends.cudnn.deterministic,
        "cudnn_bench": torch.backends.cudnn.benchmark,
        "env": os.environ.get("CUBLAS_WORKSPACE_CONFIG"),
    }


def test_use_deterministic_kernels_sets_all_knobs():
    state = _save()
    try:
        os.environ.pop("CUBLAS_WORKSPACE_CONFIG", None)
        torch.use_deterministic_algorithms(False)
        torch.backends.cudnn.deterministic = False
        torch.backends.cudnn.benchmark = True
        use_deterministic_kernels()
        assert torch.are_deterministic_algorithms_enabled()
        assert torch.backends.cudnn.deterministic
        assert not torch.backends.cudnn.benchmark
        assert os.environ["CUBLAS_WORKSPACE_CONFIG"] == ":4096:8"
    finally:
        _restore(state)


def test_existing_cublas_setting_is_kept():
    state = _save()
    try:
        os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":16:8"
        use_deterministic_kernels()
        assert os.environ["CUBLAS_WORKSPACE_CONFIG"] == ":16:8"
    finally:
        _restore(state)


def test_deterministic_defaults_off_and_is_a_cli_flag():
    assert inspect.signature(run_inference).parameters["deterministic"].default is False

    import typer

    app = typer.Typer()
    app.command("fold")(run_inference)
    result = CliRunner().invoke(app, ["--help"], env={"COLUMNS": "200"})
    assert "--deterministic" in result.output

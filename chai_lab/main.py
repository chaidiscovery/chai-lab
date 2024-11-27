# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

"""Commandline interface."""

import typer

from chai_lab.chai1 import run_inference


def cli():
    app = typer.Typer()
    app.command("fold", help="Run Chai1 to fold a complex.")(run_inference)
    app()


if __name__ == "__main__":
    cli()

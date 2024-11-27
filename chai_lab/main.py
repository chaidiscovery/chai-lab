# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

"""Commandline interface."""

import typer

from chai_lab.chai1 import app as chai1_app


def cli():
    main_app = typer.Typer()
    main_app.add_typer(chai1_app, name="fold")
    main_app()


if __name__ == "__main__":
    cli()

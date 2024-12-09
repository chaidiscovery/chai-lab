# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

"""Command line interface."""

import logging

import typer

from chai_lab.chai1 import run_inference
from chai_lab.data.parsing.msas.aligned_pqt import merge_a3m_in_directory

logging.basicConfig(level=logging.INFO)

CITATION = """
@article{Chai-1-Technical-Report,
	title        = {Chai-1: Decoding the molecular interactions of life},
	author       = {{Chai Discovery}},
	year         = 2024,
	journal      = {bioRxiv},
	publisher    = {Cold Spring Harbor Laboratory},
	doi          = {10.1101/2024.10.10.615955},
	url          = {https://www.biorxiv.org/content/early/2024/10/11/2024.10.10.615955},
	elocation-id = {2024.10.10.615955},
	eprint       = {https://www.biorxiv.org/content/early/2024/10/11/2024.10.10.615955.full.pdf}
}
""".strip()


def citation():
    """Print citation information"""
    typer.echo(CITATION)


def cli():
    app = typer.Typer()
    app.command("fold", help="Run Chai-1 to fold a complex.")(run_inference)
    app.command(
        "a3m-to-pqt",
        help="Convert all a3m files in a directory for a *single sequence* into a aligned parquet file",
    )(merge_a3m_in_directory)
    app.command("citation", help="Print citation information")(citation)
    app()


if __name__ == "__main__":
    cli()

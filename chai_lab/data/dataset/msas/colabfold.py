# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
import os
import random
import tarfile
import tempfile
import time
import typing
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

from chai_lab import __version__
from chai_lab.data.parsing.fasta import Fasta, read_fasta
from chai_lab.data.parsing.msas.aligned_pqt import expected_basename, hash_sequence
from chai_lab.data.parsing.msas.data_source import MSADataSource
from chai_lab.data.parsing.templates.m8 import parse_m8_file

logger = logging.getLogger(__name__)

TQDM_BAR_FORMAT = (
    "{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]"
)


# N.B. this function (and this function only) is copied from https://github.com/sokrypton/ColabFold
# and follows the license in that repository
# We have made modifications to how templates are returned from this function.
@typing.no_type_check  # Original ColabFold code was not well typed
def _run_mmseqs2(
    x,
    prefix,
    use_env=True,
    use_filter=True,
    use_templates=False,
    filter=None,
    use_pairing=False,
    pairing_strategy="greedy",
    host_url="https://api.colabfold.com",
    user_agent: str = "",
) -> tuple[list[str], str | None]:
    """Return a block of a3m lines and optionally template hits for each of the input sequences in x."""
    submission_endpoint = "ticket/pair" if use_pairing else "ticket/msa"

    headers = {}
    if user_agent != "":
        headers["User-Agent"] = user_agent
    else:
        logger.warning(
            "No user agent specified. Please set a user agent (e.g., 'toolname/version contact@email') to help us debug in case of problems. This warning will become an error in the future."
        )

    def submit(seqs, mode, N=101):
        n, query = N, ""
        for seq in seqs:
            query += f">{n}\n{seq}\n"
            n += 1

        while True:
            error_count = 0
            try:
                # https://requests.readthedocs.io/en/latest/user/advanced/#advanced
                # "good practice to set connect timeouts to slightly larger than a multiple of 3"
                res = requests.post(
                    f"{host_url}/{submission_endpoint}",
                    data={"q": query, "mode": mode},
                    timeout=6.02,
                    headers=headers,
                )
            except requests.exceptions.Timeout:
                logger.warning("Timeout while submitting to MSA server. Retrying...")
                continue
            except Exception as e:
                error_count += 1
                logger.warning(
                    f"Error while fetching result from MSA server. Retrying... ({error_count}/5)"
                )
                logger.warning(f"Error: {e}")
                time.sleep(5)
                if error_count > 5:
                    raise
                continue
            break

        try:
            out = res.json()
        except ValueError:
            logger.error(f"Server didn't reply with json: {res.text}")
            out = {"status": "ERROR"}
        return out

    def status(ID):
        while True:
            error_count = 0
            try:
                res = requests.get(
                    f"{host_url}/ticket/{ID}", timeout=6.02, headers=headers
                )
            except requests.exceptions.Timeout:
                logger.warning(
                    "Timeout while fetching status from MSA server. Retrying..."
                )
                continue
            except Exception as e:
                error_count += 1
                logger.warning(
                    f"Error while fetching result from MSA server. Retrying... ({error_count}/5)"
                )
                logger.warning(f"Error: {e}")
                time.sleep(5)
                if error_count > 5:
                    raise
                continue
            break
        try:
            out = res.json()
        except ValueError:
            logger.error(f"Server didn't reply with json: {res.text}")
            out = {"status": "ERROR"}
        return out

    def download(ID, path):
        error_count = 0
        while True:
            try:
                res = requests.get(
                    f"{host_url}/result/download/{ID}", timeout=6.02, headers=headers
                )
            except requests.exceptions.Timeout:
                logger.warning(
                    "Timeout while fetching result from MSA server. Retrying..."
                )
                continue
            except Exception as e:
                error_count += 1
                logger.warning(
                    f"Error while fetching result from MSA server. Retrying... ({error_count}/5)"
                )
                logger.warning(f"Error: {e}")
                time.sleep(5)
                if error_count > 5:
                    raise
                continue
            break
        with open(path, "wb") as out:
            out.write(res.content)

    # process input x
    seqs = [x] if isinstance(x, str) else x

    # compatibility to old option
    if filter is not None:
        use_filter = filter

    # setup mode
    if use_filter:
        mode = "env" if use_env else "all"
    else:
        mode = "env-nofilter" if use_env else "nofilter"

    if use_pairing:
        use_templates = False
        mode = ""
        # greedy is default, complete was the previous behavior
        if pairing_strategy == "greedy":
            mode = "pairgreedy"
        elif pairing_strategy == "complete":
            mode = "paircomplete"
        if use_env:
            mode = mode + "-env"

    # define path
    path = f"{prefix}_{mode}"
    if not os.path.isdir(path):
        os.mkdir(path)

    # call mmseqs2 api
    tar_gz_file = f"{path}/out.tar.gz"
    N, REDO = 101, True

    # deduplicate and keep track of order
    seqs_unique = []
    # TODO this might be slow for large sets
    [seqs_unique.append(x) for x in seqs if x not in seqs_unique]
    Ms = [N + seqs_unique.index(seq) for seq in seqs]
    # lets do it!
    if not os.path.isfile(tar_gz_file):
        TIME_ESTIMATE = 150 * len(seqs_unique)
        with tqdm(total=TIME_ESTIMATE, bar_format=TQDM_BAR_FORMAT) as pbar:
            while REDO:
                pbar.set_description("SUBMIT")

                # Resubmit job until it goes through
                out = submit(seqs_unique, mode, N)
                while out["status"] in ["UNKNOWN", "RATELIMIT"]:
                    sleep_time = 5 + random.randint(0, 5)
                    logger.info(f"Sleeping for {sleep_time}s. Reason: {out['status']}")
                    # resubmit
                    time.sleep(sleep_time)
                    out = submit(seqs_unique, mode, N)

                if out["status"] == "ERROR":
                    raise Exception(
                        "MMseqs2 API is giving errors. Please confirm your input is a valid protein sequence. If error persists, please try again an hour later."
                    )

                if out["status"] == "MAINTENANCE":
                    raise Exception(
                        "MMseqs2 API is undergoing maintenance. Please try again in a few minutes."
                    )

                # wait for job to finish
                ID, TIME = out["id"], 0
                pbar.set_description(out["status"])
                while out["status"] in ["UNKNOWN", "RUNNING", "PENDING"]:
                    t = 5 + random.randint(0, 5)
                    logger.info(f"Sleeping for {t}s. Reason: {out['status']}")
                    time.sleep(t)
                    out = status(ID)
                    pbar.set_description(out["status"])
                    if out["status"] == "RUNNING":
                        TIME += t
                        pbar.update(n=t)
                    # if TIME > 900 and out["status"] != "COMPLETE":
                    #  # something failed on the server side, need to resubmit
                    #  N += 1
                    #  break

                if out["status"] == "COMPLETE":
                    if TIME < TIME_ESTIMATE:
                        pbar.update(n=(TIME_ESTIMATE - TIME))
                    REDO = False

                if out["status"] == "ERROR":
                    REDO = False
                    raise Exception(
                        "MMseqs2 API is giving errors. Please confirm your input is a valid protein sequence. If error persists, please try again an hour later."
                    )

            # Download results
            download(ID, tar_gz_file)

    # prep list of a3m files
    if use_pairing:
        a3m_files = [f"{path}/pair.a3m"]
    else:
        a3m_files = [f"{path}/uniref.a3m"]
        if use_env:
            a3m_files.append(f"{path}/bfd.mgnify30.metaeuk30.smag30.a3m")

    # extract a3m files
    if any(not os.path.isfile(a3m_file) for a3m_file in a3m_files):
        with tarfile.open(tar_gz_file) as tar_gz:
            tar_gz.extractall(path)

    # templates
    template_path: str | None = None
    if use_templates:
        # print("seq\tpdb\tcid\tevalue")
        # NOTE this section has been significantly reduced to enable Chai-1 to take m8 files
        # as a common input format, while also reducing how much we ping the server.
        template_path = os.path.join(path, "pdb70.m8")
        assert os.path.isfile(template_path)

    # gather a3m lines
    a3m_lines = {}
    for a3m_file in a3m_files:
        update_M, M = True, None
        for line in open(a3m_file, "r"):
            if len(line) > 0:
                if "\x00" in line:
                    line = line.replace("\x00", "")
                    update_M = True
                if line.startswith(">") and update_M:
                    M = int(line[1:].rstrip())
                    update_M = False
                    if M not in a3m_lines:
                        a3m_lines[M] = []
                a3m_lines[M].append(line)

    a3m_lines = ["".join(a3m_lines[n]) for n in Ms]
    return a3m_lines, template_path


def _is_padding_msa_row(sequence: str) -> bool:
    """Check if the given MSA sequence is a a padding sequence."""
    seq_chars = set(sequence)
    return len(seq_chars) == 1 and seq_chars.pop() == "-"


def generate_colabfold_msas(
    protein_seqs: list[str],
    msa_dir: Path,
    msa_server_url: str,
    search_templates: bool = False,
    write_a3m_to_msa_dir: bool = False,  # Useful for manual inspection + debugging
):
    """
    Generate MSAs using the ColabFold (https://github.com/sokrypton/ColabFold)
    server. No-op if no protein sequences are given.

    N.B.:
    - the MSAs in our technical report were generated using jackhmmer, not
    ColabFold, so we would expect some difference in results.
    - this implementation relies on ColabFold's chain pairing algorithm
    rather than using Chai-1's own algorithm, which could also lead to
    differences in results.

    Places .aligned.pqt files in msa_dir; does not save intermediate a3m files.
    """
    assert msa_dir.is_dir(), "MSA directory must be a dir"
    assert not any(msa_dir.iterdir()), "MSA directory must be empty"
    if not protein_seqs:
        logger.warning("No protein sequences for MSA generation; this is a no-op.")
        return

    with tempfile.TemporaryDirectory() as tmp_dir_path:
        tmp_dir = Path(tmp_dir_path)

        mmseqs_paired_dir = tmp_dir / "mmseqs_paired"
        mmseqs_paired_dir.mkdir()

        mmseqs_dir = tmp_dir / "mmseqs"
        mmseqs_dir.mkdir()

        a3ms_dir = (tmp_dir if not write_a3m_to_msa_dir else msa_dir) / "a3ms"
        a3ms_dir.mkdir()

        # Generate MSAs for each protein chain
        logger.info(f"Running MSA generation for {len(protein_seqs)} protein sequences")

        # Identify ourselves to the ColabFold server
        user_agent = f"chai-lab/{__version__} feedback@chaidiscovery.com"

        # In paired mode, mmseqs2 returns paired a3ms where all a3ms have the same number of rows
        # and each row is already paired to have the same species. As such, we insert pairing key
        # as the i-th index of the sequence so long as it isn't a padding sequence (all -)
        paired_msas: list[str]
        if len(protein_seqs) > 1:
            paired_msas, _ = _run_mmseqs2(
                protein_seqs,
                mmseqs_paired_dir,
                use_pairing=True,
                use_templates=False,  # No templates when running paired search
                host_url=msa_server_url,
                user_agent=user_agent,
            )
        else:
            # If we only have a single protein chain, there are no paired MSAs by definition
            paired_msas = [""] * len(protein_seqs)

        # MSAs without pairing logic attached; may include sequences not contained in the paired MSA
        # Needs a second call as the colabfold server returns either paired or unpaired, not both
        per_chain_msas, template_hits_file = _run_mmseqs2(
            protein_seqs,
            mmseqs_dir,
            use_pairing=False,
            use_templates=search_templates,
            host_url=msa_server_url,
            user_agent=user_agent,
        )
        if search_templates:
            assert template_hits_file is not None and os.path.isfile(template_hits_file)
            all_templates = parse_m8_file(Path(template_hits_file))
            # query IDs are 101, 102, ... from the server; remap IDs
            query_map = {}
            for orig_query_id, orig_seq in enumerate(protein_seqs, start=101):
                h = hash_sequence(orig_seq)
                query_map[orig_query_id] = h
            all_templates["query_id"] = all_templates["query_id"].apply(query_map.get)
            assert not pd.isnull(all_templates["query_id"]).any()

            logger.info(f"Found {len(all_templates)} template hits")
            all_templates.to_csv(
                msa_dir / "all_chain_templates.m8", index=False, header=False, sep="\t"
            )

        # Process the MSAs into our internal format
        for protein_seq, pair_msa, single_msa in zip(
            protein_seqs, paired_msas, per_chain_msas, strict=True
        ):
            # Write out an A3M file for both
            hkey = hash_sequence(protein_seq.upper())
            pair_a3m_path = a3ms_dir / f"{hkey}.pair.a3m"
            pair_a3m_path.write_text(pair_msa)
            single_a3m_path = a3ms_dir / f"{hkey}.single.a3m"
            single_a3m_path.write_text(single_msa)

            ## Convert the A3M file into aligned parquet files
            # Set the pairing key as the ith-index in the sequences, skip over sequences that have
            # been inserted as padding as our internal pairing logic will match on pairing key.
            paired_fasta: list[tuple[str, str, str]] = [
                (str(pairkey), record.header, record.sequence)
                for pairkey, record in enumerate(read_fasta(pair_a3m_path))
                if not _is_padding_msa_row(record.sequence)
            ]
            pairing_key, paired_headers, paired_msa_seqs = (
                zip(*paired_fasta) if paired_fasta else ((), (), ())
            )
            unique_paired_msa_seqs = set(paired_msa_seqs)

            # Non-paired MSA sequences that weren't already covered in the paired MSA; skip header
            single_fasta: list[Fasta] = [
                record
                for i, record in enumerate(read_fasta(single_a3m_path))
                if (
                    i > 0
                    and not _is_padding_msa_row(record.sequence)
                    and record.sequence not in unique_paired_msa_seqs
                )
            ]
            single_headers = [record.header for record in single_fasta]
            single_msa_seqs = [record.sequence for record in single_fasta]
            # Create null pairing keys for each of the entries in the single MSA seq
            single_null_pair_keys = [""] * len(single_msa_seqs)

            # This shouldn't have much of an effect on the model, but we make
            # a best effort to synthesize a source database anyway
            # NOTE we already dropped the query row from the single MSAs so no need to slice
            source_databases = ["query"] + [
                (
                    MSADataSource.UNIREF90.value
                    if h.startswith("UniRef")
                    else MSADataSource.BFD_UNICLUST.value
                )
                for h in (list(paired_headers) + single_headers)[1:]
            ]

            # Combine information across paired and single hits
            all_sequences = list(paired_msa_seqs) + single_msa_seqs
            all_pairing_keys = list(pairing_key) + single_null_pair_keys
            assert (
                len(all_sequences) == len(all_pairing_keys) == len(source_databases)
            ), f"Mismatched lengths: {len(all_sequences)=} {len(all_pairing_keys)=} {len(source_databases)=}"

            # Map the MSAs to our internal format
            aligned_df = pd.DataFrame(
                data=dict(
                    sequence=all_sequences,
                    source_database=source_databases,
                    pairing_key=all_pairing_keys,
                    comment="",
                ),
            )
            msa_path = msa_dir / expected_basename(protein_seq)
            if not msa_path.exists():
                # If we have a homomer, we might see the same chain multiple
                # times. The MSAs should be identical for each.
                aligned_df.to_parquet(msa_path)

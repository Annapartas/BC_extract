#!/usr/bin/env python3
"""
BC_extract.py

JSON-driven pipeline:
1. Cutadapt finds exact 15-nt flanks in forward or reverse-complement reads.
2. Cutadapt extracts the intervening fragment and normalizes orientation.
3. Python checks the complete IUPAC fragment pattern exactly.
4. pigz compresses validated FASTQ output in parallel.
5. Per-sample count tables and reports are written.

Requirements:
    Python 3
    cutadapt
    pigz
"""

from __future__ import annotations

import argparse
import gzip
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import TextIO


IUPAC = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"},
    "S": {"C", "G"}, "W": {"A", "T"},
    "K": {"G", "T"}, "M": {"A", "C"},
    "B": {"C", "G", "T"}, "D": {"A", "G", "T"},
    "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}


def fail(message: str) -> None:
    raise ValueError(message)


def load_json(path: Path) -> dict:
    try:
        with path.open("r", encoding="utf-8") as handle:
            return json.load(handle)
    except json.JSONDecodeError as error:
        fail(
            f"Invalid JSON in {path}, line {error.lineno}, "
            f"column {error.colno}: {error.msg}"
        )


def validate_dna(sequence: str, name: str, required_length: int | None = None) -> str:
    sequence = str(sequence).strip().upper()
    if not sequence:
        fail(f"{name} is empty.")
    invalid = sorted(set(sequence) - set("ACGT"))
    if invalid:
        fail(f"{name} contains invalid characters: {', '.join(invalid)}")
    if required_length is not None and len(sequence) != required_length:
        fail(f"{name} must be {required_length} nt, but is {len(sequence)} nt.")
    return sequence


def validate_pattern(pattern: str) -> str:
    pattern = str(pattern).strip().upper()
    if not pattern:
        fail("fragment_pattern is empty.")
    invalid = sorted(set(pattern) - set(IUPAC))
    if invalid:
        fail(f"fragment_pattern contains unsupported codes: {', '.join(invalid)}")
    return pattern


def pattern_matches(sequence: str, pattern: str) -> bool:
    if len(sequence) != len(pattern):
        return False
    return all(base in IUPAC[code] for base, code in zip(sequence.upper(), pattern))


def resolve_files(patterns: list[str]) -> list[Path]:
    files: list[Path] = []
    seen: set[Path] = set()

    for pattern in patterns:
        matches = [Path(item) for item in sorted(glob.glob(pattern))]
        if not matches:
            print(f"WARNING: no input files matched: {pattern}", file=sys.stderr)
        for path in matches:
            resolved = path.resolve()
            if resolved not in seen:
                seen.add(resolved)
                files.append(path)

    if not files:
        fail("No input FASTQ files were found.")
    return files


def sample_name(path: Path) -> str:
    name = path.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return path.stem


def check_executable(name: str) -> str:
    executable = shutil.which(name)
    if executable is None:
        fail(
            f"Required program '{name}' was not found. "
            "Activate the intended Conda environment and try again."
        )
    return executable


def get_version(command: list[str]) -> str:
    result = subprocess.run(
        command,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    return result.stdout.strip().splitlines()[0]


def read_fastq(handle: TextIO, source: Path):
    record = 0
    while True:
        header = handle.readline()
        if not header:
            return

        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()
        record += 1

        if not sequence or not plus or not quality:
            fail(f"Incomplete FASTQ record {record} in {source}.")

        header = header.rstrip("\r\n")
        sequence = sequence.rstrip("\r\n")
        plus = plus.rstrip("\r\n")
        quality = quality.rstrip("\r\n")

        if not header.startswith("@"):
            fail(f"Record {record} in {source}: header does not start with '@'.")
        if not plus.startswith("+"):
            fail(f"Record {record} in {source}: third line does not start with '+'.")
        if len(sequence) != len(quality):
            fail(
                f"Record {record} in {source}: sequence and quality lengths differ."
            )

        yield header, sequence, plus, quality


def validate_candidates(
    candidate_fastq: Path,
    final_fastq_gz: Path,
    rejected_fastq_gz: Path | None,
    counts_tsv: Path,
    pattern: str,
    pigz: str,
    compression_threads: int,
    compression_level: int,
) -> dict:
    counts = Counter()
    total_candidates = 0
    validated = 0
    rejected = 0

    final_fastq_gz.parent.mkdir(parents=True, exist_ok=True)

    with final_fastq_gz.open("wb") as compressed_output:
        pigz_process = subprocess.Popen(
            [
                pigz,
                f"-{compression_level}",
                "-p", str(compression_threads),
                "-c",
            ],
            stdin=subprocess.PIPE,
            stdout=compressed_output,
        )
        if pigz_process.stdin is None:
            fail("Could not open pigz input stream.")

        rejected_handle = (
            gzip.open(rejected_fastq_gz, "wt")
            if rejected_fastq_gz is not None
            else None
        )

        try:
            with gzip.open(candidate_fastq, "rt") as candidate_handle:
                for header, sequence, plus, quality in read_fastq(
                    candidate_handle, candidate_fastq
                ):
                    total_candidates += 1

                    if pattern_matches(sequence, pattern):
                        validated += 1
                        counts[sequence] += 1
                        record = f"{header}\n{sequence}\n+\n{quality}\n"
                        pigz_process.stdin.write(record.encode("utf-8"))
                    else:
                        rejected += 1
                        if rejected_handle is not None:
                            rejected_handle.write(
                                f"{header}\n{sequence}\n{plus}\n{quality}\n"
                            )
        finally:
            pigz_process.stdin.close()
            if rejected_handle is not None:
                rejected_handle.close()

        return_code = pigz_process.wait()
        if return_code != 0:
            fail(f"pigz failed with exit code {return_code}.")

    with counts_tsv.open("w", encoding="utf-8") as handle:
        handle.write("Fragment\tCount\n")
        for sequence, count in sorted(
            counts.items(), key=lambda item: (-item[1], item[0])
        ):
            handle.write(f"{sequence}\t{count}\n")

    return {
        "cutadapt_candidates": total_candidates,
        "validated_reads": validated,
        "pattern_rejected_reads": rejected,
        "unique_validated_fragments": len(counts),
    }


def run_cutadapt(
    cutadapt: str,
    input_fastq: Path,
    candidate_fastq: Path,
    unmatched_fastq: Path | None,
    report_json: Path,
    report_text: Path,
    linked_adapter: str,
    fragment_length: int,
    cores: int,
    minimum_overlap: int,
) -> None:
    command = [
        cutadapt,
        "--cores", str(cores),
        "--revcomp",
        "--no-indels",
        "--error-rate", "0",
        "--overlap", str(minimum_overlap),
        "-Z",
        "-g", linked_adapter,
        "--minimum-length", str(fragment_length),
        "--maximum-length", str(fragment_length),
        "--json", str(report_json),
        "-o", str(candidate_fastq),
    ]

    if unmatched_fastq is not None:
        command.extend(["--untrimmed-output", str(unmatched_fastq)])
    else:
        command.append("--discard-untrimmed")

    command.append(str(input_fastq))

    with report_text.open("w", encoding="utf-8") as report_handle:
        result = subprocess.run(
            command,
            stdout=report_handle,
            stderr=subprocess.STDOUT,
            text=True,
        )

    if result.returncode != 0:
        fail(
            f"Cutadapt failed for {input_fastq} with exit code "
            f"{result.returncode}. See {report_text}."
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Extract and validate barcode fragments using Cutadapt and Python."
    )
    parser.add_argument(
        "-p", "--parameters", required=True, type=Path,
        help="JSON parameter file"
    )
    args = parser.parse_args()

    config = load_json(args.parameters)

    input_config = config.get("input", {})
    sequence_config = config.get("sequences", {})
    processing = config.get("processing", {})
    output_config = config.get("output", {})

    file_patterns = input_config.get("files", [])
    if isinstance(file_patterns, str):
        file_patterns = [file_patterns]
    if not isinstance(file_patterns, list) or not file_patterns:
        fail("input.files must be a non-empty JSON list.")

    upstream = validate_dna(
        sequence_config.get("upstream_flank", ""),
        "upstream_flank",
        required_length=15,
    )
    downstream = validate_dna(
        sequence_config.get("downstream_flank", ""),
        "downstream_flank",
        required_length=15,
    )
    pattern = validate_pattern(sequence_config.get("fragment_pattern", ""))

    cores = int(processing.get("cutadapt_cores", 4))
    compression_threads = int(processing.get("pigz_threads", 4))
    compression_level = int(processing.get("compression_level", 4))
    minimum_overlap = int(processing.get("minimum_flank_overlap", 15))
    keep_candidates = bool(processing.get("keep_cutadapt_candidates", False))

    if cores < 0:
        fail("cutadapt_cores must be 0 or a positive integer.")
    if compression_threads < 1:
        fail("pigz_threads must be at least 1.")
    if compression_level not in range(1, 10):
        fail("compression_level must be from 1 to 9.")
    if minimum_overlap != 15:
        fail("minimum_flank_overlap should be 15 for exact full-flank matching.")

    output_dir = Path(output_config.get("directory", "barcode_pipeline_output"))
    save_unmatched = bool(output_config.get("save_cutadapt_unmatched", True))
    save_pattern_rejected = bool(
        output_config.get("save_pattern_rejected", True)
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    cutadapt = check_executable("cutadapt")
    pigz = check_executable("pigz")
    files = resolve_files([str(item) for item in file_patterns])

    print(f"Cutadapt: {get_version([cutadapt, '--version'])}")
    print(f"pigz:     {get_version([pigz, '--version'])}")
    print(f"Inputs:   {len(files)}")
    print(f"Output:   {output_dir}")
    print()

    linked_adapter = f"{upstream}...{downstream}"
    combined_summary = []

    for input_fastq in files:
        sample = sample_name(input_fastq)
        sample_dir = output_dir / sample
        sample_dir.mkdir(parents=True, exist_ok=True)

        final_fastq = sample_dir / f"{sample}_validated.fastq.gz"
        counts_tsv = sample_dir / f"{sample}_fragment_counts.tsv"
        cutadapt_json = sample_dir / f"{sample}_cutadapt_report.json"
        cutadapt_text = sample_dir / f"{sample}_cutadapt_report.txt"
        unmatched_fastq = (
            sample_dir / f"{sample}_flank_unmatched.fastq.gz"
            if save_unmatched else None
        )
        rejected_fastq = (
            sample_dir / f"{sample}_pattern_rejected.fastq.gz"
            if save_pattern_rejected else None
        )

        if keep_candidates:
            candidate_fastq = sample_dir / f"{sample}_cutadapt_candidates.fastq.gz"
            temporary_directory = None
        else:
            temporary_directory = tempfile.TemporaryDirectory(
                prefix=f"{sample}_", dir=output_dir
            )
            candidate_fastq = (
                Path(temporary_directory.name) /
                f"{sample}_cutadapt_candidates.fastq.gz"
            )

        print(f"[{sample}] Cutadapt extraction")
        run_cutadapt(
            cutadapt=cutadapt,
            input_fastq=input_fastq,
            candidate_fastq=candidate_fastq,
            unmatched_fastq=unmatched_fastq,
            report_json=cutadapt_json,
            report_text=cutadapt_text,
            linked_adapter=linked_adapter,
            fragment_length=len(pattern),
            cores=cores,
            minimum_overlap=minimum_overlap,
        )

        print(f"[{sample}] Exact IUPAC validation and compression")
        validation = validate_candidates(
            candidate_fastq=candidate_fastq,
            final_fastq_gz=final_fastq,
            rejected_fastq_gz=rejected_fastq,
            counts_tsv=counts_tsv,
            pattern=pattern,
            pigz=pigz,
            compression_threads=compression_threads,
            compression_level=compression_level,
        )

        if temporary_directory is not None:
            temporary_directory.cleanup()

        summary = {
            "sample": sample,
            "input_file": str(input_fastq),
            "validated_fastq": str(final_fastq),
            "counts_tsv": str(counts_tsv),
            **validation,
        }
        combined_summary.append(summary)

        sample_summary = sample_dir / f"{sample}_pipeline_summary.json"
        sample_summary.write_text(json.dumps(summary, indent=2) + "\n")

        print(
            f"[{sample}] validated={validation['validated_reads']:,}; "
            f"pattern-rejected={validation['pattern_rejected_reads']:,}; "
            f"unique={validation['unique_validated_fragments']:,}"
        )
        print()

    summary_path = output_dir / "pipeline_summary.tsv"
    with summary_path.open("w", encoding="utf-8") as handle:
        columns = [
            "sample",
            "input_file",
            "cutadapt_candidates",
            "validated_reads",
            "pattern_rejected_reads",
            "unique_validated_fragments",
            "validated_fastq",
            "counts_tsv",
        ]
        handle.write("\t".join(columns) + "\n")
        for row in combined_summary:
            handle.write(
                "\t".join(str(row[column]) for column in columns) + "\n"
            )

    print("Pipeline complete")
    print(f"Combined summary: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

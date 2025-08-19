import sys
import re
from typing import TextIO, Iterator, Tuple


def fasta_iter(input_handle: TextIO) -> Iterator[Tuple[str, str]]:
    """
    Given a file handle, iterate over the fasta records.
    For each record, yield a tuple of (header, sequence).
    The header does NOT include the leading ">".
    The sequence is a single string with no newlines.
    """
    header = None
    sequence_parts = []

    for line in input_handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(sequence_parts)
            header = line[1:]
            sequence_parts = []
        else:
            sequence_parts.append(line)

    if header is not None:
        yield header, "".join(sequence_parts)


def read_fasta(filename: str) -> list[Tuple[str, str]]:
    """
    Read all records from a fasta file.
    """
    with open(filename, "r") as f:
        return list(fasta_iter(f))


def write_fasta(filename: str, records: list[Tuple[str, str]]):
    """
    Write a list of fasta records to a file.
    """
    with open(filename, "w") as f:
        for header, sequence in records:
            f.write(f">{header}\n")
            f.write(f"{sequence}\n")


def clstr_iter(input_handle: TextIO) -> Iterator[list[dict]]:
    """
    Given a file handle, iterate over the clusters in a .clstr file.
    For each cluster, yield a list of dictionaries, where each dictionary
    represents a sequence in the cluster.
    """
    cluster = []
    cluster_no = -1
    for line in input_handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cluster:
                yield cluster
            cluster = []
            match = re.search(r"(\d+)", line)
            if match:
                cluster_no = int(match.group(1))
        else:
            is_rep = "*" in line
            seq_id_match = re.search(r">([^\.]+)\.\.\.", line)
            seq_id = seq_id_match.group(1) if seq_id_match else ""
            length_match = re.search(r"\d+\s+(\d+)[a-z]{2}", line)
            length = int(length_match.group(1)) if length_match else 0

            if is_rep:
                identity = "100"
            else:
                identity_match = re.search(r"at\s+(.+)$", line)
                if identity_match:
                    identity = identity_match.group(1).strip()
                else:
                    # Fallback for formats without 'at'
                    identity_match = re.search(r"(\d+\.?\d*%)", line)
                    identity = identity_match.group(1) if identity_match else ""

            cluster.append({
                "cluster_no": cluster_no,
                "id": seq_id,
                "length": length,
                "is_rep": is_rep,
                "identity": identity,
            })
    if cluster:
        yield cluster

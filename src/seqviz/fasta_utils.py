"""
FASTA file parsing utilities.
Wraps BioPython's SeqIO with bioinformatics-appropriate validation.
"""

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Valid IUPAC nucleotide characters (includes ambiguity codes like N, R, Y etc.)
VALID_NUCLEOTIDES = set("ACGTURYSWKMBDHVNacgturyswkmbdhvn")


def read_fasta(filepath: str, validate: bool = True) -> str:
    """
    Parse the first sequence record from a FASTA file.

    Args:
        filepath:  Path to a .fasta or .fa file.
        validate:  If True, check the sequence contains only IUPAC nucleotides.

    Returns:
        The nucleotide sequence as a plain uppercase string.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError:        If the file is empty, not valid FASTA, or contains
                           non-IUPAC characters (when validate=True).
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"FASTA file not found: {filepath}")

    if os.path.getsize(filepath) == 0:
        raise ValueError(f"FASTA file is empty: {filepath}")

    records = list(SeqIO.parse(filepath, "fasta"))

    if not records:
        raise ValueError(
            f"No valid FASTA records found in: {filepath}\n"
            "Check that the file starts with '>' and contains a sequence."
        )

    record: SeqRecord = records[0]

    if len(records) > 1:
        print(
            f"[Warning] {filepath} contains {len(records)} records. "
            f"Using only the first: '{record.id}'"
        )

    sequence = str(record.seq).upper()

    if len(sequence) == 0:
        raise ValueError(f"Record '{record.id}' has an empty sequence body.")

    if validate:
        invalid_chars = set(sequence) - VALID_NUCLEOTIDES
        if invalid_chars:
            raise ValueError(
                f"Non-IUPAC characters found in '{record.id}': {invalid_chars}\n"
                "If this is a protein sequence, use a protein-mode parser."
            )

    return sequence


def read_fasta_all(filepath: str) -> list[tuple[str, str]]:
    """
    Parse ALL records from a multi-sequence FASTA file.

    Returns:
        List of (sequence_id, sequence) tuples.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"FASTA file not found: {filepath}")

    records = list(SeqIO.parse(filepath, "fasta"))

    if not records:
        raise ValueError(f"No valid FASTA records found in: {filepath}")

    return [(r.id, str(r.seq).upper()) for r in records]


def get_fasta_metadata(filepath: str) -> dict:
    """
    Return metadata about a FASTA file without loading the full sequence.
    Useful for logging and provenance tracking.
    """
    records = list(SeqIO.parse(filepath, "fasta"))
    if not records:
        raise ValueError(f"No records in: {filepath}")

    record = records[0]
    return {
        "id": record.id,
        "description": record.description,
        "length": len(record.seq),
        "gc_content": _gc_content(str(record.seq)),
        "num_records": len(records),
    }


def _gc_content(seq: str) -> float:
    """Calculate GC% — a basic sanity check for nucleotide sequences."""
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    return round((gc / len(seq)) * 100, 2) if seq else 0.0

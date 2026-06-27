# tests/test_fasta.py
"""
Day 10 — pytest suite for FASTA parsing utilities.
Tests both the Biopython-backed version (fasta_utils.py).

Fixture files live in tests/fixtures/ — see README for how to regenerate them.
"""

import os
import pytest
from seqviz.fasta_utils import read_fasta, read_fasta_all, get_fasta_metadata

# ─────────────────────────────────────────────────────────────────────────────
# Path helper — makes all fixture paths work regardless of where pytest is
# called from (project root, tests/ directory, or CI environment)
# ─────────────────────────────────────────────────────────────────────────────

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures")

def fixture(filename: str) -> str:
    """Return the absolute path to a fixture file."""
    return os.path.join(FIXTURES, filename)


# ─────────────────────────────────────────────────────────────────────────────
# read_fasta() — single record parsing
# ─────────────────────────────────────────────────────────────────────────────

class TestReadFasta:
    """
    Tests for read_fasta() — the primary entry point for CLI --file1/--file2.
    Every test maps to a real failure mode seen with NCBI data.
    """

    def test_reads_sequence_correctly(self):
        """
        Happy path — clean single-record FASTA returns the correct sequence.
        If this fails, something is fundamentally broken in the parser.
        """
        seq = read_fasta(fixture("valid.fasta"))
        assert seq == "ATGCTAGCTAGCTAGCTAGC"

    def test_returns_uppercase(self):
        """
        Sequences must always come back uppercase regardless of what the
        file contains. Downstream alignment matrices index by character —
        lowercase 'a' != uppercase 'A' and will silently corrupt scores.
        """
        seq = read_fasta(fixture("valid.fasta"))
        assert seq == seq.upper(), "Sequence must be uppercase"

    def test_returns_string_type(self):
        """Return type must be str, not a Biopython Seq object."""
        seq = read_fasta(fixture("valid.fasta"))
        assert isinstance(seq, str), f"Expected str, got {type(seq)}"

    def test_multiline_sequence_joined(self):
        """
        NCBI wraps FASTA sequences at 70 characters per line.
        The parser must join all lines into one continuous string.

        This is the single most common FASTA parsing bug — a naive
        line-by-line reader returns only the first 70 characters of a
        1,140 bp cytochrome b gene, and your alignment runs on a
        truncated sequence with no error message.
        """
        seq = read_fasta(fixture("multiline.fasta"))
        assert "\n" not in seq, "Newlines must be stripped from sequence"
        assert " "  not in seq, "Spaces must be stripped from sequence"
        assert len(seq) == 91,  "All 91 characters across both lines must be joined"

    def test_reads_only_first_record(self):
        """
        Multi-record files must return only the first record.
        This mirrors the --file1 CLI behaviour: one file, one sequence.
        """
        seq = read_fasta(fixture("multi_record.fasta"))
        assert seq == "ATGCATGCATGC"

    def test_does_not_return_second_record(self):
        """The second record's sequence must not appear in the output."""
        seq = read_fasta(fixture("multi_record.fasta"))
        assert "TTTTAAAACCCC" not in seq


# ─────────────────────────────────────────────────────────────────────────────
# read_fasta() — error handling
# ─────────────────────────────────────────────────────────────────────────────

class TestReadFastaErrors:
    """
    Every error test maps to a real failure mode.
    Silent failures in bioinformatics pipelines are worse than crashes —
    at least a crash tells you something went wrong.
    """

    def test_file_not_found_raises(self):
        """
        Missing file must raise FileNotFoundError, not a cryptic
        Biopython internal error.
        """
        with pytest.raises(FileNotFoundError, match="not found"):
            read_fasta(fixture("does_not_exist.fasta"))

    def test_empty_file_raises(self):
        """
        0-byte files are produced by failed NCBI downloads and broken
        shell redirects. Must raise ValueError with a clear message.
        """
        with pytest.raises(ValueError, match="empty"):
            read_fasta(fixture("empty.fasta"))

    def test_no_header_raises(self):
        """
        A file with sequence but no '>' header is not valid FASTA.
        Common when a user accidentally passes a raw .txt sequence file.
        """
        with pytest.raises(ValueError):
            read_fasta(fixture("no_header.fasta"))

    def test_invalid_characters_raises_by_default(self):
        """
        Non-IUPAC characters (digits, special chars) must raise ValueError
        when validate=True (the default).
        Catches protein sequences accidentally passed to a nucleotide tool.
        """
        with pytest.raises(ValueError, match="Non-IUPAC"):
            read_fasta(fixture("invalid_chars.fasta"))

    def test_validate_false_skips_character_check(self):
        """
        validate=False must load the sequence without raising, even if
        it contains non-standard characters.
        Needed for protein mode or custom alphabets.
        """
        seq = read_fasta(fixture("invalid_chars.fasta"), validate=False)
        assert isinstance(seq, str)
        assert len(seq) > 0

    def test_error_message_contains_filepath(self):
        """
        Error messages must include the filepath so the user knows
        which file caused the problem in a multi-file pipeline.
        """
        bad_path = fixture("does_not_exist.fasta")
        with pytest.raises(FileNotFoundError) as exc_info:
            read_fasta(bad_path)
        assert "does_not_exist.fasta" in str(exc_info.value)


# ─────────────────────────────────────────────────────────────────────────────
# read_fasta_all() — multi-record parsing
# ─────────────────────────────────────────────────────────────────────────────

class TestReadFastaAll:
    """
    Tests for read_fasta_all() — used when aligning multiple sequences
    from a single file (e.g., a multiple sequence alignment input).
    """

    def test_returns_all_records(self):
        """Must return both records from a two-record file."""
        records = read_fasta_all(fixture("multi_record.fasta"))
        assert len(records) == 2

    def test_returns_list_of_tuples(self):
        """Return type must be a list of (id, sequence) tuples."""
        records = read_fasta_all(fixture("multi_record.fasta"))
        assert isinstance(records, list)
        for item in records:
            assert isinstance(item, tuple)
            assert len(item) == 2

    def test_ids_are_correct(self):
        """Record IDs must match the headers in the file."""
        records = read_fasta_all(fixture("multi_record.fasta"))
        ids = [record_id for record_id, _ in records]
        assert "record_one" in ids
        assert "record_two" in ids

    def test_sequences_are_correct(self):
        """Sequences must match file content and be uppercase."""
        records = read_fasta_all(fixture("multi_record.fasta"))
        seq_map = {record_id: seq for record_id, seq in records}
        assert seq_map["record_one"] == "ATGCATGCATGC"
        assert seq_map["record_two"] == "TTTTAAAACCCC"

    def test_empty_file_raises(self):
        with pytest.raises(ValueError):
            read_fasta_all(fixture("empty.fasta"))

    def test_file_not_found_raises(self):
        with pytest.raises(FileNotFoundError):
            read_fasta_all(fixture("does_not_exist.fasta"))


# ─────────────────────────────────────────────────────────────────────────────
# get_fasta_metadata() — provenance and logging
# ─────────────────────────────────────────────────────────────────────────────

class TestGetFastaMetadata:
    """
    Tests for get_fasta_metadata() — used for CLI output and provenance logging.
    Provenance (knowing exactly which file and record was used) is critical
    in published bioinformatics analysis.
    """

    def test_required_keys_present(self):
        """All expected metadata fields must be present."""
        meta = get_fasta_metadata(fixture("valid.fasta"))
        for key in ["id", "description", "length", "gc_content", "num_records"]:
            assert key in meta, f"Missing metadata key: '{key}'"

    def test_id_is_correct(self):
        meta = get_fasta_metadata(fixture("valid.fasta"))
        assert meta["id"] == "test_seq_001"

    def test_length_matches_sequence(self):
        """
        Metadata length must match the actual parsed sequence length.
        A mismatch means the metadata and the sequence are out of sync.
        """
        seq  = read_fasta(fixture("valid.fasta"))
        meta = get_fasta_metadata(fixture("valid.fasta"))
        assert meta["length"] == len(seq)

    def test_gc_content_is_float(self):
        meta = get_fasta_metadata(fixture("valid.fasta"))
        assert isinstance(meta["gc_content"], float)

    def test_gc_content_in_valid_range(self):
        """GC% must be between 0.0 and 100.0 for any nucleotide sequence."""
        meta = get_fasta_metadata(fixture("valid.fasta"))
        assert 0.0 <= meta["gc_content"] <= 100.0

    def test_gc_content_known_value(self):
        """
        ATGCTAGCTAGCTAGCTAGC has 10 G/C out of 20 = 50.0% GC.
        Hardcoding a known value catches rounding or logic errors.
        """
        meta = get_fasta_metadata(fixture("valid.fasta"))
        assert meta["gc_content"] == 50.0

    def test_num_records_single_file(self):
        meta = get_fasta_metadata(fixture("valid.fasta"))
        assert meta["num_records"] == 1

    def test_num_records_multi_file(self):
        meta = get_fasta_metadata(fixture("multi_record.fasta"))
        assert meta["num_records"] == 2

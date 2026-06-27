# tests/test_alignment.py
"""
Pytest validation suite for Needleman-Wunsch and Smith-Waterman.
SeqViz: Sequence Alignment Heatmap Tool
"""

import pytest
from seqviz.algorithms.needleman_wunsch import align as needleman_wunsch
from seqviz.algorithms.smith_waterman import align as smith_waterman


# ─────────────────────────────────────────────────────────────
# FIXTURES  (reusable test data — avoids copy-pasting sequences)
# ─────────────────────────────────────────────────────────────

@pytest.fixture
def nw_seqs():
    """Classic NW benchmark pair (Durbin et al.)"""
    return "GATTACA", "GCATGCU"

@pytest.fixture
def sw_seqs():
    """Classic SW benchmark pair (Smith & Waterman 1981)"""
    return "ACACACTA", "AGCACACA"

@pytest.fixture
def standard_scoring():
    """match=1, mismatch=-1, gap=-2 — used for all tests"""
    return {"match": 1, "mismatch": -1, "gap": -2}


# ─────────────────────────────────────────────────────────────
# NEEDLEMAN-WUNSCH TESTS  (global alignment)
# ─────────────────────────────────────────────────────────────

class TestNeedlemanWunsch:
    """
    NW aligns two sequences END-TO-END.
    Every residue in both sequences must appear in the output (as a match, mismatch, or gap).
    """

    # Known score testing
    def test_nw_known_score(self, nw_seqs, standard_scoring):
        """
        GATTACA vs GCATGCU must return score = -1.
        """
        seq1, seq2 = nw_seqs
        _, _, _, score = needleman_wunsch(seq1, seq2, **standard_scoring)
        assert score == -1, f"Expected NW score -1, got {score}"

    # Known alignment testing
    def test_nw_alignment_strips_to_originals(self, nw_seqs, standard_scoring):
        """
        Remove all gaps from the aligned output — must recover the
        original sequences exactly. Validates traceback logic.
        """
        seq1, seq2 = nw_seqs
        aligned1, aligned2, _, _ = needleman_wunsch(seq1, seq2, **standard_scoring)
        assert aligned1.replace("-", "") == seq1, "seq1 residues lost in alignment"
        assert aligned2.replace("-", "") == seq2, "seq2 residues lost in alignment"

    def test_nw_aligned_lengths_equal(self, nw_seqs, standard_scoring):
        """
        Both aligned strings must be the same length.
        If they differ, your gap insertion logic is broken.
        """
        seq1, seq2 = nw_seqs
        aligned1, aligned2, _, _ = needleman_wunsch(seq1, seq2, **standard_scoring)
        assert len(aligned1) == len(aligned2), (
            f"Alignment length mismatch: {len(aligned1)} vs {len(aligned2)}"
        )

    # Identical sequences testing
    def test_nw_identical_sequences(self, standard_scoring):
        """
        Aligning a sequence against itself must score = len(seq).
        Every position is a match, zero gaps.
        """
        seq = "GATTACA"
        aligned1, aligned2, _, score = needleman_wunsch(seq, seq, **standard_scoring)
        assert score == len(seq), f"Expected {len(seq)}, got {score}"
        assert aligned1 == seq  # no gaps should be inserted
        assert aligned2 == seq

    # Completely different sequences
    def test_nw_completely_different_sequences(self, standard_scoring):
        """
        All-A vs all-T: every position is a mismatch.
        Expected score = -4 (4 mismatches × -1 each).
        This catches sign errors in your mismatch penalty.
        """
        _, _, _, score = needleman_wunsch("AAAA", "TTTT", **standard_scoring)
        assert score == -4, f"Expected -4, got {score}"

    # Empty string edge case
    def test_nw_empty_sequence_raises(self, standard_scoring):
        """
        Empty input is biologically meaningless — your code should raise
        ValueError, not silently return a garbage alignment.
        """
        with pytest.raises(ValueError):
            needleman_wunsch("", "GATTACA", **standard_scoring)

    def test_nw_invalid_characters_raises(self, standard_scoring):
        """
        Non-IUPAC characters like digits should be rejected.
        Catches missing input validation.
        """
        with pytest.raises(ValueError):
            needleman_wunsch("GATT123", "GCATGCU", **standard_scoring)


# ─────────────────────────────────────────────────────────────
# Smith-waterman tests (local alignment)
# ─────────────────────────────────────────────────────────────

class TestSmithWaterman:
    """
    SW finds the best LOCAL region of similarity.
    Key difference from NW: the matrix floors at 0 (no negative cells).
    """

    # Known score testing
    def test_sw_known_score(self, sw_seqs):
        """
        ACACACTA vs AGCACACA with match=1, mismatch=-1, gap=-2 must score 5.
        """
        seq1, seq2 = sw_seqs
        _, _, _, score = smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-2)
        assert score == 5, f"Expected SW score 5, got {score}"

    # Negative score testing
    def test_sw_score_non_negative(self, sw_seqs, standard_scoring):
        """
        SW resets any cell below 0 to 0 during matrix fill.
        A negative score means your floor logic is missing.
        """
        seq1, seq2 = sw_seqs
        _, _, _, score = smith_waterman(seq1, seq2, **standard_scoring)
        assert score >= 0, f"SW score must be >= 0, got {score}"

    # Known alignment testing
    def test_sw_local_region_is_subsequence(self, sw_seqs, standard_scoring):
        """
        The local aligned region (gaps stripped) must be a subsequence
        of the original input.
        """
        seq1, seq2 = sw_seqs
        aligned1, aligned2, _, _ = smith_waterman(seq1, seq2, **standard_scoring)
        assert aligned1.replace("-", "") in seq1, "Local region not a subsequence of seq1"
        assert aligned2.replace("-", "") in seq2, "Local region not a subsequence of seq2"

    def test_sw_aligned_lengths_equal(self, sw_seqs, standard_scoring):
        """Both local aligned strings must be the same length."""
        seq1, seq2 = sw_seqs
        aligned1, aligned2, _, _ = smith_waterman(seq1, seq2, **standard_scoring)
        assert len(aligned1) == len(aligned2)

    # Identical sequences testing
    def test_sw_identical_sequences(self, standard_scoring):
        """
        SW on identical sequences should return the full sequence as
        the local region, scoring = len(seq) × match_score.
        """
        seq = "ACACACTA"
        _, _, _, score = smith_waterman(seq, seq, **standard_scoring)
        assert score == len(seq), f"Expected {len(seq)}, got {score}"

    # Completely different sequences testing
    def test_sw_completely_different_sequences(self, standard_scoring):
        """
        All-A vs all-T: no local similarity exists.
        SW must return score = 0.
        """
        _, _, _, score = smith_waterman("AAAA", "TTTT", **standard_scoring)
        assert score == 0, f"Expected 0, got {score}"

    # Empty input edge case
    def test_sw_empty_sequence_raises(self, standard_scoring):
        """Empty input should raise ValueError."""
        with pytest.raises(ValueError):
            smith_waterman("", "AGCACACA", **standard_scoring)


# ─────────────────────────────────────────────────────────────
# PARAMETRIZED SWEEP
# ─────────────────────────────────────────────────────────────

@pytest.mark.parametrize("seq1, seq2, expected_score", [
    ("A",    "A",   1),   # single match
    ("A",    "T",  -1),   # single mismatch
    ("A",    "AA", -1),   # one gap required
    ("AA",   "A",  -1),   # gap in the other direction
])
def test_nw_parametrized_single_residues(seq1, seq2, expected_score):
    """
    Sweep tiny inputs where you can compute the answer by hand.
    If any of these fail, your scoring matrix initialisation is wrong.
    """
    _, _, _, score = needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2)
    assert score == expected_score, f"NW({seq1}, {seq2}): expected {expected_score}, got {score}"

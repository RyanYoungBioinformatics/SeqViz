import numpy as np
from seqviz.scoring import scoring, GAP, MATCH, MISMATCH


def fill_matrix(seq1: str, seq2: str, gap: int = GAP) -> np.ndarray:
    """
    Build and fill the Smith-Waterman scoring matrix.
    Scores are floored at zero — poor alignments simply reset rather than
    dragging down the score. First row and column stay zero by default
    from np.zeros(), so no initialisation loop needed.
    """
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    matrix = np.zeros((rows, cols), dtype=int)

    for i in range(1, rows):
        for j in range(1, cols):
            diagonal = matrix[i-1][j-1] + scoring(seq1[i-1], seq2[j-1], MATCH, MISMATCH)
            up       = matrix[i-1][j]   + gap
            left     = matrix[i][j-1]   + gap

            matrix[i][j] = max(0, diagonal, up, left)  # ← the only change

    return matrix


def traceback(matrix: np.ndarray, seq1: str, seq2: str, gap: int = GAP):
    """
    Walk back from the highest-scoring cell, stopping the moment we hit zero.
    That zero boundary is what makes the alignment local — we only capture
    the best-matching sub-region, not the full sequences.
    """
    aligned1 = ""
    aligned2 = ""

    # Find the highest value cell anywhere in the matrix
    i, j = np.unravel_index(np.argmax(matrix), matrix.shape)

    while matrix[i][j] > 0:
        if matrix[i][j] == matrix[i-1][j-1] + scoring(seq1[i-1], seq2[j-1], MATCH, MISMATCH):
            aligned1 = seq1[i-1] + aligned1
            aligned2 = seq2[j-1] + aligned2
            i -= 1
            j -= 1
        elif matrix[i][j] == matrix[i-1][j] + gap:
            aligned1 = seq1[i-1] + aligned1
            aligned2 = "-" + aligned2
            i -= 1
        else:
            aligned1 = "-" + aligned1
            aligned2 = seq2[j-1] + aligned2
            j -= 1

    return aligned1, aligned2


def align(
    seq1: str,
    seq2: str,
    match: int = MATCH,
    mismatch: int = MISMATCH,
    gap: int = GAP
):
    """
    Run full Smith-Waterman local alignment.
    Returns aligned subsequences, the full scoring matrix, and the best score.
    """
    matrix = fill_matrix(seq1, seq2, gap)
    aligned1, aligned2 = traceback(matrix, seq1, seq2, gap)
    best_score = int(np.max(matrix))

    return aligned1, aligned2, matrix, best_score

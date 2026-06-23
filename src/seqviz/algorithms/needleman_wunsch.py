import numpy as np
from seqviz.scoring import score, GAP, MATCH, MISMATCH


def fill_matrix(seq1: str, seq2: str, gap: int = GAP) -> np.ndarray:
    """Build and fill the Needleman-Wunsch scoring matrix."""
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    matrix = np.zeros((rows, cols), dtype=int)

    for i in range(rows):
        matrix[i][0] = i * gap
    for j in range(cols):
        matrix[0][j] = j * gap

    for i in range(1, rows):
        for j in range(1, cols):
            diagonal = matrix[i-1][j-1] + score(seq1[i-1], seq2[j-1])
            up       = matrix[i-1][j]   + gap
            left     = matrix[i][j-1]   + gap
            matrix[i][j] = max(diagonal, up, left)

    return matrix


def traceback(matrix: np.ndarray, seq1: str, seq2: str, gap: int = GAP):
    """Reconstruct the optimal alignment by walking back through the matrix."""
    aligned1 = ""
    aligned2 = ""
    i, j = len(seq1), len(seq2)

    while i > 0 or j > 0:
        if i > 0 and j > 0 and matrix[i][j] == matrix[i-1][j-1] + score(seq1[i-1], seq2[j-1]):
            aligned1 = seq1[i-1] + aligned1
            aligned2 = seq2[j-1] + aligned2
            i -= 1
            j -= 1
        elif i > 0 and matrix[i][j] == matrix[i-1][j] + gap:
            aligned1 = seq1[i-1] + aligned1
            aligned2 = "-" + aligned2
            i -= 1
        else:
            aligned1 = "-" + aligned1
            aligned2 = seq2[j-1] + aligned2
            j -= 1

    return aligned1, aligned2


def align(seq1: str, seq2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """Run full Needleman-Wunsch. Returns aligned sequences, matrix, and score."""
    matrix = fill_matrix(seq1, seq2, gap)
    aligned1, aligned2 = traceback(matrix, seq1, seq2, gap)
    final_score = int(matrix[len(seq1)][len(seq2)])
    return aligned1, aligned2, matrix, final_score

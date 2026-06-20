import numpy as np
from seqviz.scoring import scoring, GAP


def fill_matrix(seq1: str, seq2: str, gap: int = GAP) -> np.ndarray:
    """
    Build and fill the Needleman-Wunsch scoring matrix.
    Returns the completed matrix as a 2D numpy array.
    """
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    matrix = np.zeros((rows, cols), dtype=int)

    for i in range(rows):
        matrix[i][0] = i * gap

    for j in range(cols):
        matrix[0][j] = j * gap

    for i in range(1, rows):
        for j in range(1, cols):
            diagonal = matrix[i-1][j-1] + scoring(seq1[i-1], seq2[j-1])
            up       = matrix[i-1][j]   + gap
            left     = matrix[i][j-1]   + gap

            matrix[i][j] = max(diagonal, up, left)

    return matrix

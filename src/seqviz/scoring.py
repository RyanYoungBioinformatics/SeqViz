MATCH = 1
MISMATCH = -1
GAP = -2

def score(a: str, b: str, match: int = MATCH, mismatch: int = MISMATCH) -> int:
    """
    Calculate score for two sequences, a and b.
    
    DNA; Match = +1, Mismatch = -1, Gap = -2 by default.
    Protein; Possible BLOSOUM62 scoring to be added.
    """
    return match if a.upper() == b.upper() else mismatch
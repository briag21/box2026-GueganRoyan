import numpy as np
from tools import delta, rle

MARGIN = 10


def build_dp_matrix(read, extract):
    """Build the edit-distance DP matrix with row 0 set to 0 for semi-global matching."""
    matrix = np.zeros((len(read) + 1, len(extract) + 1), dtype=int)
    for i in range(len(read) + 1):
        matrix[i, 0] = i

    for i in range(1, len(read) + 1):
        for j in range(1, len(extract) + 1):
            matrix[i, j] = min(
                matrix[i - 1, j - 1] + delta(read, extract, i, j),
                matrix[i - 1, j] + 1,
                matrix[i, j - 1] + 1,
            )
    return matrix


def traceback(matrix, read, extract, end_j):
    """Backtrack from the best end column and return the start offset and RLE CIGAR."""
    operations = []
    i, j = len(read), end_j
    while i > 0:
        if j > 0 and matrix[i, j] == matrix[i - 1, j - 1] + delta(read, extract, i, j):
            operations.append("=" if read[i - 1] == extract[j - 1] else "X")
            i -= 1
            j -= 1
        elif matrix[i, j] == matrix[i - 1, j] + 1:
            operations.append("I")
            i -= 1
        else:
            operations.append("D")
            j -= 1

    operations.reverse()
    return j, rle(operations)


def extend(read, ref, start):
    """Align a read around a candidate start position and return (position, cigar, score)."""
    x = max(0, start - MARGIN)
    extract = ref[x : start + len(read) + MARGIN]

    matrix = build_dp_matrix(read, extract)
    position_min = int(np.argmin(matrix[-1]))
    score = int(matrix[-1, position_min])

    offset, cigar = traceback(matrix, read, extract, position_min)
    return x + offset, cigar, score

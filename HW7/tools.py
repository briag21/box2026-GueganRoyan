from readfa import readfq
from xopen import xopen


def rle(seq):
    """RLE a sequence (for CIGAR alignment operations)."""
    rle_seq = list()
    previous = None
    count = 0
    for element in seq:
        if element == previous:
            count += 1
        else:
            if count > 0:
                rle_seq.append(f"{count}{previous}")
            previous = element
            count = 1
    if count > 0:  # dernier run
        rle_seq.append(f"{count}{previous}")
    return "".join(rle_seq)


def delta(s1, s2, i, j):
    """Return 0 if characters at indices i and j match, else 1."""
    if s1[i - 1] != s2[j - 1]:
        return 1
    else:
        return 0


def kmers(seq, k):
    """Yield all k-mers of seq with their start indices."""
    for i in range(len(seq) - k + 1):
        yield (seq[i : i + k], i)


def revComp(seq):
    """Return the reverse complement of a DNA sequence."""
    revComp = list()
    for char in reversed(seq):
        match char:
            case "A":
                revComp.append("T")
            case "T":
                revComp.append("A")
            case "C":
                revComp.append("G")
            case "G":
                revComp.append("C")
            case _:
                revComp.append("N")
    return "".join(revComp)


def read_seqs(path):
    """Yield (name, sequence) pairs from a FASTA or FASTQ file."""
    with xopen(path) as fasta:
        for name, seq, _ in readfq(fasta):
            yield name, seq

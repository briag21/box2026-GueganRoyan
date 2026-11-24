from tools import kmers


def build_index(ref, k=11):
    """Index all k-mers in the reference sequence with their start positions."""
    index: dict[str, list[int]] = dict()
    for kmer, i in kmers(ref, k):
        if kmer in index:
            index[kmer].append(i)
        else:
            index[kmer] = [i]
    return index


def find_candidates(read, index, k=11):
    """Find candidate start positions via offset voting and return the top-voted one(s) and vote count."""
    start_positions = dict()
    for kmer, i in kmers(read, k):
        if kmer in index:
            for position in index[kmer]:
                start = position - i
                if start not in start_positions:
                    start_positions[start] = 1
                else:
                    start_positions[start] += 1
    maxs = list()
    max_val = 0
    for element in start_positions:
        val = start_positions[element]
        if val > max_val:
            max_val = val
            maxs = [element]
        elif val == max_val:
            maxs.append(element)
    return maxs, max_val

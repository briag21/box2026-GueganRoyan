from tools import read_seqs, revComp
from seeding import build_index, find_candidates
from extension import extend


def align_read(name, read, ref_name, ref, index, k=11):
    """Find best-scoring alignments for a single read across both strands."""
    seeds = []
    for strand, seq in (("F", read), ("R", revComp(read))):
        cands, votes = find_candidates(seq, index, k)
        for start in cands:
            seeds.append((votes, strand, seq, start))

    # No alignment candidate found
    if not seeds:
        return [(name, read, ref_name, -1, "*", "*", -1)]

    # Prefilter
    best_votes = max(s[0] for s in seeds)
    seeds = [s for s in seeds if s[0] >= 0.8 * best_votes]

    results = set()
    for _, strand, seq, start in seeds:
        pos, cigar, score = extend(seq, ref, start)
        results.add((score, pos, cigar, strand))

    # Best score is the lowest cost (edit distance)
    best_score = min(r[0] for r in results)
    rows = []

    # Sorted(results) is necessary for a deterministic output
    for score, pos, cigar, strand in sorted(results):
        if score == best_score:
            rows.append((name, read, ref_name, pos, cigar, strand, score))
    return rows


def align(ref_name, ref, reads, k=11):
    """Index reference k-mers and yield best alignment records for each read."""
    index = build_index(ref, k)
    for read_name, read in reads:
        yield from align_read(read_name, read, ref_name, ref, index, k)


def main(
    ref_path="./inputs/ref.fa",
    reads_path="./inputs/reads.fq",
    out_path="./outfile/out.tsv",
):
    ref_name, ref = next(read_seqs(ref_path))

    with open(out_path, "w") as out:
        for row in align(ref_name, ref, read_seqs(reads_path)):
            out.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()

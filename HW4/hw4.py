import matplotlib.pyplot as plt
import networkx as nx
from readfa import readfq
from xopen import xopen
from kmer_treatment import make_canonical, reverse_complement
from collections import Counter

FILES = [
    "./fastas/ecoli_sample_perfect_reads_forward.fasta.gz",
    "./fastas/ecoli_sample_perfect_reads.fasta.gz",
    "./fastas/ecoli_sample_reads_001.fasta.gz",
    "./fastas/ecoli_sample_reads_01.fasta.gz",
]
K = 31

def read_seqs(path):
    with xopen(path) as fasta:
        for _, seq, _ in readfq(fasta):
            yield seq


def count_kmers(seqs, k):
    counts = Counter()
    for seq in seqs:
        for i in range(len(seq) - k + 1):
            counts[make_canonical(seq[i:i + k])] += 1
    return counts


def create_dbg(path, k, t=1):
    return create_dbg_from_seqs(read_seqs(path), k, t)


def create_dbg_from_seqs(seqs, k, t=1):
    return {x: c for x, c in count_kmers(seqs, k).items() if c >= t}


def in_dbg(dbg, x):
    return make_canonical(x) in dbg


def successors(dbg, x):
    return [x[1:] + b for b in "ACGT" if in_dbg(dbg, x[1:] + b)]


def to_nx(dbg):
    G = nx.DiGraph()
    for x in dbg:
        for y in (x, reverse_complement(x)):
            G.add_node(y)
            for s in successors(dbg, y):
                G.add_edge(y, s)
    return G


def plot_dbg(seqs, k, name):
    G = to_nx(create_dbg_from_seqs(seqs, k))
    nx.draw(G, nx.spring_layout(G, seed=0), with_labels=True, font_size=6)
    plt.savefig(name, dpi=200, bbox_inches="tight")
    plt.close()


def kmer_stats(path, k=31):
    counts = count_kmers(read_seqs(path), k)
    mult = Counter(counts.values())
    print(path)
    print("k-mers processed:", sum(counts.values()))
    print("distinct k-mers (in the graph):", len(counts))
    return mult


def plot_mult(mults, names, out="mult.png"):
    for m, name in zip(mults, names):
        xs = sorted(m)
        plt.plot(xs, [m[x] for x in xs], marker=".", label=name)
    plt.yscale("log")
    plt.xlabel("multiplicity")
    plt.ylabel("number of k-mers")
    plt.legend()
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()

def evaluate_t(path, truth):
    for t in range(1, 8):
        kept = create_dbg(path, 31, t).keys()
        print(t, len(kept), len(kept - truth), len(truth - kept))

def predecessors(dbg, x):
    return [b + x[:-1] for b in "ACGT" if in_dbg(dbg, b + x[:-1])]


def extend_right(dbg, s, visited):
    bases = [] # list of added bps
    cur = s
    while True:
        succ = successors(dbg, cur)
        if len(succ) != 1: #exactly one successor to continue the unitig
            break
        nxt = succ[0]
        if len(predecessors(dbg, nxt)) != 1: #exactly one predecessor for the successor to continue the unitig
            break
        if make_canonical(nxt) in visited:  #cycle detection
            break
        visited.add(make_canonical(nxt))
        bases.append(nxt[-1])
        cur = nxt
    return "".join(bases)


def unitig_from(dbg, s, visited=None):
    if not in_dbg(dbg, s):
        raise ValueError("k-mer not in the graph")
    if visited is None:
        visited = set()
    visited.add(make_canonical(s))
    right = extend_right(dbg, s, visited)
    left = extend_right(dbg, reverse_complement(s), visited)
    return reverse_complement(left) + s + right

def unitig_sizes(path, start):
    for t in (1,2,3,4):
        dbg = create_dbg(path, 31, t)
        if in_dbg(dbg, start):
            print(t, len(unitig_from(dbg, start)))
        else:
            print(t, "k-mer absent")

def create_unitigs(dbg):
    visited = set()
    unitigs = []
    for kmer in dbg:
        if kmer not in visited:
            unitigs.append(unitig_from(dbg, kmer, visited))
    return unitigs

def unitig_totals(path):
    for t in (1,2,3,4):
        unitigs = create_unitigs(create_dbg(path, 31, t))
        print(t, sum(len(u) for u in unitigs))

def main():
    # Q3: toy graphs
    plot_dbg(["TAGGCATCCGTTAGCAGTCAAGCT"], 5, "./toygraphs/simple.png")
    plot_dbg(["ACGTTGCAGTTGCA"], 5, "./toygraphs/repeat.png")

    print("\nQ4: build one de Bruijn graph per sample")
    for f in FILES:
        print(f, len(create_dbg(f, K)))

    print("\nQ5: statistics and multiplicity distribution")
    mults = [kmer_stats(f) for f in FILES]
    plot_mult(mults, FILES)

    print("\nQ7: evaluating threshold values on the files with errors")
    truth = create_dbg("ecoli_sample_perfect_reads.fasta.gz", K).keys()
    for f in FILES[2:]:
        print(f)
        evaluate_t(f, truth)

    print("\nQ9: recovering a unitig")
    kmer = "CGCTCTGTGTGACAAGCCGGAAACCGCCCAG"
    for f in FILES:
        print(f)
        unitig_sizes(f, kmer)

    print("\nQ11: sum of the size of the recovered unitigs")
    for f in FILES:
        print(f)
        unitig_totals(f)

if __name__ == "__main__":
    main()
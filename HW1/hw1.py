import sys
import random
import matplotlib.pyplot as plt

K_MIN = 2
K_MAX = 30


def read_fa(file):
    """Parse .fa text and return the list of contig sequences."""
    contigs = []
    for block in file.split(">"):
        lines = block.splitlines()
        if len(lines) > 1:
            seq = "".join(lines[1:])
            contigs.append(seq)
    return contigs


def compute_subword_complexity(seqs, k):
    """Compute the subword-complexity F(k) for a string or a list of contigs."""
    if isinstance(seqs, str):
        seqs = [seqs]

    distinct_substring = set()
    for seq in seqs:
        n = len(seq)
        for i in range(n - k + 1):
            distinct_substring.add(seq[i : i + k])
    return len(distinct_substring)


def compute_subword_complexities(dna_seq):
    """Compute F(k) for k ranging from K_MIN to K_MAX (inclusive)."""
    return [compute_subword_complexity(dna_seq, k) for k in range(K_MIN, K_MAX + 1)]


def seq_size(seqs):
    if isinstance(seqs, str):
        seqs = [seqs]

    return sum(len(seq) for seq in seqs)


def random_seq(size):
    return "".join(random.choices("ACGT", k=size))


def fibo_inc(word):
    neword = ""
    for letter in word:
        if letter == "A":
            neword += "AB"
        if letter == "B":
            neword += "A"
    return neword


def fibo_word(size):
    """Generate the Fibonacci word of the given length (size in params)."""
    word = "A"
    while len(word) < size:
        word = fibo_inc(word)
    return word[:size]


def main():
    if len(sys.argv) <2 :
        genome_path = "./genome_hw1.fa"
    else:
        genome_path = sys.argv[1]

    file = open(genome_path)
    genome = read_fa(file.read())
    Fk_genome = compute_subword_complexities(genome)
    genome_size = seq_size(genome)
    print("Genome size is : " + str(genome_size)) 

    rand_seq = random_seq(genome_size)
    Fk_random = compute_subword_complexities(rand_seq)

    fibonacci_seq = fibo_word(genome_size)
    Fk_fibo = compute_subword_complexities(fibonacci_seq)

    ks = list(range(K_MIN, K_MAX + 1))
    plt.plot(ks, Fk_genome, label="Streptococcus pneumoniae")
    plt.plot(ks, Fk_random, label="Random DNA sequence")
    plt.plot(ks, Fk_fibo, label="Fibonacci word")
    plt.yscale("log")
    plt.legend()
    plt.xlabel("k")
    plt.ylabel("F(k)")
    plt.title("Subword-complexity comparison")
    plt.savefig("./analysis.png")


if __name__ == "__main__":
    main()

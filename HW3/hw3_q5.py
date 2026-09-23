from hw3_q4 import tspToSCS
from hw3_q1 import greedySCS
from readfa import readfq
from simreads import simulate_reads_of_fixed_size
import time
import matplotlib.pyplot as plt


def main():
    path = "./genome.fa"
    file = open(path)
    seqs = []
    for _, seq, _ in readfq(file):
        seqs.append(seq)
    file.close()

    genome = "".join(seqs)

    try:
        n = int(input("Give a prefix length (by default : 1000) : "))
    except ValueError:
        n = 1000

    print("Choosen prefix size : {} bps".format(n))
    genome = genome[:n]

    S = simulate_reads_of_fixed_size(genome, 100, 42)

    s_Greedy = greedySCS(S)
    print(
        "Reconstructed sequence (length={}bps) using greedy algorithm :".format(
            len(s_Greedy)
        )
    )
    print(s_Greedy)

    s_TSPh = tspToSCS(S, 0)
    print(
        "Reconstructed sequence (length={}bps) using openTSP (heuristic):".format(
            len(s_TSPh)
        )
    )
    print(s_TSPh)

    # s_TSPe = tspToSCS(S, 1)
    # print("Reconstructed sequence (length={}bps) using openTSP (exact):".format(len(s_TSPe)))
    # print(s_TSPe)

    timeTabG = []
    timeTabH = []

    lengthTabG = []
    lengthTabH = []

    # Various reads size
    for r in range(50, 1000, 10):
        S = simulate_reads_of_fixed_size(genome, r, 42)

        startG = time.perf_counter()
        s_Greedy = greedySCS(S)
        endG = time.perf_counter()

        timeTabG.append(endG - startG)
        lengthTabG.append(len(s_Greedy))

        startH = time.perf_counter()
        s_TSPh = tspToSCS(S, 0)
        endH = time.perf_counter()

        timeTabH.append(endH - startH)
        lengthTabH.append(len(s_TSPh))

    x = list(range(50, 1000, 10))
    plt.plot(x, lengthTabH, label="Heuristic TSP")
    plt.plot(x, lengthTabG, label="Greedy")
    plt.legend()
    plt.xlabel("Reads length")
    plt.ylabel("Sequence length")
    plt.savefig("./variousReads_length.png")

    plt.clf()
    plt.plot(x, timeTabH, label="Heuristic TSP")
    plt.plot(x, timeTabG, label="Greedy")
    plt.legend()
    plt.xlabel("Reads length")
    plt.ylabel("Computation time")
    plt.savefig("./variousReads_time.png")


if __name__ == "__main__":
    main()

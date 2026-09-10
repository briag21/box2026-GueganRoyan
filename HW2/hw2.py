import matplotlib.pyplot as plt
import sys
from readfa import readfq

NUC_CODE = {"A" : 0, "T" : 1, "C" : 2, "G" : 3}

def make_canonical(kmer):
     """Compute the reverse-complement and return the lexicographically smaller sequence."""
     reverse_kmer = ""
     for letter in kmer:
          match letter:
               case "A": reverse_kmer += "T"
               case "T": reverse_kmer += "A"
               case "C": reverse_kmer += "G"
               case "G": reverse_kmer += "C"
     reverse_kmer = reverse_kmer[::-1]
     if kmer < reverse_kmer:
          canonical_kmer = kmer
     else:
          canonical_kmer = reverse_kmer
     return canonical_kmer

def compute_kmers(canonical:bool):
     """Extract 20-mer sets from all files, optionally converting them to canonical form."""
     sets = [set() for _ in range(6)]
     for i in range (6):
               path = "file{}.fa".format(i+1)
               file = open(path)
               nb_seq = 0
               bp_size = 0
               for _, seq, _ in readfq(file):
                    k=0
                    while len(seq) - (k+19) > 0:
                         kmer = seq[k:k+20].upper()
                         if not is_ambiguous(kmer):
                              if canonical:
                                   sets[i].add(make_canonical(kmer))
                              else:
                                   sets[i].add(kmer)
                         k += 1
               file.close()
     return sets

def jaccard(setA, setB):
     """Compute the Jaccard index between two sets of k-mers."""
     union = setA.union(setB)
     intersec = setA.intersection(setB)
     jaccard_index = len(intersec)/len(union)
     return jaccard_index

def kmer_to_numbers(nuc_string):
     """Encode a nucleotide string into a 2-bit per base integer."""
     new_kmer = 0
     for letter in nuc_string:
          new_kmer = (new_kmer << 2) | NUC_CODE[letter]
     return new_kmer

def sets_to_numbers(sets):
     """Convert k-mers in each set from strings to 2-bit encoded integers."""
     new_sets = []
     for i in range(len(sets)):
          new_set = set()
          for kmer in sets[i]:
               new_kmer = kmer_to_numbers(kmer)
               new_set.add(new_kmer)
          new_sets.append(new_set)

     return new_sets
                    
def is_ambiguous_lett(letter):
     """Return True if the nucleotide is ambiguous, False otherwise."""
     if (letter != "A") & (letter != "T") & (letter != "C") & (letter != "G"):
                    return True
     return False

def is_ambiguous(kmer):
    """Return True if the k-mer contains any ambiguous nucleotide, False otherwise."""
    for letter in kmer:
         if is_ambiguous_lett(letter):
              return True
    return False

def counting(kmers, file):
     """Count occurrences of canonical k-mers in the given file."""
     abundance = dict()
     for kmer in kmers:
          abundance[kmer] = 0
     for _, seq, _ in readfq(file):
          k = 0
          while len(seq) - (k +19) > 0:
               current = seq[k:k+20].upper()
               if not is_ambiguous(current):
                    current = kmer_to_numbers(make_canonical(current))
                    abundance[current] += 1
               k += 1
     return abundance

def main():
    # Measures
    for i in range (6):
        path = "file{}.fa".format(i+1)
        file = open(path)
        nb_seq = 0
        bp_size = 0
        ambiguous_char = 0
        for _, seq, _ in readfq(file):
            nb_seq += 1
            for elem in seq.upper():
                if is_ambiguous_lett(elem):
                     ambiguous_char += 1
                bp_size += 1
        print("File {} : {} sequence(s), cumulative length : {}. (Ambiguous chars : {})".format(i+1, nb_seq, bp_size, ambiguous_char))
        file.close()
    
    print("")

    #k-mers (k=20)
    sets = compute_kmers(False)

    with open("./kmerssets.txt", "w") as file:
         for elem in sets:
          for kmer in elem:
               file.write(kmer + " ")
          file.write("\n")
          
    ##Jaccard indexes
    
    print("Jaccard indexes : ")
    for i in range(5):
         for j in range (i+1, 6):
              jaccard_index = jaccard(sets[i], sets[j])
              print("J({},{}) = {}".format(i+1, j+1, jaccard_index))

    ## Encoding nucleotides in numbers

    sets = sets_to_numbers(sets)

    with open("./kmerssets_encoded.txt", "w") as file:
         for elem in sets:
              for kmer in elem:
                   file.write(str(kmer) + " ")
              file.write("\n")

     # Canonical encoded k-mers
    canonical_sets = compute_kmers(True)
    canonical_sets = sets_to_numbers(canonical_sets)

    with open("./canonical_kmerssets_encoded.txt", "w") as file:
          for elem in canonical_sets:
               for kmer in elem:
                    file.write(str(kmer) + " ")
               file.write("\n")

    print("")


     # Jaccard indexes with canonical k-mers
    print("Jaccard indexes (with canonical k-mers): ")
    for i in range(5):
          for j in range (i+1, 6):
               jaccard_index = jaccard(canonical_sets[i], canonical_sets[j])
               print("J({},{}) = {}".format(i+1, j+1, jaccard_index))


    # Abundance histogram of file 1
    if(len(sys.argv)>1):
          file_index = int(sys.argv[1])
    else:
         file_index = 1

    kmers1 = canonical_sets[file_index-1]
    file = open("./file{}.fa".format(file_index))
    abundance = counting(kmers1, file)
    file.close()    
    hist = dict()
    for count in abundance.values():
     hist[count] = hist.get(count, 0) + 1

    x = list(hist.keys())
    y = list(hist.values())

    plt.bar(x, y)
    plt.yscale('log')
    plt.xlabel("Abundance")
    plt.ylabel("Number of distinct k-mers")
    plt.title("20-mers abundance histogram of file{}.fa".format(file_index))
    plt.savefig("./abundance_hist.png")

if __name__ == "__main__":
    main()
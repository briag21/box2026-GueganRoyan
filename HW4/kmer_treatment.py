
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
 
 
def reverse_complement(kmer):
    return "".join(COMPLEMENT[letter] for letter in reversed(kmer))
 
 
def make_canonical(kmer):
    return min(kmer, reverse_complement(kmer))
from sympy import randprime

mu = {"A" : 1, "T" : 2, "C" : 3, "G" : 4}
SIZE = len(mu)

def pick_prime(pattern_len, text_len):
    return randprime(2, pattern_len * (text_len**2))

def hash_kr(seq, prime):
    hash = 0
    exp = len(seq) - 1
    for i in range(len(seq)):
        char = mu[seq[i]]
        hash += (char * pow(SIZE, exp, prime)) % prime
        exp -= 1 
    return hash % prime

def kr(seq, pattern):
    """LocateAll variant for Karp-Rabin pattern matching"""
    seqLen = len(seq)
    patLen = len(pattern)
    prime = pick_prime(patLen, seqLen)

    trueHash = hash_kr(pattern, prime)
    startingPositions = []

    for i in range(0, seqLen-patLen + 1):
        if(i==0):
            hash = hash_kr(seq[0:patLen], prime)
        else:
            hash = (SIZE * (hash - mu[seq[i-1]] * pow(SIZE, patLen-1, prime)) + mu[seq[i+patLen -1]])% prime

        if hash == trueHash:
            if seq[i:i + patLen] == pattern :
                startingPositions.append(i)

    return startingPositions

def compute_lps(pattern):
    patternLen = len(pattern)
    lps = len(pattern) * [0]

    i=1
    lpsLen=0
    while i<patternLen:
        if pattern[i] == pattern[lpsLen]:
            lpsLen += 1
            lps[i]=lpsLen
            i += 1
        else:
            if lpsLen != 0:
                lpsLen = lps[lpsLen - 1]
            else:
                lps[i] = 0
                i += 1
    return lps
        
        
def kmp(seq, pattern):
    """LocateAll variant for Knuth-Morris-Prath pattern matching"""
    seqLen = len(seq)
    patLen = len(pattern)
    lps = compute_lps(pattern)

    occurencies = []
    i = 0
    j = 0
    while(i < seqLen):
        if seq[i] == pattern[j]:
            i += 1
            j += 1
            if j == patLen:
                occurencies.append(i - patLen + 1)
                j = lps[j-1]
        else:
            if j!= 0:
                j = lps[j -1]
            else:
                i += 1

def main():
    nothing = 0
if __name__ == '__main__':
    main()

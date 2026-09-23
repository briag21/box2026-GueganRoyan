from readfa import readfq


def overlap(u, v):
    for i in range(min(len(u), len(v)), 0, -1):
        if u[-i:] == v[:i]:
            return i
    return 0


def merge(u, v):
    return u + v[overlap(u, v) :]


def getStringsWithMaximalOverlap(S):
    u, v = None, None
    best = -1
    for i in range(len(S)):
        for j in range(len(S)):
            if i != j:
                ov = overlap(S[i], S[j])
                if ov > best:
                    best = ov
                    u, v = S[i], S[j]
    return u, v


def greedySCS(S):
    S = list(S)
    while len(S) > 1:
        u, v = getStringsWithMaximalOverlap(S)
        w = merge(u, v)
        S.remove(u)
        S.remove(v)
        S.append(w)
    return S[0]


def main():
    tests = {"abc", "bcde"}

    word = greedySCS(tests)
    print(word)


if __name__ == "__main__":
    main()

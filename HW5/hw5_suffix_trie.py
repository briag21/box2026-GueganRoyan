from readfa import readfq
from hw5_pattern_matching import kr
from collections import deque
def main():
    seq = read("./sHDAg.fa")
    suffixTrieHepatitis = SuffixTrie(seq)
    lrf = suffixTrieHepatitis.longest_repeated_factor()
    print("Longest repeated factor :" + lrf)
    position = kr(seq, lrf)
    print("Position of the longest repeated factor: " + position)

    seq = read("./hiv1.fa")
    suffixTrieHIV = SuffixTrie(seq)
    answer = shortest_discriminant(suffixTrieHepatitis, suffixTrieHIV)
    print("Shortest discriminant for HIV and Hepatitis: " + answer)

def read(path):
    with open(path, "r") as fasta:
        _, seq, _ = next(readfq(fasta))
    return seq

class SuffixTrie:
    class Node:
        def __init__(self, letter: str = None):
            self.__letter = letter
            self.__childs = []
            self.__marked = False

        def gotChild(self, letter: str) -> bool:
            for node in self.__childs:
                if node.__letter == letter:
                    return True
            return False

        def child(self, letter: str):
            for child_node in self.__childs:
                if child_node.__letter == letter:
                    return child_node
            return None

        def mark(self):
            self.__marked = True

        @property
        def childs(self):
            return self.__childs

        @property
        def letter(self):
            return self.__letter

    def __init__(self, S: str):
        S = S + "$"
        self.__root = self.Node(None)
        
        suffixes = [S[i:] for i in range(len(S))]

        for suffix in suffixes:
            node = self.__root
            for i in range(len(suffix)):
                if not node.gotChild(suffix[i]):
                    newNode = self.Node(suffix[i])
                    node.childs.append(newNode)
                node = node.child(suffix[i])
            node.mark()

    @property
    def root(self):
        return self.__root

    def longest_repeated_factor(self) -> str:
        longest_factor = ""

        def dfs(node, current_prefix : str):
            nonlocal longest_factor
            if len(node.childs) > 1:
                clean_prefix = current_prefix.replace("$", "")
                if len(clean_prefix) > len(longest_factor):
                    longest_factor = clean_prefix

            for child in node.childs:
                dfs(child, current_prefix + child.letter)

        dfs(self.__root, "")
        return longest_factor

def shortest_discriminant(ST1 : SuffixTrie, ST2 : SuffixTrie):
    def bfs(root1 : SuffixTrie.Node, root2 : SuffixTrie.Node):

        queue = deque([(root1, root2, "")])

        while queue:
            n1, n2, prefix = queue.popleft()
            children1 = {child.letter: child for child in n1.childs}
            children2 = {child.letter: child for child in n2.childs}

            for char, child1 in children1.items():
                if char == "$":
                    continue
                new_prefix = prefix + char
                if char not in children2:
                    return new_prefix
                else:
                    queue.append((child1, children2[char], new_prefix))

            for char, _ in children2.items():
                if char == "$":
                    continue
                new_prefix = prefix + char
                if char not in children1:
                    return new_prefix
    return bfs(ST1.root, ST2.root)

            




   

    
                    

if __name__ == "__main__":
    main()


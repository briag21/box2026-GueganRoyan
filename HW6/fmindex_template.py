# >>===================================================<<
# ||+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+||
# ||--------------------|B|I|F|-----------------------|||
# |||B|u|r|r|o|w|s|-|W|h|e|e|l|e|r|-|T|r|a|n|s|f|o|r|m|||
# ||+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+||
# >>===================================================<<

# Template by Leo Ackermann (2025)

from ks import simple_kark_sort

# >>==========================================================================<<
# >>                        Textbook FM-index class                           <<
# >>==========================================================================<<

class TextbookFMindex:

    # >>===[ Class constructor ]================================================

    def __init__(self, seq, verbose=False):
        '''
        Create a new FMindex (textbook version)
        '''

        # Fields declaration
        self.bwt      = None
        self.sa       = None
        self.countmap = None
        self.ranks    = None
        self.next_smallest_letter = None
        self.total_count = None

        # Fields initialization
        # TODO
        augmentedSeq = seq + "$"
        self.sa = self._compute_sa(augmentedSeq)
        self.bwt = self._compute_bwt(augmentedSeq)
        self.ranks = self._compute_ranks()
        self.countmap = self._compute_countmap()
        self.next_smallest_letter = self._compute_next_smallest_letter()




    # >>===[ Pretty printer ]===================================================

    def pretty_print(self):
        '''
        Pretty printer for the FM index
        '''

        # This function take a list ((name_i, array_i))_i as input, and print
        # them vertically, side by side, using name_i's as headers
        # NOTE. We assume here that {header}\t and {content}\t reach the same
        # column
        def print_arrays_vertically(entitled_arrays, tabsize=4):

            # Prettier printing of None
            def xstr(s):
                return '·' if s is None else str(s)

            if entitled_arrays == []:
                raise Exception("<entitled_arrays> cannot be empty")

            header = "#\t|\t"
            separator = '-' * len(f"#\t".expandtabs(tabsize)) \
                + '+' + '-' * (len(f"+\t".expandtabs(tabsize))-1)
            for (name, array) in entitled_arrays:
                header += f"{name}\t"
                separator += '-' * len(f"{name}\t".expandtabs(tabsize))
            print(header)
            print(separator)

            len_arrays = (len(entitled_arrays[0][1]))
            for i in range(len_arrays):
                row = f"{i}\t|\t"
                for (name, array) in entitled_arrays:
                    row += f"{xstr(array[i])}\t"
                print(row)

        print()
        print_arrays_vertically([])
        print()



    # >>===[ Attribute initialisation functions ]===============================

    def _compute_sa(self, seq):
        '''
        Compute the suffix array of the sequence (with termination symbol), in
        linear time
        '''

        return simple_kark_sort(seq)



    def _compute_bwt(self, seq):
        '''
        Compute the Burrows-Wheeler transform using the SA-based definition, in
        linear time
        '''

        n = len(self.sa)
        bwt = [""]*n

        i = 0
        for index in self.sa:
            bwt[i] = seq[index - 1]
            i += 1

        return "".join(bwt)



    def _compute_ranks(self):
        '''
        Compute the ranks arrays while scanning the Burrows-Wheeler transform
        '''

        ranks = {}
        last = {}
        for i in range(len(self.bwt)):
            element = self.bwt[i]
            if element not in ranks:
                ranks[element] = [0]*len(self.bwt)
                last[element] = 0
            ranks[element][i] = last[element] + 1
            last[element] += 1

        self.total_count = last
        return ranks



    def _compute_countmap(self):
        '''
        Compute the countmap array from the ranks arrays arrays
        '''
        characters = sorted(self.ranks.keys())
        countmap = {}
        count = 0
        for element in characters:
            countmap[element] = count
            count += self.total_count[element]

        return countmap



    def _compute_next_smallest_letter(self):
        '''
        Compute the map of <Sigma> -> <Sigma>+{None} that maps a letter to the
        next one in the lexicographic order, using countmap
        '''

        output = {}
        characters = sorted(self.ranks.keys())

        for i in range(len(characters)):
            if i == (len(characters) - 1):
                output[characters[i]] = None
            else:
                output[characters[i]] = characters[i+1]

        return output



    # >>===[ LF-mapping ]=======================================================

    def _lf_mapping(self, i):
        '''
        The LF-mapping function, with O(1) lookup time
        '''
        char = self.bwt[i]
        occurencies = self.ranks[char][i]
        charCount = self.countmap[char]
            
        return charCount + occurencies - 1




    # >>===[ Burrows-Wheeler transform inversion ]==============================

    def get_string(self):
        '''
        Retrieve the string encoded within the BWT, in linear time
        '''
        output = ""
        index = 0
        for i in range(len(self.bwt)):
            output = self.bwt[index] + output
            index = self._lf_mapping(index)
        return output[1:]

    def get_string__naive(self):
        '''
        Retrieve the string encoded within the BWT, naive version in quadratic time
        '''
        n  = len(self.sa)
        current = [""] * n
        for i in range (n):
            for j in range (n):
                current[j] = self.bwt[j] + current[j]
            current.sort()

        for row in current:
            if row.endswith("$"):
                original = row
    
        return original[:-1]

    # >>===[ Pattern matching functions ]=======================================

    def _get_pattern_interval(self, pattern):
        '''
        Retrieve the F-interval (i_min, i_max) where a pattern lies. If no such
        interval exists, returns (None, None).
        '''

        i_min = 0
        i_max = len(self.bwt) - 1

        for char in reversed(pattern):
            if char not in self.ranks:
                return (None , None)
            else:
                if i_min > 0:
                    i_min = self.countmap[char] + self.ranks[char][i_min - 1]
                else:
                    i_min = self.countmap[char]
                i_max = self.countmap[char] + self.ranks[char][i_max] - 1
                if i_min > i_max:
                    return (None, None)

        return (i_min, i_max)




    def membership(self, pattern):
        '''
        FM-index based membership pattern-matching
        '''

        (i_min, i_max) = self._get_pattern_interval(pattern)

        if (i_min, i_max) == (None, None):
            return False
        else:
            return True



    def count(self, pattern):
        '''
        FM-index based count pattern-matching
        '''

        (i_min, i_max) = self._get_pattern_interval(pattern)

        if (i_min, i_max) == (None, None):
            return 0
        
        return i_max - i_min + 1



    def locate(self, pattern):
        '''
        FM-index based locate pattern-matching
        '''

        (i_min, i_max) = self._get_pattern_interval(pattern)

        if (i_min, i_max) == (None, None):
            return []

        positions = []
        for i in range (i_min, i_max + 1):
            positions.append(self.sa[i])

        return positions




# >>==========================================================================<<
# >>                              FM-index class                              <<
# >>==========================================================================<<

class FMindex(TextbookFMindex):

    # >>===[ Class constructor ]================================================

    def __init__(self, seq, verbose=False):
        '''
        Create a new FMindex
        '''
        if verbose:
            print("-> FM-index: Init started")

        # Inherit fields and initialisation from the Textbook FMindex class
        super().__init__(seq, verbose)

        # Extra fields declaration
        self.fm_count = None
        self.fm_rank = None
        self.fm_ranks = None
        self.next_smallest_letter = None

        # Extra fields initialization
        self.fm_count = self.countmap
        self.next_smallest_letter = self.next_smallest_letter

        if verbose:
            print("-> FM-index: Subsampling ranks & SA")

        self.fm_ranks = self._subsample_ranks()
        self.sa = self._subsample_sa()

        if verbose:
            print("-> FM-index: Done")




    # >>===[ Subsampling and on-the-fly retrieval functions ]===================

    def _subsample_ranks(self):
        '''
        Subsample the ranks arrays, setting most entries to <None>
        '''
        subranks = {}

        for char in self.ranks.keys():
            tab = self.ranks[char]
            i=0
            subranks[char] = []
            for element in tab:
                if i % 32 == 0:
                    subranks[char].append(element)
                else:
                    subranks[char].append(None)
                i+=1

        return subranks


    def _get_ranks_entry(self, letter, i):
        '''
        Equivalent to TextbookFMindex's <self.ranks[letter][i]>
        '''

        if letter not in self.fm_ranks:
            return 0

        idx = i
        while idx >= 0 and self.fm_ranks[letter][idx] is None:
            idx -= 1

        if idx >= 0:
            rank = 0
            start_scan = 0
        else:
            rank = self.fm_ranks[letter][idx]
            start_scan = idx + 1

        for j in range(start_scan, i + 1):
            if(self.bwt[j] == letter):
                rank += 1

        return rank



    def _subsample_sa(self):
        '''
        Subsample the suffix array, setting most entries to <None>
        '''
        sub_sa = []
        i = 0

        for val in self.sa:
            if i % 32 == 0:
                sub_sa.append(val)
            else:
                sub_sa.append(None)
            i += 1
        return sub_sa



    def _get_sa_entry(self, i):
        '''
        Equivalent to TextbookFMindex's <self.sa[i]>
        '''

        count = 0
        pass #TODO !


    # >>===[ OVERRIDE functions calling <self.sa> or <self.ranks> ]=============

    # TODO


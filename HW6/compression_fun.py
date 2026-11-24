# >>==========================================================================<<
# >>                       (De)compression functions                          <<
# >>==========================================================================<<

def run_length_encode(seq):
    '''
    TODO Docstring
    '''
    output = []
    for item in seq:
        if output and output[-1][0] == item:
            output[-1] = (item, output[-1][1] + 1)
        else:
            output.append((item, 1))
    return output

def run_length_decode(pairs):
    '''
    TODO Docstring
    '''
    output = []
    for item, count in pairs:
        output.extend([item] * count)
    return output

MTF_ALPHABET = [chr(i) for i in range(256)]

def move_to_front_encode(seq):
    '''
    TODO Docstring
    '''
    alphabet = MTF_ALPHABET.copy()
    output = []

    for element in seq:
        index = alphabet.index(element)
        output.append(index)
        alphabet.pop(index)
        alphabet.insert(0, element)

    return output

def move_to_front_decode(encoded_seq):
    '''
    TODO Docstring
    '''
    alphabet = MTF_ALPHABET.copy()
    output = []

    if encoded_seq == "":
        return ""

    for code in encoded_seq:
        char = alphabet.pop(code)
        alphabet.insert(0, char)
        output.append(char)

    return output

def combined_encode(seq):
    '''
    TODO Docstring
    '''
    return run_length_encode(move_to_front_encode(seq))

def combined_decode(seq):
    '''
    TODO Docstring
    '''
    return "".join(move_to_front_decode(run_length_decode(seq)))
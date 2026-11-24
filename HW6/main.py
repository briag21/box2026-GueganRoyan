from fmindex_template import *
from compression_fun import *
from sys import getsizeof
import matplotlib.pyplot as plt

PATHS = ["./dickens.txt", "dudh.txt", "ecoli.fa", "random.txt"]

def openFiles():
    strings= {}
    for path in PATHS:
        strings[path] = ""
        with open(path) as f:
            for line in f:
                if not line.startswith(">"):
                    strings[path] += line
    return strings

def bwtEach(strings):
    bwts = {}
    for path in PATHS:
        newClass = TextbookFMindex(strings[path])
        bwts[path] = newClass.bwt
    return bwts
def compressEach(bwts):
    compressedString = {}
    for path in PATHS:
        compressedString[f"rl_{path}_encoded"] = run_length_encode(bwts[path])
        compressedString[f"mtf_{path}_encoded"] = move_to_front_encode(bwts[path])
        compressedString[f"combined_{path}_encoded"] = combined_encode(bwts[path])
    return compressedString

def evaluate_compression(strings, compressed_strings):
  paths = list(strings.keys())
  metrics = ['Original', 'RLE', 'MTF', 'Combined']

  bpc_data = {m: [] for m in metrics}

  for path, orig_text in strings.items():
    orig_len = len(orig_text)
    if orig_len == 0:
      for m in metrics:
        bpc_data[m].append(0)
      continue

    bpc_data['Original'].append(8.0)
    bpc_data['RLE'].append(
        (len(compressed_strings[f'rl_{path}_encoded']) * 3 * 8) / orig_len
    )
    bpc_data['MTF'].append(
        (len(compressed_strings[f'mtf_{path}_encoded']) * 1 * 8) / orig_len
    )
    bpc_data['Combined'].append(
        (len(compressed_strings[f'combined_{path}_encoded']) * 3 * 8) / orig_len
    )

  x = range(len(paths))
  width = 0.2

  fig, ax = plt.subplots(figsize=(10, 5))
  for i, method in enumerate(metrics):
    ax.bar(
        [pos + i * width for pos in x],
        bpc_data[method],
        width,
        label=method,
    )

  ax.set_ylabel('Bits per Character (BPC)')
  ax.set_xticks([pos + width * 1.5 for pos in x])
  ax.set_xticklabels(paths, rotation=15)
  ax.legend()
  ax.axhline(8, color='red', linestyle='--', linewidth=1, label='Baseline (8 BPC)')

  plt.tight_layout()
  plt.show()

def main():
    textIndex = TextbookFMindex("BANANA")
    print(textIndex.bwt) # should output ANNB$AA 
    print(textIndex.get_string__naive())
    print(textIndex.get_string())

    text = "banaaaaaana"
    print(run_length_encode(text))
    print(move_to_front_encode(text))
    print(combined_encode(text))

    print(combined_decode(combined_encode(text)))
    assert(text == combined_decode(combined_encode(text)))
    strings = openFiles()
    bwts = bwtEach(strings)
    compressedStrings = compressEach(bwts)
    evaluate_compression(strings, compressedStrings)


if __name__ == "__main__":
    main()
#Convolution with one spectrum

import sys
from collections import Counter

def spectral_convolution(spectrum):
    differences = []
    for i in range(len(spectrum)):
        for j in range(len(spectrum)):
            diff = spectrum[j] - spectrum[i]
            if diff > 0:
                differences.append(diff)
    return differences

if __name__ == "__main__":
    input_file = sys.argv[1] if len(sys.argv) > 1 else "bioinfo_genome_sequencing/datasets/dataset_21.txt"
    spectrum = list(map(int, open(input_file).read().split()))
    result = spectral_convolution(spectrum)
    print(" ".join(map(str, result)))
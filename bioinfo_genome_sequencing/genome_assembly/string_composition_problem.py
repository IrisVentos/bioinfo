#reconstructing strings from k mers
#Given a string Text, its k-mer composition Compositionk(Text) is the collection 
#of all k-mer substrings of Text (including repeated k-mers). For example,
# Composition3(TATGGGGTGC) = {ATG, GGG, GGG, GGT, GTG, TAT, TGC, TGG}.

#hypothesis = full coverage of genome (ideal situation) + order of k-mer doesn't matter yet

def composition_k(k, text):
    """
    Find all k-mers in the given text.
    
    Args:
        k: Length of each k-mer
        text: The DNA string to analyze
    
    Returns:
        List of k-mers
    """
    kmers = []
    
    # Slide a window of size k across the text
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        kmers.append(kmer)
    
    return kmers


def main():
    # Read input from file
    with open('datasets/dataset_1.txt', 'r') as f:
        lines = f.read().strip().split('\n')
    
    # Parse input
    k = int(lines[0])
    text = lines[1].strip()
    
    # Get k-mers
    kmers = composition_k(k, text)
    
    # Output result (space-separated)
    result = ' '.join(kmers)
    print(result)
    
    # Optionally write to output file
    with open('datasets/output_1.txt', 'w') as f:
        f.write(result)


if __name__ == "__main__":
    main()
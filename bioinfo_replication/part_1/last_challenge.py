def FindClumps(Genome, k, L, t):
    """
    Find all k-mers forming (L, t)-clumps in Genome
    
    Args:
        Genome (str): Input genome string
        k (int): Length of k-mers
        L (int): Window length
        t (int): Minimum frequency threshold
    
    Returns:
        list: All distinct k-mers forming (L, t)-clumps
    """
    clumps = set()
    
    # Slide window of length L across the genome
    for i in range(len(Genome) - L + 1):
        window = Genome[i:i + L]
        
        # Count k-mers in this window
        kmer_count = {}
        for j in range(len(window) - k + 1):
            kmer = window[j:j + k]
            kmer_count[kmer] = kmer_count.get(kmer, 0) + 1
        
        # Find k-mers that appear at least t times
        for kmer, count in kmer_count.items():
            if count >= t:
                clumps.add(kmer)
    
    return list(clumps)

# Solve from file
with open('data/last_sample.txt', 'r') as file:
    lines = file.read().strip().split('\n')
    
    genome = lines[0].strip()
    k = int(lines[1].strip())
    L = int(lines[2].strip())
    t = int(lines[3].strip())
    
    result = FindClumps(genome, k, L, t)
    result_str = " ".join(result)
    
    with open('output_5.txt', 'w') as output:
        output.write(result_str)
    
    print(f"Found {len(result)} clumps")
    print(f"Result: {result_str}")
    print("Result saved to output_5.txt")
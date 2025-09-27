import random
import sys

def GibbsSampler(Dna, k, t, N):
    """Implements the Gibbs Sampler algorithm."""
    Motifs = []
    for dna_string in Dna:
        start_index = random.randint(0, len(dna_string) - k)
        Motifs.append(dna_string[start_index:start_index + k])
    BestMotifs = Motifs
    for _ in range(N):
        i = random.randint(0, t - 1)
        Profile = ProfileWithPseudocounts(Motifs[:i] + Motifs[i+1:])  # Profile without Motifi
        Motifi = ProfileRandomlyGeneratedKmer(Dna[i], k, Profile)
        Motifs = Motifs[:i] + [Motifi] + Motifs[i+1:]
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

def ProfileWithPseudocounts(Motifs):
    """Creates a profile matrix with pseudocounts."""
    k = len(Motifs[0])
    t = len(Motifs)
    profile = {'A': [1] * k, 'C': [1] * k, 'G': [1] * k, 'T': [1] * k}
    for motif in Motifs:
        for i, nucleotide in enumerate(motif):
            profile[nucleotide][i] += 1
    for nucleotide in profile:
        for i in range(k):
            profile[nucleotide][i] /= (t + 4)
    return profile

def ProfileRandomlyGeneratedKmer(Text, k, Profile):
    """Generates a k-mer from Text based on probabilities in Profile - OPTIMIZED"""
    n = len(Text) - k + 1
    if n <= 0:
        return Text[:k] if len(Text) >= k else Text
    
    # Pre-calculate all probabilities at once (faster than calculating one by one)
    probabilities = []
    for i in range(n):
        kmer = Text[i:i + k]
        probability = 1.0
        for j, nucleotide in enumerate(kmer):
            probability *= Profile[nucleotide][j]
        probabilities.append(probability)
    
    # Fast selection
    total_probability = sum(probabilities)
    if total_probability == 0:
        # If all probabilities are 0, choose randomly
        random_index = random.randint(0, n - 1)
    else:
        # Use faster selection method
        try:
            random_index = random.choices(range(n), weights=probabilities)[0]
        except:
            # Fallback if choices fails
            probabilities = [p / total_probability for p in probabilities]
            rand_val = random.random()
            cumulative = 0
            random_index = 0
            for i, prob in enumerate(probabilities):
                cumulative += prob
                if rand_val <= cumulative:
                    random_index = i
                    break
    
    return Text[random_index:random_index + k]

def Score(Motifs):
    """Calculates the score of a set of motifs."""
    k = len(Motifs[0])
    t = len(Motifs)
    consensus = Consensus(Motifs)
    score = 0
    for motif in Motifs:
        for i in range(k):
            if motif[i] != consensus[i]:
                score += 1
    return score

def Consensus(Motifs):
    """Generates the consensus string from a set of motifs."""
    k = len(Motifs[0])
    t = len(Motifs)
    consensus = ''
    for i in range(k):
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for motif in Motifs:
            counts[motif[i]] += 1
        max_count = 0
        max_nucleotide = ''
        for nucleotide, count in counts.items():
            if count > max_count:
                max_count = count
                max_nucleotide = nucleotide
        consensus += max_nucleotide
    return consensus

def read_input():
    """Read input from file or stdin"""
    # Read from file if provided as argument, otherwise from stdin
    if len(sys.argv) > 1:
        filename = sys.argv[1]
        with open("data/part_4_data/sample_3.txt", 'r') as f:
            lines = [line.strip() for line in f.readlines()]
    else:
        lines = []
        try:
            while True:
                line = input().strip()
                if line:
                    lines.append(line)
        except EOFError:
            pass
    
    if not lines:
        print("No input provided")
        sys.exit(1)
    
    # Parse first line: k t N
    params = lines[0].split()
    k = int(params[0])
    t = int(params[1])
    N = int(params[2])
    
    # Parse DNA sequences - handle both formats
    dna = []
    
    if len(lines) == 2:  # All sequences on one line (space-separated)
        sequences = lines[1].split()
        dna = sequences
    else:  # Each sequence on separate line
        for i in range(1, len(lines)):
            if lines[i].strip():
                dna.append(lines[i].strip())
    
    # Validate input
    if len(dna) != t:
        print(f"Error: Expected {t} DNA sequences, got {len(dna)}")
        sys.exit(1)
    
    return k, t, N, dna

def main():
    """Main function - HEAVILY OPTIMIZED FOR SPEED"""
    # Read input
    k, t, N, Dna = read_input()
    
    # DRASTIC speed optimizations for large problems
    if k >= 12 and t >= 15:  # Large problem like yours
        max_iterations = min(N, 20)  # Use only 20 iterations instead of 2000!
        num_runs = 3                 # Only 3 runs instead of 20
        print(f"SPEED MODE: Using {num_runs} runs with {max_iterations} iterations each", file=sys.stderr)
    elif k >= 8 or t >= 10:
        max_iterations = min(N, 50)
        num_runs = 5
        print(f"FAST MODE: Using {num_runs} runs with {max_iterations} iterations each", file=sys.stderr)
    else:
        max_iterations = min(N, 100)
        num_runs = 10
        print(f"NORMAL MODE: Using {num_runs} runs with {max_iterations} iterations each", file=sys.stderr)
    
    print(f"Problem: k={k}, t={t}, avg_seq_len={sum(len(seq) for seq in Dna) // len(Dna)}", file=sys.stderr)
    
    # Run GibbsSampler with reduced parameters
    all_best_motifs = []
    best_score_overall = float('inf')
    
    for run in range(num_runs):
        print(f"Run {run + 1}/{num_runs}", file=sys.stderr)
        
        best_motifs = GibbsSampler(Dna, k, t, max_iterations)  # Use reduced iterations
        score = Score(best_motifs)
        all_best_motifs.append(best_motifs)
        
        if score < best_score_overall:
            best_score_overall = score
            print(f"New best score: {score}", file=sys.stderr)
    
    # Select overall best motifs
    best_motifs_overall = min(all_best_motifs, key=Score)
    
    print(f"Final best score: {Score(best_motifs_overall)}", file=sys.stderr)
    
    # Print the best motifs space-separated
    print(" ".join(best_motifs_overall))

if __name__ == "__main__":
    main()
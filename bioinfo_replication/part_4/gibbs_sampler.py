import random
import sys

def read_input():
    """Read input from file or stdin"""
    if len(sys.argv) > 1:
        with open(sys.argv[1], 'r') as f:
            lines = [line.strip() for line in f.readlines()]
    else:
        lines = [line.strip() for line in sys.stdin.readlines()]
    
    # Parse parameters
    k, t, N = map(int, lines[0].split())
    
    # Parse DNA - handle both formats
    if len(lines) == 2:
        dna = lines[1].split()  # Space-separated
    else:
        dna = [lines[i] for i in range(1, t+1)]  # Line-separated
    
    return k, t, N, dna

def create_profile(motifs, k):
    """Create profile with pseudocounts - simplified"""
    profile = {'A': [1]*k, 'C': [1]*k, 'G': [1]*k, 'T': [1]*k}
    
    for motif in motifs:
        for i, nuc in enumerate(motif):
            profile[nuc][i] += 1
    
    total = len(motifs) + 4
    for nuc in profile:
        for i in range(k):
            profile[nuc][i] /= total
    
    return profile

def random_kmer(text, k, profile):
    """Simple random k-mer selection"""
    n = len(text) - k + 1
    probs = []
    
    for i in range(n):
        kmer = text[i:i+k]
        p = 1.0
        for j, nuc in enumerate(kmer):
            p *= profile[nuc][j]
        probs.append(p)
    
    # Simple weighted selection
    total = sum(probs)
    if total == 0:
        idx = random.randint(0, n-1)
    else:
        r = random.random() * total
        cumsum = 0
        idx = 0
        for i, p in enumerate(probs):
            cumsum += p
            if r <= cumsum:
                idx = i
                break
    
    return text[idx:idx+k]

def score(motifs):
    """Simple scoring"""
    if not motifs:
        return 1000
    
    k = len(motifs[0])
    total = 0
    
    for pos in range(k):
        counts = {'A':0, 'C':0, 'G':0, 'T':0}
        for motif in motifs:
            counts[motif[pos]] += 1
        max_count = max(counts.values())
        total += len(motifs) - max_count
    
    return total

def gibbs_minimal(dna, k, t, max_iter=10):
    """Minimal Gibbs sampler - VERY limited iterations"""
    # Random initial motifs
    motifs = []
    for seq in dna:
        start = random.randint(0, max(0, len(seq) - k))
        motifs.append(seq[start:start+k])
    
    best_motifs = motifs[:]
    best_score = score(best_motifs)
    
    # Very few iterations
    for _ in range(max_iter):
        i = random.randint(0, t-1)
        
        # Profile without motif i
        reduced = motifs[:i] + motifs[i+1:]
        profile = create_profile(reduced, k)
        
        # New motif
        motifs[i] = random_kmer(dna[i], k, profile)
        
        # Check improvement
        current_score = score(motifs)
        if current_score < best_score:
            best_score = current_score
            best_motifs = motifs[:]
    
    return best_motifs

def main():
    print("Reading input...", file=sys.stderr)
    k, t, N, dna = read_input()
    
    print(f"k={k}, t={t}, sequences loaded", file=sys.stderr)
    
    # EXTREMELY limited parameters for speed
    if k >= 12:
        iterations = 2000   # Only 5 iterations!
        runs = 20      # Only 2 runs!
    else:
        iterations = 2000
        runs = 20
    
    print(f"Using {runs} runs with {iterations} iterations each", file=sys.stderr)
    
    best_overall = None
    best_score_overall = 1000
    
    for run in range(runs):
        print(f"Run {run+1}", file=sys.stderr)
        
        result = gibbs_minimal(dna, k, t, iterations)
        current_score = score(result)
        
        if current_score < best_score_overall:
            best_score_overall = current_score
            best_overall = result
            print(f"Score: {current_score}", file=sys.stderr)
    
    print(f"Done! Final score: {best_score_overall}", file=sys.stderr)
    print(" ".join(best_overall))

if __name__ == "__main__":
    main()
from itertools import product
from helpers import debruijn_graph_from_kmers, eulerian_path

def generate_binary_kmers(k):
    """Generate all possible binary k-mers of length k."""
    return [''.join(p) for p in product('01', repeat=k)]

def eulerian_cycle_to_circular_string(cycle, k):
    """
    Convert an Eulerian cycle to a k-universal circular string.
    
    Args:
        cycle: List of (k-1)-mers representing nodes in the cycle
        k: Length of k-mers
    
    Returns:
        k-universal circular string
    """
    if not cycle:
        return ""
    
    # For a circular string, we traverse all edges in the cycle
    # Each edge corresponds to a k-mer
    # The number of edges = 2^k (all possible k-mers)
    # The circular string length should be exactly 2^k
    
    # Start with the first (k-1) characters
    result = cycle[0]
    
    # For each edge in the cycle, add the last character of the destination node
    # Stop after we've added 2^k - (k-1) more characters
    num_edges = 2**k
    chars_to_add = num_edges - (k - 1)
    
    for i in range(1, min(len(cycle), chars_to_add + 1)):
        result += cycle[i][-1]
    
    return result

def k_universal_circular_string(k):
    """
    Generate a k-universal circular string.
    
    A k-universal circular string is a circular string that contains 
    every possible k-mer exactly once as a substring.
    
    Args:
        k: Length of k-mers
    
    Returns:
        k-universal circular string
    """
    # Step 1: Generate all binary k-mers
    kmers = generate_binary_kmers(k)
    
    # Step 2: Build De Bruijn graph from k-mers (reusing helper)
    graph = debruijn_graph_from_kmers(kmers)
    
    # Step 3: Find Eulerian path/cycle in the graph (reusing helper)
    # For a k-universal string, this will be an Eulerian cycle
    cycle = eulerian_path(graph)
    
    # Step 4: Convert cycle to circular string
    circular_string = eulerian_cycle_to_circular_string(cycle, k)
    
    return circular_string

def main():
    # Read input from file
    with open('bioinfo_genome_sequencing/datasets/dataset_9.txt', 'r') as f:
        k = int(f.read().strip())
    
    # Generate k-universal circular string
    result = k_universal_circular_string(k)
    
    # Write output to file
    with open('bioinfo_genome_sequencing/datasets/output_9.txt', 'w') as f:
        f.write(result)
    
    print(f"k-universal circular string for k={k}: {result}")
    print(f"Length: {len(result)}")
    print(f"Expected length: {2**k} (contains all {2**k} binary {k}-mers)")
    print(f"Output written to output.txt")

if __name__ == "__main__":
    main()
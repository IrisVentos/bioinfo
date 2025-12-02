from collections import defaultdict

def build_paired_debruijn_graph(paired_reads):
    """
    Build a De Bruijn graph from paired k-mers.
    
    Args:
        paired_reads: List of tuples (first_kmer, second_kmer)
    
    Returns:
        Dictionary representing adjacency list where nodes are paired (k-1)-mers
    """
    graph = defaultdict(list)
    
    for first, second in paired_reads:
        # Prefix: first k-1 chars of both k-mers
        prefix = (first[:-1], second[:-1])
        # Suffix: last k-1 chars of both k-mers
        suffix = (first[1:], second[1:])
        
        graph[prefix].append(suffix)
    
    return graph

def find_eulerian_path_paired(graph):
    """Find Eulerian path in paired De Bruijn graph."""
    # Calculate degrees
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    
    for node in graph:
        out_degree[node] = len(graph[node])
        for neighbor in graph[node]:
            in_degree[neighbor] += 1
    
    # Find start node
    start_node = None
    all_nodes = set(graph.keys()) | set(in_degree.keys())
    
    for node in all_nodes:
        if out_degree[node] - in_degree[node] == 1:
            start_node = node
            break
    
    if start_node is None:
        start_node = next(iter(graph.keys()))
    
    # Build path using Hierholzer's algorithm
    edges = defaultdict(list)
    for node in graph:
        edges[node] = graph[node][:]
    
    stack = [start_node]
    path = []
    
    while stack:
        current = stack[-1]
        if edges[current]:
            next_node = edges[current].pop(0)
            stack.append(next_node)
        else:
            path.append(stack.pop())
    
    path.reverse()
    return path

def string_spelled_by_gapped_patterns(gapped_patterns, k, d):
    """
    Reconstruct string from ordered gapped patterns.
    
    Args:
        gapped_patterns: Ordered list of (k,d)-mer tuples
        k: Length of k-mers
        d: Gap distance
    
    Returns:
        Reconstructed string
    """
    first_patterns = [pattern[0] for pattern in gapped_patterns]
    second_patterns = [pattern[1] for pattern in gapped_patterns]
    
    # Build prefix string from first patterns
    prefix_string = first_patterns[0]
    for i in range(1, len(first_patterns)):
        prefix_string += first_patterns[i][-1]
    
    # Build suffix string from second patterns
    suffix_string = second_patterns[0]
    for i in range(1, len(second_patterns)):
        suffix_string += second_patterns[i][-1]
    
    # They should overlap at position k+d
    # Verify overlap and combine
    for i in range(k + d, len(prefix_string)):
        if prefix_string[i] != suffix_string[i - k - d]:
            return None
    
    # Return the combined string
    return prefix_string + suffix_string[len(prefix_string) - k - d:]

def string_reconstruction_from_read_pairs(k, d, paired_reads):
    """
    Reconstruct string from collection of paired k-mers.
    
    Args:
        k: Length of k-mers
        d: Gap distance
        paired_reads: Collection of (k,d)-mer tuples
    
    Returns:
        Reconstructed string
    """
    # Build paired De Bruijn graph
    graph = build_paired_debruijn_graph(paired_reads)
    
    # Find Eulerian path
    path = find_eulerian_path_paired(graph)
    
    # Convert path nodes to gapped patterns
    # Each edge in the path represents a (k,d)-mer
    gapped_patterns = []
    for i in range(len(path) - 1):
        first_prefix, second_prefix = path[i]
        first_suffix, second_suffix = path[i + 1]
        
        # Reconstruct the k-mer from prefix and suffix
        first_kmer = first_prefix + first_suffix[-1]
        second_kmer = second_prefix + second_suffix[-1]
        
        gapped_patterns.append((first_kmer, second_kmer))
    
    # Spell string from gapped patterns
    result = string_spelled_by_gapped_patterns(gapped_patterns, k, d)
    
    return result

def main():
    # Read input
    with open('bioinfo_genome_sequencing/datasets/dataset_10.txt', 'r') as f:
        lines = f.readlines()
    
    # Parse k and d
    k, d = map(int, lines[0].strip().split())
    
    # Parse paired reads
    paired_reads = []
    read_strings = lines[1].strip().split()
    
    for read_str in read_strings:
        parts = read_str.split('|')
        paired_reads.append((parts[0], parts[1]))
    
    # Reconstruct string
    result = string_reconstruction_from_read_pairs(k, d, paired_reads)
    
    if result is None:
        print("Error: Cannot reconstruct string")
        result = "Error"
    
    # Write output
    with open('bioinfo_genome_sequencing/datasets/output_10.txt', 'w') as f:
        f.write(result)
    
    print(f"Reconstructed string: {result}")
    print(f"Output written to output.txt")

if __name__ == "__main__":
    main()
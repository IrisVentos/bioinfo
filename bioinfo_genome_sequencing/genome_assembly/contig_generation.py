from collections import defaultdict
from helpers import debruijn_graph_from_kmers

def calculate_degrees(graph):
    """Calculate in-degree and out-degree for each node."""
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    
    for node in graph:
        out_degree[node] = len(graph[node])
        for neighbor in graph[node]:
            in_degree[neighbor] += 1
    
    # Ensure all nodes are in both dictionaries
    all_nodes = set(graph.keys()) | set(in_degree.keys())
    for node in all_nodes:
        if node not in in_degree:
            in_degree[node] = 0
        if node not in out_degree:
            out_degree[node] = 0
    
    return in_degree, out_degree

def is_one_in_one_out(node, in_degree, out_degree):
    """Check if a node has exactly one incoming and one outgoing edge."""
    return in_degree[node] == 1 and out_degree[node] == 1

def find_maximal_non_branching_paths(graph):
    """
    Find all maximal non-branching paths in the graph.
    These paths form the contigs.
    """
    paths = []
    in_degree, out_degree = calculate_degrees(graph)
    
    # Track which nodes have been used as path starts
    used = set()
    
    # Find all maximal non-branching paths starting from branching nodes
    for node in graph:
        # If node is not 1-in-1-out, it's a branching node
        if not is_one_in_one_out(node, in_degree, out_degree):
            # Explore all outgoing edges
            if out_degree[node] > 0:
                for neighbor in graph[node]:
                    # Start a new path
                    path = [node, neighbor]
                    current = neighbor
                    
                    # Extend path while current node is 1-in-1-out
                    while is_one_in_one_out(current, in_degree, out_degree):
                        used.add(current)
                        # Get the single outgoing neighbor
                        current = graph[current][0]
                        path.append(current)
                    
                    paths.append(path)
    
    # Find isolated cycles (all nodes are 1-in-1-out)
    for node in graph:
        if node not in used and is_one_in_one_out(node, in_degree, out_degree):
            # Start a cycle from this node
            path = [node]
            current = graph[node][0]
            
            while current != node:
                path.append(current)
                used.add(current)
                if current not in graph or not graph[current]:
                    break
                current = graph[current][0]
            
            if current == node:
                path.append(node)
                paths.append(path)
    
    return paths

def path_to_contig(path):
    """Convert a path of nodes to a contig string."""
    if not path:
        return ""
    
    # Start with first node
    contig = path[0]
    
    # Add last character of each subsequent node
    for i in range(1, len(path)):
        contig += path[i][-1]
    
    return contig

def generate_contigs(patterns):
    """
    Generate all contigs from a collection of k-mers.
    
    Args:
        patterns: List of k-mers
    
    Returns:
        List of contig strings
    """
    # Build De Bruijn graph
    graph = debruijn_graph_from_kmers(patterns)
    
    # Find maximal non-branching paths
    paths = find_maximal_non_branching_paths(graph)
    
    # Convert paths to contigs
    contigs = [path_to_contig(path) for path in paths]
    
    # Sort contigs lexicographically
    contigs.sort()
    
    return contigs

def main():
    # Read input
    with open('bioinfo_genome_sequencing/datasets/dataset_11.txt', 'r') as f:
        line = f.read().strip()
    
    # Parse k-mers
    patterns = line.split()
    
    # Generate contigs
    contigs = generate_contigs(patterns)
    
    # Write output (space-separated)
    with open('bioinfo_genome_sequencing/datasets/output_11.txt', 'w') as f:
        f.write(' '.join(contigs))
    
    print(f"Generated {len(contigs)} contigs:")
    print(' '.join(contigs))
    print(f"Output written to output.txt")

if __name__ == "__main__":
    main()
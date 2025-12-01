def debruijn_graph_from_kmers(kmers):
    """
    Construct the De Bruijn graph from a collection of k-mers.
    
    Args:
        kmers: List of k-mers (may contain duplicates)
    
    Returns:
        Dictionary representing adjacency list
    """
    graph = {}
    
    # For each k-mer, create an edge from prefix to suffix
    for kmer in kmers:
        prefix = kmer[:-1]  # First k-1 characters
        suffix = kmer[1:]   # Last k-1 characters
        
        # Add edge (duplicates are preserved)
        if prefix not in graph:
            graph[prefix] = []
        graph[prefix].append(suffix)
    
    return graph

from collections import defaultdict

def read_graph(filename):
    """Read adjacency list from input file."""
    graph = defaultdict(list)
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split(':')
            node = int(parts[0].strip())
            
            if len(parts) > 1 and parts[1].strip():
                neighbors = [int(x.strip()) for x in parts[1].strip().split()]
                graph[node] = neighbors
            else:
                graph[node] = []
    
    return graph

def create_edge_list(graph):
    """Convert adjacency list to a list of available edges."""
    edges = defaultdict(list)
    for node in graph:
        for neighbor in graph[node]:
            edges[node].append(neighbor)
    return edges

def calculate_degrees(graph):
    """Calculate in-degree and out-degree for each node."""
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    
    for node in graph:
        out_degree[node] += len(graph[node])
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

def find_start_node(graph):
    """Find the starting node for Eulerian path."""
    in_degree, out_degree = calculate_degrees(graph)
    
    # Look for node with out-degree - in-degree = 1 (start of path)
    start_node = None
    for node in set(list(in_degree.keys()) + list(out_degree.keys())):
        if out_degree[node] - in_degree[node] == 1:
            start_node = node
            break
    
    # If no such node exists, start from any node with edges (Eulerian cycle case)
    if start_node is None:
        for node in graph:
            if graph[node]:
                start_node = node
                break
    
    return start_node

def find_cycle(edges, start):
    """Find a cycle starting from a given node."""
    cycle = [start]
    current = start
    
    while edges[current]:
        next_node = edges[current].pop(0)
        cycle.append(next_node)
        current = next_node
    
    return cycle

def has_unexplored_edges(edges):
    """Check if there are any unexplored edges remaining."""
    for node in edges:
        if edges[node]:
            return True
    return False

def find_node_with_unexplored_edges(cycle, edges):
    """Find a node in the cycle that still has unexplored edges."""
    for node in cycle:
        if edges[node]:
            return node
    return None

def eulerian_path(graph):
    """Find an Eulerian path in the graph."""
    edges = create_edge_list(graph)
    
    # Find the starting node
    start_node = find_start_node(graph)
    
    if start_node is None:
        return []
    
    # Form initial cycle/path
    path = find_cycle(edges, start_node)
    
    # While there are unexplored edges
    while has_unexplored_edges(edges):
        # Find a node in the current path with unexplored edges
        new_start = find_node_with_unexplored_edges(path, edges)
        
        if new_start is None:
            break
        
        # Form a new cycle starting from new_start
        new_cycle = find_cycle(edges, new_start)
        
        # Merge: insert the new cycle at the position of new_start
        idx = path.index(new_start)
        path = path[:idx] + new_cycle + path[idx+1:]
    
    return path

def path_to_genome(kmers):
    """
    Reconstruct a string from its genome path.
    
    Args:
        kmers: List of k-mers where consecutive k-mers overlap by k-1 symbols
    
    Returns:
        Reconstructed genome string
    """
    if not kmers:
        return ""
    
    # Start with the first k-mer
    genome = kmers[0]
    
    # Add the last symbol of each subsequent k-mer
    for i in range(1, len(kmers)):
        genome += kmers[i][-1]
    
    return genome


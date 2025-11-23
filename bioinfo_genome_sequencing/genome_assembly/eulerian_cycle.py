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

def find_cycle(edges, start):
    """Find a cycle starting from a given node."""
    cycle = [start]
    current = start
    
    while edges[current]:
        # Randomly pick next edge (or just take first available)
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

def rotate_cycle(cycle, new_start):
    """Rotate cycle to start at new_start node."""
    if new_start not in cycle:
        return cycle
    
    idx = cycle.index(new_start)
    # Rotate: keep last element (which equals first) at the end
    return cycle[idx:-1] + cycle[:idx+1]

def eulerian_cycle(graph):
    """Find an Eulerian cycle in the graph."""
    edges = create_edge_list(graph)
    
    # Start from any node with edges
    start_node = None
    for node in edges:
        if edges[node]:
            start_node = node
            break
    
    if start_node is None:
        return []
    
    # Form initial cycle
    cycle = find_cycle(edges, start_node)
    
    # While there are unexplored edges
    while has_unexplored_edges(edges):
        # Find a node in the current cycle with unexplored edges
        new_start = find_node_with_unexplored_edges(cycle, edges)
        
        if new_start is None:
            break
        
        # Form a new cycle starting from new_start
        new_cycle = find_cycle(edges, new_start)
        
        # Merge cycles: rotate original cycle to start at new_start,
        # then insert the new cycle
        idx = cycle.index(new_start)
        cycle = cycle[:idx] + new_cycle + cycle[idx+1:]
    
    return cycle

def write_output(filename, cycle):
    """Write the Eulerian cycle to output file."""
    with open(filename, 'w') as f:
        f.write(' '.join(map(str, cycle)) + '\n')

def main():
    # Input and output file names
    input_file = "datasets/dataset_6.txt"
    output_file = "datasets/output_6.txt"
    
    # Read graph from input file
    graph = read_graph(input_file)
    
    # Find Eulerian cycle
    cycle = eulerian_cycle(graph)
    
    # Write output
    write_output(output_file, cycle)
    
    print(f"Eulerian cycle found: {' '.join(map(str, cycle))}")
    print(f"Output written to {output_file}")

if __name__ == "__main__":
    main()
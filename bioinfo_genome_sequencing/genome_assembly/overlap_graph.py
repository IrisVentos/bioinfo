def overlap_graph(patterns):
    """
    Construct the overlap graph from a collection of k-mers.
    
    Args:
        patterns: List of k-mers
    
    Returns:
        Dictionary representing adjacency list (node -> list of neighbors)
    """
    graph = {}
    
    # For each pair of k-mers, check if they overlap
    for pattern1 in patterns:
        neighbors = []
        suffix = pattern1[1:]  # Last k-1 symbols of pattern1
        
        for pattern2 in patterns:
            if pattern1 != pattern2:  # Don't connect a k-mer to itself
                prefix = pattern2[:-1]  # First k-1 symbols of pattern2
                
                if suffix == prefix:
                    neighbors.append(pattern2)
        
        # Only add to graph if there are outgoing edges
        if neighbors:
            graph[pattern1] = neighbors
    
    return graph

def main():
    # Read input from file
    with open('datasets/dataset_3.txt', 'r') as f:
        content = f.read().strip()
    
    # Parse k-mers (space-separated)
    patterns = content.split()
    
    # Build overlap graph
    graph = overlap_graph(patterns)
    
    # Build output and print
    output_lines = []
    for node in sorted(graph.keys()):
        neighbors = ' '.join(graph[node])
        line = f"{node}: {neighbors}"
        print(line)  # Print to console
        output_lines.append(line)  # Save for file
    
    # Save to output file
    with open('datasets/output_3.txt', 'w') as f:
        f.write('\n'.join(output_lines))

if __name__ == "__main__":
    main()
def PatternMatching(Pattern, Genome):
    """
    Find all starting positions where Pattern appears as a substring of Genome
    
    Args:
        Pattern (str): The pattern to search for
        Genome (str): The genome string to search in
    
    Returns:
        list: List of starting positions (0-indexed)
    """
    positions = []
    pattern_length = len(Pattern)
    genome_length = len(Genome)
    
    # Search through the genome
    for i in range(genome_length - pattern_length + 1):
        # Check if pattern matches at position i
        if Genome[i:i + pattern_length] == Pattern:
            positions.append(i)
    
    return positions

# Test with sample data
def test_sample():
    print("=== Testing with sample data ===")
    pattern = "ATAT"
    genome = "GATATATGCATATACTT"
    
    positions = PatternMatching(pattern, genome)
    result_str = " ".join(map(str, positions))
    
    print(f"Pattern: {pattern}")
    print(f"Genome:  {genome}")
    print(f"Positions: {result_str}")
    print(f"Expected:  1 3 9")
    print(f"Correct:   {result_str == '1 3 9'}")
    print()
    
    # Show visual representation
    print("Visual representation:")
    print(f"Genome:   {genome}")
    print("          ", end="")
    for i, char in enumerate(genome):
        if i in positions:
            print("^", end="")
        else:
            print(" ", end="")
    print()
    print()
    
    return result_str

# Function to solve from file
def solve_from_file(filename):
    """
    Read Pattern and Genome from file and solve the Pattern Matching problem
    Expected file format:
    Line 1: Pattern
    Line 2: Genome
    """
    try:
        with open(filename, 'r') as file:
            lines = file.read().strip().split('\n')
            
            if len(lines) >= 2:
                pattern = lines[0].strip()
                genome = lines[1].strip()
                
                print(f"Pattern: {pattern}")
                print(f"Genome:  {genome}")
                
                positions = PatternMatching(pattern, genome)
                result_str = " ".join(map(str, positions))
                
                print(f"Starting positions: {result_str}")
                
                # Save result to output file
                with open('output_bonus_v2.txt', 'w') as output:
                    output.write(result_str)
                print("Result saved to output_bonus_v2.txt")
                
                return result_str
            else:
                print("Error: File should have at least 2 lines (Pattern and Genome)")
                
    except FileNotFoundError:
        print(f"File '{filename}' not found")
    except Exception as e:
        print(f"Error: {e}")

# Debug function to show all matches with context
def debug_matches(Pattern, Genome):
    """
    Show detailed information about each match
    """
    positions = PatternMatching(Pattern, Genome)
    
    print(f"=== Debug: Finding '{Pattern}' in genome ===")
    print(f"Genome length: {len(Genome)}")
    print(f"Pattern length: {len(Pattern)}")
    print(f"Number of matches: {len(positions)}")
    print()
    
    if positions:
        print("Matches found at positions:")
        for pos in positions:
            # Show context around each match
            start_context = max(0, pos - 5)
            end_context = min(len(Genome), pos + len(Pattern) + 5)
            
            context = Genome[start_context:end_context]
            match_start = pos - start_context
            match_end = match_start + len(Pattern)
            
            print(f"Position {pos}: ...{context}...")
            print(f"           {' ' * (len(str(pos)) + 2 + len('Position ') + 3 + match_start)}{'^' * len(Pattern)}")
            print()
    else:
        print("No matches found.")

# Main execution
if __name__ == "__main__":
    print("Pattern Matching Problem Solver")
    print("===============================")
    
    # Test with sample data
    test_sample()
    
    # Debug the sample
    debug_matches("ATAT", "GATATATGCATATACTT")
    
    # Solve from your file
    print("="*50)
    solve_from_file('data/part_1_data/bonus_sample.txt')  
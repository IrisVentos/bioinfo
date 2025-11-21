from helpers import approximate_pattern_matching

def approximate_pattern_count(pattern, text, d):
    """
    Count the number of approximate occurrences of a pattern in a string.
    
    Args:
        pattern (str): The pattern to search for
        text (str): The text to search in
        d (int): Maximum number of mismatches allowed
        
    Returns:
        int: Count of approximate pattern occurrences
    """
    # Reuse the approximate_pattern_matching function
    positions = approximate_pattern_matching(pattern, text, d)
    
    # Return the count of positions found
    return len(positions)

# Read input from file
with open('data/part_2_data/sample_5.txt', 'r') as f:
    lines = f.read().strip().split('\n')
    pattern = lines[0]
    text = lines[1]
    d = int(lines[2])

# Count approximate pattern matches
result = approximate_pattern_count(pattern, text, d)

# Print the count
print(result)
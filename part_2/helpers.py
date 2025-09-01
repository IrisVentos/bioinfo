#DNA analysis functions

def calculate_skew(genome):
    """
    Calculate the skew values for a DNA string.
    
    Skew_i(Genome) = difference between total G's and total C's 
    in the first i nucleotides of Genome.
    
    Args:
        genome (str): DNA string containing nucleotides A, T, G, C
        
    Returns:
        list: List of skew values from position 0 to len(genome)
    """
    skew_values = [0]  # Skew_0 is always 0
    current_skew = 0
    
    for nucleotide in genome:
        if nucleotide == 'G':
            current_skew += 1
        elif nucleotide == 'C':
            current_skew -= 1
        # For A and T, skew remains the same
        
        skew_values.append(current_skew)
    
    return skew_values

def find_minimum_skew(genome):
    """
    Find all positions where the skew diagram attains a minimum.
    
    Args:
        genome (str): DNA string containing nucleotides A, T, G, C
        
    Returns:
        list: All positions (indices) where skew is minimum
    """
    # Reuse the calculate_skew function
    skew_values = calculate_skew(genome)
    
    # Find the minimum skew value
    min_skew = min(skew_values)
    
    # Find all positions where this minimum occurs
    min_positions = []
    for i, skew in enumerate(skew_values):
        if skew == min_skew:
            min_positions.append(i)
    
    return min_positions

def solve_minimum_skew_problem(genome):
    """
    Complete solution to the Minimum Skew Problem.
    
    Args:
        genome (str): DNA string
        
    Returns:
        dict: Dictionary containing minimum positions, skew values, and min skew value
    """
    min_positions = find_minimum_skew(genome)
    skew_values = calculate_skew(genome)
    min_skew_value = min(skew_values)
    
    return {
        'min_positions': min_positions,
        'skew_values': skew_values,
        'min_skew_value': min_skew_value
    }

def format_minimum_skew_result(genome):
    """
    Format the minimum skew result for display
    
    Args:
        genome (str): DNA string
        
    Returns:
        str: Formatted result string
    """
    result = solve_minimum_skew_problem(genome)
    return ' '.join(map(str, result['min_positions']))

from helpers import profile_probability, profile_most_probable_kmer


def parse_profile_matrix(profile_lines):
    """
    Parse the profile matrix from input lines.
    Each line contains space-separated probabilities for one nucleotide.
    """
    profile = []
    for line in profile_lines:
        # Convert string numbers to floats
        row = [float(x) for x in line.strip().split()]
        profile.append(row)
    
    return profile


def normalize_profile_matrix(profile):
    """
    Normalize the profile matrix so that each column sums exactly to 1.0.
    This fixes small floating-point precision errors.
    """
    if len(profile) == 0:
        return profile
    
    k = len(profile[0])  # number of positions
    normalized_profile = [row[:] for row in profile]  # deep copy
    
    for col in range(k):
        # Calculate current column sum
        column_sum = sum(profile[row][col] for row in range(4))
        
        # Normalize each entry in this column
        if column_sum > 0:  # avoid division by zero
            for row in range(4):
                normalized_profile[row][col] = profile[row][col] / column_sum
    
    return normalized_profile


def validate_profile_matrix(profile):
    """
    Validate that the profile matrix is properly formatted.
    Each column should sum to 1.0 (allowing for floating point precision errors).
    """
    if len(profile) != 4:
        return False, "Profile must have exactly 4 rows (A, C, G, T)"
    
    k = len(profile[0])
    
    # Check that all rows have the same length
    for i, row in enumerate(profile):
        if len(row) != k:
            return False, f"Row {i} has {len(row)} elements, expected {k}"
    
    # Check that each column sums to approximately 1.0
    # Use a more generous tolerance for floating point errors
    tolerance = 0.01  # Allow up to 1% deviation
    problematic_columns = []
    
    for col in range(k):
        column_sum = sum(profile[row][col] for row in range(4))
        deviation = abs(column_sum - 1.0)
        
        if deviation > tolerance:
            problematic_columns.append((col, column_sum, deviation))
    
    if problematic_columns:
        error_details = []
        for col, sum_val, deviation in problematic_columns:
            error_details.append(f"Column {col}: sum={sum_val:.6f} (deviation: {deviation:.6f})")
        
        return False, f"Columns don't sum to 1.0:\n" + "\n".join(error_details)
    
    # Show any minor deviations as warnings but still proceed
    minor_warnings = []
    for col in range(k):
        column_sum = sum(profile[row][col] for row in range(4))
        deviation = abs(column_sum - 1.0)
        
        if deviation > 0.001:  # Show warning for deviations > 0.001
            minor_warnings.append(f"Column {col}: sum={column_sum:.6f} (small deviation: {deviation:.6f})")
    
    if minor_warnings:
        warning_msg = "Profile matrix valid with minor floating-point deviations:\n" + "\n".join(minor_warnings)
        return True, warning_msg
    
    return True, "Profile matrix is valid"


def display_profile_analysis(text, k, profile, result_kmer, position, probability):
    """
    Display detailed analysis of the profile and results.
    """
    print("=== PROFILE MATRIX ANALYSIS ===")
    print("Profile matrix (A, C, G, T rows):")
    nucleotides = ['A', 'C', 'G', 'T']
    
    for i, nucleotide in enumerate(nucleotides):
        print(f"{nucleotide}: {' '.join(f'{prob:.3f}' for prob in profile[i])}")
    
    print(f"\nPosition:  {' '.join(str(i).rjust(5) for i in range(1, k+1))}")
    
    # Show most likely nucleotide at each position
    consensus = ""
    for col in range(k):
        max_prob = max(profile[row][col] for row in range(4))
        best_nucleotides = [nucleotides[row] for row in range(4) if profile[row][col] == max_prob]
        consensus += best_nucleotides[0]  # Take first in case of tie
    
    print(f"Consensus: {' '.join(nuc.rjust(5) for nuc in consensus)}")
    
    print(f"\n=== SEARCH RESULTS ===")
    print(f"Text: {text}")
    print(f"k-mer length: {k}")
    print(f"Most probable k-mer: {result_kmer}")
    print(f"Position in text: {position} (0-indexed)")
    print(f"Probability: {probability:.10f}")
    
    # Show the k-mer highlighted in the text
    highlighted_text = text[:position] + "[" + result_kmer + "]" + text[position + k:]
    print(f"Highlighted: {highlighted_text}")
    
    # Show probability calculation breakdown
    print(f"\n=== PROBABILITY BREAKDOWN FOR '{result_kmer}' ===")
    nucleotide_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    prob_parts = []
    
    for i, nuc in enumerate(result_kmer):
        row = nucleotide_to_index[nuc]
        pos_prob = profile[row][i]
        prob_parts.append(f"P({nuc}|pos{i+1})={pos_prob:.3f}")
    
    print(" × ".join(prob_parts))
    print(f"= {probability:.10f}")


def main():
    # Read input from file
    try:
        with open('data/part_3_data.sample_3.txt', 'r') as file: 
            lines = file.read().strip().split('\n')
    except FileNotFoundError:
        print("Error: input.txt file not found!")
        print("\nExpected format:")
        print("Line 1: DNA text string")
        print("Line 2: k (k-mer length)")
        print("Lines 3-6: Profile matrix (4 rows for A,C,G,T)")
        return
    
    # Parse input
    if len(lines) < 6:
        print("Error: Input file must have at least 6 lines")
        return
    
    text = lines[0].strip()
    k = int(lines[1].strip())
    
    # Parse profile matrix (4 rows)
    profile_lines = lines[2:6]
    profile = parse_profile_matrix(profile_lines)
    
    # Validate profile matrix
    is_valid, message = validate_profile_matrix(profile)
    
    if not is_valid:
        print(f"Profile validation failed: {message}")
        print("\nAttempting to fix by normalizing columns...")
        
        # Try to fix by normalizing
        normalized_profile = normalize_profile_matrix(profile)
        is_valid_after_fix, fix_message = validate_profile_matrix(normalized_profile)
        
        if is_valid_after_fix:
            print("✓ Successfully normalized the profile matrix")
            profile = normalized_profile
            print(f"Status: {fix_message}")
        else:
            print(f"✗ Could not fix profile matrix: {fix_message}")
            return
    else:
        print(f"Profile matrix validation: {message}")
    
    # Find the most probable k-mer using imported function from helpers.py
    # Note: helpers.py function returns (kmer, position, probability)
    result = profile_most_probable_kmer(text, k, profile)
    
    # Handle different return formats from helpers.py
    if isinstance(result, tuple) and len(result) == 3:
        result_kmer, position, probability = result
    else:
        # If helpers.py only returns the k-mer string
        result_kmer = result
        # Calculate position and probability manually
        position = -1
        probability = 0.0
        max_prob = -1.0
        for i in range(len(text) - k + 1):
            kmer = text[i:i + k]
            prob = profile_probability(kmer, profile)
            if prob > max_prob:
                max_prob = prob
                position = i
                probability = prob
    
    if result_kmer:
        # Display detailed analysis
        display_profile_analysis(text, k, profile, result_kmer, position, probability)
        
        # Output the answer (matching expected format)
        print(f"\n=== FINAL ANSWER ===")
        print(result_kmer)
        
    else:
        print("No valid k-mer found in the text")


if __name__ == "__main__":
    main()
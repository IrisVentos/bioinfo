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

# Test with the sample input
sample_input = "CATGGGCATCGGCCATACGCC"
sample_skew = calculate_skew(sample_input)
print("Sample Input:", sample_input)
print("Sample Output:", ' '.join(map(str, sample_skew)))

# Exercise: Calculate skew for GAGCCACCGCGATA
exercise_input = "GAGCCACCGCGATA"
exercise_skew = calculate_skew(exercise_input)
print("\nExercise Input:", exercise_input)
print("Exercise Output:", ' '.join(map(str, exercise_skew)))

# Let's also show step-by-step calculation for the exercise
print("\nStep-by-step calculation for", exercise_input)
print("Position | Nucleotide | Skew")
print("---------|------------|-----")
current_skew = 0
print(f"    0    |     -      |  {current_skew}")

for i, nucleotide in enumerate(exercise_input):
    if nucleotide == 'G':
        current_skew += 1
        change = "+1"
    elif nucleotide == 'C':
        current_skew -= 1
        change = "-1"
    else:
        change = " 0"
    
    print(f"    {i+1}    |     {nucleotide}      |  {current_skew} ({change})")
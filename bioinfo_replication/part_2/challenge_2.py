from helpers import find_minimum_skew

# Read genome from file
with open('data/part_2_data/sample_2.txt', 'r') as f:
    genome = f.read().strip()

# Find and print minimum skew positions
result = find_minimum_skew(genome)
print(' '.join(map(str, result)))
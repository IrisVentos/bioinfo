from helpers import hamming_distance,approximate_pattern_matching


# Read input from file
with open('data/part_2_data/test.txt', 'r') as f:
    lines = f.read().strip().split('\n')
    pattern = lines[0]
    text = lines[1]
    d = int(lines[2])

# Find approximate pattern matches
result = approximate_pattern_matching(pattern, text, d)

# Print result as space-separated integers
print(' '.join(map(str, result)))
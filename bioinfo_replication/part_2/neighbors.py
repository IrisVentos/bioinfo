from helpers import neighbors

# Read input from file
with open('data/part_2_data/final_sample.txt', 'r') as f:
    lines = f.read().strip().split('\n')
    pattern = lines[0]
    d = int(lines[1])

# Find d-neighborhood
result = neighbors(pattern, d)

# Print each string on a separate line
for neighbor in sorted(result):  # sorted for consistent output
    print(neighbor)


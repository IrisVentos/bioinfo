from helpers import hamming_distance

# Read input from file
with open('data/part_2_data/exam.txt', 'r') as f:
    lines = f.read().strip().split('\n')
    string1 = lines[0]
    string2 = lines[1]

# Calculate and print Hamming distance
result = hamming_distance(string1, string2)
print(result)
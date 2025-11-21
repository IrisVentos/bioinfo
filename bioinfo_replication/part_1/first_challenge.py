# Simple PatternCount Implementation
def PatternCount(Text, Pattern):
    count = 0
    pattern_length = len(Pattern)
    
    for i in range(len(Text) - pattern_length + 1):
        if Text[i:i + pattern_length] == Pattern:
            count += 1
    
    return count

# Step 1: Let's test with simple data first
print("Step 1: Testing the function")
test_text = "GCGCG"
test_pattern = "GCG"
result = PatternCount(test_text, test_pattern)
print(f"Text: {test_text}")
print(f"Pattern: {test_pattern}")
print(f"Result: {result}")
print("-" * 30)

# Step 2: Actual file 
print("Step 2")
try:
    with open('data/first_sample.txt', 'r') as file:
        content = file.read()
        print("File content:")
        print(repr(content))  
        print("-" * 30)
        
        lines = content.strip().split('\n')
        print(f"Number of lines found: {len(lines)}")
        
        if len(lines) >= 2:
            text = lines[0]
            pattern = lines[1]
            
            print(f"Text (line 1): {text}")
            print(f"Pattern (line 2): {pattern}")
            
            result = PatternCount(text, pattern)
            print(f"Final Result: {result}")
            
            # Save to output file
            with open('output_v1.txt', 'w') as output_file:
                output_file.write(str(result))
            print("Result saved to output.txt")
        else:
            print("Error: Need at least 2 lines in the file")
            
except FileNotFoundError:
    print("File not found in current directory")
    print("Make sure the file is in the same folder as this Python script")
    
except Exception as e:
    print(f"Error reading file: {e}")
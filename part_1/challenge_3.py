#Thinking of the necessary steps first:
#1. Mapping complementary nucleotides (A&T, C&G)
#2. Apply complementaty nucleotide to each input nucleotide
#3. Reverse str
#4. Print output

def ReverseComplement(Pattern):
    """
    Find the reverse complement of a DNA string
    
    Args:
        Pattern (str): Input DNA string
    
    Returns:
        str: Reverse complement of the input pattern
    """
    # Define complement mapping
    complement_map = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    
    # Step 1: Take complement of each nucleotide
    complement = ""
    for nucleotide in Pattern:
        complement += complement_map[nucleotide]
    
    # Step 2: Reverse the complement string
    reverse_complement = complement[::-1]
    
    return reverse_complement

# Test with sample data
def test_sample():
    print("=== Testing with sample data ===")
    pattern = "AAAACCCGGT"
    
    result = ReverseComplement(pattern)
    
    print(f"Input:    {pattern}")
    print(f"Output:   {result}")
    print(f"Expected: ACCGGGTTTT")
    print(f"Correct:  {result == 'ACCGGGTTTT'}")
    print()

# Function to solve from file
def solve_from_file(filename):
    """
    Read DNA pattern from file and find its reverse complement
    Expected file format:
    Line 1: DNA pattern
    """
    try:
        with open(filename, 'r') as file:
            lines = file.read().strip().split('\n')
            
            if len(lines) >= 1:
                pattern = lines[0].strip()
                
                print(f"Input pattern: {pattern}")
                
                result = ReverseComplement(pattern)
                
                print(f"Reverse complement: {result}")
                
                # Save result to output file
                with open('output_3.txt', 'w') as output:
                    output.write(result)
                print("Result saved to output_3.txt")
                
                return result
            else:
                print("Error: File should contain at least 1 line with DNA pattern")
                
    except FileNotFoundError:
        print(f"File '{filename}' not found")
    except KeyError as e:
        print(f"Error: Invalid nucleotide found: {e}")
    except Exception as e:
        print(f"Error: {e}")

# Main execution
if __name__ == "__main__":
    print("Reverse Complement Problem Solver")
    print("=================================")
    
    # Test with sample data
    test_sample()
    
    # Solve from file
    solve_from_file('data/sample_3.txt')  
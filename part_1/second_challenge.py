#second_challenge

def FrequentWords(Text, k):
    """
    Find all most frequent k-mers in Text
    
    Args:
        Text (str): Input DNA string
        k (int): Length of k-mers to find
    
    Returns:
        list: All most frequent k-mers
    """
    # Dictionary to count frequency of each k-mer
    frequency_map = {}
    
    # Generate all k-mers and count their frequencies
    for i in range(len(Text) - k + 1):
        kmer = Text[i:i + k]
        if kmer in frequency_map:
            frequency_map[kmer] += 1
        else:
            frequency_map[kmer] = 1
    
    # Find the maximum frequency
    max_frequency = max(frequency_map.values())
    
    # Find all k-mers with maximum frequency
    frequent_kmers = []
    for kmer, frequency in frequency_map.items():
        if frequency == max_frequency:
            frequent_kmers.append(kmer)
    
    return frequent_kmers

# Test with the sample data
def test_sample():
    print("=== Testing with sample data ===")
    text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
    k = 4
    
    result = FrequentWords(text, k)
    result_str = " ".join(result)  # Format as required
    
    print(f"Text: {text}")
    print(f"k: {k}")
    print(f"Most frequent {k}-mers: {result_str}")
    print(f"Expected: CATG GCAT")
    print()
    
    return result_str

# Function to solve from file
def solve_from_file(filename):
    """
    Read Text and k from file and solve the Frequent Words problem
    Expected file format:
    Line 1: Text
    Line 2: k (integer)
    """
    try:
        with open(filename, 'r') as file:
            lines = file.read().strip().split('\n')
            
            if len(lines) >= 2:
                text = lines[0].strip()
                k = int(lines[1].strip())
                
                print(f"Text: {text}")
                print(f"k: {k}")
                
                result = FrequentWords(text, k)
                result_str = " ".join(result)
                
                print(f"Most frequent {k}-mers: {result_str}")
                
                # Save result to output file
                with open('output.txt', 'w') as output:
                    output.write(result_str)
                print("Result saved to output.txt")
                
                return result_str
            else:
                print("Error: File should have at least 2 lines (Text and k)")
                
    except FileNotFoundError:
        print(f"File '{filename}' not found")
    except ValueError:
        print("Error: Second line should be an integer (k value)")
    except Exception as e:
        print(f"Error: {e}")

# Interactive mode
def solve_interactive():
    """
    Interactive input for the Frequent Words problem
    """
    print("=== Interactive mode ===")
    
    text = input("Enter the DNA text: ").strip()
    k = int(input("Enter k (length of k-mers): ").strip())
    
    result = FrequentWords(text, k)
    result_str = " ".join(result)
    
    print(f"\nText: {text}")
    print(f"k: {k}")
    print(f"Most frequent {k}-mers: {result_str}")
    
    return result_str

# Debug function to show all k-mers and their frequencies
def debug_kmers(Text, k):
    """
    Show all k-mers and their frequencies for debugging
    """
    print(f"=== Debug: All {k}-mers in {Text} ===")
    
    frequency_map = {}
    
    for i in range(len(Text) - k + 1):
        kmer = Text[i:i + k]
        if kmer in frequency_map:
            frequency_map[kmer] += 1
        else:
            frequency_map[kmer] = 1
    
    # Sort by frequency (descending) then by k-mer name
    sorted_kmers = sorted(frequency_map.items(), key=lambda x: (-x[1], x[0]))
    
    print("K-mer\tFrequency")
    print("-" * 15)
    for kmer, freq in sorted_kmers:
        print(f"{kmer}\t{freq}")
    
    max_freq = max(frequency_map.values())
    print(f"\nMaximum frequency: {max_freq}")

# Main execution
if __name__ == "__main__":
    print("Frequent Words Problem Solver")
    print("============================")
    
    # Test with sample data
    test_sample()
    
    # Show debug information for sample
    debug_kmers("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4)
    
    # Try to solve from file (adjust filename as needed)
    print("\n" + "="*50)
    solve_from_file('data/second_sample.txt')  # Change this to your file path
    
    # Uncomment for interactive mode
    # solve_interactive()
# Goal is to count the number of linear subpeptides of a peptide of length n.
# A linear peptide of length n has n*(n+1)/2 subpeptides, including the empty peptide and the full peptide itself.

def count_linear_subpeptides(n):
    return n * (n + 1) // 2 + 1

n = int(input("Enter peptide length n: "))
print(count_linear_subpeptides(n))
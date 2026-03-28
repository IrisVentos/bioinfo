# We say that a DNA string Pattern encodes an amino acid string Peptide if the RNA string transcribed from either Pattern or its reverse complement Pattern translates into Peptide.
# For example, the DNA string GAAACT is transcribed into GAAACU and translated into ET. 
# The reverse complement of this DNA string, AGTTTC, is transcribed into AGUUUC and translated into SF. Thus, GAAACT encodes both ET and SF.

# Let us find all substrings of a DNA string Text that encode a given amino acid sequence 

from protein_translation import translate_rna_to_aa_string, parse_codon_table

def transcribe_dna_to_rna(dna_string):
    """Transcribe a DNA string into an RNA string."""
    return dna_string.replace('T', 'U')

def reverse_complement(dna_string):
    """Return the reverse complement of a DNA string."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(dna_string))

def peptide_encoding_from_file(filepath, codon_table):
    """Read DNA and peptide from file and find all substrings encoding the peptide."""
    with open(filepath, 'r') as f:
        text = f.readline().strip()
        peptide = f.readline().strip()

    substring_length = len(peptide) * 3
    results = []

    for i in range(len(text) - substring_length + 1):
        substring = text[i:i + substring_length]

        if translate_rna_to_aa_string(transcribe_dna_to_rna(substring), codon_table) == peptide:
            results.append(substring)

        rev_comp = reverse_complement(substring)
        if translate_rna_to_aa_string(transcribe_dna_to_rna(rev_comp), codon_table) == peptide:
            results.append(substring)

    return results


codon_table = parse_codon_table("bioinfo_genome_sequencing/datasets/RNA_codon_table_1.txt")
results = peptide_encoding_from_file("bioinfo_genome_sequencing/datasets/dataset_13.txt", codon_table)
for r in results:
    print(r)


# After solving the Peptide Encoding Problem for Tyrocidine B1, we should be able to find a 30-mer in the Bacillus brevis genome encoding Tyrocidine B1, and yet no such 30-mer exists!
# But this bacteria produces that antibiotic
# enzymatically produced instead of genetically encoded
# Antibiotics

A substance that kills bacteria. It occurs naturally because of millions of years of evolutionary warfare, like the Penicillin discovered by Fleming in 1920's. It's produced by fungi (e.g molds) and bacteria.

On the molecular level
We will study Tyrocidine B1, an antibiotic produced by bacteria Bacillus Brevis. B1 is a mini-protein called a peptide, short string of amino acids.

## How are antibiotics produced?

Quick reminder of the central "dogma" of molecular biology : DNA is transcribed into RNA by protein RNA polymerase, RNA is translated into proteins by enzyme ribosome.
First laid out by Francis Crick.

Codon = triplet (3-mer) of nucleotides.
Genetic code : assignment of codons to amino acids to make proteins.

Goal of the week :
Find a 30-mer in the Bacillus brevis genome that transcribes and translates into Tyrocidine B1 (peptide of length 10). 
- Unfortunately, thousands candidates could translate into said antibiotic.
- Also translation can start anywhere in the genome : 6 different reading frames (3 on each strand, twice in different directions)
- And Tyrocidine B1 is cyclic, so ten different linear representations depending on where we start on the string of amino acids

Comes dodging the dogma... In the 60's, Edward Tatum inhibits the ribosome in Bacillus Brevis. But instead of stopping all translation, some peptides including tyrocidines are still being produced.
--> Tyrocidines are non-ribosomal peptides (NRPs), produced by protein NRP synthetase. 10 modules, each responsible for one amino-acid, and when ready the peptide circularizes. 


## How do we sequence antibiotics?

In this case, the clue is not hidden in the genome, but in NRPs/PKs enzymes.

Genome mining : by sequencing the genome to identify NRPs/PKs gene clusters, we can predict the structure of the antibiotic it will produce. With tools like antiSMASH, we scan bacterial genomes for biosynthetic gene clusters (BGCs) and predict what natural products they might encode - including antibiotics never yet characterized. 

Molecular weight = 1 dalton (Da aka mass of proton/neutron) x nb of protons in atoms/molecules = approx integer mass

Each amino acid has an integer mass, some have the same so we move from an alphabet of 20 AA to 18 integer masses.
Working with the Mass Spectrometer, building a theoretical spectrum = the list of mass of every possible subpeptide, plus the mass of the peptide plus 0.

Towards a computational problem : going from spectrum to peptide.

Note : Some antibiotics are RiPPs (ribosomally synthesized and post translationally modified peptides) are a large and growing class of antibiotics that are encoded by DNA, transcribed into mRNA, and translated by ribosomes like any normal protein. The resulting peptide is then heavily modified by dedicated enzymes. Examples include:
- Nisin (an antibiotic used in food preservation)
- Microcin antibiotics in bacteria
- Thiopeptides like GE2270A

## Cyclopeptide Sequencing Problem 

Reconstruct a cyclic peptide from its theoretical spectrum
Experimental spectrum → comes from the lab (mass spectrometer breaks the real peptide and measures the fragments)
Theoretical spectrum → comes from math (you take a candidate sequence and compute what its fragments should weigh)
The goal → find the candidate whose theoretical spectrum best matches the experimental one

### Brute force algorithm

The mass of the entire peptide is usually known.
1. Generate all peptides with g iven mass (in this case, trillions)
2. Form their theoretical spectra
3. Look for matches with the given spectrum, try all candidates

Spectrum = sum along the chain of peptides (growing and growing). If same, potential candidates can remain.

### Cyclopeptide sequencing with Branch-and-Bound

It is a method for solving optimization problems by breaking them down into smaller subproblems and using a bounding function to eliminate subproblems that cannot contain the optimal solution.
The algorithm explores branches of this tree, which represent subsets of the solution set. 
Before enumerating the candidate solutions of a branch, the branch is checked against upper and lower estimated bounds on the optimal solution, and is discarded if it cannot produce a better solution than the best one found so far by the algorithm.

Bounds = potential solutions
then branches again to remove hopeless candidates
until final bounds are consistent with Spectrum.
The goal is to trim the initial list

B&B for cyclopeptide sequencing :
1. Find all amino acids whose masses occur in Spectrum. Add to List.
2. Extend each peptide in List by each of 18 different AA masses.
3. Trim inconsistent peptides from List.
4. Return any peptides in List whose theoretical spectra match Spectrum.
5. Iterate Steps 2-4 until List is empty (final bounds)

Note!

The number of DNA strings that encode a given amino acid sequence is the product of the number of codons for each amino acid.

## From theory to experience 

The experimental spectra (as opposed to theoretical) is flawed. The spectrum is noisy in practice.
It produces false masses or hides some.
A false mass is present in the experimental spectrum but absent from the theoretical spectrum; a missing mass is present in the theoretical spectrum but absent from the experimental spectrum.
We need a new approach as the current algo matches exactly a peptide's theoretical spectrum with the experimental spectrum.

--> scoring a peptide on how many masses its spectrum shares with the experimental spectrum and find best candidates, including ties.

### Leaderboard Cyclopeptide Sequencing

The gist is :

1. Add "O peptide" to Leaderboard as LeaderPeptide
2. Extend each peptide in Leaderboard by each 18 different amino-acid masses
3. Cut low-scoring peptides (keep top N with ties)
4. Update LeaderPeptide with candidates with mass(Leaderboard) = parent mass
5. Eliminate others
6. Iterate
7. Return LeaderPeptide

It may generate errors as it is heuristic (trial & error, aka eliminates good candidates in the beginning)

### NRPs exception

As they fall outside of the central dogma, NRPs contain more non-standard amino acids.
From 18 to 100+ amino masses.
It complicates things even more for noisy spectra.

### The spectral convolution saves the day

Goal : reduce the nb of aa that we need to consider
Spectral convolution : positive difference between every pair of masses in spectrum.
 
The convolution of two spectra S₁ and S₂ is defined as the multiset of all pairwise differences between peaks:

Conv(S₁, S₂) = { s₁ − s₂ | s₁ ∈ S₁, s₂ ∈ S₂ }

When you convolve a theoretical spectrum against an experimental spectrum, peaks (of intensity aka # of ions of peptide fragments) that appear frequently in the result correspond to mass shifts that are consistent across many fragment pairs — which reveals likely amino acid masses in the sequence.

The Sequencing Pipeline : 

1. Acquire MS/MS spectrum of the antibiotic peptide
2. Generate candidate amino acid alphabet using spectral convolution (the most frequent differences ≈ likely residue masses)
3. Score candidate sequences by comparing their theoretical spectra against the experimental one
4. Search cyclically for cyclic peptides using a leaderboard algorithm
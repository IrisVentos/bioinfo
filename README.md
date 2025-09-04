# bioinfo
Bioinformatic course - Coursera

----------------------------------------------------------------------
First Code Challenge: Implement PatternCount
    To find surprinsingly frequent substrings/patterns in DNA strings 
    Goal : find ori
    Input: Strings Text and Pattern.
     Output: Count(Text, Pattern).

PatternCount(Text, Pattern)
  count ← 0
  for i ← 0 to |Text| − |Pattern|
    if Text(i, |Pattern|) = Pattern
      count ← count + 1
  return count
Some notes on how code challenges work:

------------------------------------------------------------------------
A straightforward algorithm for finding the most frequent k-mers in a string Text checks all k-mers appearing in this string (there are |Text| − k + 1 such k-mers) and then computes how many times each k-mer appears in Text. To implement this algorithm, called FrequentWords, we will need to generate an array Count, where Count(i) stores Count(Text, Pattern) for Pattern = Text(i, k) (see figure below).



Figure: The array Count for Text = ACTGACTCCCACCCC and k = 3. For example, Count(0) = Count(4) = 2 because ACT (shown in boldface) appears twice in Text.

The pseudocode for FrequentWords is shown below.

FrequentWords(Text, k)
    FrequentPatterns ← an empty set
    for i ← 0 to |Text| − k
        Pattern ← the k-mer Text(i, k)
        Count(i) ← PatternCount(Text, Pattern)
    maxCount ← maximum value in array Count
    for i ← 0 to |Text| − k
        if Count(i) = maxCount
            add Text(i, k) to FrequentPatterns
    remove duplicates from FrequentPatterns
    return FrequentPatterns

    Although FrequentWords finds most frequent k-mers, it is not very efficient. Each call to PatternCount(Text, Pattern) checks whether the k-mer Pattern appears in position 0 of Text, position 1 of Text, and so on. Since each k-mer requires |Text| − k + 1 such checks, each one requiring as many as k comparisons, the overall number of steps of PatternCount(Text, Pattern) is (|Text| − k + 1) · k. Furthermore, FrequentWords must call PatternCount |Text| − k + 1 times (once for each k-mer of Text), so that its overall number of steps is (|Text| − k + 1) · (|Text| − k + 1) · k. To simplify the matter, computer scientists often say that the runtime of FrequentWords has an upper bound of |Text|^2 · k steps and refer to the complexity of this algorithm as O(|Text|2 · k). For more details, see "DETOUR: Big-O Notation" in the print companion.


FrequencyTable(Text, k)
    freqMap ← empty map
    n ← |Text|
    for i ← 0 to n − k
        Pattern ← Text(i, k)
        if freqMap[Pattern] doesn't exist
            freqMap[Pattern]← 1
        else
           freqMap[pattern] ←freqMap[pattern]+1 
    return freqMap

    BetterFrequentWords(Text, k)
    FrequentPatterns ← an array of strings of length 0
    freqMap ← FrequencyTable(Text, k)
    max ← MaxMap(freqMap)
    for all strings Pattern in freqMap
        if freqMap[pattern] = max
            append Pattern to frequentPatterns
    return frequentPatterns

    However, before concluding that we have found the DnaA box of Vibrio cholerae, the careful bioinformatician should check if there are other short regions in the Vibrio cholerae genome exhibiting multiple occurrences of ATGATCAAG (or CTTGATCAT). After all, maybe these strings occur as repeats throughout the entire Vibrio cholerae genome, rather than just in the ori region. To this end, we need to solve the following problem.

Pattern Matching Problem: Find all occurrences of a pattern in a string.

Input: Strings Pattern and Genome.
Output: All starting positions in Genome where Pattern appears as a substring.

-----------------------------------------------------------------------
Looking for hidden messages in multiple genomes
We should not jump to the conclusion that ATGATCAAG/CTTGATCAT is a hidden message for all bacterial genomes without first checking whether it even appears in known ori regions from other bacteria. After all, maybe the clumping effect of ATGATCAAG/CTTGATCAT in the ori region of Vibrio cholerae is simply a statistical fluke that has nothing to do with replication. Or maybe different bacteria have different DnaA boxes…

Before we lose all hope, let’s change our computational focus: instead of finding clumps of a specific k-mer, let’s try to find every k-mer that forms a clump in the genome. Hopefully, the locations of these clumps will shed light on the location of ori.

Our plan is to slide a window of fixed length L along the genome, looking for a region where a k-mer appears several times in short succession. The parameter value L = 500 reflects the typical length of ori in bacterial genomes.

We defined a k-mer as a "clump" if it appears many times within a short interval of the genome. More formally, given integers L and t, a k-mer Pattern forms an (L, t)-clump inside a (longer) string Genome if there is an interval of Genome of length L in which this k-mer appears at least t times. (This definition assumes that the k-mer completely fits within the interval. This also does not take reverse complements into account yet.) For example, TGCA forms a (25,3)-clump in the following Genome:

gatcagcataagggtccCTGCAATGCATGACAAGCCTGCAGTtgttttac

From our previous examples of ori regions, ATGATCAAG forms a (500,3)-clump in the Vibrio cholerae genome, and CCTACCACC forms a (500,3)-clump in the Thermotoga petrophila genome. We are now ready to formulate the following problem.

Clump Finding Problem: Find patterns forming clumps in a string.
     Input: A string Genome, and integers k, L, and t.
     Output: All distinct k-mers forming (L, t)-clumps in Genome.


----------------------------------------------------------------------
We do not know what purpose — if any — these other 9-mers serve in the E. coli genome, but we do know that there are many different types of hidden messages in genomes; these hidden messages have a tendency to cluster within a genome, and most of them have nothing to do with replication. One example is the regulatory DNA motifs responsible for gene expression that we will study soon. The important lesson is that existing approaches to ori prediction remain imperfect and are sometimes inconclusive. However, even providing biologists with a small collection of 9-mers as candidate DnaA boxes is a great aid as long as one of these 9-mers is correct.

Thus, the moral of this chapter is that even though computational predictions can be powerful, bioinformaticians should collaborate with biologists to verify their computational predictions. Or improve these predictions: the next question hints at how ori predictions can be carried out using comparative genomics, a bioinformatics approach that uses evolutionary similarities to answer difficult questions about genomes.

We have considered three genomes and found three different hypothesized 9-mers encoding ATGATCAAG in Vibrio cholerae, CCTACCACC in Thermotoga petrophila, and TTATCCACA in E. coli. We must warn you that finding ori is often more complex than in the three examples we considered. Some bacteria have even fewer DnaA boxes than E. coli, making it difficult to identify them. The ter region is often located not directly opposite to ori but may be significantly shifted, resulting in reverse and forward half-strands having substantially different lengths. The position of the skew minimum is often only a rough indicator of ori position, which forces researchers to expand their windows when searching for DnaA boxes, bringing in extraneous repeated substrings. Finally, skew diagrams do not always look as nice as that of E. coli; for example, the skew diagram for Thermotoga petrophila is shown below.


Bibliography Notes
Using the skew to find replication origins was first proposed by Lobry, 1996 and also described in Grigoriev, 1998. Grigoriev, 2011 provides an excellent introduction to the skew approach, and Sernova and Gelfand, 2008 gave a review of algorithms and software tools for finding replication origins in bacteria. Lundgren et al., 2004 demonstrated that archaea may have multiple ori. Wang et al., 2011 inserted an artificial ori into the E. coli genome and showed that it triggers replication. Xia, 2012 was the first to conjecture that bacteria may have multiple replication origins. Gao and Zhang, 2008 developed the Ori-Finder software program for finding bacterial replication origins.

Liachko et al., 2013 provided the most comprehensive description of the replication originsofyeast. Solov’ev, 1966 was the first to derive accurate formulas for approximating the probabilities of patterns in a string. Gardner, 1974 wrote an excellent introductory article about the Best Bet for Simpletons paradox. Guibas and Odlyzko, 1981 provided an excellent coverage of the overlapping words paradox that illustrates the complexity of computing the probabilities of patterns in a random text. They also derived a rather complicated proof of Conway’s formula for Best Bet for Simpletons. Sedgewick and Flajolet, 2013 gave an overview of various approaches for computing the probabilities of patterns in a string.
"The availability of hundreds of complete bacterial genomes has created new challenges and simultaneously opportunities for bioinformatics. In the area of statistical analysis of genomic sequences, the studies of nucleotide compositional bias and gene bias between strands and replichores paved way to the development of tools for prediction of bacterial replication origins. Only a few (about 20) origin regions for eubacteria and archaea have been proven experimentally. One reason for that may be that this is now considered as an essentially bioinformatics problem, where predictions are sufficiently reliable not to run labor-intensive experiments, unless specifically needed. Here we describe the main existing approaches to the identification of replication origin (oriC) and termination (terC) loci in prokaryotic chromosomes and characterize a number of computational tools based on various skew types and other types of evidence. We also classify the eubacterial and archaeal chromosomes by predictability of their replication origins using skew plots. Finally, we discuss possible combined approaches to the identification of the oriC sites that may be used to improve the prediction tools, in particular, the analysis of DnaA binding sites using the comparative genomic methods."
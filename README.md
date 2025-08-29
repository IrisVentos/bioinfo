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

You should write your code implementing PatternCount first; you can choose any programming language you like.
When you click "Download Dataset", you will receive a randomized dataset. In this problem, the dataset will contain two lines: the first line contains Text, and the second line contains Pattern. In general, the "Sample Input" section shows how 
Run your program (in the programming language of your choice) on the dataset, and then return the output of your program in the text field below. (Please do not enter your code in the browser.)
There is a time limit on each problem to ensure that your code is sufficiently efficient to return the correct output quickly.
You can see how you should format your answer by looking at the sample output.
You have unlimited attempts to answer the question, but each time you click "Try Again", you will need to download a new dataset.
We also provide additional small datasets (see link below) to help you debug your code.


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

If |Text| and k are small, as is the case when looking for DnaA boxes in the typical bacterial ori, then an algorithm with running time of O(|Text|2 · k) is perfectly acceptable. But once we find some new biological application requiring us to solve the Frequent Words Problem for a very long Text, we will quickly run into trouble. What can we do instead?

We know that an array of length n is an ordered table of values, where we access the values using the integer indices 0 through n-1. The frequency table is a generalized version of an array called a map or dictionary for which the indices are allowed to be arbitrary values (in this case, they are strings). More precisely, the indices of a map are called keys.

Given a map dict, we can access the value associated with a key key using the notation dict[key]. In the case of a frequency table called freq, we can access the value associated with some key string pattern using the notation freq[pattern]. The following pseudocode function takes a string text and an integer k as input and returns their frequency table as a map of string keys to integer values.

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

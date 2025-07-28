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

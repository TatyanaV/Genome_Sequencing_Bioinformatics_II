'''
CODE CHALLENGE: Solve the Overlap Graph Problem (restated below).
     Input: A collection Patterns of k-mers.
     Output: The overlap graph Overlap(Patterns), in the form of an adjacency list. (You may return the edges in any order.)
      Sample Input:

ATGCG
GCATG
CATGC
AGGCA
GGCAT

Sample Output:

CATGC -> ATGCG
GCATG -> CATGC
GGCAT -> GCATG
AGGCA -> GGCAT
'''

patterns = [s for s in"""
ATGCG
GCATG
CATGC
AGGCA
GGCAT
""".split() if s]

for iterator1, pattern1 in enumerate(patterns):
    for iterator2, pattern2 in enumerate(patterns):
        if pattern1[1:] == pattern2[:len(pattern2) - 1] and iterator1 != iterator2:
            print (pattern1 + ' -> ' + pattern2)
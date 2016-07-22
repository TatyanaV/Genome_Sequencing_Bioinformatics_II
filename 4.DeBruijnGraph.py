'''
CODE CHALLENGE: Solve the De Bruijn Graph from a String Problem.
     Input: An integer k and a string Text.
     Output: DeBruijnk(Text), in the form of an adjacency list.
      Sample Input:

4
AAGATTCTCTAAGA

Sample Output:

AAG -> AGA,AGA
AGA -> GAT
ATT -> TTC
CTA -> TAA
CTC -> TCT
GAT -> ATT
TAA -> AAG
TCT -> CTA,CTC
TTC -> TCT
'''

from collections import defaultdict
k = 4
text = """
AAGATTCTCTAC
"""
print (text)
tuples = []
for i in range(1, len(text) - k):
    tuples.append((text[i: i+k-1], text[i+1: i+k]))
dd = defaultdict(set)
for t in tuples:
    dd[t[0]].add(t[1])

print (*(sorted([key + ' -> ' + ','.join([v for v in value])
                       for key, value in dd.items()])),sep = '\n')


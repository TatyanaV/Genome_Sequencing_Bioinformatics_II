'''
DeBruijn Graph from k-mers Problem: Construct the de Bruijn graph from a set of k-mers.
     Input: A collection of k-mers Patterns.
     Output: The adjacency list of the de Bruijn graph DeBruijn(Patterns).

CODE CHALLENGE: Solve the de Bruijn Graph from k-mers Problem.
 Sample Input:

GAGG
CAGG
GGGG
GGGA
CAGG
AGGG
GGAG

Sample Output:

AGG -> GGG
CAG -> AGG,AGG
GAG -> AGG
GGA -> GAG
GGG -> GGA,GGG
'''
#https://github.com/AnnaUfliand/Bioinformatics/blob/18c9c385c8f62c2895dc029ec5a07c416a6363fa/HW4/ReconstructString.py
from collections import defaultdict
kmers = """
GAGG
CAGG
GGGG
GGGA
CAGG
AGGG
GGAG
""".split()
'''
tuples = [(kmer[:len(kmer) - 1], kmer[1:]) for kmer in kmers]
dd = defaultdict(set)
for t in tuples:
    dd[t[0]].add(t[1])
print (*(sorted([key + ' -> ' + ','.join(sorted([v for v in value]))
                       for key, value in dd.items()])),sep ='\n')

'''
print(kmers)
composition = kmers

def deBruijnGraph(kmers):
    k = len(kmers[0])
    graph = {}
    for i in range(len(kmers)):
        try:
            graph[kmers[i][:-1]].append(kmers[i][1:])
        except:
            graph[kmers[i][:-1]] = [kmers[i][1:]]
    return graph

graph = deBruijnGraph(composition)
print(graph)
print (*(sorted([key + ' -> ' + ','.join(sorted([v for v in value]))
                       for key, value in graph.items()])),sep ='\n')
#output = out_string.split('->')
'''
CODE CHALLENGE: Solve the String Reconstruction Problem.
     Input: An integer k followed by a list of k-mers Patterns.
     Output: A string Text with k-mer composition equal to Patterns. (If multiple answers exist, you may
     return any one.)

Sample Input:
     4
     CTTA
     ACCA
     TACC
     GGCT
     GCTT
     TTAC

Sample Output:
     GGCTTACCA
'''
from collections import defaultdict

from collections import defaultdict

kmers = """
CTTA
ACCA
TACC
GGCT
GCTT
TTAC
""".split()
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
'''
kmers = """
CTTA
ACCA
TACC
GGCT
GCTT
TTAC
""".split()
tuples = [(kmer[:len(kmer) - 1], kmer[1:]) for kmer in kmers]
dd = defaultdict(set)
for t in tuples:
    dd[t[0]].add(t[1])
print (*(sorted([key + ' -> ' + ','.join(sorted([v for v in value]))
                       for key, value in dd.items()])),sep ='\n')
sequence =  ([key + ',' + ','.join(sorted([v for v in value]))
                       for key, value in dd.items()])
print (sequence)

#edges = [tuple(edge.split(' -> ')) for edge in input if edge]
'''
input = """
CTT -> TTA
ACC -> CCA
TAC -> ACC
GGC -> GCT
GCT -> CTT
TTA -> TAC
""".split('\n')

#with open('/Users/dbarabanov/Downloads/dataset_57_6.txt') as f:
#    input = [s.strip() for s in f.readlines()]
#with open('string_reconstruction_input.txt') as f:
    #input = [s.strip() for s in f.readlines()]
'''
edges = [tuple(edge.split(' -> ')) for edge in input if edge]
print("from the new code")
print (edges)
edges = [(t[0], [i for i in t[1].split(',')]) for t in edges]
print(edges)
graph = {x: y for x, y in edges}
print (graph)
'''
degrees = defaultdict(int)
for k in graph:
    for v in graph[k]:
        degrees[k] += 1
        degrees[v] -= 1
source = [k for k, v in degrees.items() if v == 1][0]
sinc = [k for k, v in degrees.items() if v == -1][0]
list(graph)
start = list(graph)[0]
# not sure what to do with this start = graph.keys()[0]
#print 'source: %s, sinc: %s' % (source, sinc)

if sinc in graph.keys():
    graph[sinc].append(source)
else:
    graph[sinc] = [source]

cycles = {}
while graph:
    current = next(iter(graph))
    cycle = [current]
    cycles[current] = cycle
    while current in graph:
        _next = graph[current][0]
        del graph[current][0]
        if len(graph[current]) == 0:
            del graph[current]
        current = _next
        cycle.append(_next)


def traverse(tree, root):
    out = []
    for r in tree[root]:
        if r != root and r in tree:
            out += traverse(tree, r)
        else:
            out.append(r)
    return out

cycle = traverse(cycles, start)
for i in range(1, len(cycle)):
    if cycle[i-1] == sinc and cycle[i] == source:
        boarder = i
path = cycle[boarder:]+cycle[1:boarder]
print (*((([s[0] for s in list(path)]) + list(path[-1][1:]))),sep = '')
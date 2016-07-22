'''
An Euler path is a path that uses every edge of a graph exactly once.
CODE CHALLENGE: Solve the Eulerian Path Problem.
     Input: The adjacency list of a directed graph that has an Eulerian path.
     Output: An Eulerian path in this graph.

Sample Input:
     0 -> 2
     1 -> 3
     2 -> 1
     3 -> 0,4
     6 -> 3,7
     7 -> 8
     8 -> 9
     9 -> 6

Sample Output:
     6->7->8->9->6->3->0->2->1->3->4
'''

from collections import defaultdict
input = """
     0 -> 2
     1 -> 3
     2 -> 1
     3 -> 0,4
     6 -> 3,7
     7 -> 8
     8 -> 9
     9 -> 6
""".split('\n')

#with open('/Users/dbarabanov/Downloads/dataset_57_5 (3).txt') as f:
#    input = f.readlines()
#with open('eulerian_path.txt') as f:
    #input = f.readlines()
edges = [tuple(edge.split(' -> ')) for edge in input if edge]
print(edges)
edges = [(int(t[0]), [int(i) for i in t[1].split(',')]) for t in edges]
print(edges)

graph = {x: y for x, y in edges}
print(graph)
degrees = defaultdict(int)
for k in graph:
    for v in graph[k]:
        degrees[k] += 1
        degrees[v] -= 1
source = [k for k, v in degrees.items() if v == 1][0]
sinc = [k for k, v in degrees.items() if v == -1][0]
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
        next_ = graph[current][0]
        del graph[current][0]
        if len(graph[current]) == 0:
            del graph[current]
        current = next_
        cycle.append(next_)


def traverse(tree, root):
    out = []
    for r in tree[root]:
        if r != root and r in tree:
            out += traverse(tree, r)
        else:
            out.append(r)
    return out

cycle = traverse(cycles, 0)
for i in range(1, len(cycle)):
    if cycle[i-1] == sinc and cycle[i] == source:
        boarder = i
path = cycle[boarder:]+cycle[1:boarder]
print ('->'.join([str(i) for i in path]))
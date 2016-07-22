'''
 MaximalNonBranchingPaths(Graph)
        Paths ? empty list
        for each node v in Graph
            if v is not a 1-in-1-out node
                if out(v) > 0
                    for each outgoing edge (v, w) from v
                        NonBranchingPath ? the path consisting of the single edge (v, w)
                        while w is a 1-in-1-out node
                            extend NonBranchingPath by the outgoing edge (w, u) from w
                            w ? u
                        add NonBranchingPath to the set Paths
        for each isolated cycle Cycle in Graph
            add Cycle to Paths
        return Paths

CODE CHALLENGE: Implement MaximalNonBranchingPaths.
     Input: The adjacency list of a graph whose nodes are integers.
     Output: The collection of all maximal nonbranching paths in this graph.
 Sample Input:

1 -> 2
2 -> 3
3 -> 4,5
6 -> 7
7 -> 6

Sample Output:

1 -> 2 -> 3
3 -> 4
3 -> 5
7 -> 6 -> 7
'''

import sys
import random
from copy import deepcopy
from collections import Counter

input = """
1 -> 2
2 -> 3
3 -> 4,5
6 -> 7
7 -> 6
"""
edges = [tuple(edge.split(' -> ')) for edge in input.split('\n') if edge]
print (edges)
edges = [(int(t[0]), [int(i) for i in t[1].split(',')]) for t in edges]
print(edges)
graph = {x: y for x, y in edges}
print (graph)

def maximalNonBranchingPathsInTrie(graph, start_node=-1, end_node=-1):

    paths = []
    dict = {}
    np = []

    in_keys = graph.keys ()
    in_keys.sort ()

    #print "GRAPH", graph
    #print "LEN OF KEEYS", len (graph.keys ())
    for node in graph.keys ():
        if len (graph[node]) > 1:

            nbp = []
            for n in graph[node]:
                nbp = [node]
                nbp.append(n)
                w = n

                while graph.has_key (w) and (len (graph[w]) == 1):
                    u = graph[w][0]
                    nbp.append(u)
                    w = u
                paths.append(nbp)
                np.append(nbp[1:])
            dict[node] = np[:]
            np = []
            continue


    #print "PATSH", paths
    #for i in paths:
    #    print ' -> '.join (str (x) for x in i)
    return  dict

print(maximalNonBranchingPathsInTrie(graph, start_node=-1, end_node=-1))
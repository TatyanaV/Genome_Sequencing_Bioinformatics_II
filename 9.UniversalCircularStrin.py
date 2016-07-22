'''
https://github.com/AnnaUfliand/Bioinformatics/blob/18c9c385c8f62c2895dc029ec5a07c416a6363fa/HW4/UniversalString.py
CODE CHALLENGE: Solve the k-Universal Circular String Problem.
     Input: An integer k.
     Output: A k-universal circular string.

Sample Input:
     4

Sample Output:
     0000110010111101
'''
ALPHABET = ['1', '0']

import random
def deBruijnGraph(kmers):
    k = len(kmers[0])
    graph = {}
    for i in range(len(kmers)):
        try:
            graph[kmers[i][:-1]].append(kmers[i][1:])
        except:
            graph[kmers[i][:-1]] = [kmers[i][1:]]
    return graph
def findCycle(graph, node, previous_path):
    path = []

    if len(previous_path) != 0:
        node_index = 0
        for i in range(len(previous_path)):
            if previous_path[i][0] == node:
                node_index = i
                break
        for i in reversed(range(len(previous_path))):
            path.append(previous_path[node_index - 1 - i])

    edge = (node, graph[node][random.randint(0, len(graph[node]) - 1)])

    while edge not in path and len(graph) != 0:
        path.append(edge)
        if len(graph[edge[0]]) > 1:
            graph[edge[0]].remove(edge[1])
        else:
            del graph[edge[0]]
        current_point = edge[1]
        if graph.get(current_point):
            edge = (current_point, graph[current_point][random.randint(0, len(graph[current_point]) - 1)])
    return len(graph), path


def findWholeCycle(graph):
    l = len(graph)
    node = sorted(graph.keys())[random.randint(0, len(graph) - 1)]
    path = []
    while l != 0:
        l, path = findCycle(graph, node, path)
        nodes = []
        if l != 0 and len(path):
            for edge in path:
                if edge[0] in graph:
                    nodes.append(edge[0])
            node = nodes[random.randint(0, len(nodes) - 1)]
    return path

def constructGenome(output):
    genome = output[0]
    for i in range(1, len(output)):
        genome += output[i][-1]
    return genome

def constructAllNumbers(k):
    combinations = []
    if k == 1:
        return ALPHABET
    for a in ALPHABET:
        suffixes = constructAllNumbers(k - 1)
        for suffix in suffixes:
            combinations.append(a + suffix)
    return combinations

k = 8
composition = constructAllNumbers(k)
graph = deBruijnGraph(composition)
output = findWholeCycle(graph)
print (constructGenome([i[0] for i in output])[:-k + 2])
'''
https://github.com/stevehaigh/Bioinformatics/blob/229e4f68efb75c67d14bc914bc63e5e4490c2c6a/Chapter%204/OverlapGraph.py
'''

from io import open
import sys
from builtins import len, dict, range


class Node:
    def __init__(self, label, edges):
        self.label = label
        self.edges = edges
        self.in_degree = 0


def calc_indegrees(graph):
    new_nodes = []

    for node in graph.values():
        for edge in node.edges:
            # might find a node with no out going edges this way...
            if edge in graph:
                graph[edge].in_degree += 1
            else:
                new_nodes.append(Node(edge, []))

    for new_node in new_nodes:
        if new_node.label in graph:
            graph[new_node.label].in_degree += 1
        else:
            new_node.in_degree = 1
            graph[new_node.label] = new_node

    return graph


def build_graph(lines):
    """
    convert list of lines of the form 1 -> 2,3
    into a graph in the form of a dict with a list.
    :param lines:
    :return:
    """
    result = dict()

    for line in lines:
        pair = line.split("->")
        if len(pair) == 2:
            from_node = pair[0].strip(' ')
            to_nodes = pair[1].strip(' ').split(",")

            if from_node in result:
                result[from_node].edges.extend(to_nodes)
            else:
                result[from_node] = Node(from_node, to_nodes)

    result = calc_indegrees(result)
    return result


def find_start_node(graph):
    """
    start node has more edges than in-degree.
    :param graph:
    :return:
    """
    more_out_than_in = []

    for node in graph.values():
        if node.in_degree < len(node.edges):
            more_out_than_in.append(node)

    min_in = len(graph)
    min_node = None

    for node in more_out_than_in:
        if node.in_degree < min_in:
            min_in = node.in_degree
            min_node = node

    return min_node


def get_path_from_node(current_node, graph, cycle):
    while True:
        if len(current_node.edges) == 0:
            break
        next_node = graph[current_node.edges[0]]
        current_node.edges = current_node.edges[1:]
        cycle.append(next_node)
        current_node = next_node

    return cycle


def find_cycle(graph):
    # select any node to start cycle, may as well be the first
    cycle = []
    current_node = find_start_node(graph)
    # list(graph.values())[0]
    cycle.append(current_node)

    # do initial cycle
    cycle = get_path_from_node(current_node, graph, cycle)
    done = False

    while not done:
        c = cycle[:]
        for i in range(0, len(c)):
            if len(c[i].edges) > 0:
                cycle_prime = []
                cycle_prime.append(c[i])
                cycle_prime = get_path_from_node(c[i], graph, cycle_prime)
                cycle = cycle[:i] + cycle_prime + cycle[i + 1:]
                break

                # if we get here an there are no more nodes then we're done
            if i == len(c) - 1:
                done = True

    return cycle


def main(argv=None):
    """
    :param argv: the command line args
    :return: nothing
    """
    if argv is None:
        argv = sys.argv

    with open("cycles.txt") as contents:
        lines = [line.rstrip('\n') for line in contents]

    graph = build_graph(lines)

    cycle = find_cycle(graph)

    temp = []

    for node in cycle:
        temp.append(node.label)

    with open("eulerian_path.txt", "w") as text_file:
        text_file.write("->".join(temp))


if __name__ == "__main__":
    sys.exit(main())
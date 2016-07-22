'''
https://github.com/hkpcmit/BioInfo/tree/6b245435515102f1ddc3c74003d614f4296c2563
'''
import collections


class Error(Exception):
    """Base error class."""


def CircularStringOld(k):
    formatter = '{{:0{}b}}'.format(k)
    kmers = [formatter.format(i) for i in xrange(2**k)]
    graph = KmerDeBruijnGraph(kmers)
    out_degrees = {key: set(values) for key, values in graph.iteritems()}
    cycle = EulerCycleHelper(out_degrees, start='0'*(k-1))
    string = cycle[0]
    for node in cycle[1:(1-k)]:
        string += node[-1]
    return string


def CircularString(k):
    formatter = '{{:0{}b}}'.format(k)
    kmers = [formatter.format(i) for i in xrange(2**k)]
    graph = KmerDeBruijnGraph(kmers)
    cycle = EulerCycleHelper(graph, start='0'*(k-1))
    string = cycle[0]
    for node in cycle[1:(1-k)]:
        string += node[-1]
    return string


def Contigs(kmers):
    adj, one_in_one_outs, out_nodes = ContigsAdjacencies(kmers)
    paths = []
    for node in out_nodes:
        node_list = [node]
        while True:
            next_node = adj[node][-1]
            node_list.append(next_node)
            if len(adj[node]) == 1:
                del adj[node]
            else:
                adj[node] = adj[node][:-1]
            if next_node not in one_in_one_outs:
                break
            node = next_node
        path = node_list[0]
        for node in node_list[1:]:
            path += node[-1]
        paths.append(path)
    if adj:
        raise Error('Unused adjacencies: {}'.format(adj))
    return sorted(paths)


def ContigsAdjacencies(kmers):
    adj = collections.defaultdict(list)
    in_degrees, out_degrees = collections.Counter(), collections.Counter()
    out_nodes, one_in_one_outs = [], set()
    for kmer in kmers:
        source, dest = kmer[:-1], kmer[1:]
        out_degrees[source] += 1
        in_degrees[dest] += 1
        adj[source].append(dest)
    for node in adj:
        if (in_degrees[node] == 1) and (out_degrees[node] == 1):
            one_in_one_outs.add(node)
        else:
            out_nodes.extend([node] * len(adj[node]))
    return adj, one_in_one_outs, out_nodes


def DeBruijnGraph(k, text):
    if not k:
        raise Error('Zero k.')
    text_length = len(text)
    if k >= text_length:
        raise Error('Text length is too short.')
    graph = collections.defaultdict(list)
    for i in xrange(len(text)-k+1):
        node1, node2 = text[i:i+k-1], text[i+1:i+k]
        graph[node1].append(node2)
    return graph


def EulerCycle(segments):
    out_degrees = EulerCycleOutDegrees(segments)
    return EulerCycleHelper(out_degrees)


def EulerCycleGetCycleOld(out_degrees, start):
    cycle = [start]
    this_node = start
    while this_node in out_degrees:
        neighbors = out_degrees[this_node]
        next_node = next(node for node in neighbors)
        cycle.append(next_node)
        if len(neighbors) == 1:
            del out_degrees[this_node]
        else:
            out_degrees[this_node] = neighbors - set([next_node])
        if next_node == start:
            return cycle
        this_node = next_node
    raise Error('GetCycle: Invalid out_degrees: {}'.format(out_degrees))


def EulerCycleGetCycle(out_degrees, start):
    cycle = [start]
    this_node = start
    while this_node in out_degrees:
        neighbors = out_degrees[this_node]
        next_node = neighbors[-1]
        cycle.append(next_node)
        if len(neighbors) == 1:
            del out_degrees[this_node]
        else:
            out_degrees[this_node] = neighbors[:-1]
        if next_node == start:
            return cycle
        this_node = next_node
    raise Error('GetCycle: Invalid out_degrees: {}'.format(out_degrees))


def EulerCycleHelper(out_degrees, start=None):
    if not start:
        start = next(node for node, neighbors in out_degrees.iteritems()
                     if len(neighbors) > 1)
    old_cycle = None
    while out_degrees:
        cycle = EulerCycleGetCycle(out_degrees, start)
        if old_cycle:
            cycle = EulerCycleJoinCycles([old_cycle, cycle])
            if len(cycle) <= len(old_cycle):
                raise Error('Bad join!')
        if not out_degrees:
            break
        start = EulerCyclePickNextNode(out_degrees, cycle)
        old_cycle = cycle
    return cycle


def EulerCycleJoinCycles(cycle_list):
    if len(cycle_list) == 1:
        return cycle_list[0]
    remain_cycle = EulerCycleJoinCycles(cycle_list[1:])
    i = cycle_list[0].index(remain_cycle[0])
    return cycle_list[0][:i] + remain_cycle + cycle_list[0][i+1:]


def EulerCycleOutDegreesOld(segments):
    out_degrees = collections.defaultdict(set)
    for segment in segments:
        node, neighbors = segment.split(' -> ')
        out_degrees[node] = set(neighbors.split(','))
    return out_degrees


def EulerCycleOutDegrees(segments):
    out_degrees = collections.defaultdict(set)
    for segment in segments:
        node, neighbors = segment.split(' -> ')
        out_degrees[node] = neighbors.split(',')
    return out_degrees


def EulerCyclePickNextNode(out_degrees, cycle):
    for node in cycle[:-1]:
        if node in out_degrees:
            return node
    raise Error('PickNextNode: Invalid out_degrees: {}'.format(out_degrees))


def EulerPathOld(segments):
    out_degrees, extra_in, extra_out = EulerPathProcessSegments(segments)
    out_degrees[extra_in[0]].add(extra_out[0])
    cycle = EulerCycleHelper(out_degrees, start=extra_out[0])
    if cycle[-2] != extra_in[0]:
        i = cycle.index(extra_in[0])
        return cycle[i+1:] + cycle[1:i+1]
    return cycle[:-1]


def EulerPath(segments):
    out_degrees, extra_in, extra_out = EulerPathProcessSegments(segments)
    out_degrees[extra_in[0]].append(extra_out[0])
    cycle = EulerCycleHelper(out_degrees, start=extra_out[0])
    if cycle[-2] != extra_in[0]:
        i = cycle.index(extra_in[0])
        return cycle[i+1:] + cycle[1:i+1]
    return cycle[:-1]


def EulerPathProcessSegmentsOld(segments):
    out_degrees = collections.defaultdict(set)
    degree_counter = collections.Counter()
    for segment in segments:
        node, neighbors = segment.split(' -> ')
        out_degrees[node] = set(neighbors.split(','))
        degree_counter[node] += len(out_degrees[node])
        for neighbor in out_degrees[node]:
            degree_counter[neighbor] -= 1
    extra_in, extra_out = [], []
    for node, ctr in degree_counter.iteritems():
        if ctr > 0:
            extra_out.append(node)
        elif ctr < 0:
            extra_in.append(node)
    return out_degrees, extra_in, extra_out


def EulerPathProcessSegments(segments):
    out_degrees = collections.defaultdict(list)
    degree_counter = collections.Counter()
    for segment in segments:
        node, neighbors = segment.split(' -> ')
        out_degrees[node] = neighbors.split(',')
        degree_counter[node] += len(out_degrees[node])
        for neighbor in out_degrees[node]:
            degree_counter[neighbor] -= 1
    extra_in, extra_out = [], []
    for node, ctr in degree_counter.iteritems():
        if ctr > 0:
            extra_out.append(node)
        elif ctr < 0:
            extra_in.append(node)
    return out_degrees, extra_in, extra_out


def EulerReconstructString(k, kmers):
    graph = KmerDeBruijnGraph(kmers)
    segments = [' -> '.join([node, ','.join(neighbors)])
                for node, neighbors in graph.iteritems()]
    path = EulerPath(segments)
    string = path[0]
    for node in path[1:]:
        string += node[-1]
    return string


def GenomePathString(genomes):
    if not genomes:
        raise Error('Empty genomes input.')
    result = genomes[0]
    genome_length = len(result)
    for genome in genomes[1:]:
        if not result.endswith(genome[:genome_length-1]):
            return result
        result += genome[-1]
    return result


def KmerDeBruijnGraph(kmers):
    if not kmers:
        raise Error('Empty kmers input.')
    graph = collections.defaultdict(list)
    for kmer in kmers:
        graph[kmer[:-1]].append(kmer[1:])
    return graph


def OverlapGraph(patterns):
    if not patterns:
        raise Error('Empty patterns input.')
    graph = collections.defaultdict(list)
    for i, p1 in enumerate(patterns):
        for j in xrange(i+1, len(patterns)):
            p2 = patterns[j]
            if p1.endswith(p2[:-1]):
                graph[p1].append(p2)
            elif p2.endswith(p1[:-1]):
                graph[p2].append(p1)
    return graph


def ReadPairsAdjacenciesOld(pairs):
    adj = collections.defaultdict(set)
    degree_counter = collections.Counter()
    for pair in pairs:
        nodes = tuple(pair.split('|'))
        orig_keys = degree_counter.keys()
        degree_counter[nodes] = 0
        for orig_nodes in orig_keys:
            if (orig_nodes[0].endswith(nodes[0][:-1]) and
                orig_nodes[1].endswith(nodes[1][:-1])):
                adj[orig_nodes].add(nodes)
                degree_counter[orig_nodes] += 1
                degree_counter[nodes] -= 1
            if (nodes[0].endswith(orig_nodes[0][:-1]) and
                nodes[1].endswith(orig_nodes[1][:-1])):
                adj[nodes].add(orig_nodes)
                degree_counter[nodes] += 1
                degree_counter[orig_nodes] -= 1
    extra_in, extra_out = [], []
    for nodes, freq in degree_counter.iteritems():
        if freq < 0:
            extra_in.append(nodes)
        if freq > 0:
            extra_out.append(nodes)
    return adj, extra_in, extra_out


def ReadPairsAdjacencies(pairs):
    adj = collections.defaultdict(list)
    degree_counter = collections.Counter()
    for pair in pairs:
        nodes = tuple(pair.split('|'))
        source = (nodes[0][:-1], nodes[1][:-1])
        dest = (nodes[0][1:], nodes[1][1:])
        adj[source].append(dest)
        degree_counter[source] += 1
        degree_counter[dest] -= 1
    extra_in, extra_out = [], []
    for nodes, freq in degree_counter.iteritems():
        if freq < 0:
            extra_in.append(nodes)
        if freq > 0:
            extra_out.append(nodes)
    return adj, extra_in, extra_out


def ReconstructString(k, text):
    return sorted(text[i:i+k] for i in xrange(len(text)-k+1))


def ReconstructStringReadPairsOld(k, d, pairs):
    adj, extra_in, extra_out = ReadPairsAdjacencies(pairs)
    adj[extra_in[0]].add(extra_out[0])
    cycle = EulerCycleHelper(adj, start=extra_out[0])
    string1, string2 = cycle[0]
    for nodes in cycle[1:-1]:
        string1 += nodes[0][-1]
        string2 += nodes[1][-1]
    common = string1[k+d:]
    if not string2.startswith(common):
        raise Error('Invalid string1: {}, string2: {}'.format(string1, string2))
    return string1 + string2[len(common):]


def ReconstructStringReadPairs(k, d, pairs):
    adj, extra_in, extra_out = ReadPairsAdjacencies(pairs)
    adj[extra_in[0]].append(extra_out[0])
    cycle = EulerCycleHelper(adj, start=extra_out[0])
    string1, string2 = cycle[0]
    for nodes in cycle[1:-1]:
        string1 += nodes[0][-1]
        string2 += nodes[1][-1]
    common = string1[k+d:]
    if not string2.startswith(common):
        raise Error('Invalid string1: {}, string2: {}'.format(string1, string2))
    return string1 + string2[len(common):]
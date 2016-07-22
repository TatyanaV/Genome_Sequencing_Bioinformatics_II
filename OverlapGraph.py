from builtins import dict, open
import sys


def build_prefix_graph(kmers):
    graph = dict()

    for kmer in kmers:
        prefix = kmer[:-1]
        if prefix in graph:
            graph[prefix].append(kmer)
        else:
            graph[prefix] = [kmer]

    return graph


def find_overlaps(kmers):
    """
    :param kmers: list of kmers
    :return: a dict of kmer -> list-of-overlaps-for-kmer
    """

    # build dict of suffixes to kmers
    graph = build_prefix_graph(kmers)

    result = dict()

    for kmer in kmers:
        suffix = kmer[1:]
        if kmer in result:
            if suffix in graph:
                result[kmer].extend(graph[suffix])
        else:
            if suffix in graph:
                result[kmer] = graph[suffix]

    return result


def main(argv=None):
    """
    :param argv: the command line args
    :return: nothing
    """
    if argv is None:
        argv = sys.argv

    with open("kmers-2") as contents:
        kmers = [line.rstrip('\n') for line in contents]

    result = find_overlaps(kmers)

    # for key, value in result.items():
    # print(key + " -> " + " ".join(value))

    with open("overlaps.txt", "w") as text_file:
        for key, value in result.items():
            text_file.write(key + " -> " + " ".join(value) + "\n")


if __name__ == "__main__":
    sys.exit(main())
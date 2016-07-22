#!/usr/bin/python

'''
https://github.com/hkpcmit/BioInfo/blob/6b245435515102f1ddc3c74003d614f4296c2563/reconstruct_string1_input.txt
AWESOME WEBSITE
'''

from AssembleGenomes import ReconstructString, GenomePathString, OverlapGraph
from AssembleGenomes import DeBruijnGraph, KmerDeBruijnGraph, EulerCycle, EulerPath
from AssembleGenomes import EulerReconstructString, CircularString
from AssembleGenomes import ReconstructStringReadPairs, Contigs
from SequenceAntiBiotics_test import ReadFile, WriteFile
import collections
import unittest


class ReconstructStringTest(unittest.TestCase):

    def ReadInput(self, filename):
        input = ReadFile(filename)
        return int(input[0]), input[1]

    def testData1(self):
        k = 5
        text = 'CAATCCAAC'
        output = ReconstructString(k, text)
        expect = ['AATCC', 'ATCCA', 'CAATC', 'CCAAC', 'TCCAA']
        self.assertEqual(expect, output)

    def testData2(self):
        k, text = self.ReadInput('string_com_input.txt')
        output = ReconstructString(k, text)
        expect = ReadFile('string_com_output.txt')
        self.assertEqual(expect, output)

    def testData3(self):
        k, text = self.ReadInput('dataset_197_3.txt')
        output = ReconstructString(k, text)
        expect = ReadFile('dataset_197_3_output.txt')
        self.assertEqual(expect, output)


class GenomePathStringTest(unittest.TestCase):

    def testData1(self):
        genomes = ['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT']
        output = GenomePathString(genomes)
        expect = 'ACCGAAGCT'
        self.assertEqual(expect, output)

    def testData2(self):
        genomes = ReadFile('GenomePathString_input.txt')
        output = GenomePathString(genomes)
        expect = ReadFile('GenomePathString_output.txt')[0]
        self.assertEqual(expect, output)

    def testData3(self):
        genomes = ReadFile('dataset_198_3.txt')
        output = GenomePathString(genomes)
        expect = ReadFile('dataset_198_3_output.txt')[0]
        self.assertEqual(expect, output)


class OverlapGraphTest(unittest.TestCase):

    def FormatGraph(self, graph):
        return ['{} -> {}'.format(key, graph[key][0])
                for key in sorted(graph)]

    def testData1(self):
        patterns = ['ATGCG', 'GCATG', 'CATGC', 'AGGCA', 'GGCAT']
        graph = OverlapGraph(patterns)
        expect = ReadFile('overlap_graph1_output.txt')
        self.assertEqual(expect, self.FormatGraph(graph))

    def testData2(self):
        patterns = ReadFile('overlap_graph_1_input.txt')
        graph = OverlapGraph(patterns)
        expect = ReadFile('overlap_graph_1_output.txt')
        self.assertEqual(expect, self.FormatGraph(graph))

    def testData3(self):
        patterns = ReadFile('dataset_198_9.txt')
        graph = OverlapGraph(patterns)
        expect = ReadFile('dataset_198_9_output.txt')
        output = self.FormatGraph(graph)
        self.assertEqual(expect, output)


class DeBruijnGraphTest(unittest.TestCase):

    def ReadInput(self, filename):
        input = ReadFile(filename)
        return int(input[0]), input[1]

    def FormatGraph(self, graph):
        return ['{} -> {}'.format(key,
                                  ','.join(sorted(graph[key])))
                for key in sorted(graph)]

    def testData1(self):
        k = 4
        text = 'AAGATTCTCTAAGA'
        graph = DeBruijnGraph(k, text)
        expect = ReadFile('debruijn_graph1_output.txt')
        output = self.FormatGraph(graph)
        self.assertEqual(expect, output)

    def testData2(self):
        k, text = self.ReadInput('debruijn_graph_string_input.txt')
        graph = DeBruijnGraph(k, text)
        expect = ReadFile('debruijn_graph_string_output.txt')
        output = self.FormatGraph(graph)
        self.assertEqual(expect, output)

    def testData3(self):
        k, text = self.ReadInput('dataset_199_6.txt')
        graph = DeBruijnGraph(k, text)
        expect = ReadFile('dataset_199_6_output.txt')
        output = self.FormatGraph(graph)
        self.assertEqual(expect, output)

    def testData4(self):
        kmers = ['GAGG', 'GGGG', 'GGGA', 'CAGG', 'AGGG', 'GGAG']
        graph = KmerDeBruijnGraph(kmers)
        expect = ReadFile('debruijn_graph2_output.txt')
        output = self.FormatGraph(graph)
        self.assertEqual(expect, output)

    def testData5(self):
        kmers = ReadFile('debruijn_graph_kmers_input.txt')
        graph = KmerDeBruijnGraph(kmers)
        expect = ReadFile('debruijn_graph_kmers_output.txt')
        output = self.FormatGraph(graph)
        self.assertEqual(expect, output)

    def testData6(self):
        kmers = ReadFile('dataset_200_7.txt')
        graph = KmerDeBruijnGraph(kmers)
        expect = ReadFile('dataset_200_7_output.txt')
        output = self.FormatGraph(graph)
        self.assertEqual(expect, output)


class EulerCycleTest(unittest.TestCase):

    def CheckOutput(self, filename, output):
        expect = ReadFile(filename)[0]
        i = output[1:].index(output[0])
        subset = output[1:i+2]
        self.assertIn('->'.join(subset), expect)

    def testData1(self):
        segments = ReadFile('euler_cycle_input.txt')
        output = EulerCycle(segments)
        self.CheckOutput('euler_cycle_output.txt', output)

    def testData2(self):
        segments = ReadFile('eulerian_cycle_input.txt')
        output = EulerCycle(segments)
        self.CheckOutput('eulerian_cycle_output.txt', output)

    def testData3(self):
        segments = ReadFile('dataset_203_2.txt')
        output = EulerCycle(segments)
        self.CheckOutput('dataset_203_2_output.txt', output)


class EulerPathTest(unittest.TestCase):

    def CheckOutput(self, filename, output, num=7):
        expect = ReadFile(filename)[0]
        for i in xrange(len(output)):
            try:
                j = output[i+1:].index(output[1])
                subset = output[i+1:j+2]
            except ValueError:
                continue
        self.assertIn('->'.join(subset), expect)

    def testData1(self):
        segments = ReadFile('euler_path_input.txt')
        output = EulerPath(segments)
        self.CheckOutput('euler_path_output.txt', output)

    def testData2(self):
        segments = ReadFile('eulerian_path_input.txt')
        output = EulerPath(segments)
        self.CheckOutput('eulerian_path_output.txt', output)

    def testData3(self):
        segments = ReadFile('dataset_203_5.txt')
        output = EulerPath(segments)
        self.CheckOutput('dataset_203_5_output.txt', output)


class EulerReconstructStringTest(unittest.TestCase):

    def ReadInput(self, filename):
        input = ReadFile(filename)
        return int(input[0]), input[1:]

    def testData1(self):
        k, kmers = self.ReadInput('reconstruct_string_input.txt')
        output = EulerReconstructString(k, kmers)
        self.assertEqual('GGCTTACCA', output)

    def testData2(self):
        k, kmers = self.ReadInput('StringReconstructionProblem_input.txt')
        output = EulerReconstructString(k, kmers)
        expect = ReadFile('StringReconstructionProblem_output.txt')[0]
        self.assertEqual(expect, output)

    def testData3(self):
        k, kmers = self.ReadInput('dataset_203_6.txt')
        output = EulerReconstructString(k, kmers)
        expect = ReadFile('dataset_203_6_output.txt')[0]
        self.assertEqual(expect, output)

    def testData4(self):
        k, kmers = self.ReadInput('reconstruct_string1_input.txt')
        output = EulerReconstructString(k, kmers)
        self.assertEqual('CAAATGCATCATACGCTCACCCAG', output)


class CircularStringTest(unittest.TestCase):

    def CheckOutput(self, expect, output):
        c_expect = collections.Counter(expect)
        c_output = collections.Counter(output)
        self.assertEqual(c_expect, c_output)

    def testData1(self):
        k = 4
        output = CircularString(k)
        expect = '0000110010111101'
        self.CheckOutput(expect, output)

    def testData2(self):
        k = 14
        output = CircularString(k)
        expect = ReadFile('uni_str_output.txt')[0]
        self.CheckOutput(expect, output)

    def testData3(self):
        k = 8
        output = CircularString(k)
        expect = '0000000010101011011001110111111011100101110001111000100100010001110100111001100000110001100110110100100101101010100101011101101110101000101001100100110101100101000011011110010000100111101100001011000101111010111110100011010000010000001110000111110011111111'
        self.CheckOutput(expect, output)


class ReconstructStringReadPairsTest(unittest.TestCase):

    def ReadInput(self, filename):
        input = ReadFile(filename)
        k, d = [int(i) for i in input[0].split()]
        return k, d, input[1:]

    def testData1(self):
        k, d, pairs = self.ReadInput('reconstruct_string_read_pair_input.txt')
        output = ReconstructStringReadPairs(k, d, pairs)
        expect = 'GTGGTCGTGAGATGTTGA'
        self.assertEqual(expect, output)

    def testData2(self):
        k, d, pairs = self.ReadInput('StringReconstructionFromReadPairs_input.txt')
        output = ReconstructStringReadPairs(k, d, pairs)
        expect = ReadFile('StringReconstructionFromReadPairs_output.txt')[0]
        self.assertEqual(expect, output)

    def testData3(self):
        k, d, pairs = self.ReadInput('dataset_204_14.txt')
        output = ReconstructStringReadPairs(k, d, pairs)
        expect = ReadFile('dataset_204_14_output.txt')[0]
        self.assertEqual(expect, output)

    def testData4(self):
        k, d, pairs = self.ReadInput('quiz2.txt')
        output = ReconstructStringReadPairs(k, d, pairs)
        expect = ''
        self.assertEqual(expect, output)


class ContigsTest(unittest.TestCase):

    def testData1(self):
        kmers = ReadFile('contigs_input.txt')
        output = Contigs(kmers)
        expect = ['AGA', 'ATG', 'ATG', 'CAT', 'GAT', 'TGGA', 'TGT']
        self.assertEqual(expect, output)

    def testData2(self):
        kmers = ReadFile('contig_generation_3_input.txt')
        output = Contigs(kmers)
        expect = ReadFile('contig_generation_3_output.txt')
        self.assertEqual(expect, output)

    def testData3(self):
        kmers = ReadFile('dataset_205_5.txt')
        output = Contigs(kmers)
        expect = ReadFile('dataset_205_5_output.txt')
        self.assertEqual(expect, output)


if __name__ == '__main__':
    unittest.main()
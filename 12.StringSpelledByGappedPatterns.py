'''
    StringSpelledByGappedPatterns(GappedPatterns, k, d)
        FirstPatterns ? the sequence of initial k-mers from GappedPatterns
        SecondPatterns ? the sequence of terminal k-mers from GappedPatterns
        PrefixString ? StringSpelledByPatterns(FirstPatterns, k)
        SuffixString ? StringSpelledByPatterns(SecondPatterns, k)
        for i = k + d + 1 to |PrefixString|
            if the i-th symbol in PrefixString does not equal the (i - k - d)-th symbol in SuffixString
                return "there is no string spelled by the gapped patterns"
        return PrefixString concatenated with the last k + d symbols of SuffixString


CODE CHALLENGE: Implement StringSpelledByGappedPatterns.

https://github.com/hkpcmit/BioInfo/blob/master/AssembleGenomes.py

Sample Input:

4 2
GACC|GCGC
ACCG|CGCC
CCGA|GCCG
CGAG|CCGG
GAGC|CGGA

Sample Output:

GACCGAGCGCCGGA

https://github.com/AnnaUfliand/Bioinformatics/blob/18c9c385c8f62c2895dc029ec5a07c416a6363fa/HW4/ConstructGappedString.py
'''

k = 50
d = 200
input = """
GACC|GCGC
ACCG|CGCC
CCGA|GCCG
CGAG|CCGG
GAGC|CGGA
""".split('\n')


pairs = [tuple(edge.split('|')) for edge in input if edge]
print(pairs)
import codecs


def recounstruct(pairs, k, d):
    out = pairs[0][0]
    for i in range(1, len(pairs)):
        out += pairs[i][0][-1]
    s = 0
    if len(out) - (2 * k + d) < 0:
        s = abs(len(out) - (2 * k + d))
    else:
        s = k
    j = max(len(out) - (2 * k + d), 0)

    out += pairs[j][1][:-s]
    for i in range(j + 1, len(pairs)):
        out += pairs[i][1][-1]
    return out

print (recounstruct(pairs, k, d))
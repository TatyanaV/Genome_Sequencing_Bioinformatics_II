'''
String Spelled by a Genome Path Problem. Reconstruct a string from its genome path.
     Input: A sequence of k-mers Pattern1, Patternn such that the last k - 1 symbols of Patterni are

                equal to the first k-1 symbols of Patterni+1 for 1 ? i ? n-1.
     Output: A string Text of length k+n-1 such that the i-th k-mer in Text is equal to Patterni  (for 1 ? i ? n)
      Sample Input:

ACCGA
CCGAA
CGAAG
GAAGC
AAGCT

Sample Output:

ACCGAAGCT
'''

dna = [s for s in"""
AAAT
AATG
ACCC
ACGC
ATAC
ATCA
ATGC
CAAA
CACC
CATA
CATC
CCAG
CCCA
CGCT
CTCA
GCAT
GCTC
TACG
TCAC
TCAT
TGCA
""".split() if s]

def reconstruct(kmers):
    genomePath = kmers[0]
    for i in range(1, len(kmers)):
        genomePath += kmers[i][-1]
    return genomePath

print(reconstruct(dna))
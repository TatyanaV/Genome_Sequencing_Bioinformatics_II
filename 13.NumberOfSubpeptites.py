'''
For now, we will assume for simplicity that the mass spectrometer breaks the copies of a cyclic peptide at every possible two bonds, so that the resulting experimental spectrum contains the masses of all possible linear fragments of the peptide, which are called subpeptides. For example, the cyclic peptide NQEL has 12 subpeptides: N, Q, E, L, NQ, QE, EL, LN, NQE, QEL, ELN, and LNQ. We will also assume that subpeptides may occur more than once if an amino acid occurs multiple times in the peptide (e.g., ELEL also has 12 subpeptides: E, L, E, L, EL, LE, EL, LE, ELE, LEL, ELE, and LEL.

EXERCISE BREAK: How many subpeptides does a cyclic peptide of length n have?

Sample Input:
     31315

Sample Output:
     980597910
'''
#https://quizlet.com/54070358/ch-2-bioinformatics-flash-cards/

def subpep(n):
	return n*(n-1)

print(subpep(13952))

'''
Generating Theoretical Spectrum Problem: Generate the theoretical spectrum of a cyclic peptide.
     Input: An amino acid string Peptide.
     Output: Cyclospectrum(Peptide).

CODE CHALLENGE: Solve the Generating Theoretical Spectrum Problem.

Sample Input:
     LEQN

Sample Output:
     0 113 114 128 129 227 242 242 257 355 356 370 371 484
'''

input = 'CTVTARVHMAKKE'
masses = {}
with open('integer_mass_table.txt') as f:
    for line in f:
        pair = line.strip().split()
        masses[pair[0]] = int(pair[1])
#print masses


def get_mass(peptide):
    return sum([masses[p] for p in peptide])


def subpeptide(peptide, pos, length):
    if pos+length <= len(peptide):
        return peptide[pos: pos+length]
    else:
        return peptide[pos:] + peptide[:length + pos - len(peptide)]

combos = [subpeptide(input, p, l) for p in range(len(input))
          for l in range(1, len(input))]

print (*([str(j) for j in
                sorted([get_mass(i) for i in combos] + [0, get_mass(input)])]), sep = ' ')
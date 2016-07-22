'''
    CYCLOPEPTIDESEQUENCING(Spectrum)
        Peptides ? a set containing only the empty peptide
        while Peptides is nonempty
            Peptides ? Expand(Peptides)
            for each peptide Peptide in Peptides
                if Mass(Peptide) = ParentMass(Spectrum)
                    if Cyclospectrum(Peptide) = Spectrum
                        output Peptide
                    remove Peptide from Peptides
                else if Peptide is not consistent with Spectrum
                    remove Peptide from Peptides

Sample Input:
     0 113 128 186 241 299 314 427

Sample Output:
     186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186
'''

input = '0 71 101 113 131 184 202 214 232 285 303 315 345 416'
with open('integer_mass_table.txt') as f:
    masses = list(set([int(line.strip().split()[1]) for line in f]))

spectrum = [int(i) for i in input.strip().split()]
#print spectrum


def expand_list(peptides, masses):
    if len(peptides) == 0:
        return [([m], [0, m]) for m in masses]

    def combine_spectrum(peptide, mass):
        def extend_spectrum(masses, m):
            return masses + [m] + [(sum(peptide[0][i:]) + m)
                                   for i in range(len(peptide[0]))]
        return (peptide[0]+[mass], extend_spectrum(peptide[1], mass))
    return [combine_spectrum(p, m) for p in peptides for m in masses]


def is_consistent(linear, spectrum):
    for e in linear:
        if linear.count(e) > spectrum.count(e):
            return False
    return True

consistent = []

for i in range(20):
    consistent = [cand for cand in expand_list(consistent, masses)
                  if is_consistent(cand[1], spectrum)]
    print ("i: %s, len(consistent): %s" % (i, len(consistent)))
    if len(consistent) == 0:
        break
    print (*(sorted(['-'.join([str(i) for i in c[0]])
                   for c in consistent])), sep = ' ')
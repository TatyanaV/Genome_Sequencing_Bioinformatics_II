'''
LinearSpectrum(Peptide, AminoAcid, AminoAcidMass)
    PrefixMass(0) ? 0
    for i ? 1 to |Peptide|
        for j ? 1 to 20
            if AminoAcid(j) =  i-th amino acid in Peptide
                PrefixMass(i) ? PrefixMass(i ? 1) + AminoAcidMass(j)
    LinearSpectrum ? a list consisting of the single integer 0
    for i ? 0 to |Peptide| ? 1
        for j ? i + 1 to |Peptide|
            add PrefixMass(j) ? PrefixMass(i) to LinearSpectrum
    return sorted list LinearSpectrum

CODE CHALLENGE: Implement LinearSpectrum.
     Input: An amino acid string Peptide.
     Output: The linear spectrum of Peptide.

 Sample Input:

NQEL

Sample Output:

0 113 114 128 129 242 242 257 370 371 484
https://github.com/AnnaUfliand/Bioinformatics/blob/1a38fc077eaef5cf176fecf97153ad7f78f3deab/HW5/TheoreticalSpectrumOfLinearPeptide.py
'''

masses = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115,
          'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}


def linearSpectrum(peptide):
    prefixMass = [0]
    for i in range(0, len(peptide) - 1):
        prefixMass.append(prefixMass[i] + masses[peptide[i]])
    prefixMass.append(prefixMass[-1] + masses[peptide[-1]])
    spectrum = [0]
    for i in range(len(peptide) + 1):
        for j in range(i + 1, len(peptide) + 1):
            spectrum.append(prefixMass[j] - prefixMass[i])
    return sorted(spectrum)


peptide = 'VAQ'
print (*linearSpectrum(peptide), sep = ' ')
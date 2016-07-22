'''

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


def countOccurrances(item, list):
    count = 0
    for i in list:
        if i == item:
            count += 1
    return count


def score(spectrum, masses):
    done = set([])
    count = 0
    for mass in spectrum:
        if mass not in done and mass in masses:
            occurrances = countOccurrances(mass, spectrum)
            occurrancesPeptide = countOccurrances(mass, masses)
            count += min(occurrances, occurrancesPeptide)
            for i in range(occurrances):
                done.add(mass)
    return count


def linearScore(peptide, spectrum):
    linear_spectrum = linearSpectrum(peptide)
    return score(spectrum, linear_spectrum)


peptide = 'NQEL'
spectrum = [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]
print (linearScore(peptide, spectrum))
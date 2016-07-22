import collections
import itertools
from ReverseComplement import ReverseComplement


B1_SPECTRUM = [0, 97, 99, 113, 114, 128, 128, 147, 147, 163, 186, 227, 241, 242, 244,
               260, 261, 262, 283, 291, 333, 340, 357, 388, 389, 390, 390, 405, 430, 430,
               447, 485, 487, 503, 504, 518, 543, 544, 552, 575, 577, 584, 631, 632, 650,
               651, 671, 672, 690, 691, 738, 745, 747, 770, 778, 779, 804, 818, 819, 835,
               837, 875, 892, 892, 917, 932, 932, 933, 934, 965, 982, 989, 1031, 1039, 1060,
               1061, 1062, 1078, 1080, 1081, 1095, 1136, 1159, 1175, 1175, 1194, 1194, 1208, 1209, 1223,
               1225, 1322]
B1_MASS_COUNTER = collections.Counter(B1_SPECTRUM)
CODE_MAP = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K',
            'AAU': 'N', 'ACA': 'T', 'ACC': 'T',
            'ACG': 'T', 'ACU': 'T', 'AGA': 'R',
            'AGC': 'S', 'AGG': 'R', 'AGU': 'S',
            'AUA': 'I', 'AUC': 'I', 'AUG': 'M',
            'AUU': 'I', 'CAA': 'Q', 'CAC': 'H',
            'CAG': 'Q', 'CAU': 'H', 'CCA': 'P',
            'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R',
            'CGU': 'R', 'CUA': 'L', 'CUC': 'L',
            'CUG': 'L', 'CUU': 'L', 'GAA': 'E',
            'GAC': 'D', 'GAG': 'E', 'GAU': 'D',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A',
            'GCU': 'A', 'GGA': 'G', 'GGC': 'G',
            'GGG': 'G', 'GGU': 'G', 'GUA': 'V',
            'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
            'UAC': 'Y', 'UAU': 'Y', 'UCA': 'S',
            'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
            'UGC': 'C', 'UGG': 'W', 'UGU': 'C',
            'UUA': 'L', 'UUC': 'F', 'UUG': 'L',
            'UUU': 'F'}
INT_MASS_MAP = {'G': 57, 'A': 71, 'S': 87,
                'P': 97, 'V': 99, 'T': 101,
                'C': 103, 'I': 113, 'L': 113,
                'N': 114, 'D': 115, 'K': 128,
                'Q': 128, 'E': 129, 'M': 131,
                'H': 137, 'F': 147, 'R': 156,
                'Y': 163, 'W': 186}
INT_MASS_SET = {m for m in INT_MASS_MAP.itervalues()}


def CheckConsistentSpectrum(peptide, spectrum_counter):
    peptide_spectrum = LinearSpectrumByMass(peptide)
    peptide_counter = collections.Counter(peptide_spectrum)
    for mass, freq in peptide_counter.iteritems():
        if mass not in spectrum_counter:
            return False
        if freq > spectrum_counter[mass]:
            return False
    return True


def Convolution(spectrum, mass_map=None):
    m_map = collections.defaultdict(set)
    sorted_spectrum = sorted(spectrum)
    output_counter = collections.Counter()
    for m1 in sorted_spectrum:
        if not m1:
            continue
        for m2 in sorted_spectrum:
            if m1 <= m2:
                continue
            mass_diff = m1 - m2
            output_counter[mass_diff] += 1
            m_map[m1].add(mass_diff)
    if mass_map is not None:
        mass_map.update(m_map)
    return list(itertools.chain.from_iterable(
            [m] * freq for m, freq in output_counter.iteritems()))


def ConvolutionCycloPeptideSequencing(M, N, spectrum):
    mass_map = {}
    conv = (m
            for m in Convolution(spectrum, mass_map=mass_map)
            if 57 <= m <= 200)
    conv_counter = collections.Counter(conv)
    candidates = sorted(
        ((freq, m) for m, freq in conv_counter.iteritems()),
        reverse=True)
    min_length = min(M, len(candidates))
    min_freq = candidates[min_length - 1][0]
    freq_mass = [m for freq, m in candidates if freq >= min_freq]
    return LeaderboardConvolutionCycloPeptideSequencing(
        sorted(spectrum), N, freq_mass)


def CycloPeptideScoreByMass(peptide_mass, spectrum):
    theoretic_spectrum = TheoreticalSpectrumByMass(peptide_mass)
    theoretic_counter = collections.Counter(theoretic_spectrum)
    observe_counter = collections.Counter(spectrum)
    return sum(min(freq, theoretic_counter[m])
               for m, freq in observe_counter.iteritems()
               if m in theoretic_counter)


def CycloPeptideScoring(peptide, spectrum):
    theoretic_spectrum = TheoreticalSpectrum(peptide)
    theoretic_counter = collections.Counter(theoretic_spectrum)
    observe_counter = collections.Counter(spectrum)
    return sum(min(freq, theoretic_counter[m])
               for m, freq in observe_counter.iteritems()
               if m in theoretic_counter)


def CycloPeptideSequencing(spectrum):
    result = []
    spectrum_counter = collections.Counter(spectrum)
    peptides = collections.deque([[]])
    while peptides:
        peptide = peptides.popleft()
        candidates = [peptide + [m] for m in INT_MASS_SET]
        for p in candidates:
            if sum(p) == spectrum[-1]:
                found = '-'.join(str(m) for m in p)
                result.append(found)
            elif CheckConsistentSpectrum(p, spectrum_counter):
                peptides.append(p)
    return result


def FindSubstringsPeptideEncoding(dna, peptide):
    unit_length = len(peptide) * 3
    result = []
    for i in xrange(len(dna) - unit_length + 1):
        substring = dna[i:i+unit_length]
        peptides = GetPeptides(substring)
        if peptide in peptides:
            result.append(substring)
    return result


def EncodeDna(dna):
    # Transcribe and translate this DNA.
    rna = ''.join(('U' if c == 'T' else c) for c in dna)
    return Translation(rna)
    return [Translation(TranscribeDna(d))
            for d in (dna, ReverseComplement(dna))]


def ExpandPeptide(peptide):
    return [peptide + c for c in INT_MASS_MAP]


def ExpandPeptideByAlphabet(peptide, alphabet_list):
    return [peptide + [m] for m in alphabet_list]


def GetPeptides(dna):
    return [EncodeDna(d) for d in (dna, ReverseComplement(dna))]


def LeaderboardCycloPeptideSequencing(spectrum, N):
    leader_peptide = ''
    leader_score = 1
    leaderboard = ['']
    spectrum_counter = collections.Counter(spectrum)
    while leaderboard:
        candidates = [p
                      for peptide in leaderboard
                      for p in ExpandPeptide(peptide)]
        leaderboard = []
        for p in candidates:
            peptide_mass = sum(INT_MASS_MAP[c] for c in p)
            if peptide_mass == spectrum[-1]:
                peptide_score = LinearScore(p, spectrum)
                if peptide_score > leader_score:
                    leader_peptide = p
                    leader_score = peptide_score
                leaderboard.append(p)
            elif peptide_mass < spectrum[-1]:
                leaderboard.append(p)
        if len(leaderboard) > N:
            leaderboard = Trim(leaderboard, spectrum, N)
    return '-'.join(str(INT_MASS_MAP[c]) for c in leader_peptide)


def LeaderboardConvolutionCycloPeptideSequencing(spectrum, N, alphabet_list):
    leader_peptide, leaders = [], []
    leader_score = 1
    leaderboard = [[]]
    spectrum_counter = collections.Counter(spectrum)
    while leaderboard:
        candidates = [
            p
            for peptide in leaderboard
            for p in ExpandPeptideByAlphabet(peptide, alphabet_list)]
        leaderboard = []
        for p in candidates:
            peptide_mass = sum(p)
            if peptide_mass == spectrum[-1]:
                peptide_score = CycloPeptideScoreByMass(p, spectrum)
                if peptide_score > leader_score:
                    leader_peptide = p
                    leader_score = peptide_score
                    leaders = [p]
                elif peptide_score == leader_score:
                    leaders.append(p)
                leaderboard.append(p)
            elif peptide_mass < spectrum[-1]:
                leaderboard.append(p)
        if len(leaderboard) > N:
            leaderboard = CyclicTrim(leaderboard, spectrum, N)
    return {'-'.join(str(m) for m in l)
            for l in leaders}


def LinearScore(peptide, spectrum):
    input_counter = collections.Counter(spectrum)
    linear_spectrum = LinearSpectrum(peptide)
    linear_counter = collections.Counter(linear_spectrum)
    return sum(min(freq, input_counter[m])
               for m, freq in linear_counter.iteritems() if m in input_counter)


def LinearSpectrum(peptide):
    prefix_mass = [0]
    for c in peptide:
        prefix_mass.append(prefix_mass[-1] + INT_MASS_MAP[c])
    spectrum = [0]
    peptide_length = len(peptide)
    for i in xrange(peptide_length):
        for j in xrange(i+1, peptide_length+1):
            spectrum.append(prefix_mass[j] - prefix_mass[i])
    return sorted(spectrum)


def LinearSpectrumByMass(peptide_mass):
    prefix_mass = [0]
    for m in peptide_mass:
        prefix_mass.append(prefix_mass[-1] + m)
    spectrum = [0]
    peptide_length = len(peptide_mass)
    for i in xrange(peptide_length):
        for j in xrange(i+1, peptide_length+1):
            spectrum.append(prefix_mass[j] - prefix_mass[i])
    return sorted(spectrum)


def TheoreticalSpectrum(peptide):
    # Add 0 & total mass.
    result = [0,
              sum(INT_MASS_MAP[c] for c in peptide)]
    peptide_length = len(peptide)
    for i in xrange(peptide_length):
        for length in xrange(1, peptide_length):
            max_i = i + length
            if max_i <= peptide_length:
                sub_peptide = peptide[i:max_i]
            else:
                sub_peptide = peptide[:max_i - peptide_length] + peptide[i:]
            result.append(sum(INT_MASS_MAP[c] for c in sub_peptide))
    return sorted(result)


def TheoreticalSpectrumByMass(peptide_mass):
    # Add 0 & total mass.
    result = [0] + peptide_mass
    peptide_length = len(peptide_mass)
    for i in xrange(peptide_length):
        for length in xrange(1, peptide_length):
            max_i = i + length
            if max_i <= peptide_length:
                sub_peptide = peptide_mass[i:max_i]
            else:
                sub_peptide = peptide_mass[:max_i - peptide_length] + peptide_mass[i:]
            result.append(sum(sub_peptide))
    return sorted(result)


def Translation(pattern):
    result = ''
    for i in xrange(0, len(pattern), 3):
        code = pattern[i:i+3]
        if CODE_MAP.get(code):
            result += CODE_MAP[code]
    return result


def Trim(leaderboard, spectrum, N):
    sorted_list = sorted(((LinearScore(p, spectrum), p)
                          for p in leaderboard),
                         reverse=True)
    min_score = sorted_list[N-1][0]
    return [tu[1] for tu in sorted_list if tu[0] >= min_score]


def CyclicTrim(leaderboard, spectrum, N):
    sorted_list = sorted(((CycloPeptideScoreByMass(p, spectrum), p)
                          for p in leaderboard),
                         reverse=True)
    min_score = sorted_list[N-1][0]
    return [tu[1] for tu in sorted_list if tu[0] >= min_score]
'''
Protein Translation Problem: Translate an RNA string into an amino acid string.
     Input: An RNA string Pattern and the array GeneticCode.
     Output: The translation of Pattern into an amino acid string Peptide.

CODE CHALLENGE: Solve the Protein Translation Problem.

'''
'''
with open('rna.txt', 'r') as myfile:
    rna=myfile.read().replace('\n', '')
'''
rna = 'CCCAGUACCGAGAUGAAU'


def rna_codon_map():
    codon_map = {}
    with open('RNA_codon_table_1.txt') as f:
        for line in f:
            pair = line.strip().split()
            if len(pair) == 2:
                codon_map[pair[0]] = pair[1]
            else:
                codon_map[pair[0]] = ''
    return codon_map


print (*([rna_codon_map()[rna[i:i+3]] for i in range(0, len(rna), 3)]), sep = '')
#print (*(sorted([text[i:i+k] for i in range(len(text) - k + 1)])), sep = '\n')

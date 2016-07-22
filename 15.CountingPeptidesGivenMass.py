'''
Counting Peptides with Given Mass Problem: Compute the number of peptides of given mass.
     Input: An integer m.
     Output: The number of linear peptides having integer mass m.

     Sample Input:
     1024

Sample Output:
     14712706211
'''
m = 1323 + 1
with open('integer_mass_table.txt') as f:
    masses = [int(line.strip().split()[1]) for line in f]
masses = [v for i, v in enumerate(masses)
          if i == 0 or i > 0 and masses[i-1] != v]

peptides = [0 for i in range(m)]
for i in range(m):
    for k in masses:
        if i-k > 0:
            peptides[i] += peptides[i-k]
        elif i-k == 0:
            peptides[i] += 1
print (peptides[-1])



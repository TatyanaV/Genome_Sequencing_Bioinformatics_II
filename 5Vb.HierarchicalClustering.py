'''

CODE CHALLENGE: Implement HierarchicalClustering.
     Input: An integer n, followed by an n x n distance matrix.
     Output: The result of applying HierarchicalClustering to this distance matrix (using Davg),
     with each newly created cluster listed on each line.

 Sample Input:

4
0.00 20.0 9.00 11.0
20.0 0.00 17.0 11.0
9.00 17.0 0.00 8.00
11.0 11.0 8.00 0.00

7
0.00 0.74 0.85 0.54 0.83 0.92 0.89
0.74 0.00 1.59 1.35 1.20 1.48 1.55
0.85 1.59 0.00 0.63 1.13 0.69 0.73
0.54 1.35 0.63 0.00 0.66 0.43 0.88
0.83 1.20 1.13 0.66 0.00 0.72 0.55
0.92 1.48 0.69 0.43 0.72 0.00 0.80
0.89 1.55 0.73 0.88 0.55 0.80 0.00
https://github.com/chrisKingsley/Coursera_Bioinformatics_Algorithms2/blob/eb4c2e72e6b8e014baf5a213385edda089baf727/3_clustering/yeast.py

 Sample Output:

4 6
5 7
3 4 6
1 2
5 7 3 4 6
1 2 5 7 3 4 6

https://github.com/chrisKingsley/Coursera_Bioinformatics_Algorithms2/blob/eb4c2e72e6b8e014baf5a213385edda089baf727/3_clustering/yeast.py
https://github.com/johnmerm/bioinfo/blob/e6a255a461ca83c58e56c83d4717bc13f8880c5d/src/main/java/bioinfo/yeast/hier_clustering.py
'''


from timeit import itertools

def dist(clust1,clust2,Data):
    return sum([ sum([ Data[i][j] for j in clust2]) for i in clust1])/(len(clust1)*len(clust2))

def hier_clustering(Data,n):
    gen = n
    
    Clusters = {i:[i] for i in range(n)}
    all_Clusters = {i:[i] for i in range(n)}
    
    T = {i:[] for i in range(n)}
    while len(Clusters)>1:
        permuts = itertools.permutations(Clusters.keys(),2)
        closest = min(permuts,key = lambda x:dist(Clusters[x[0]],Clusters[x[1]],Data))
        Clusters[gen] = Clusters[closest[0]]+Clusters[closest[1]]
        T[gen] = [closest[0],closest[1]]
        all_Clusters[gen] = Clusters[closest[0]]+Clusters[closest[1]]
        print (' '.join([str(c+1) for c in Clusters[closest[0]]])+" "+' '.join([str(c+1) for c in Clusters[closest[1]]]))
        del Clusters[closest[0]]
        del Clusters[closest[1]]
        gen = gen+1
    return T,all_Clusters
    
def test_hier_clustering():
    n=7
    Data = [[0.00, 0.74, 0.85, 0.54, 0.83, 0.92, 0.89],
            [0.74, 0.00, 1.59, 1.35, 1.20, 1.48, 1.55],
            [0.85, 1.59, 0.00, 0.63, 1.13, 0.69, 0.73],
            [0.54, 1.35, 0.63, 0.00, 0.66, 0.43, 0.88],
            [0.83, 1.20, 1.13, 0.66, 0.00, 0.72, 0.55],
            [0.92, 1.48, 0.69, 0.43, 0.72, 0.00, 0.80],
            [0.89, 1.55, 0.73, 0.88, 0.55, 0.80, 0.00]]
    
    hier_clustering(Data, n)

def exam_hier_clustering():
    al = list(open('distMat.txt'))
    n = int( al[0].strip())
    Data = [ [float(x) for x in a.strip().split(' ')] for a in al[1:] ]
    
    hier_clustering(Data, n)
    
if __name__ == '__main__':
    exam_hier_clustering()
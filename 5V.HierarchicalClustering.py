'''


CODE CHALLENGE: Implement HierarchicalClustering.
     Input: An integer n, followed by an n x n distance matrix.
     Output: The result of applying HierarchicalClustering to this distance matrix (using Davg),
     with each newly created cluster listed on each line.

 Sample Input:

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


'''


#!/usr/bin/env python

import math, re, sys

# reads a matrix of points.  Reads centers if the
# passed distortion flag is true
def readMatrixFile(fileName, readNumClusters=True,
                   distortion=False, betaParam=False):
    matrix, centers, readCenters = [], [], True

    infile = open(fileName, 'r')
    if readNumClusters:
        k, m = [ int(x) for x in infile.readline().rstrip().split() ]
    else:
        m = int(infile.readline().rstrip())
    if betaParam:
        beta = float(infile.readline().rstrip())

    for line in infile:
        if distortion and re.match('-+', line):
            readCenters = False
            continue

        vals = [ float(x) for x in line.rstrip().split() ]
        if readCenters and distortion:
            centers.append(vals)
        else:
            matrix.append(vals)
    infile.close()

    if distortion:
        return  k, m, matrix, centers
    elif betaParam:
        return k, m, beta, matrix
    elif readNumClusters:
        return  k, m, matrix
    else:
        return m, matrix

# returns hash key for two numerical indices
def getKey(i, j):
    return '%s %s' % (i, j)


# returns the index of the two closest clusters based on the distance matrix
def closestClusters(clusters, dists, m):
    c1, c2, minDist = 0, 0, sys.maxsize

    for i in clusters:
        for j in clusters:
            if i!=j and dists[ getKey(i,j) ] < minDist:
                c1, c2, minDist = i, j, dists[ getKey(i,j) ]

    return c1, c2


# prints the members of the newly formed cluster
def printNewCluster(clusters, c1, c2):
    clust1 = ' '.join(str(x+1) for x in clusters[c1] )
    clust2 = ' '.join(str(x+1) for x in clusters[c2] )
    print ('%s %s' % (clust1, clust2))

# calculates cluster to cluster distance based on average distance
# of the members of each cluster to each other
def clusterAveDist(clusters, dists, k, c1):
    n, dist = 0.0, 0.0
    for i in clusters[k]:
        for j in clusters[c1]:
            dist += dists[ getKey(i, j) ]
            n += 1

    return dist/n

# merges the two passed clusters (indexed by c1/c2) into one cluster
# deletes one of the old clusters and updates the distance matrix
def joinClusters(clusters, dists, c1, c2):
    clusters[c1] = clusters[c1] + clusters[c2]
    del clusters[c2]

    for k in clusters:
        dist = clusterAveDist(clusters, dists, k, c1)
        dists[ getKey(k, c1) ] = dist
        dists[ getKey(c1, k) ] = dist


# performs heirarchical clustering on the passed distance matrix of
# dimension m.  Prints cluster members for each new cluster generated
# in each round of joining.
def hierarchicalCluster(m, distMat):
    dists = { getKey(i,j): distMat[i][j]
               for i in range(m) for j in range(m) }
    clusters = { i:[i] for i in range(m) }

    while len(clusters) > 1:
        c1, c2 = closestClusters(clusters, dists, m)
        printNewCluster(clusters, c1, c2)
        joinClusters(clusters, dists, c1, c2)

m, distMat = readMatrixFile('distMat.txt', readNumClusters=False)
hierarchicalCluster(m, distMat)

# prints the members of the newly formed cluster
def printNewCluster(clusters, c1, c2):
    clust1 = ' '.join(str(x+1) for x in clusters[c1] )
    clust2 = ' '.join(str(x+1) for x in clusters[c2] )
    print ('%s %s' % (clust1, clust2))


# returns the index of the two closest clusters based on the distance matrix
def closestClusters(clusters, dists, m):
    c1, c2, minDist = 0, 0, sys.maxint

    for i in clusters:
        for j in clusters:
            if i!=j and dists[ getKey(i,j) ] < minDist:
                c1, c2, minDist = i, j, dists[ getKey(i,j) ]

    return c1, c2

# merges the two passed clusters (indexed by c1/c2) into one cluster
# deletes one of the old clusters and updates the distance matrix
def joinClusters(clusters, dists, c1, c2):
    clusters[c1] = clusters[c1] + clusters[c2]
    del clusters[c2]

    for k in clusters:
        dist = clusterAveDist(clusters, dists, k, c1)
        dists[ getKey(k, c1) ] = dist
        dists[ getKey(c1, k) ] = dist


# calculates cluster to cluster distance based on average distance
# of the members of each cluster to each other
def clusterAveDist(clusters, dists, k, c1):
    n, dist = 0.0, 0.0
    for i in clusters[k]:
        for j in clusters[c1]:
            dist += dists[ getKey(i, j) ]
            n += 1

    return dist/n
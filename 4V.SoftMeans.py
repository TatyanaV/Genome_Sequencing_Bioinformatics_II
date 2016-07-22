#!/usr/bin/env python
'''


CODE CHALLENGE: Implement the expectation maximization algorithm for soft k-means clustering.
     Input: Integers k and m, followed by a stiffness parameter ?, followed by a set of points
     Data in m-dimensional space.
     Output: A set Centers consisting of k points (centers) resulting from applying the
     expectation maximization algorithm for soft k-means clustering. Select the first k points
     from Data as the first centers for the algorithm and run the algorithm for 100 E-steps
     and 100 M-steps. Results should be accurate up to three decimal places.
      Sample Input:

2 2
2.7
1.3 1.1
1.3 0.2
0.6 2.8
3.0 3.2
1.2 0.7
1.4 1.6
1.2 1.0
1.2 1.1
0.6 1.5
1.8 2.6
1.2 1.3
1.2 1.0
0.0 1.9

Sample Output:

1.662 2.623
1.075 1.148
https://github.com/chrisKingsley/Coursera_Bioinformatics_Algorithms2/blob/eb4c2e72e6b8e014baf5a213385edda089baf727/3_clustering/yeast.py

https://github.com/chrisKingsley/Coursera_Bioinformatics_Algorithms2/blob/eb4c2e72e6b8e014baf5a213385edda089baf727/3_clustering/yeast.py
'''
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
    
    
# returns the Euclidean distance between two lists
def euclideanDist(list1, list2):
    dist = 0.0
    for i in range(len(list1)):
        dist += (list1[i]-list2[i])*(list1[i]-list2[i])
    
    return math.sqrt(dist)

    
# returns the dot product between two lists
def dotProd(list1, list2):
    prod = 0.0
    for i in range(len(list1)):
        prod += (list1[i]*list2[i])
    
    return prod
    

# returns the largest distance between the data points and their
# closest center
def maxDistance(data, centers):
    dists = []
    
    for i in range(len(data)):
        minDist = sys.maxsize
        for j in range(len(centers)):
            dist = euclideanDist(data[i], centers[j])
            if dist < minDist:
                minDist = dist
        dists.append(minDist)
        
    return max(dists)
               
    
# for the given set of points and value of k, returns the
# center points that have the largest distance to the
# closest center
def farthestFirstTraversal(k, matrix):
    centers = [ matrix[0] ]
    
    while len(centers) < k:
        maxDist, maxIdx = -sys.maxint, 0
        for i in range(1, len(matrix)):
            dists = []
            for j in range(len(centers)):
                dists.append(euclideanDist(matrix[i], centers[j]))
            dist = min(dists)
            if dist > maxDist:
                maxDist, maxIdx = dist, i
        centers.append( matrix[maxIdx] )
    
    return centers

# k, m, data = readMatrixFile('geMatrix.txt')
# centers = farthestFirstTraversal(k, data)
# for x in centers:
    # print ' '.join([ str(y) for y in x ])
    

# for the given set of data of points and centers, returns the 
# squared error distortion (mean squared minimum distance between
# each point and all centers)
def squaredErrorDistortion(data, centers):
    distortionDist = 0.0
    
    for i in range(len(data)):
        closestDist = sys.maxsize
        for j in range(len(centers)):
            dist = euclideanDist(data[i],centers[j])
            if dist < closestDist:
                closestDist = dist
        distortionDist += closestDist*closestDist
        
    return distortionDist/len(data)

# k, m, data, centers = readMatrixFile('distortionMat.txt', distortion=True)
# print squaredErrorDistortion(data, centers)


# for the passed data points and centers, assigns each point
# to its closest center.  Returns a dict where center number
# points to the indices of the points in data that are in 
# the cluster with that center
def assignClusters(data, centers, k):
    clusters = { x:[] for x in range(k) }
    
    for i in range(len(data)):
        closestDist, centerIdx = sys.maxint, 0
        for j in range(k):
            dist = euclideanDist(data[i],centers[j])
            if dist < closestDist:
                closestDist, centerIdx  = dist, j
                
        clusters[ centerIdx ].append( i )
    
    return clusters

    
# for the given set of points and their dimension m, returns
# the center of gravity (mean of the sum of each component)
def centerOfGravity(points, m):
    center = [0.0]*m
    
    for i in range(len(points)):
        for j in range(m):
            center[j] += float(points[i][j])/len(points)
        
    return center


# for the given set of data, clusters, k (# of clusters), and
# m (dimension of points), determines the center of each cluster
# and returns in a 2D array
def getCenters(data, clusters, k, m):
    centers = []
    
    for i in clusters:
        points = []
        for j in clusters[ i ]:
            points.append( data[j] )
        center = centerOfGravity(points, m)
        centers.append( center )
            
    return centers
    
    
# clusters the set of m dimensional points in data
# into k clusters using Lloyd's algorithm
def kMeansClustering(data, k, m):
    centers = []
    for i in range(k):
        centers.append( data[i] )
        
    while True:
        clusters = assignClusters(data, centers, k)
        newCenters = getCenters(data, clusters, k, m)
        if newCenters==centers:
            break
        centers = newCenters
        
    return centers
    
# k, m, data = readMatrixFile('kMeansMatrix.txt')
# centers = kMeansClustering(data, k, m)
# for center in centers:
    # print ' '.join([ '%0.3f' % x for x in center ])

    
# expectation step in EM algorithm for soft k-means clustering.
# estimates new responsibility matrix based on cluster centers
def eStep(data, beta, centers):
    resp = [ [0]*len(data) for x in range(len(centers)) ]
    
    # calculate responsibility of each cluster for each data point
    for i in range(len(data)):
        sum = 0.0
        for j in range(len(centers)):
            dist = euclideanDist(centers[j], data[i])
            resp[j][i] = math.exp(-beta * dist)
            sum += resp[j][i]
        for j in range(len(centers)):
            resp[j][i] /= sum
            
    return resp
    
    
# maximum likelihood step for EM algorithm for soft k-means clustering.
# estimates new cluster centers based on responsibility matrix
def mStep(data, respMatrix, m):
    centers = []
    
    # calculate new centers based on data and responsibility matrix
    for i in range(len(respMatrix)):
        center = []
        for j in range(m):
            dotProd, denom = 0.0, 0.0
            for k in range(len(data)):
                dotProd += data[k][j]*respMatrix[i][k]
                denom += respMatrix[i][k]
            center.append(dotProd/denom)
        centers.append(center)
        
    return centers

    
# soft k-means algorithm that uses EM
def soft_kMeans(data, k, m, beta, numIter=100):
    centers = [ data[x] for x in range(k) ]
    
    for i in range(numIter):
        # print centers
        respMatrix = eStep(data, beta, centers)
        centers = mStep(data, respMatrix, m)
        # print respMatrix
        
    return centers
'''
k, m, beta, data = readMatrixFile('soft_kMeansMatrix.txt', betaParam=True)
centers = soft_kMeans(data, k, m, beta, numIter=100)
for center in centers:
    print (' '.join([ str(x) for x in center ]))
'''
    
# returns hash key for two numerical indices
def getKey(i, j):
    return '%s %s' % (i, j)
    
    
# returns the index of the two closest clusters based on the distance matrix
def closestClusters(clusters, dists, m):
    c1, c2, minDist = 0, 0, sys.maxint
    
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

# m, distMat = readMatrixFile('distMat.txt', readNumClusters=False)
# hierarchicalCluster(m, distMat)


# Quiz questions
data = [[2,6], [4,9], [5,7], [6,5], [8,3]]
centers = [[4,5], [7,4]]
print (maxDistance(data, centers))

data = [[2,8], [2,5], [6,9], [7,5], [5,2]]
centers = [[3,5], [5,4]]
print (squaredErrorDistortion(data, centers))

data = [[1, 3, -1], [9, 8, 14], [6, 2, 10], [4, 3, 1]]
print (centerOfGravity(data, 3))


# version of the expectation step for clustering, but using
# Newtonian distance
def eStep_newton(data, centers):
    resp = [ [0]*len(data) for x in range(len(centers)) ]
    
    # calculate responsibility of each cluster for each data point
    for i in range(len(data)):
        sum = 0.0
        for j in range(len(centers)):
            dist = euclideanDist(centers[j], data[i])
            resp[j][i] = 1/(dist*dist)
            sum += resp[j][i]
        for j in range(len(centers)):
            resp[j][i] /= sum
            
    return resp

'''
Data: (2,8), (2,5), (6,9), (7,5), (5,2)
Centers: (3,5), (5,4)
'''
print ("quiz question")
data = [[2,8], [2,5], [6,9], [7,5], [5,2]]
centers = [[3,5], [5,4]]
resp = eStep(data, 1, centers)
print (resp)

print("HERE")
data = [[2,6], [4,9], [5,7], [6,5], [8,3]]
hiddenMat = [[0.6,0.1,0.8,0.5,0.7],[0.4,0.9,0.2,0.5,0.3]]
center = mStep(data, hiddenMat, 2)
print (center)
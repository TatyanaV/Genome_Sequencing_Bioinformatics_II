'''
CODE CHALLENGE: Implement the Lloyd algorithm for k-means clustering.
     Input: Integers k and m followed by a set of points Data in m-dimensional space.
     Output: A set Centers consisting of k points (centers) resulting from applying the
     Lloyd algorithm to Data and Centers, where the first k points from Data are selected
     as the first k centers.

 Sample Input:

2 2
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

1.800 2.867
1.060 1.140
'''
def pointDist(a,b):
    return sum([(float(a[i])-float(b[i]))**2 for i in range(len(a))])

def centerOffGravity(points):
    m = len(next(iter(points)))
    return tuple([sum([float(d[i]) for d in points])/len(points) for i in range(m)])

def lloyd(data,k):
    centers = {c for c in data[:k]}
    moves,dists = lloydStep(data, centers)
    while(sum(dists.values())>0.000001):
        centers = moves.values()
        moves,dists = lloydStep(data, centers)
    return moves.values()

def lloydStep(Data,centers):
    clusters = {}
    for dataPoint in Data:
        cluster = min(centers,key = lambda center:pointDist(dataPoint, center))
        if cluster in clusters:
            clusters[cluster].add(dataPoint)
        else:
            clusters[cluster] = {dataPoint}

    moves = {}
    dists = {}
    for c,dc in clusters.items():
        cn = centerOffGravity(dc)
        moves[c] = cn
        dists[c] = pointDist(c, cn)
    return moves,dists
'''
def exam_lloyd():
    al = list(open('dataset_10928_3.txt'))
    k,m = [int(x) for x in al[0].strip().split(' ')]
    Data = [ tuple([float(x) for x in a.strip().split(' ')]) for a in al[1:] ]

    centers = lloyd(Data,k)
    for c in centers:
        print ' '.join([str(e) for e in c])
'''

k = 2
data ="""
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
""".split('\n')

print(data)

Data = [tuple(edge.split(' ')) for edge in data if edge]
#Data = [ tuple([float(x) for x in a.strip().split(' ')]) for a in data ]

print(Data)

centers = lloyd(Data, k)
for c in centers:
    print (' '.join([str(e) for e in c]))
'''
CODE CHALLENGE: Implement the FarthestFirstTraversal clustering heuristic.
     Input: Integers k and m followed by a set of points Data in m-dimensional space.
     Output: A set Centers consisting of k points (centers) resulting from applying
     FarthestFirstTraversal(Data, k), where the first point from Data is chosen as the
     first center to initialize the algorithm

 FarthestFirstTraversal(Data, k)
?Centers ? the set consisting of a single randomly chosen point from Data
??while |Centers| < k
???DataPoint ? the point in Data maximizing d(DataPoint, Centers)
???add DataPoint to Centers
?return Centers
 Sample Input:

3 2
0.0 0.0
5.0 5.0
0.0 5.0
1.0 1.0
2.0 2.0
3.0 3.0
1.0 2.0

Sample Output:

0.0 0.0
5.0 5.0
0.0 5.0

https://github.com/johnmerm/bioinfo/blob/e6a255a461ca83c58e56c83d4717bc13f8880c5d/src/main/java/bioinfo/yeast/kmeans.py
https://github.com/ilap/Bioinformatics_Part2/blob/8c4fcec1df79e45b80818a281683578876060a6b/CH_3_How_Did_Yeast_Become_a_Wine-Maker/7_3_IntroToClustering.py
'''

k = 5
data = """
0.0 0.0
5.0 5.0
0.0 5.0
1.0 1.0
2.0 2.0
3.0 3.0
1.0 2.0
""".split('\n')

print(data)

Data = [tuple(edge.split(' ')) for edge in data if edge]

print(Data)

def pointDist(a,b):
    return sum([(float(a[i])-float(b[i]))**2 for i in range(len(a))])

def dist(a,centers):
    d = {b:pointDist(a,b) for b in centers}
    dm = min(d.items(),key = lambda x:x[1])
    return dm[1]

def farthestFirstTraversal(Data, k):
    dataPoint =  Data[0]
    centers = {dataPoint}
    while len(centers)<k:
        dataPoint = max(Data,key = lambda x:dist(x,centers))
        centers.add(dataPoint)
    return centers


centers = farthestFirstTraversal(Data, k)
print (centers)

for arr in centers:
    print (' '.join ([str ( val) for val in arr]))

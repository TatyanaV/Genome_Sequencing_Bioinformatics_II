'''
Squared Error Distortion Problem:? Compute the squared error distortion of a set of data points with respect to a set of centers.
? Input: A set of points Data and a set of centers Centers.?
? Output: The squared error distortion Distortion(Data, Centers).

CODE CHALLENGE: Solve the Squared Error Distortion Problem.
     Input: Integers k and m, followed by a set of centers Centers and a set of points Data.
     Output: The squared error distortion Distortion(Data, Centers).

 Sample Input:

2 2
2.31 4.55
5.96 9.08
--------
3.42 6.03
6.23 8.25
4.76 1.64
4.47 4.33
3.95 7.61
8.93 2.97
9.74 4.03
1.73 1.28
9.72 5.01
7.27 3.77

Sample Output:

18.246
'''
k = 2
center ="""
2.31 4.55
5.96 9.08
""".split('\n')

print(center)

centers = {tuple(edge.split(' ')) for edge in center if edge}

print(centers)

data = """
3.42 6.03
6.23 8.25
4.76 1.64
4.47 4.33
3.95 7.61
8.93 2.97
9.74 4.03
1.73 1.28
9.72 5.01
7.27 3.77
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

def distortion(Data,k,centers):
    n = len(Data)
    dm =[dist(a,centers) for a in Data]
    return sum(dm)/n

print (distortion(Data, k, centers))



import numpy as np
import pandas as pd
import argparse
import myspkmeans

parser = argparse.ArgumentParser(description='Calculate K_meanks of DataPoints')
parser.add_argument('FirstParam', type=int, help='number of Clusters')
parser.add_argument('SecondParam', type=str, help='goal')
parser.add_argument('filename', type=str, help='DataPoints1')
args = parser.parse_args()



def k_means_pp(points,clusterss,N,K):

    if K==0 or K==1:
        print ("An Error Has Occurred")
        exit()

    size = N
    indexes = np.array([i for i in range(size)])
    arr = np.array([0 for i in range(K)])
    clussind = np.zeros(K)
    np.random.seed(0)
    firststrandom = np.random.choice(indexes)
    arr[0] = firststrandom
    clussind[0] = firststrandom
    clusterss[0] = points[firststrandom]


    j = 1
    while  K > j:

        d = np.zeros(size)
        for i in range(size):
            d[i] = nearest(clusterss,points[i],j)


        sum = np.sum(d)
        
            
        probarr = np.zeros(size)
        for i in range(size):

            probarr[i] = d[i]/sum

        othertstrandom = np.random.choice(indexes,p=probarr)
        clusterss[j] = points[othertstrandom]

        arr[j] = othertstrandom
        clussind[j] = othertstrandom

        j = j +1

    for i in range(K-1):
        print(str(int(clussind[i]))+",", end="")
    print(int(clussind[-1]))
    
   
def nearest(clus, p12, r):
    min = calc(p12, clus[0])

    for i in range(r):

        c1 = calc(p12, clus[i])
        if  min >= c1:

            min = c1

    return min


def calc(x, y):
    c2 = 0.0
    for i in range(len(x)):
        c2 = c2 +  ((x[i] - y[i])*(x[i] - y[i]))
    return c2



K = args.FirstParam
goal = args.SecondParam

df = pd.read_csv(args.filename, header=None)
N = len(df)
d = len(df.columns)

DataPoints = df.to_numpy()
DataPointslist = DataPoints.tolist()

DataPointslist = myspkmeans.fit1(K, DataPointslist, goal)

if goal == "spk":
    
    DataPoints = np.asarray(DataPointslist)
    K = int(DataPoints.size / N)

    if K >= N:
        print("Invalid Input!")
        exit()
    
    indexes = np.zeros(K, dtype=int)
    
    clusterss = np.zeros((K, len(DataPoints[0])))
    k_means_pp(DataPoints,  clusterss , N, K)

    DataPointslist = DataPoints.tolist()
    centroidslist = clusterss.tolist()
    myspkmeans.fit2(DataPointslist, centroidslist)

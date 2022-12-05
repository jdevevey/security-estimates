import numpy as np
from math import sqrt

#Computes median value for S

n=256

#############################
######## Reference parameters
#############################

#These values are given in the documentation of GDilithium and wan be used to compare with this simulation
# We take k+l here as we consider only rot(s1) now and not (rot(s1)^top rot(s2)^top)^top

#Weak
#eta = 7
#l = 4
# Median is approximately 230

#Medium
#eta = 6
#l = 6
# Median is approximately 225

#High
#eta = 5
#l = 8
# Median is approximately 210

#Very High
#eta = 3
#l = 10
# Median is approximately 145

#############################
######### Parameters ########
#############################

eta = 1
k = 4
l = 5
tau = 49

def S_estimate(N,n,l,eta):
    """
    Runs N computation of secret key and returns the median value of the largest singular value of (rot(s1) rot(s2)).
    """
    res = []
    
    rots1 = [[0 for i in range(n)] for j in range(n*l)]
    
    for loop in range(N):
        s1 = np.random.randint(-eta,eta+1,size=n*l)
        #Build the rotational of s1, which is such that rot(s1)*vec(c) = vec(s1*c)
        for k in range(l):
            for i in range(n):
                for j in range(n):
                    rots1[i+k*n][j] = (-1)**((i+j)//n)*s1[((i+j)%n)+k*n]
        u,s,vh = np.linalg.svd(rots1)
        res.append(max(s))
    
    print("First quartile: "+str(np.nanquantile(res, 0.25)))
    print("Median: "+str(np.nanquantile(res, 0.5)))
    print("Average: "+str(np.average(res)))
    print("sigma = "+str(6.85*sqrt(tau)*np.nanquantile(res,0.5)))
    print("old sigma = "+str(11*sqrt(tau)*np.nanquantile(res,0.5)))

N = 10
S_estimate(N,n,l+k,eta)

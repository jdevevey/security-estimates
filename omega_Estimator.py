from numpy.polynomial import Polynomial
import numpy as np
from math import sqrt

n = 256
q = 8380417

############################
####### REFERENCE PARAMETERS
############################

distribution = "uniform"

#Medium
k = 2
l = 3
eta = 10
d=8
S=295
tau = 39
sigma = 5.7*S*sqrt(tau)
gamma2 = (q-1)//2112

#Recommended
#k = 5
#l = 5
#eta = 4
#d = 13
#S = 147
#tau = 49
#sigma = 5.7*S*sqrt(tau)
#gamma2 = (q-1)//32

#Very High
#k = 8
#l = 7
#eta = 2
#d = 13
#S = 90
#tau = 60
#sigma = 5.7*S*sqrt(tau)
#gamma2 = (q-1)//32

N = 100
res = []

def SampleInBall(n,tau):
    c = [0 for i in range(n)]
    for i in range(n-tau,n):
        j = np.random.randint(0,i+1)
        s = np.random.randint(0,2)
        c[i] = c[j]
        c[j] = (-1)**s
    return Polynomial(c)

def modq(r):
    return np.array([Polynomial([(r[i].coef[j])%q for j in range(r[i].degree())]) for i in range(k)])

def centeredmod(r,v):
    return np.array([Polynomial([int(r[i].coef[j])%v if int(r[i].coef[j])%v<=(v//2) else (int(r[i].coef[j])%v)-v for j in range(r[i].degree())]) for i in range(k)])

def Power2Round(r,d):
    v = 2**(d)
    r = np.array([Polynomial([c%q for c in r[i].coef]) for i in range(k)])
    r0 = np.array([Polynomial([int(c)%v if int(c)%v<=(v//2) else (int(c)%v)-v for c in r[i].coef]) for i in range(k)])
    temp = r-r0
    return (np.array([Polynomial([int(c)//v for c in temp[i].coef]) for i in range(k)]),r0)

def MakeHint(z,r,alpha):
    def Highbits(truc, gamma):
        truc = modq(truc)
        r0 = centeredmod(truc,gamma)
        return np.array([Polynomial([int(truc[i].coef[j]-r0[i].coef[j])//gamma if int(truc[i].coef[j]-r0[i].coef[j])!=q-1 else 0 for j in range(n)]) for i in range(k)])
    h1 = Highbits(r-z, alpha)
    h2 = Highbits(r,alpha)
    return [[0 if int(h1[i].coef[j])==int(h2[i].coef[j]) else 1 for j in range(n)] for i in range(k)]

for loop in range(N):
    #Sample signing key
    A =  np.array([[Polynomial([np.random.randint(0,q) for ell in range(n)]) for j in range(l)] for i in range(k)])
    s1 = np.array([Polynomial([np.random.randint(-eta,eta+1) for i in range(n)]) for j in range(l)])
    s2 = np.array([Polynomial([np.random.randint(-eta,eta+1) for i in range(n)]) for j in range(k)])
    t = np.matmul(A,s1)+s2
    (t1,t0) = Power2Round(t,d)
    #Compute a (simulated) signature
    if distribution == "uniform":
        y = np.array([Polynomial([np.random.randint(-gamma1+beta,gamma1-beta+1) for i in range(n)]) for j in range(l)])
    w = np.matmul(A,y)
    c = SampleInBall(n,tau)
    z = np.array([c*p for p in t0])
    h = MakeHint(-z,w-np.array([c*p for p in s2]),2*gamma2)
    res.append(sum([sum(i) for i in h]))

print("Estimated value for omega (0.99 quantile): "+str(np.nanquantile(res, 0.99)))
print("Higher value for omega (max): "+str(max(res)))

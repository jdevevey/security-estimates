from numpy.polynomial import Polynomial
import numpy as np
from math import sqrt
from S_Estimator import S_estimate

n = 256
q = 8380417

############################
####### REFERENCE PARAMETERS
############################

#Intuition: d has the most impact. The smaller d, the smaller the lower bound on S.
#tau can be used to fine tune -> the smaller tau, the smaller the lower bound on S.
#eta has almost no impact on gamma2.

#Medium
k = 2
l = 2
eta = 7
tau = 60
d = 11
current_gamma = (q-1)//(512)

#Recommended
#k = 5
#l = 5
#eta = 4
#d = 13
#tau = 49

#Very High
#k = 8
#l = 7
#eta = 2
#d = 13
#tau = 60

N = 100
res = []

print("Lower bound on S with current gamma2: "+str(sqrt(k/(l*tau))*(2*current_gamma+1+2**(d-1)*tau)/5.985))
print("Current Acceptance proba: "+str(np.exp(-256*k*tau*eta/current_gamma)))

def SampleInBall(n,tau):
    c = [0 for i in range(n)]
    for i in range(n-tau,n):
        j = np.random.randint(0,i+1)
        s = np.random.randint(0,2)
        c[i] = c[j]
        c[j] = (-1)**s
    return Polynomial(c)

def Power2Round(r,d):
    v = 2**(d)
    r = np.array([Polynomial([c%q for c in r[i].coef]) for i in range(k)])
    r0 = np.array([Polynomial([int(c)%v if int(c)%v<=(v//2) else (int(c)%v)-v for c in r[i].coef]) for i in range(k)])
    temp = r-r0
    return (np.array([Polynomial([int(c)//v for c in temp[i].coef]) for i in range(k)]),r0)

for loop in range(N):
    #Sample signing key
    A =  np.array([[Polynomial([np.random.randint(0,q) for ell in range(n)]) for j in range(l)] for i in range(k)])
    s1 = np.array([Polynomial([np.random.randint(-eta,eta+1) for i in range(n)]) for j in range(l)])
    s2 = np.array([Polynomial([np.random.randint(-eta,eta+1) for i in range(n)]) for j in range(k)])
    t = np.matmul(A,s1)+s2
    (t1,t0) = Power2Round(t,d)
    #Compute a c*t_0
    c = SampleInBall(n,tau)
    z = np.array([c*p for p in t0])
    res.append(max([max([abs(i) for i in j.coef])for j in z]))

print("Estimated lower bound for gamma2 (0.99 quantile): "+str(np.nanquantile(res, 0.99)))
print("Implying lower bound on S: "+str(sqrt(k/(l*tau))*(2*np.nanquantile(res,0.99)+1+2**(d-1)*tau)/5.985))
print("Higher lower bound for gamma2 (max): "+str(max(res)))
print("Implying lower bound on S: "+str(sqrt(k/(l*tau))*(2*max(res)+1+2**(d-1)*tau)/5.985))

print("Running estimations for S")
S_estimate(10,n,l,eta)

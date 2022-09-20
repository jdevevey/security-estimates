import scipy.special as sc
import math
import numpy as np


# Consider a vector uniformly chosen in the hyperball of dimension m and radius r
# Returns a Pr[ t1 <= X1 <= t2 ] for some t1 <= t2 in interval [-r, r]
def prob_slab(m,r,t1,t2):
    # want 1 - Pr[ X > t2] - Pr[X < t1] = Pr[ X > t1] - Pr[ X > t2 ]
    # note that if t >= 0, then Pr[ X > t] = 1/2*sc.betainc((m+1)/2, 1/2, 1-t**2)
    #           else, it is  Pr[ X > t] = 1 - 1/2*sc.betainc((m+1)/2, 1/2, 1-t**2)
    if t1 > t2 or abs(t1) > r or abs(t2) > r: raise NameError('Invalid interval')
    t1 = t1/r;
    t2 = t2/r;
    if t1 > 0: p1 = 1/2*sc.betainc((m+1)/2, 1/2, 1-t1**2);
    else: p1 = 1 - 1/2*sc.betainc((m+1)/2, 1/2, 1-t1**2);
    if t2 > 0: p2 = 1/2*sc.betainc((m+1)/2, 1/2, 1-t2**2);
    else: p2 = 1 - 1/2*sc.betainc((m+1)/2, 1/2, 1-t2**2);
    return p1 - p2;

# Returns x*log2(x)
def xlog2x(x):
    return (x+1e-16)*np.log((x+1e-16))/np.log(2);


# Returns the Shannon entropy of a coordinate of a vector chosen uniformly
# in the hyperball of dimension m and radius r
def compute_entropy_coordinate(m,r):
    prob = np.zeros(2*r+1);
    count = 0;
    prob[count] = prob_slab(m,r,-r,-r+1/2);
    count = count + 1;
    for i in range(-r+1,r):
        prob[count] = prob_slab(m,r,i-1/2,i+1/2);
        count = count + 1;
    prob[count] = prob_slab(m,r,r-1/2,r);
    return -np.sum(xlog2x(prob));

#k    = 4
#l    = 3
#S    = 91
#tau  = 39
#beta = math.sqrt(tau)*S
#m=256*(k+l)
#cut = 1-sc.betaincinv((m+1)/2,1/2,2**(-63))
#r = math.floor((math.sqrt(cut)+math.sqrt(cut+16**(1/m)-1))/(16**(1/m)-1)*beta )

#entropy = compute_entropy_coordinate(m,r)
#print('Entropy m='+str(m)+' and r='+str(r)+' is '+ str(entropy))
#print('Expected signature size is '+str((entropy+1)*32*l+32+32*2.5*k))
#print('Estimation is '+str((math.log(2*np.pi*r**2/m,2)+2/m)*32*l/2+32+32*2.5*k))

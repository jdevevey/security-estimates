from MSIS_security import MSIS_summarize_attacks, MSISParameterSet
from MLWE_security import MLWE_summarize_attacks, MLWEParameterSet
from math import sqrt, log, pi, floor
from scipy.special import betaincinv, gammaln
from entropy_coordinate_hyperball import compute_entropy_coordinate

# This script supports the experimental section of AC22 submission titled
# "On Rejection Sampling in Lyubashevsky's Signature Scheme" 
# 
# This script is a modification of the one for Dilithium, available at
# https://github.com/pq-crystals/security-estimates
#
# For comparisons, we rely on the following, referred to as Dilithium-NIST:
# https://pq-crystals.org/dilithium/data/dilithium-specification-round3-20210208.pdf
# 
# For comparisons, we also rely on the following, referred to as Dilithium-G:
# https://eprint.iacr.org/eprint-bin/getfile.pl?entry=2017/633&version=20170627:201152&file=633.pdf



#Select which estimators should be ran
#The key recovery estimator is much slower than the others and should be toggled off when k,l and eta are unchanged
# weak_uf by default as this is the default in Table 1 of Dilithium-NIST
weak_uf = True
strong_uf = weak_uf
key_recovery = True
size = True



#Select variants to run estimator on. Supported distributions are uniform, gaussian and continuous
#"uniform" means uniform in a hypercube intersected with integer vectors
#"gaussian" means discrete integer Gaussian
#"hyperball" means (continuous) uniform in a hyperball
variants = ["uniform",
            "gaussian-eprint",
            "gaussian",
            "hyperball"
           ]



###########################################################################
######################## PRELIMINARIES ####################################
###########################################################################

class UniformDilithiumParameterSet(object):
    def __init__(self, n, k, l, gamma1, gamma2, tau, q, eta, omega, pkdrop=0):
        self.n = n
        self.k = k
        self.l = l
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.q = q
        self.eta = eta
        self.beta = tau*eta #ell_oo bound on the shift
        # SIS ell_oo bound for weak unforgeability
        # Defined in Section 6.2 from Round 3 specification of Dilithium
        self.zeta = max(gamma1-self.beta, 2*gamma2 + 1 + 2**(pkdrop-1)*tau)
        # SIS ell_oo bound for strong unforgeability
        # Defined in Section 6.2 from Round 3 specification of Dilithium
        self.zeta_prime = max(2*(gamma1-self.beta), 4*gamma2 + 1)                
        self.pkdrop = pkdrop
        self.omega = omega


class GaussianDilithiumParameterSet(object):
    def __init__(self, n, k, l, S, gamma2, tau, q, eta, pkdrop=0, old = False):
        self.n = n
        self.k = k
        self.l = l
        self.tau = tau
        self.S = S   
        #S is a prob. upper bound (approx. the median) on the largest singular value of the secret key                                                      
        #Can be estimated using S_Estimator.py
        #beta is an ell_2 upper bound on the shift
        self.beta = S*sqrt(tau)                                          
        if old:
            self.sigma = 11*self.beta                                     
            #Old value for sigma (see Dilithium-G)
        else:
            self.sigma = 6.85*self.beta                                       
            #6.85 is computed with M=4 and eps = 1/qs = 2**(-64) using the values computed in the submission (second part of Lemma C.2)
            #1/(sqrt(2)*log(4))*(sqrt(log(2**64))+sqrt(log(2**64)+log(4)))
        self.q = q
        self.eta = eta
        self.pkdrop = pkdrop
        self.gamma2 = gamma2
        #B is from Eq.(1) of Dilithium-G
        self.B = sqrt(1.05**2*self.sigma**2*((self.l+self.k)*n)+4**(self.pkdrop-1)*(self.tau*self.n*self.k))
        # SIS ell_2 bound for unforgeability, using selftargetMSIS
        self.zeta = self.B
        # SIS ell_2 bound for strong unforgeability
        self.zeta_prime = 2*self.B
        self.delta = log(self.sigma,2)


class HyperballDilithiumParameterSet(object): 
    def __init__(self, n, k, l, S, gamma2, tau, q, eta, pkdrop=0):
        self.n = n
        self.k = k
        self.l = l
        self.tau = tau
        self.S = S                             
	  #S is a prob. upper bound (approx. the median) on the largest singular value of the secret key                                                      
        #Can be estimated using S_Estimator.py
        #beta is an ell_2 upper bound on the shift
        self.beta = S*sqrt(tau)                                          
        #ell_2 bound on the shift
        self.m = self.n*(self.k+self.l)
        # cut is the 1/eta^2 from Lemma 5.1 of the submission
        self.cut = 1-betaincinv((self.m+1)/2,1/2,2**(-63))
        self.gamma1 = (sqrt(self.cut)+sqrt(self.cut+16**(1/self.m)-1))/(16**(1/self.m)-1)*self.beta 
        #Using Lemma 4.1 from the submission
        self.q = q
        self.eta = eta
        self.pkdrop = pkdrop
        self.gamma2 = gamma2
        self.B = sqrt(self.gamma1**2+4**(self.pkdrop-1)*(self.tau*self.n*self.k))
        # SIS ell_2 bound for unforgeability
        self.zeta = self.B
        # SIS ell_2 bound for strong unforgeability
        self.zeta_prime = 2*self.B


n = 256

# Schemes parameters
####################


q = 8380417
# Maple code: q := 8380417: n:=256: l := numelems((Factors(x^n + 1) mod q)[2]): isprime(q), l, ifactor(q-1);
# Maple answer is "true, 256, 2^13*3*11*31"

# Parameters taken from Table 2 of Dilithium-NIST
UnifMediumDilithium           = UniformDilithiumParameterSet(n, 4, 4, 2**17, (q-1)/88, 39, q, 2, 80, pkdrop=13)
UnifRecommendedDilithium      = UniformDilithiumParameterSet(n, 6, 5, 2**19, (q-1)/32, 49, q, 4, 55, pkdrop=13)
UnifVeryHighDilithium         = UniformDilithiumParameterSet(n, 8, 7, 2**19, (q-1)/32, 60, q, 2, 75, pkdrop=13)


# Parameters taken from Table 2 of Dilithium-G
# old=true means the 11 constant is chosen for sigma, rather than its improved value (6.85)
GaussianDilithiumEprint1           = GaussianDilithiumParameterSet(n, 2, 2, 230,    (q-1)/512, 60, q, 7, pkdrop=11, old=True)
GaussianDilithiumEprint2           = GaussianDilithiumParameterSet(n, 3, 3, 225,    (q-1)/512, 60, q, 6, pkdrop=11, old=True)
GaussianDilithiumEprint3           = GaussianDilithiumParameterSet(n, 4, 4, 210,    (q-1)/512, 60, q, 5, pkdrop=11, old=True)
GaussianDilithiumEprint4           = GaussianDilithiumParameterSet(n, 5, 5, 145,    (q-1)/512, 60, q, 3, pkdrop=11, old=True)


q = 1038337
# Maple code: q := 1038337: n:=256: l := numelems((Factors(x^n + 1) mod q)[2]): isprime(q), l, ifactor(q-1);
# Maple answer is "true, 256, 2^11*3*13^2"

#q = 520193#238081 #202753 #125441 #520193 #254977

q = 238081

GaussianMediumDilithium            = GaussianDilithiumParameterSet(n, 4, 3, 91,    (q-1)/156, 39, q, 2, pkdrop=13)
GaussianMediumDilithium_old        = GaussianDilithiumParameterSet(n, 4, 3, 91,    (q-1)/78,  39, q, 2, pkdrop=13, old=True)

HyperballMediumDilithium          = HyperballDilithiumParameterSet(n, 4, 3, 91,    (q-1)/156, 39, q, 2, pkdrop=10)

q = 254977

GaussianRecommendedDilithium       = GaussianDilithiumParameterSet(n, 5, 4, 134,    (q-1)/96,  49, q, 3, pkdrop=11)
GaussianRecommendedDilithium_old   = GaussianDilithiumParameterSet(n, 5, 4, 134,    (q-1)/48,  49, q, 3, pkdrop=10,old=True)

HyperballRecommendedDilithium     = HyperballDilithiumParameterSet(n, 6, 4, 140,    (q-1)/96,  49, q, 3, pkdrop=12)

q = 254977

GaussianVeryHighDilithium          = GaussianDilithiumParameterSet(n, 7, 6, 111,    (q-1)/156, 60, q, 2, pkdrop=12)
GaussianVeryHighDilithium_old      = GaussianDilithiumParameterSet(n, 7, 6, 111,    (q-1)/78,  60, q, 2, pkdrop=11,old=True)

HyperballVeryHighDilithium        = HyperballDilithiumParameterSet(n, 8, 6, 115,    (q-1)/156, 60, q, 2, pkdrop=11)


all_params_uniform = [("Uniform Dilithium Medium", UnifMediumDilithium),
                   ("Uniform Dilithium Recommended", UnifRecommendedDilithium),
                   ("Uniform Dilithium Very High", UnifVeryHighDilithium)]

all_params_gaussian_eprint = [("G-Dilithium-eprint Weak", GaussianDilithiumEprint1), ("G-Dilithium-eprint Medium", GaussianDilithiumEprint2), ("G-Dilithium-eprint Recommended", GaussianDilithiumEprint3), ("G-Dilithium-eprint Very High", GaussianDilithiumEprint4)]

all_params_gaussian = [("Gaussian Dilithium Medium", GaussianMediumDilithium),
                       ("Gaussian Dilithium Medium (old)", GaussianMediumDilithium_old),
                       ("Gaussian Dilithium Recommended", GaussianRecommendedDilithium),
                       ("Gaussian Dilithium Recommended (old)", GaussianRecommendedDilithium_old),
                       ("Gaussian Dilithium Very High", GaussianVeryHighDilithium),
                       ("Gaussian Dilithium Very High (old)", GaussianVeryHighDilithium_old),
                      ]

all_params_hyperball = [("Hyperball Dilithium Medium", HyperballMediumDilithium),
                         ("Hyperball Dilithium Recommended", HyperballRecommendedDilithium),
                         ("Hyperball Dilithium Very High", HyperballVeryHighDilithium),
                        ]

all_params = { "uniform" : all_params_uniform,
               "gaussian-eprint" : all_params_gaussian_eprint, 
               "gaussian" : all_params_gaussian,
               "hyperball" : all_params_hyperball}


#########################
# Conversion to MSIS/MLWE
#########################

def Dilithium_to_MSIS(dps, strong_uf = False):
    if type(dps)==GaussianDilithiumParameterSet or type(dps)==HyperballDilithiumParameterSet:
        if strong_uf:
            return MSISParameterSet(dps.n,dps.k + dps.l, dps.k, dps.zeta_prime, dps.q, norm="l2")
        return MSISParameterSet(dps.n, dps.k + dps.l + 1, dps.k, dps.zeta, dps.q, norm="l2")
    if strong_uf:
        return MSISParameterSet(dps.n, dps.k + dps.l, dps.k, dps.zeta_prime, dps.q, norm="linf")
    else:
        return MSISParameterSet(dps.n, dps.k + dps.l + 1, dps.k, dps.zeta, dps.q, norm="linf")


def Dilithium_to_MLWE(dps):
    return MLWEParameterSet(dps.n, dps.l, dps.k, dps.eta, dps.q, distr="uniform")

text_SIS = ["BKZ block-size $b$ to break SIS","Best Known Classical bit-cost","Best Known Quantum bit-cost","Best Plausible bit-cost"]
text_LWE = ["BKZ block-size $b$ to break LWE","Best Known Classical bit-cost","Best Known Quantum bit-cost","Best Plausible bit-cost"]


##################
# Size Computation
##################

def Dilithium_Signature_Size(dps):
    """
    Computes the expected size of a signature depending on the type of distribution used.
    """
    #Size of tilde{c} is n bits
    size_c = dps.n 
    if type(dps)==UniformDilithiumParameterSet:
        # size of h is "enough bits to write all omega positions". According to the reference implem, it is omega+k bytes
        size_c_h = size_c+(dps.omega+dps.k)*8
        #Size of z = l*n*(number of bits to write one coeff)
        return size_c_h+dps.l*dps.n*(int(log(2*(dps.gamma1-dps.beta),2))+1) 
    if type(dps)==GaussianDilithiumParameterSet:
        # See p.17 of Dilithium-G
        size_c_h = size_c + 2.5*dps.n*dps.k
        return size_c_h+(2.25+dps.delta)*dps.n*dps.l 
    else:
        m = dps.n*(dps.l+dps.k)
        #Keep same approximation as Gaussians for shape of h
        size_c_h = size_c + 2.5*dps.n*dps.k 
        return size_c_h + (compute_entropy_coordinate(m,floor(dps.gamma1))+1)*dps.n*dps.l #Huffman coding of each coordinate.

def Dilithium_Entropy(dps):
    """
    Computes the optimal expected size of a signature by replacing the encoding of z by its entropy
    """
    #Size of tilde{c} is n bits
    size_c = dps.n 
    if type(dps)==UniformDilithiumParameterSet:
        # size of h is "enough bits to write all omega positions". According to the reference implem, it is omega+k bytes
        size_c_h = size_c+(dps.omega+dps.k)*8
        #Size of z = l*n*(number of bits to write one coeff)
        return size_c_h+(dps.l+dps.k)*dps.n*(int(log(2*(dps.gamma1-dps.beta),2))+1) 
    if type(dps)==GaussianDilithiumParameterSet:
        # See p.17 of Dilithium-G
        size_c_h = size_c + 2.5*dps.n*dps.k
        return size_c_h+(1.8257+dps.delta)*dps.n*dps.l
    else:
        m = dps.n*(dps.l+dps.k)
        #Keep same approximation as Gaussians for shape of h
        size_c_h = size_c + 2.5*dps.n*dps.k 
        return size_c_h + (compute_entropy_coordinate(m,floor(dps.gamma1)))*dps.n*dps.l #rANS coding of each coordinate.


def Dilithium_PK_Size(dps):
    """
    Computes the expected size of a verification key. This does not depend on the distribution used.
    """
    return (256 + dps.k*dps.n*(int(log(dps.q,2)+1-dps.pkdrop)))


# rest of script is just formatting
#############################################
######################### ANALYSIS AND REPORT 
#############################################



table_weak_SIS   = { "uniform" : [len(all_params_uniform)*[0] for i in range(4)],
                     "gaussian-eprint" : [len(all_params_gaussian_eprint)*[0] for i in range(4)],
                     "gaussian" : [len(all_params_gaussian)*[0] for i in range(4)],
                     "hyperball" : [len(all_params_hyperball)*[0] for i in range(4)]}
table_strong_SIS = { "uniform" : [len(all_params_uniform)*[0] for i in range(4)],
                     "gaussian-eprint" : [len(all_params_gaussian_eprint)*[0] for i in range(4)],
                     "gaussian" : [len(all_params_gaussian)*[0] for i in range(4)],
                     "hyperball" : [len(all_params_hyperball)*[0] for i in range(4)]}
table_LWE        = { "uniform" : [len(all_params_uniform)*[0] for i in range(4)],
                     "gaussian-eprint" : [len(all_params_gaussian_eprint)*[0] for i in range(4)],
                     "gaussian" : [len(all_params_gaussian)*[0] for i in range(4)],
                     "hyperball" : [len(all_params_hyperball)*[0] for i in range(4)]}
table_size       = { "uniform" : [0 for i in range(len(all_params_uniform))],
                     "gaussian-eprint" : [0 for i in range(len(all_params_gaussian_eprint))],
                     "gaussian" : [0 for i in range(len(all_params_gaussian))],
                     "hyperball" : [0 for i in range(len(all_params_hyperball))]}
table_entropy    = { "uniform" : [0 for i in range(len(all_params_uniform))],
                     "gaussian-eprint" : [0 for i in range(len(all_params_gaussian_eprint))],
                     "gaussian" : [0 for i in range(len(all_params_gaussian))],
                     "hyperball" : [0 for i in range(len(all_params_hyperball))]}
table_pk         = { "uniform" : [0 for i in range(len(all_params_uniform))],
                     "gaussian-eprint" : [0 for i in range(len(all_params_gaussian_eprint))],
                     "gaussian" : [0 for i in range(len(all_params_gaussian))],
                     "hyperball" : [0 for i in range(len(all_params_hyperball))]}


#For each selected scheme, build the estimate cost of selected attacks
for distribution in variants:
    j = 0
    for (scheme, param) in all_params[distribution]:
        print("\n"+scheme)
        print(param.__dict__)
        print("")
        if weak_uf:
            print("=== WEAK UF")
            v = MSIS_summarize_attacks(Dilithium_to_MSIS(param))
            for i in range(4):
                table_weak_SIS[distribution][i][j] = v[i]
        if strong_uf:
            print("=== STRONG UF")
            v = MSIS_summarize_attacks(Dilithium_to_MSIS(param, strong_uf=True))
            for i in range(4):
                table_strong_SIS[distribution][i][j] = v[i]
        if key_recovery:
            print("=== SECRET KEY RECOVERY")
            v = MLWE_summarize_attacks(Dilithium_to_MLWE(param))
            for i in range(4):
                table_LWE[distribution][i][j] = v[i]
        if size:
            print("=== SIGNATURE SIZE")
            table_size[distribution][j] = Dilithium_Signature_Size(param)
            print(table_size[distribution][j])
            print("=== SIGNATURE ENTROPY")
            table_entropy[distribution][j] = Dilithium_Entropy(param)
            print(table_entropy[distribution][j])
            print("=== PK SIZE")
            table_pk[distribution][j] = Dilithium_PK_Size(param)
            print(table_pk[distribution][j])
        j+=1


for distribution in variants:
    print(distribution.upper()+" DILITHIUM TABLE")
    print("========================")
    print("\\hline")
    for j in range(4):
        print(text_SIS[j]+"".join([" & "+str(table_weak_SIS[distribution][j][i])+" ("+str(table_strong_SIS[distribution][j][i])+")" for i in range(len(all_params[distribution]))]))
        print("\\\\")
    print("\\hline")
    for j in range(4):
        print(text_LWE[j]+"".join([" & "+str(table_LWE[distribution][j][i]) for i in range(len(all_params[distribution]))]))
        print("\\\\")
    print("\\hline")
    print("Expected signature size"+"".join([" & "+str(int(table_size[distribution][i]/8)) for i in range(len(all_params[distribution]))]))
    print("\\\\")
    print("Best signature size"+"".join([" & "+str(int(table_entropy[distribution][i]/8)) for i in range(len(all_params[distribution]))]))
    print("\\\\")
    print("Expected public key size"+"".join([" & "+str(int(table_pk[distribution][i]/8)) for i in range(len(all_params[distribution]))]))
    print("\\\\")
    print("\\hline")
    print("========================")



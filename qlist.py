import sympy

k = 0
while(512*k+1 <= 760000):
    if sympy.isprime(512*k+1):
        print(512*k+1)
        print(sympy.factorint(512*k))
    k+=1

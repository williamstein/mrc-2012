from psage.modform.hilbert.sqrt5.hmf import *
from sage.all import *

def prime_powers_of_norm(p,e):
    """
    Returns ideals in F of norm p^e.  
    If none, returns [].
    """
    if e == 0: return [F.ideal(1)]
    ideals = []
    if F.prime_above(p).residue_class_degree() == 2:
        if e % 2 == 0:
            ideals.append(F.ideal(p**(e/2)))
        else: 
            return []
    if p == 5:
        ideals.append(F.prime_above(5)**e)
    elif F.prime_above(p).norm() == p: 
        for i in range(e+1):
            ideals.append(F.primes_above(p)[0]**i*F.primes_above(p)[1]**(e-i))    
    return ideals

def ideals_of_norm(N):
    """
    Returns a list of ideals in F of norm N
    """
    if N < 1 : return []
    if N == 1 : return [F.ideal(1)]
    NN = Integer(N)
    fac = NN.factor()
    PP = []
    ideals = []
    for factor in fac:
        ppideals = prime_powers_of_norm(factor[0],factor[1])
        if ppideals == []:
            return []
        else:
            PP.append(ppideals)
    if len(fac) == 1: return PP[0]
    for L in cartesian_product_iterator(PP):
        P = F.ideal(1)
        for I in L:
            P = P * I
        ideals.append(P)
    return ideals

def ideals_norm_range(N1,N2):
    ideals = []
    for N in range(N1,N2):
        idsN = ideals_of_norm(N)
        if idsN != []:
            for ideal in idsN:
                ideals.append(ideal)
    return ideals

def ideals(N1,N2):
    for ideal in ideals_norm_range(N1,N2):
        yield ideal


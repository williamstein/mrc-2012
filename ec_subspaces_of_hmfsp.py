from psage.modform.hilbert.sqrt5.hmf import *

def quick_eig(s,V):
    v = s.basis_matrix().row(0)
    n = v*V
    b = v.list()[v.nonzero_positions()[0]]
    t = n.list()[n.nonzero_positions()[0]]
    return t/b

def small_prime(N,prime_list):
    for p in prime_list:
        if not p.divides(N):
            return p
    raise ValueError, "prime list needs more primes."
    


def old_form_dims(N,set_of_levels):
    """
    INPUT:
        `N` - the level
        set_of_leves - this should be a set
    OUTPUT:
        All possible dimensions `d` of oldforms in the space of Hilbert modular forms of level N
        We don't need to compute all divisors of `N`, just the ones which correspond to levels with
        elliptic curves.
    """
    dims = []
    levels = []
    divs = divisors(N)
    for M in divs:
        if M in set_of_levels:
            sig = sigma_0(N,M)
            dims.append(sig)
            levels.append(M)
    return dims,dims_and_levels
    
    

def foo(L,x,i,F):
    n = len(F)
    if i == n:
        L.append(x)
    else:
        p, e = F[i]
        q = 1
        for j in range(e+1):
            foo(L,x*q,i+1,F)
            q = q*p

def divisors(N):
    """
    Returns a list of divisors of N
    """
    fac = N.factor()
    L = []
    foo(L,1,0,fac)
    return L
    
def sigma_0(N,M):
    """
    INPUT:
        `N,M` - ideals of F, M divides N
    OUTPUT:
        Number of divisors of N/M
    """
    I = N/M
    return prod([ expt+1 for p, expt in factor(I) ])
    
def next_small_prime(N,p,prime_list):
    i = prime_list.index(p)
    for j in range(i+1,len(prime_list)):
        q = prime_list[j]
        if not q.divides(N):
            return q
    raise ValueError, 'prime list needs more primes'
    
    
def find_ecs(Ecs,N,H,V,p,prime_list,oldform_dims):
    """
    Goal: find elliptic curves and cut out oldforms and abelian varieties
    INPUT:
        `V` - a subspace of H = Hilbert modular forms of level N
        prime_list - a list of small primes
        
    OUTPUT:
        Appends to Ecs a list of dimension 1 subspaces 
        (i.e., eigenvalues) of H of multiplicity one with rational 
        eigenvalues satisfying the Hasse bound.
    """
    #Initialize list of elliptic curve eigenvalues
    #compute next Tp
    #first find the next small prime coprime to N
    p = next_small_prime(N,p,prime_list)
    #and Hasse Bound
    b = 2*p.norm().sqrt().n().floor()
    #remove it from the list to create a new list
    prime_list.remove(p)
    #compute Tp
    Tp = H.hecke_matrix(p)
    #decompose with respect to the space V
    Vps = Tp.decomposition_of_subspace(V)
    
    #now run through the subspaces and 
    #check, throw out, or cut down
    for S in Vps:
        S = S[0]
        if S.dimension() == 1:
            e = quick_eig(S,V)
            if e.abs() <= b:
                #found one!
                ECs.append(S)
        if S.dimension() > 1:
            #check for oldforms
            if not oldform_check(S,N):
                #repeat for the next prime and S = V
                q = p
                q = next_small_prime(N,q,prime_list)
                find_ecs(Ecs,N,H,S,p,prime_list,oldform_dims)
            
            
    
def oldform_check(ap_list,prime_list,S,N):
    """
    INPUT:
        `ap_list` - the list of a_p values for the subspace we wish to test
        `prime_list` - the corresponding list of primes in F=Q(sqrt(5))
        `S` - subspace
        `N` - Level
    OUTPUT:
        Whether or not the subspace comes from oldforms.
        Returns `True` if it can be discarded and `False` otherwise
    """
    dims,levels = old_form_dims(N)
    d = S.dimension()
    s = session() #??? Is this the correct place for this
    IRO = is_rational_old(s, ap_list, prime_list, N)
    #IRO,M = is_rational_old(s,ap_list,prime_list,N)
    dM = dims[levels.index(M)]
    if dM = d and IRO:
        return True
    #??? Do we need to still check the dimensions?  
    #Does William's code do this? Nope - so we need it to return the level.
    return False
            

def find_ecs_from_N(N,prime_list):
    H = HilbertModularForms(N)
    p = small_prime(N,prime_list)
    prime_list.remove(p)
    b = 2*p.norm().sqrt().n().floor()
    Tp = H.hecke_matrix(p)
    ECs = []
    for a in range(-b,b+1):
        #print "in for loop for a = ", a
        V = (Tp-a)
        K = V.kernel()
        d = K.dimension()
        if d == 1:
            print "appending in d = 1"
            print K
            ECs.append(K)
        elif d > 1:
            print " d = ", d, K
            S = Tp.decomposition_of_subspace(K)
            for s in S:
                find_ecs(Ecs,N,H,S,p,prime_list,oldform_dims)
    return ECs
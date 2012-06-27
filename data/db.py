########################################################
# SQLalchemy Schema
########################################################

from sqlalchemy import create_engine
engine = create_engine('postgresql://mrc@geom.math.washington.edu:6432/mrc', echo=True)
#engine = create_engine('sqlite:///a.sqlite3')

from sqlalchemy.ext.declarative import declarative_base
Base = declarative_base()

from sqlalchemy import (Boolean, Column, DateTime, Float, Integer, LargeBinary,
                        String, ForeignKey, UniqueConstraint)
from sqlalchemy.orm import relationship, backref

class Space(Base):
    """
    We represent the level space by (x,y,z), which means the ideal (x+y*a, x), with Hermite normal
    form [x,y;0,z]
    """
    __tablename__ = "spaces"
    id = Column(Integer, primary_key=True, autoincrement=True)
    norm = Column(Integer)
    dimension = Column(Integer)
    x = Column(Integer)
    y = Column(Integer)
    z = Column(Integer)
    nf_bound = Column(Integer)      # upper bound on number of rational newforms
    __table_args__ = (    
        UniqueConstraint(x,y,z),
    )

    def hmf(self):
        from psage.modform.hilbert.sqrt5.hmf import HilbertModularForms
        return HilbertModularForms(tuple_to_ideal((self.x,self.y,self.z)))


    def __repr__(self):
        return "<HMF: level=[%s,%s;0,%s], norm=%s, dimension=%s>"%(
            self.x,self.y,self.z, self.norm,  self.dimension)

class RationalNewform(Base):
    __tablename__ = "rational_newforms"
    id = Column(Integer, primary_key=True)
    space_id = Column(Integer, ForeignKey("spaces.id"))
    vector = Column(String)
    dual_vector = Column(String)
    space = relationship("Space", backref=backref("rational_newforms", order_by=id))
    __table_args__ = (
        UniqueConstraint(space_id, vector),
    )

    def store_eigenvalue(self, P, value):
        self.eigenvalues.append(RationalEigenvalue(P, value))

    def __repr__(self):
        return "<Rational newform in %s given by the vector %s>"%(self.space, self.vector)

class RationalEigenvalue(Base):
    __tablename__ = "rational_eigenvalues"
    p = Column(Integer, primary_key=True)   # residue characteristic
    r = Column(Integer, primary_key=True)   # image of a under reduction mod p, or 0 if p is inert
    norm = Column(Integer)
    value = Column(Integer)                 # the actual eigenvalue

    newform_id = Column(Integer, ForeignKey("rational_newforms.id"), primary_key=True)
    newform = relationship("RationalNewform", backref=backref("eigenvalues", order_by=norm))

    def __repr__(self):
        return "<P=(a-%s,%s), a_P=%s>"%(self.r, self.p, self.value)

    def __init__(self, P, value):
        P = fast_ideal(P)
        self.p = P.p
        self.r = P.r
        self.norm = P.norm()
        self.value = value

def create():
    Base.metadata.create_all(engine)

def session():
    from sqlalchemy.orm import sessionmaker
    Session = sessionmaker(bind=engine)
    return Session()

########################################################
# Convenience functions to use the database
########################################################
from sage.all import QQ, NumberField, polygen, dumps, gcd, parallel, divisors
from sage.rings.all import is_Ideal

x = polygen(QQ, 'x')
F = NumberField(x**2 - x - 1, 'a')
a = F.gen()

def ideal_to_tuple(N):
    v = N.free_module().echelonized_basis_matrix().list()
    return int(v[0]), int(v[1]), int(v[3])

def tuple_to_ideal(t):
    return F.ideal([t[0] + a*t[1], t[2] + a*t[2]])

def fast_ideal(P):
    if is_Ideal(P):
        import psage.number_fields.sqrt5.prime
        P = psage.number_fields.sqrt5.prime.Prime(P)
    return P
    
def store_space(s, H):
    """
    s = session
    H = hilbert modular forms space
    """
    V = Space()
    I = H.level()
    V.x, V.y, V.z = ideal_to_tuple(I)
    V.norm = int(I.norm())
    V.dimension = H.dimension()
    s.add(V)
    s.commit()

def canonically_scale(v):
    v = v.denominator() * v
    v = v / gcd(v)
    if v[v.nonzero_positions()[0]] < 0:
        v *= -1
    return v

def ns_str(t):
    return ''.join(str(t).split())

def store_rational_newform(s, M, v, vdual):
    """
    M = ambient Hilbert modular forms space
    v = rational eigenvector
    vdual = dual eigenvector
    s = session
    """
    x,y,z = ideal_to_tuple(M.level())
    V = s.query(Space).filter(Space.x==x).filter(Space.y==y).filter(Space.z==z).one()
    f = RationalNewform()
    f.vector = ns_str(canonically_scale(v))
    f.dual_vector = ns_str(canonically_scale(vdual))
    V.rational_newforms.append(f)
    return f

def get_space(s, N):
    """
    Get space with a given level N.
    
    s = session
    N = integral ideal of F
    """
    x,y,z = ideal_to_tuple(N)
    return s.query(Space).filter(Space.x==x).filter(Space.y==y).filter(Space.z==z).one()
    

def know_all_rational_eigenvectors(s, N):
    """
    Return True if we know all rational eigenvectors of level N.

    This is by definition that nf_bound equals the number of recorded
    rational newforms.
    """
    H = get_space(s, N)
    return H.nf_bound == len(H.rational_newforms)
    
########################################################
# Computing Data
########################################################
def compute_spaces(s, B1, B2):
    """
    Compute basic data (dimension) about all spaces in the given range of levels.
    
    s = session
    B1, B2 = integers
    """
    B1 = max(2,B1)
    from psage.modform.hilbert.sqrt5.hmf import HilbertModularForms
    v = F.ideals_of_bdd_norm(B2)
    for I in sum([z for _, z in v.iteritems()],[]):
        if I.norm() >= B1:
            store_space(s, HilbertModularForms(I))

def proper_divisors(N):
    return [I for I in divisors(N) if I!=1 and I!=N]

def is_rational_old(s, v, primes, N):
    """
    s = session
    v = a_p's (or ?) as output by E.aplist(...)
    primes = list of psage primes in F=Q(sqrt(5))
    N = level
    """
    w = [(P.r, P.p, int(v[i])) for i, P in enumerate(primes) if P.sage_ideal().is_coprime(N)]
    for M in proper_divisors(N):
        if not know_all_rational_eigenvectors(s, M):
            raise RuntimeError, "all newforms of level %s (norm=%s) not known so we can't tell if this is a rational newform or not"%(M, M.norm())
        # query the eigenvalues determined by primes above for rational newforms of level M.
        for f in get_space(s, M).rational_newforms:
            all_eigs_same = True
            for r, p, ap in w:
                if s.query(RationalEigenvalue).filter(
                         RationalEigenvalue.newform_id==f.id).filter(
                         RationalEigenvalue.p==p).filter(
                         RationalEigenvalue.r==r).one().value != ap:
                    all_eigs_same = False
                    break
            if all_eigs_same:
                return True
    return False

def compute_rational_eigenvectors(s, N, bound=100):
    """
    N = ideal, the level
    bound = compute this many a_p and use this many for ruling out oldforms
    
    This functions calls the "inefficient" elliptic_curve_factors method on the Hilbert
    modular forms space of level N, which has some problems:
        (1) some of the returned factors may actually be *old*
        (2) it is potentially very slow, since it computes the factors of all dimensions

    WARNING!: Just *loading* from disk the pre-computed elts of norm pi_p up to
    norm 4096 takes nearly *20 minutes*, but afterwards computing particular
    aplists for any form takes way under 20 *seconds*.
    """
    # First do a query to see if we already know *all* of the rational
    # eigenvectors.
    if know_all_rational_eigenvectors(s, N):
        print "We already know all rational eigenvectors at level %s (of norm %s)"%(N, N.norm())
        return 
    H = get_space(s, N)
    M = H.hmf()
    from psage.modform.hilbert.sqrt5.hmf import primes_of_bounded_norm
    primes = primes_of_bounded_norm(bound)
    num_forms = 0
    for E in M.elliptic_curve_factors():
        if not N.is_prime():
            # worry about E possibly actually being old
            v = E.aplist(bound)
            if is_rational_old(s, v, primes, N):
                print "*"*80
                print "skipping a form that is old"
                print "*"*80                
        f = store_rational_newform(s, M, E._S.basis()[0], E.dual_eigenvector())
        num_forms += 1
        print "*"*80
        print "form %s"%num_forms
        print "*"*80
        # store eigenvalues
        v = E.aplist(bound)
        for i in range(len(v)):
            if v[i] != '?':
                f.store_eigenvalue(primes[i], int(v[i]))
    H.nf_bound = num_forms
    s.commit()
              
def compute_rational_eigenvectors_range(s, B1, B2, bound=100):
    B1 = max(2,B1)
    v = F.ideals_of_bdd_norm(B2)
    for N in sum([z for _, z in v.iteritems()],[]):
        if N.norm() >= B1:
            print "*"*80
            print N.norm()
            print "*"*80
            compute_rational_eigenvectors(s, N, bound=bound)

def compute_rational_eigenvectors_in_parallel(B1, B2, bound=100, ncpus=8):
    @parallel(ncpus)
    def f(N):
        s = session()
        compute_rational_eigenvectors(s, N, bound=bound)
        s.commit()
    B1 = max(2,B1)
    v = F.ideals_of_bdd_norm(B2)
    for X in f([N for N in sum([z for _, z in v.iteritems()],[]) if N.norm() >= B1]):
        print X


########################################################
# SQLalchemy Schema
########################################################

from sqlalchemy import create_engine
engine = create_engine('postgresql://mrc@geom.math.washington.edu:6432/mrc2', echo=False)
#engine = create_engine('sqlite:///')

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
    number_of_rational_newforms = Column(Integer)
    rational_oldform_dimension = Column(Integer)
    time_to_compute_newforms = Column(Float)
    __table_args__ = (    
        UniqueConstraint(x,y,z),
    )

    def hmf(self):
        from  sage.modular.hilbert.sqrt5_hmf import QuaternionicModule
        return QuaternionicModule(tuple_to_ideal((self.x,self.y,self.z)))
    
    def level(self):
        return tuple_to_ideal((self.x,self.y,self.z))

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
    root_number = Column(Integer)
    root_number_prec = Column(Integer)    
    rank = Column(Integer)
    rank_prec = Column(Integer)

    __table_args__ = (UniqueConstraint(space_id, vector),)

    def store_eigenvalue(self, P, value):
        self.eigenvalues.append(RationalEigenvalue(P, int(ZZ(value))))

    def level(self):
        return self.space.level()

    def biggest_good_ap_normp_known(self):
        v = self.eigenvalues
        if len(v) == 0:
            return 0
        else:
            return v[-1].norm

    def known_aplist(self):
        """
        Output is a 2-tuple (aplist, primes,), where:
            aplist = list of integers of None
            primes = list of psage fast primes
        """
        if len(self.eigenvalues) == 0:
            return [], [], []
        B = self.eigenvalues[-1].norm
        primes = primes_of_bounded_norm(B+1)
        aplist = [None]*len(primes)
        last_p = None
        for ap in self.eigenvalues:
            P = Prime(ap.p, ap.r, first = (last_p != ap.p))
            last_p = ap.p
            i = primes.index(P)  # potentially "slow" but nothing compared to DB accesses...
            aplist[i] = ap.value
        return aplist, primes

    def __repr__(self):
        return "<Rational newform in %s given by the vector %s>"%(self.space, self.vector)

class RationalEigenvalue(Base):
    __tablename__ = "rational_eigenvalues"
    p = Column(Integer, primary_key=True)   # residue characteristic
    r = Column(Integer, primary_key=True)   # image of a under reduction mod p, or 0 if p is inert
    norm = Column(Integer)
    value = Column(Integer)                 # the actual eigenvalue

    newform_id = Column(Integer, ForeignKey("rational_newforms.id"), primary_key=True)
    newform = relationship("RationalNewform", backref=backref("eigenvalues", order_by=(norm, r)))

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
from sage.all import (QQ, ZZ, NumberField, polygen, dumps, gcd, parallel, divisors,
                      cartesian_product_iterator)
from sage.rings.all import is_Ideal
from psage.modform.hilbert.sqrt5.hmf import primes_of_bounded_norm
from psage.number_fields.sqrt5.prime import Prime

x = polygen(QQ, 'x')
F = NumberField(x**2 - x - 1, 'a')
a = F.gen()

def ideal_to_tuple(N):
    v = N.free_module().echelonized_basis_matrix().list()
    return int(v[0]), int(v[1]), int(v[3])

def tuple_to_ideal(t):
    return F.ideal([t[0] + a*t[1], t[2] + a*t[2]])

def fast_ideal(P):
    return Prime(P) if is_Ideal(P) else P
    
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
    return V

def canonically_scale(v):
    v = v.denominator() * v
    v = v / gcd(v)
    if v[v.nonzero_positions()[0]] < 0:
        v *= -1
    return v

def ns_str(t):
    return ''.join(str(t).split())

def store_rational_newform(s, M, v=None, vdual=None):
    """
    M = ambient Hilbert modular forms space
    v = rational eigenvector
    vdual = dual eigenvector
    s = session
    """
    x,y,z = ideal_to_tuple(M.level())
    V = s.query(Space).filter(Space.x==x).filter(Space.y==y).filter(Space.z==z).one()
    f = RationalNewform()
    if v is not None:
        f.vector = ns_str(canonically_scale(v))
    if vdual is not None:
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
    

def know_all_rational_newforms(s, N):
    """
    Return True if we know all rational eigenvectors of level N.
    """
    x,y,z = ideal_to_tuple(N)
    q = s.query(Space).filter(Space.x==x).filter(Space.y==y).filter(Space.z==z)
    if q.count() == 0:
        return False
    H = q.one()
    return H.number_of_rational_newforms == len(H.rational_newforms)
    
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
    from  sage.modular.hilbert.sqrt5_hmf import QuaternionicModule
    v = F.ideals_of_bdd_norm(B2)
    for I in sum([z for _, z in v.iteritems()],[]):
        if I.norm() >= B1:
            store_space(s, QuaternionicModule(I))

def proper_divisors(N):
    return [I for I in divisors(N) if I!=1 and I!=N]

def is_rational_old(s, v, primes, N):
    """
    s = session
    v = a_p's (or ?) as output by E.aplist(...)
    primes = list of psage primes in F=Q(sqrt(5))
    N = level
    """
    print 'in is_rational_oldform'
    w = [(P.r, P.p, int(v[i])) for i, P in enumerate(primes) if P.sage_ideal().is_coprime(N)]
    print 'w', w
    for M in proper_divisors(N):
        #print 'in for loop on M = ', M
        if not know_all_rational_newforms(s, M):
            #print 'no enough info, trying again'
            return False, None #I want to keep the space and try again later.
            #raise RuntimeError, "all newforms of level %s (norm=%s) not known so we can't tell if this is a rational newform or not"%(M, M.norm())
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
                print 'found oldform!',M
                return True, M
    return False, None

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
    if know_all_rational_newforms(s, N):
        print "We already know all rational eigenvectors at level %s (of norm %s)"%(N, N.norm())
        return
    H = get_space(s, N)
    M = H.hmf()
    primes = primes_of_bounded_norm(bound)
    num_forms = 0
    for E in M.elliptic_curve_factors():
        if not N.is_prime():
            # worry about E possibly actually being old
            v = E.aplist(bound)
            if is_rational_old(s, v, primes, N)[0]:
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
    H.number_of_rational_newforms = num_forms
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

def compute_more_rational_eigenvalues(s, N, bound, verb=False):
    """
    Computes rational eigenvalues for all rational eigenvectors in the
    given space up to the given bound.  Commits result to database only
    if it succeeds at computing all the (good) eigenvalues.
    """
    if not know_all_rational_newforms(s, N):
        if verb: print "Can't compute more rational eigenvalues until we know all rational eigenvectors at level %s (of norm %s)"%(N, N.norm())
        return
    else:
        print "Computing more rational eigenvalues for eigenvectors of level %s (of norm %s)"%(N, N.norm())
    NrmN = N.norm()
    H = get_space(s, N)
    M = H.hmf()
    V = M.vector_space()
    I = M._icosians_mod_p1
    import sys
    for f in H.rational_newforms:
        vector = V(eval(f.vector)) if f.vector is not None else None
        dual_vector = V(eval(f.dual_vector))
        i = dual_vector.nonzero_positions()[0]
        c = dual_vector[i]
        for P in primes_of_bounded_norm(bound):
            print P,; sys.stdout.flush()
            Ps = P.sage_ideal()
            # 1. Do we already know this eigenvalue or not?
            if s.query(RationalEigenvalue).filter(
                         RationalEigenvalue.newform_id == f.id).filter(
                         RationalEigenvalue.p == P.p).filter(
                         RationalEigenvalue.r == P.r).count() > 0:
                continue
                
            # 2. We do not know it, so compute it -- two cases:
            #   2a. if it has residue char coprime to the level, use dual eigenvector
            ap = None
            if NrmN % P.p != 0:
                ap = I.hecke_operator_on_basis_element(Ps, i).dot_product(dual_vector)/c
            
            #   2b. if it has residue char not coprime to level, use slower direct approaches
            else:
                if (Ps*Ps).divides(N):
                    ap = 0
                elif not Ps.divides(N):
                    if vector is None:
                        if verb: print "You need to compute subspace vector to compute a_p at present."
                    else:
                        j = vector.nonzero_positions()[0]
                        ap = (vector*M.hecke_matrix(Ps))[j] / vector[j]
                else:
                    # we have no algorithm directly at present to decide if this is 1 or -1.
                    if verb: print "Can't compute ap for p=%s"%P
            
            # 3. Store eigenvalue ap in the database.
            if ap is not None:
                f.store_eigenvalue(P, ap)

def compute_more_rational_eigenvalues_in_parallel(B1, B2, bound, ncpus=8):
    print "Preloading Hecke operator sets before forking."
    from sage.modular.hilbert.sqrt5 import hecke_elements
    from ideals_of_norm import ideals_of_norm
    
    for P in primes_of_bounded_norm(bound):
        print P, len(hecke_elements(P.sage_ideal()))
    
    @parallel(ncpus)
    def f(N):
        for I in ideals_of_norm(N): 
            s = session()
            compute_more_rational_eigenvalues(s, I, bound=bound)
            s.commit()
        
    for X in f(range(max(2,B1), B2+1)):
        print X
    

def compute_lseries(s, f, prec, T=1.05):
    """
    s = session
    f = rational newform object (got using the session s!)
    prec = bits of precision

    This function computes the L-series attached to the given newform
    and determine the a_p at the primes of bad reduction (up to the
    biggest p such that a_p is known) if they are not known, and saves
    those a_p in the database.  It uses prec bits of precision, and if
    not enough a_p are known in the database, then it will fail.

    DOES NOT COMMIT.  Call s.commit() after calling this to save.
    """
    
    #
    # 1. Get list of primes and the corresponding known a_P (for all P
    #    up to some bound):
    #        - query for biggest Norm(P) so we know a_P (separate function)
    #        - get list of all primes q with Norm(q) <= Norm(P).
    #        - make another list of primes q such that q exactly divides level
    #        - query database and get corresponding good a_q, or raise
    #          error if some missing.

    # Compute aplist:
    #     - aplist = list of integers of None
    #     - primes = list of psage fast primes
    #     - unknown = list of ints i such that ap[i] = None
    aplist, primes = f.known_aplist()
    unknown = [i for i in range(len(aplist)) if aplist[i] is None]

    # Check that the unknown primes all exactly divide the level.
    level = f.level()
    for i in unknown:
        P = primes[i].sage_ideal()
        if not P.divides(level):
            raise RuntimeError, "gap in list of eigenvalues a_P; missing good P=%s"%P
        if (P*P).divides(level):
            raise RuntimeError, "gap: additive a_P for P=%s not set"%P
    
    # 2. For each possibility for bad a_P, construct the L-series,
    #    making a list of the L-series that actually work.
    #        - use cartesian product iterator over [[-1,1]]*n
    #        - we use a custom L-series class deriving from what
    #          is in psage defined above
    lseries_that_work = []
    for bad_ap in cartesian_product_iterator([[-1,1]]*len(unknown)):
        print "bad_ap = ", bad_ap
        aplist1 = list(aplist)
        for i, j in enumerate(unknown):
            aplist1[j] = bad_ap[i]
        for eps in [-1,1]:
            print "trying eps = ", eps
            L = LSeries(level=level, aplist=aplist1, primes=primes, prec=prec, root_number=eps)
            try:
                L._function(prec=prec, T=T)
            except RuntimeError:
                print "Definitely does not satisfy functional equation..."
            else:
                print "It seems to satisfy functional equation..."
                lseries_that_work.append((L, bad_ap))
    if len(lseries_that_work) > 1:
        print "WARNING: %s choices of bad a_p's seem to work -- functional equation doesn't nail down one choice.\nWe will check numerically that all choices are consistent, and only save a_p that we are certain about."%len(lseries_that_work)
    if len(lseries_that_work) == 0:
        raise RuntimeError, "no choices of bad a_p's seem to work -- please increase precision!"
    
    # 3. If *exactly one* L-series works, save to the database the
    #    corresponding bad a_p, and return the L-series.  Otherwise,
    #    raise an error.
    # Only save bad_ap's that are the same for all L-series.
    
    
    # save missing a_p, which we now know, to the database
    print "saving missing eigenvalues, which we just determined, to the database..."
    for i, j in enumerate(unknown):
        # only store ones such that multiple distinct choices didn't work.
        if len([bad_ap[i] for _, bad_ap in lseries_that_work]) == 1:
            f.store_eigenvalue(primes[j], lseries_that_work[0][1][i])

    # return one of the L-series, after doing a double check that they are
    # all basically the same numerically.
    if len(lseries_that_work) > 1:
        ts = lseries_that_work[0][0].taylor_series(prec=prec)
        for i in range(1, len(lseries_that_work)):
            ts2 = lseries_that_work[i][0].taylor_series(prec=prec)
            if ts != ts2:
                raise RuntimeError, "ts=%s, ts2=%s"%(ts, ts2)

    L = lseries_that_work[0][0]
    
    # store epsilon factor (sign of f.e.) to the database
    if f.root_number is None:
        root_number = L.epsilon()
        print "saving root number (=%s) to database..."%root_number
        f.root_number = root_number
        f.root_number_prec = root_number_prec

    if f.rank is None:
        rank = L.analytic_rank(prec=prec)
        print "saving rank (=%s) to database..."%rank
        f.rank = rank
        f.rank_prec = prec
        
        
    return L


from psage.lseries.eulerprod import LSeriesAbstract, prime_below
class LSeries(LSeriesAbstract):
    def __init__(self, level, aplist, primes, root_number, prec):
        self._level = level
        self._aplist = aplist
        self._primes = primes
        self._prec = prec
        self._root_number = root_number
        LSeriesAbstract.__init__(self, conductor = level.norm() * 25,
                                 hodge_numbers = [0]*2+[1]*2, weight = 2, epsilon = root_number,
                                 poles = [], residues=[], base_field = F, prec=prec)

    def _primes_above(self, p):
        """
        Return the primes above p.  This function returns a special
        optimized prime of the ring of integers of Q(sqrt(5)).
        """
        from psage.number_fields.sqrt5.prime import primes_above
        return primes_above(p)

    def _local_factor(self, P, prec):
        T = ZZ['T'].gen()
        try:
            i = self._primes.index(P)
        except ValueError:
            msg = "Constructing L-series of form of norm level %s to precision %s: not enough precision -- don't know a_P for P (=%s) of norm %s"%(self._level.norm(), self._prec, P, P.norm())
            print msg
            raise Exception, msg
        ap = self._aplist[i]
        q = P.norm()
        Ps = P.sage_ideal()
        p = prime_below(Ps)
        f = ZZ(q).ord(p)
        if Ps.divides(self._level):
            return 1 - ap*(T**f)
        else:
            return 1 - ap*(T**f) + q*(T**(2*f))
                           
    
def find_all_curves(B1, B2, verb=False, ncpu=20, maxtime=60*5):
    from  sage.modular.hilbert.sqrt5_hmf import F, QuaternionicModule
    from sage.all import cputime
    import ideals_of_norm

    v = range(max(2,B1), B2+1)
        
    def g(N):
        from sage.misc.misc import alarm, cancel_alarm
        alarm(maxtime)
        from sqlalchemy import create_engine
        from sqlalchemy.orm import sessionmaker
        engine = create_engine('postgresql://mrc@geom.math.washington.edu:6432/mrc2', echo=False)
        s = sessionmaker(bind=engine)()
        t = cputime()
        x,y,z = ideal_to_tuple(N)
        q = s.query(Space).filter(Space.x==x).filter(Space.y==y).filter(Space.z==z)
        if q.count() > 0:
            M = q.one()
            if M.number_of_rational_newforms == len(M.rational_newforms):
                return "Already done with level %s (norm = %s)"%(N, N.norm())
            H = M.hmf()
        else:
            H = QuaternionicModule(N)
            M = store_space(s, H)
        V, rational_oldform_dimension = H.dual_rational_newforms(verb=verb)
        for vdual, aplist in V:
            f = RationalNewform()
            for p, ap in aplist:
                f.store_eigenvalue(p, ap)
            M.rational_newforms.append(f)
            f.dual_vector = ns_str(canonically_scale(vdual))
        M.number_of_rational_newforms = len(V)
        M.rational_oldform_dimension = int(rational_oldform_dimension)
        M.time_to_compute_newforms = cputime(t)
        s.commit()

        cancel_alarm()
        return N.norm(), cputime(t), len(V)

    @parallel(ncpu)
    def f(N):
        for I in ideals_of_norm.ideals_of_norm(N):
            return g(I)
             
    if ncpu > 1:
        for X in f(v):
            print X
    else:
        for X in v:
            print X, f(X)

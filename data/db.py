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
# Use the database
########################################################
from sage.all import QQ, NumberField, polygen, dumps, gcd
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

def store_spaces(s, B1, B2):
    """
    Compute basic data about all spaces in the given range of levels.
    
    s = session
    B1, B2 = integers
    """
    B1 = max(2,B1)
    from psage.modform.hilbert.sqrt5.hmf import HilbertModularForms
    v = F.ideals_of_bdd_norm(B2)
    for I in sum([z for _, z in v.iteritems()],[]):
        if I.norm() >= B1:
            store_space(s, HilbertModularForms(I))
    
def canonically_scale(v):
    v = v.denominator() * v
    v = v / gcd(v)
    if v[v.nonzero_positions()[0]] < 0:
        v *= -1
    return v

def ns_str(t):
    return ''.join(str(t).split())

def store_rational_newform(s, M, v):
    """
    M = ambient Hilbert modular forms space
    v = rational eigenvector
    s = session
    """
    x,y,z = ideal_to_tuple(M.level())
    V = s.query(Space).filter(Space.x==x).filter(Space.y==y).filter(Space.z==z).one()
    f = RationalNewform()
    f.vector = ns_str(canonically_scale(v))
    V.rational_newforms.append(f)

def get_space(s, N):
    """
    Get space with a given level N.
    
    s = session
    N = integral ideal of F
    """
    x,y,z = ideal_to_tuple(N)
    return s.query(Space).filter(Space.x==x).filter(Space.y==y).filter(Space.z==z).one()
    

    
    

########################################################
# SQLalchemy Schema
########################################################

from sqlalchemy import create_engine
engine = create_engine('postgresql://mrc@geom.math.washington.edu:6432/mrc', echo=True)

from sqlalchemy.ext.declarative import declarative_base
Base = declarative_base()

from sqlalchemy import (Boolean, Column, DateTime, Float, Integer, LargeBinary, String, ForeignKey)
from sqlalchemy.orm import relationship, backref

class Space(Base):
    __tablename__ = "spaces"
    id = Column(Integer, primary_key=True)
    level = Column(String, primary_key=True)
    dim = Column(Integer)
    sobj = Column(LargeBinary)

def create():
    Base.metadata.create_all(engine)

def session():
    from sqlalchemy.orm import sessionmaker
    Session = sessionmaker(bind=engine)
    return Session()


########################################################
# Use the database
########################################################
from sage.all import QQ, NumberField, polygen

x = polygen(QQ, 'x')
F = NumberField(x**2 - x - 1, 'a')
a = F.gen()

def ideal_to_tuple(N):
    return tuple(N.free_module().echelonized_basis_matrix().list())

def tuple_to_ideal(t):
    return F.ideal([t[0] + a*t[1], t[2] + a*t[3]])

def store_space(M):
    s = session()
    A = Space()
    A.level = ideal_to_tuple(M.level())

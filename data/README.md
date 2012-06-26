We will put summary data here.  
It just must not take up too much space, and must be in plain text format.

Representation
--------------

We make one choice, which is to represent Q(sqrt(5)) as F=Q(a)=Q[x]/(x^2-x-1), 
so a=root of x^2-x-1=(1+sqrt(5))/2.

We represent elements in terms of a.

We represent an ideal by a 3-tuple of integers, using the upper
triangular Hermite normal form of the free module it defines, in terms
of the basis 1,a.

Thus:

sage: F.<a> = NumberField(x^2 - x - 1)
sage: I = F.primes_above(11)[0]; I
Fractional ideal (-3*a + 2)
sage: M = I.free_module()
sage: M.echelonized_basis_matrix()
[ 1  4]
[ 0 11]
sage: M.echelonized_basis_matrix().list()
[1, 4, 0, 11]
sage: F.ideal([1+4*a, 11*a])
Fractional ideal (-3*a + 2)
sage: import db
sage: db.ideal_to_tuple(I)
(1, 4, 11)
sage: db.tuple_to_ideal((1,4,11))
Fractional ideal (-3*a + 2)


Database
--------

PostreSQL on port 6432 on geom.math.washington.edu.  The username is
"mrc", and one can connect from Snowbird or on geom, but not from
anywhere else.






from sympy import *
# from numpy import *

a = Symbol('a')
b = Symbol('b')
u = Symbol('u')
v = Symbol('v')
w = Symbol('w')
c = Symbol('a')
d = Symbol('d')


u = array([a,b])
v = array([a,-b])

# Axiom 1
print "The sum of u and v, denoted by u + v, is in V."
print u+v
print " "

# Axiom 2
print "u + v = v + u"
print u+v==v+u
print " "
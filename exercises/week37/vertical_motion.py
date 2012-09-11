# Exercise 24

from math import pi
import sys

def Stokes_const(d, mu, rho_b, V):
    return (-3*pi*d*mu)/(rho_b*V)


def quadratic_const(C_D, rho, A, rho_b, V):
    return (-C_D*rho*A)/(2*rho_b*V)

rho_b = float(sys.argv[1])
m = float(sys.argv[2])
V = float(m/rho_b)
d = float(sys.argv[3])
A = pi*(d/2)**2
rho = float(sys.argv[4])
C_D = float(sys.argv[5])
v0 = float(sys.argv[6])


"""def skrivUt(a, b, c, d, e, f, g, h):
    print "rho_b = %g" % a
    print "m = %g" % b
    print "V = %g" % c
    print "d = %g" % d
    print "A = %g" % e
    print "rho = %g" % f
    print "C_D = %g" % g
    print "v0 = %g" % h
    return True

test = skrivUt(rho_b, m, V, d, A, rho, C_D, v0)
print ""
print test
"""    
    


from numpy import *
I = 0.1
a = 1.0
dt = 0.01
x = a*dt

C2 = I*(1 + sqrt(x**2 +1))/(2*sqrt(x**2 + 2))
C1 = I - C2

A1 = -x - sqrt(x**2 + 1)
A2 = -x + sqrt(x**2 + 1)


T = 4
N = int(round(T/dt))
t = linspace(0, T, N+1)
u_e = I*exp(-a*t)

for n in range(0, N+1):
    print "n = %d" % n
    print ("%.12f        %.12f        %.12f") % (C1*(A1**n), C2*(A2**n), u_e[n])


# Exercise 24

from math import pi

"""
def Stokes_const(d, mu, rho_b, V) :
    return (3*pi*d*mu)/(rho_b*V)
"""

"""
def quadratic_const(C_D, rho, A, rho_b, V) :
    return (C_D*rho*A)/(2*rho_b*V)
"""
def solver(dt, rho_b, m, V, d, A, rho, C_D, v_0, g=9.81, v):
    v[0] = v_0
    for i in range(len(v)-1):
        Re = (rho*d*v[i])/mu
        if Re < 1:
            a = (3*pi*d*mu)/(rho_b*V)
        else:
            a = (C_D*rho*A)/(2*rho_b*V)

        v[i+1] = (v[i] + dt*b)/(1 + dt*a*abs(v[i]))
        
        return v







dt = float(sys.argv[1])
T = 20.0 # seconds of simulation
N = T/dt

t = linspace(0,T,N)
v = zeros(N)

rho_b = 1003  # kg/m^3
m = 80  # kg
V = float(m/rho_b) #m^3
d = 0.5  # m
A = pi*(d/2)**2
rho = 0.79  # kg/m^3
C_D = 1.2
v_0 = 0
g = 9.81  # m/s^2

b = g*((rho/rho_b) - 1)


soluton = solver(dt, rho_b, m, V, d, A, rho, C_D, v_0, g, v)

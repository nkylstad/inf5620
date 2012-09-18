# Exercise 15

from numpy import *
from matplotlib.pyplot import *

def leapfrogSolver(I, a, b, T, dt):
    """ Solve u'(t) = -a(t)u(t) + b(t), u(0) = I.
    Uses Forward Euler to compute u(1)."""
    dt = float(dt)
    N = int(round(T/dt))
    T = N*dt
    u = zeros(N+1)
    t = linspace(0, T, N+1)

    u[0] = I

    u[1] = u[0] + dt*(-a[0]*u[0] + b[0])

    for n in range(1, N):
        u[n+1] = u[n-1] + 2*dt*(-a[n]*u[n] + b[n])
    return u



def test_linear_solution():
    """
    Test problem where u = c*t + I is the exact solution, to be reproduced (to machine precision) by any relevant method.
    """

    def exact_solution(t):
        return c*t + I

    def a(t):
        return t**0.5

    def b(t):
        return c + a(t)*exact_solution(t)


    I = 0.1; dt = 0.1; c = -0.5
    T = 4
    N = int(T/dt)
    t = linspace(0,T, N+1)

    a1 = a(t)
    b1 = b(t)
    
    u = leapfrogSolver(I=I, a=a1, b=b1, T=N*dt, dt=dt)
    u_e = exact_solution(t)

    difference = abs(u_e - u).max()
    print difference

    figure()
    plot(t, u)
    plot(t, u_e)
    legend(['numerical', 'exact'])
    savefig('linear_test.png')
    show()



test_linear_solution()

def test_case():
    """ Test case when u'(t) = -u(t) + 1, u(0) = 0"""

    def exact_solution(t):
        return -exp(-t) + 1

    I = 0
    dt = 0.1
    T = 4
    N = int(T/dt)
    t = linspace(0,T, N+1)
    a = ones(N+1)
    b = ones(N+1)
    
    u = leapfrogSolver(I=I, a=a, b=b, T=N*dt, dt=dt)
    u_e = exact_solution(t)

    difference_case = abs(u_e - u).max()
    print difference_case

    figure()
    plot(t, u)
    plot(t, u_e)
    legend(['numerical', 'exact'])
    savefig('case_test.png')
    show()


test_case()


    

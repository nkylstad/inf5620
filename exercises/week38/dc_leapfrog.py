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

    u[0] = I
    u[1] = u[0] + dt*(-a[0]*u[0] + b[0])

    for n in range(1, N):
        u[n+1] = u[n-1] + 2*dt*(-a[n]*u[n] + b[n])
    return u


def conv_rate(t, I, a, T, u, u_e):
    E_values = []
    dt_values=[0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001]
    e = u_e - u
    for dt in dt_values:
        E = sqrt(dt*sum(e**2))
        E_values.append(E)
        #print "e = %g" % sum(e**2)
        #print "E = %g" % E

    #print E_values

    m = len(dt_values)
    r = [log(E_values[i-1]/E_values[i])/log(dt_values[i-1]/dt_values[i]) for i in range(1, m, 1)]

    print "Convergence rates:"
    print ' '.join(['%.2f' % i for i in r])
    return r
    
    


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

    conv_rate(t, I, a, T, u_e, u)



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
    conv_rate(t, I, a, T, u_e, u)
    


test_linear_solution()
test_case()



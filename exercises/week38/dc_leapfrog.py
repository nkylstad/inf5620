﻿# Exercise 15

from numpy import *
from matplotlib.pyplot import *


def leapfrogSolver(I, a, b, T, dt):
    """ Solve u'(t) = -a(t)u(t) + b(t), u(0) = I.
    Uses Forward Euler to compute u(1)."""
    dt = float(dt)
    N = int(round(T/dt))
    T = N*dt
    t = linspace(0, T, N+1)
    u = zeros(N+1)
    u[0] = I
    u[1] = u[0] + dt*(-a[0]*u[0] + b[0])
    for n in range(1, N):
        u[n+1] = u[n-1] + 2*dt*(-a[n]*u[n] + b[n])
    return u, t


def exact_solution(t, I, a, b):
    u_e = zeros(N+1)
    for i in range(N+1):
        u_e[i] = b[i]/a[i] + (I - (b[i]/a[i]))*exp(-a[i]*t[i])
    return u_e


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
    


def explore(I, a, b, T, dt, makeplot=True): 
    """ 
    Run a case with the solver, compute error measure, 
    and plot the numerical and exact solutions (if makeplot=True). 
    """ 
    u, t = leapfrogSolver(I, a, b, T, dt) # Numerical solution
    u_e = exact_solution(t, I, a) 
    e = u_e - u 
    E = sqrt(dt*sum(e**2)) 
    if makeplot: 
        figure() # create new plot 
        t_e = linspace(0, T, 1001) # fine mesh for u_e 
        u_e = exact_solution(t_e, I, a) 
        plot(t, u, ’r--o’) # red dashes w/circles 
        plot(t_e, u_e, ’b-’) # blue line for exact sol. 
        legend([’numerical’, ’exact’]) 
        xlabel(’t’) 
        ylabel(’u’) 
        title(’theta=%g, dt=%g’ % (theta, dt)) 
        theta2name = {0: ’FE’, 1: ’BE’, 0.5: ’CN’} 
        savefig(’%s_%g.png’ % (theta2name[theta], dt)) 
        savefig(’%s_%g.pdf’ % (theta2name[theta], dt)) 
        savefig(’%s_%g.eps’ % (theta2name[theta], dt)) 
        show() 
    return E 


def main(): 
    #I, a, T, makeplot, dt_values = read_command_line()
    
    I = 0.1; dt = 0.1; T = 4; N = int(T/dt)
    a = ones(N+1)
    n = ones(N+1)

    dt_values = [0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001]
    
    E_values = [] 
    for dt in dt_values: 
        E = explore(I, a, b, T, dt, makeplot=False) 
        E_values.append(E)
    # Compute convergence rates 
    m = len(dt_values) 
    r = [log(E_values[i-1]/E_values[i])/log(dt_values[i-1]/dt_values[i]) for i in range(1, m, 1)] 
    print 'Convergence rates for Leapfrog'
    print ' '.join(['%.2f' % i for i in r]) 
    return r 

main()

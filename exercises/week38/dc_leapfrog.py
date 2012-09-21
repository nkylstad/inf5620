from numpy import *
from matplotlib.pyplot import*


def solver(I, a, b, T, dt):
    """
    Solve u'(t) = -a(t)u(t) + b(t), u(0) = I,
    for T in (0,T] with steps dt.
    a and b are python functions of t.
    """
    dt = float(dt)
    N = int(round(T/dt))
    T = N*dt
    u = zeros(N+1)
    t = linspace(0, T, N+1)

    u[0] = I
    u[1] = u[0] + dt*(b(t[0]) - a(t[0])*u[0])

    for n in range(1, N):
        u[n+1] = u[n-1] + 2*dt*(b(t[n]) - a(t[n])*u[n])
    return u, t


def test_linear_solution():
    """
    Test problem where u=c*t + I is the exact solution, to be
    reproduced (to machine precision) by any relevant method.
    """
    def exact_solution(t):
        return c*t + I

    def a(t):
        return t**0.5  # can be arbitrary

    def b(t):
        return c + a(t)*exact_solution(t)
    I = 0.1; dt = 0.1; c = -0.5
    T = 4
    N = int(round(T/dt))
    u,t = solver(I=I, a=a, b=b, T=N*dt, dt=dt)
    u_e = exact_solution(t)
    difference = abs(u_e - u).max()
    print difference


def test_case():
    """
    Test problem where a = b = 1, such that
    u'(t)=-u(t)+1, u(0)=0.
    """

    def exact_solution(t):
        return 1 - exp(-t)

    def a(t):
        return 1

    def b(t):
        return 1

    I = 0; dt = 0.1;
    T = 4
    N = int(round(T/dt))
    u, t = solver(I=I, a=a, b=b, T=N*dt, dt=dt)
    u_e = exact_solution(t)

    plot(t, u, 'r--')
    plot (t, u_e, 'b-')
    legend(['numerical', 'exact'], 'best')
    xlabel('t')
    ylabel('u(t)')
    title('Leapfrog scheme for a=b=1')
    savefig('LF_test_case.png')
    show()


def exact_solution(t, I, a, b):
    u_e = zeros(len(t))
    for i in range(len(t)):
        u_e[i] = b(t[i])/a(t[i]) + (I- b(t[i])/a(t[i]))*exp(-a(t[i])*t[i])
    return u_e


def find_error(I, a, b, T, dt):
    u, t = solver(I, a, b, T, dt)
    u_e = exact_solution(t, I, a, b)
    e = u_e - u
    E = sqrt(dt*sum(e**2))
    return E

def convergence_rates():

    def a(t):
        return 1

    def b(t):
        return 1
    
    I = 0.1;
    T = 4
    dt_list = [0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001]
    E_values = []
    
    for dt in dt_list:
        E = find_error(I, a, b, T, dt)
        E_values.append(E)

    m = len(dt_list)
    r = [log(E_values[i-1]/E_values[i])/
         log(dt_list[i-1]/dt_list[i])
         for i in range(1, m, 1)]

    print r

convergence_rates()
test_case()



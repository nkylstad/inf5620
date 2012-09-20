from numpy import *
from matplotlib.pyplot import *

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


def exact_solution(t, I, a, b):
    u_e = zeros(len(t))
    for i in range(len(t)):
        u_e[i] = b(t[i])/a(t[i]) + (I- b(t[i])/a(t[i]))*exp(-a(t[i])*t[i])
    return u_e


def plot_solution(u, t, u_e, t_e, dt, T):
    figure()
    plot(t, u, 'r--')
    plot(t_e, u_e, 'b-')
    legend(['numerical', 'exact'], 'best')
    xlabel('t')
    ylabel('u(t)')
    title('Leapfrog scheme for T = %g, dt = %g' % (T, dt))
    savefig('LF_T%g_dt%g.png' % (T, dt))
    show()


def dt_experiments(I, a, b, T, dt_values):

    M = 1000
    t_e = linspace(0, T, M+1)
    u_e = exact_solution(t_e, I, a, b)
    
    for dt in dt_values:
        N = int(round(T/dt))
        u, t = solver(I, a, b, T, dt)

        plot_solution(u, t, u_e, t_e, dt, T)

def T_experiments(I, a, b, T_values, dt):
    M = 1000
    
    
    for T in T_values:
        N = int(round(T/dt))
        u, t = solver(I, a, b, T, dt)
        t_e = linspace(0, T, M+1)
        u_e = exact_solution(t_e, I, a, b)
        plot_solution(u, t, u_e, t_e, dt, T)


def main():
    I = 0.1
    T = 4
    dt = 0.001
    dt_values = [0.5, 0.1, 0.05, 0.01, 0.005, 0.001]
    T_values = [2, 4, 6, 8, 10]

    def a(t):
        return 1

    def b(t):
        return 1

    #dt_experiments(I, a, b, T, dt_values)
    T_experiments(I, a, b, T_values, dt)

main()
        

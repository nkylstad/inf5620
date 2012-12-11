from dolfin import *
import numpy as np
import time


def run_simulation(version):
	if version=="test-case":

		# Create mesh
		mesh = UnitCube(8,8,8)

		# Create function space
		V = VectorFunctionSpace(mesh, "Lagrange", 1)

		# Elasticity parameters
		E, nu = 1.0, 0.1
		mu, lmbda = E/(2.0*(1.0 + nu)), E*nu/((1.0 + nu)*(1.0 - 2.0*nu))
		rho = 1.0
		alpha = 1.0
		gamma = -lmbda*alpha/(2*(mu + lmbda))
		dt = 0.1

		# Source term
		b = Expression(("2.0*alpha*x[0]", "2.0*gamma*x[1]", "2.0*gamma*x[2]"),\
		 alpha=alpha, gamma=gamma)

		# Initial conditions
		u0 = Expression(("0.0", "0.0", "0.0"))
		v0 = Expression(("0.0", "0.0", "0.0"))

		exact = Expression(("t*t*alpha*x[0]", "t*t*gamma*x[1]", "t*t*gamma*x[2]"),\
		 alpha=alpha, gamma=gamma, t=0.0)

		solver(mesh, V, u0, v0, b, mu, lmbda, rho, alpha, gamma, dt, exact=exact)



def solver(mesh, V, u0, v0, f, mu, lmbda, rho, alpha, gamma, dt, version="lvp", exact=None):
	
	# Create test and trial functions
	u = TrialFunction(V)
	v = TestFunction(V)
	n = FacetNormal(mesh)

	def eps(u):
		return (1.0/2.0)*(nabla_grad(u) + transpose(nabla_grad(u)))

	# Stress
	def sigma(u):
	 	return  2*mu*eps(u) + lmbda*tr(eps(u))*Identity(v.cell().d)


	# Boundary conditions
	def left_boundary(x, on_boundary):
		tol = 1E-14
		return on_boundary and abs(x[0]) < tol

	def right_boundary(x, on_boundary):
		tol = 1E-14
		return on_boundary and abs(x[0] - 1) < tol

	def update_boundary(t=0.0):
		r = Expression(("0.0","t*t*gamma*x[1]","t*t*gamma*x[2]"), gamma=gamma, t=t)
		l = Expression(("t*t*alpha","t*t*gamma*x[1]","t*t*gamma*x[2]"), \
			alpha=alpha, gamma=gamma, t=t)

		bc_left = DirichletBC(V, r, left_boundary)
		bc_right = DirichletBC(V, l, right_boundary)
		bcs = [bc_left, bc_right]
		return bcs



	#------------------- Special case for first time step ----------------------#
	u2 = interpolate(u0, V)

	F = (dot(u,v) - dot(u2,v) - dt*dot(v0,v) - (dt**2/2)*dot(f,v))*dx + \
	dt**2/(2.0*rho)*inner(sigma(u2),nabla_grad(v))*dx - dt**2/(2.0*rho)*dot(dot(sigma(u2),n),v)*ds

	a, L = lhs(F), rhs(F)

	A = assemble(a)
	b = assemble(L)

	bcs = update_boundary(t=dt)
	for bc in bcs:
		bc.apply(A,b)

	u1 = Function(V)
	solve(A, u1.vector(), b)

	if exact:
		exact.t = dt
		u_e = interpolate(exact, V)
		max_error = np.max(u_e.vector().array() - u1.vector().array())
		print "first step: ",max_error
	#----------------------------------------------------------------------------#


	# Governing balance equation for the remaining time steps

	F = (dot(u,v) - 2*dot(u1,v) + dot(u2,v) - dt**2*dot(f,v))*dx + \
		(dt**2/rho)*inner(sigma(u1), grad(v))*dx - dot(dot(sigma(u),n),v)*ds
	#a = lhs(F)
	#A = assemble(a)
	u = Function(V)
	T = 1.0
	t = 2*dt
	counter = 2
	t1 = time.clock()
	while t <= T:
		print t
		# Extract linear and bilinear forms from F
		L = rhs(F)
		bcs = update_boundary(t=t)
		b = assemble(L)
		for bc in bcs:
			bc.apply(A,b)
		solve(A, u.vector(), b)
		u2.assign(u1)
		u1.assign(u)
		#u = TrialFunction(V)

		if exact:
			exact.t = t
			u_e = interpolate(exact, V)
			max_error = np.max(u_e.vector().array() - u1.vector().array())
			print max_error

		t += dt
		counter += 1

		# plot(u1, mode="displacement", axes=True, title="t=%g"%(t-dt),\
		#  label="t=%g"%(t-dt), basename="time_dependent_elast")
		# #title("t="+(t-dt))
		# # plot(mesh)
		# interactive()
	t2 = time.clock()
	print "Time taken: ", t2-t1


run_simulation("test-case")
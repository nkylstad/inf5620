from dolfin import *
import numpy as np

# Create mesh
mesh = UnitCube(8,8,8)

# Create function space
V = VectorFunctionSpace(mesh, "Lagrange", 1)

# Create test and trial functions
u = TrialFunction(V)
v = TestFunction(V)
n = FacetNormal(mesh)



# Elasticity parameters
E, nu = 1.0, 0.1
mu, lmbda = E/(2.0*(1.0 + nu)), E*nu/((1.0 + nu)*(1.0 - 2.0*nu))
rho = 1.0
alpha = 1.0
gamma = -lmbda*alpha/(2*(mu + lmbda))

# Source term
b = Expression(("0.0", "0.0", "0.0"))
def eps(u):
	return (1.0/2.0)*(nabla_grad(u) + transpose(nabla_grad(u)))

# Stress
def sigma(u):
 	return  2*mu*eps(u) + lmbda*tr(eps(u))*Identity(v.cell().d)


# Governing balance equation
F = inner(sigma(u), grad(v))*dx - rho*dot(b,v)*dx #- dot(dot(sigma(u),n),v)*ds


# Extract linear and bilinear forms from F
a, L = lhs(F), rhs(F)

# Boundary conditions
def left_boundary(x, on_boundary):
	tol = 1E-14
	return on_boundary and abs(x[0]) < tol

def right_boundary(x, on_boundary):
	tol = 1E-14
	return on_boundary and abs(x[0] - 1) < tol

#gamma = 0.0
c = Expression(("0.0","gamma*x[1]","gamma*x[2]"), gamma=gamma)
r = Expression(("alpha","gamma*x[1]","gamma*x[2]"), alpha=alpha, gamma=gamma)
bc_left = DirichletBC(V, c, left_boundary)
bc_right = DirichletBC(V, r, right_boundary)
bcs = [bc_left, bc_right]


u = Function(V)
problem = LinearVariationalProblem(a, L, u, bcs=bcs)
solver = LinearVariationalSolver(problem)
solver.solve()

plot(u, mode="displacement", axes=True)
plot(mesh)
interactive()


V_ = TensorFunctionSpace(mesh, "Lagrange", 1)
sigma_ = project(sigma(u), V_)
size = len(sigma_.vector().array())
print size
tol = 1E-14
counter = 0
for i in range(size-1):
	if sigma_.vector().array()[i] <= tol:
		counter += 1
print counter

test = np.zeros((3,3))
test[0,0] = sigma_.vector().array()[0]
test[0,1] = sigma_.vector().array()[729]
test[0,2] = sigma_.vector().array()[1458]
test[1,0] = sigma_.vector().array()[2187]
test[1,1] = sigma_.vector().array()[2916]
test[1,2] = sigma_.vector().array()[3655]
test[2,0] = sigma_.vector().array()[4374]
test[2,1] = sigma_.vector().array()[5103]
test[2,2] = sigma_.vector().array()[5832]

print test
print sigma_.vector().array()[728]
print sigma_.vector().array()[729]
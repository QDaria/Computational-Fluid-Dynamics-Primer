from __future__ import print_function
from dolfin import *

# Print log messages only from the root process in parallel
parameters["std_out_all_processes"] = False;

# Load mesh from file
mesh = Mesh("Bi2_01.xml")

# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)

# Define trial and test functions
u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

# Set parameter values
dt = 0.01
T = 3
nu = 0.01

# Define time-dependent pressure boundary condition
#p_in = Expression("sin(3.0*t)", t=0.0)
p_in = Constant(0.54)
# p_in = Expression("x[0]*x[0]")

# Define boundary conditions
noslip  = DirichletBC(V, (0, 0),
                      "on_boundary && \
                       (x[0] < 2.0 - DOLFIN_EPS | x[1] < DOLFIN_EPS | \
                       (x[0] < - 2.0 + DOLFIN_EPS && x[1] > 2.0 - DOLFIN_EPS))")
inflow  = DirichletBC(Q, p_in, "x[0] < DOLFIN_EPS")
outflow = DirichletBC(Q, 0, " x[1] <  23.2 - x[0] && \
                       x[1] <  - 8.6 + x[0]*tan(pi/18)")
# outflow2 = DirichletBC(Q, 0, " x[1] > -6 + x[0]*tan(pi/18) - DOLFIN_EPS")
bcu = [noslip]
# bcp = [inflow, outflow1, outflow2 ]
bcp = [inflow, outflow]

# Create functions
u0 = Function(V)
u1 = Function(V)
p1 = Function(Q)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))

# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx + inner(grad(u0)*u0, v)*dx + \
     nu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure update
a2 = inner(grad(p), grad(q))*dx
L2 = -(1/k)*div(u1)*q*dx

# Velocity update
a3 = inner(u, v)*dx
L3 = inner(u1, v)*dx - k*inner(grad(p1), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Use amg preconditioner if available
prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

# Create files for storing solution
ufile = File("results/velocity.pvd")
pfile = File("results/pressure.pvd")

# Time-stepping
t = dt
while t < T + DOLFIN_EPS:

    # Update pressure boundary condition
    p_in.t = t

    # Compute tentative velocity step
    begin("Computing tentative velocity")
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u1.vector(), b1, "gmres", "default")
    end()

    # Pressure correction
    begin("Computing pressure correction")
    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in bcp]
    solve(A2, p1.vector(), b2, "cg", prec)
    end()

    # Velocity correction
    begin("Computing velocity correction")
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "gmres", "default")
    end()

    # Plot solution
    plot(p1, title="Pressure", rescale=True)
    plot(u1, title="Velocity", rescale=True)

    # Save to file
    ufile << u1
    pfile << p1

    # Move to next time step
    u0.assign(u1)
    t += dt
    print("t =", t)

# Hold plot
interactive()

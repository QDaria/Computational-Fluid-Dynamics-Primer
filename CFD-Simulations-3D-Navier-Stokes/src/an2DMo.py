from __future__ import print_function
from dolfin import *

# Print log messages only from the root process in parallel
parameters["std_out_all_processes"] = False;

# Load mesh from file
mesh = Mesh("AnMicro.xml")

# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "Lagrange", 1)
Q = FunctionSpace(mesh, "Lagrange", 1)

# Define trial and test functions
u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

# Set parameter values
T = 4000 	;dt = 20.
nu = 1.004	;Re = 84
a = 400.; b = 185.

Dh = (2.*a*b)/(a + b)
Umax = Re/(nu*Dh)  
Umean = Umax*2./3
#print Umax, Umean


# Define time-dependent pressure boundary condition
# p_in = Expression("sin(3.0*t)", t=0.0)
#p_in = Constant(0.5)

#----------------------------------------------------------------------------
# inflow 
class inlet(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[0] + DOLFIN_EPS, 0.) and on_boundary

# outflows
class outlet1(SubDomain):
	def inside(self, x, on_boundary):
		return x[1] > 1924.788 - x[0] - 1e-3 and on_boundary

class outlet2(SubDomain):
	def inside(self, x, on_boundary):
		return x[1] < - 862.154 + x[0]*tan(pi/18) + 1e-3 and on_boundary

# No-slip BC"s
class top(SubDomain):
	def inside(self, x, on_boundary):
		return (x[1] > 200 - 1e-3) and between(x[0], (0., 597.11)) and on_boundary
		   
class bottom(SubDomain):
	def inside(self, x, on_boundary):
		return (x[1] < -200 + 1e-3) and between(x[0], (0., 300.)) and on_boundary
		   
class out1L(SubDomain):
	def inside(self, x, on_boundary):
		return x[1] > -397.109 + x[0] - 1e-3 and \
		   	   between(x[0], (597., 1161.)) and \
		       between(x[1], (200., 764.)) and on_boundary

class out1R(SubDomain):
	def inside(self, x, on_boundary):
		return x[1] < -538.53 + x[0] + 1e-3 and \
		   	   between(x[0], (651.062788375, 1232.)) and \
		       between(x[1], (293.128838748, 695.)) and on_boundary

class out2L(SubDomain):		   
	def inside(self, x, on_boundary):
		return x[1] < 1501.384 - x[0]*tan(4*pi/9) + 1e-3 and \
		   	   between(x[0], (300., 404.2)) and \
		       between(x[1], (-790.9, -200.)) and on_boundary	

class out2R(SubDomain):
	def inside(self, x, on_boundary):
		return x[1] > 3804.89273914 - x[0]*tan(4*pi/9) - 1e-3 and \
		       between(x[0], (651.062788375, 798.113)) and \
		   	   between(x[1], (-721.43, 112.5322 )) and on_boundary 
class an(SubDomain):
	def inside(self, x, on_boundary):
		return x[1]**2  > - (x[0] - 767.819926237)**2 + 300.**2 - 130 and \
		   	   between(x[0], (723.215, 1067.82)) and \
		       between(x[1], (-300., 293.13)) and on_boundary
#----------------------------------------------------------------------------
ff = FacetFunction("size_t", mesh, 0)
inlet = inlet() 	; inlet.mark(ff, 1)
outlet1 = outlet1() ; outlet1.mark(ff, 2)
outlet2 = outlet2()	; outlet2.mark(ff, 3)
top = top()			; top.mark(ff, 4)
bottom = bottom()	; bottom.mark(ff, 5)
out1L = out1L()		; out1L.mark(ff, 6)
out1R = out1R()		; out1R.mark(ff, 7)
out2L = out2L()		; out2L.mark(ff, 8)
out2R = out2R()		; out2R.mark(ff, 9)
an = an()			; an.mark(ff, 10)
#----------------------------------------------------------------------------

u_in = Expression(("Umax*(1 -  pow(x[1],2)/pow(200,2))","0"), Umax=Umax) 

pc21 = DirichletBC(Q, 0.05 , outlet1)
pc22 = DirichletBC(Q, 0.005 , outlet2)
# velocity
bc11 = DirichletBC(V, u_in, inlet)
bc12 = DirichletBC(V, (0,0), inlet)
bc21 = DirichletBC(V, (0,0), outlet1)
bc22 = DirichletBC(V, (0,0), outlet2)
# noslip
bc01 = DirichletBC(V, (0,0), top)
bc02 = DirichletBC(V, (0,0), bottom)
bc03 = DirichletBC(V, (0,0), out1L)
bc04 = DirichletBC(V, (0,0), out1R)
bc05 = DirichletBC(V, (0,0), out2L)
bc06 = DirichletBC(V, (0,0), out2R)
bc07 = DirichletBC(V, (0,0), an)
#bcs = []
bcu = [bc01, bc02, bc03, bc04, bc05, bc06, bc07, bc11] #, bc21, bc22
#bcs= [bc01, bc02, bc03, bc04, bc05, bc06, bc07, bc12] #, bc21, bc22
#collect pressure bcs:
bcp = [pc21, pc22]

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
ufile = File("Mo005_0005/Mo005_0005.pvd")
#pfile = File("Mo1_6results/pressure.pvd")

# Time-stepping
t = dt
while t < T + DOLFIN_EPS:

    # Update pressure boundary condition
    #p_in.t = t

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
    #pfile << p1

    # Move to next time step
    u0.assign(u1)
    t += dt
    print("t =", t)

# Hold plot
interactive()

"""
Constant lines (ordered), and circle:

c0 = 5.97109248121
c1 = 23.218975998
c2 = - 8.62154061343
c3 = -3.97109248121
c4 = -5.38530604356
c5 = 15.3674649314
c6 = 38.7471849205
x0 = 7.67819926237


 (x[0] - 7.67819926237)**2 + x[1]**2 > 9 |

noslip  = DirichletBC(V, (0, 0),
                      "on_boundary && pow(x[0] - 7.67819926237,2) + pow(x[1],2) > (9-1e-7) | \
                       (x[0] < 5.97109248121 - DOLFIN_EPS | \
						x[1] < -3.97109248121 + x[0] + DOLFIN_EPS | \
						x[1] < -5.38530604356 + x[0] - DOLFIN_EPS | \
						x[1] < 15.3674649314 - x[0]*tan(4*pi/9) + DOLFIN_EPS | \
						x[1] < 38.7471849205 - x[0]*tan(4*pi/9) - DOLFIN_EPS)")
inflow  = DirichletBC(Q, p_in, "x[0] < DOLFIN_EPS")
outflow = DirichletBC(Q, 0, " x[1] <  23.218975998 - x[0] && \
                       x[1] < - 8.62154061343 + x[0]*tan(pi/18)")


# outflow2 = DirichletBC(Q, 0, " x[1] > -6 + x[0]*tan(pi/18) - DOLFIN_EPS")
bcu = [noslip]
# bcp = [inflow, outflow1, outflow2 ]
bcp = [inflow, outflow]
"""

						









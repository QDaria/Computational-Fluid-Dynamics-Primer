from ..NSfracStep import *

# Mesh
mesh = Mesh("mesh/an2Df.xml")
Re = 84 ; nu = 1.004
a = 400.; b = 185.
Dh = (2.*a*b)/(a + b)
Umax = Re/(nu*Dh) ; Umean = Umax*2./3

NS_parameters.update(dict(
    nu = nu,
    T = 100, 
    dt = 10,
    Re = Re,
    Dh = Dh,
    Umean = Umean,
    Umax = Umax,
    folder = "an2D",
    max_iter = 1,
    plot_interval = 1,
    save_step = 1000000000000,      
    checkpoint = 1000000000000,
    velocity_degree = 1,
    use_krylov_solvers = True
  )
)

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

u_in = Expression("Umax*(1 -  pow(x[1],2)/pow(200,2))", Umax=Umax) 

def create_bcs(V, Q, sys_comp, **NS_namespace):
	bcs = dict((ui, []) for ui in sys_comp)

	#pc11 = DirichletBC(Q, 0., inlet)
	pc21 = DirichletBC(Q, 0. , outlet1)
	pc22 = DirichletBC(Q, 0. , outlet2)
	# velocity
	bc11 = DirichletBC(V, u_in, inlet)
	bc12 = DirichletBC(V, 0., inlet)

	bc21 = DirichletBC(V, 0., outlet1)
	bc22 = DirichletBC(V, 0., outlet2)
	# noslip
	bc01 = DirichletBC(V, 0., top)
	bc02 = DirichletBC(V, 0., bottom)
	bc03 = DirichletBC(V, 0., out1L)
	bc04 = DirichletBC(V, 0., out1R)
	bc05 = DirichletBC(V, 0., out2L)
	bc06 = DirichletBC(V, 0., out2R)
	bc07 = DirichletBC(V, 0., an)
	# collect no slip:
	bcs['u0'] = [bc01, bc02, bc03, bc04, bc05, bc06, bc07, bc11] #, bc21, bc22
	bcs['u1'] = [bc01, bc02, bc03, bc04, bc05, bc06, bc07, bc12] #, bc21, bc22
    #collect pressure bcs:
	bcs['p'] = [pc21, pc22] 
	return bcs
#----------------------------------------------------------------------------

def pressure_hook(u_, dot, ds, bcs, p0bc, p1bc, **NS_namespace):
    n = FacetNormal(mesh)
    flux2 = dot(u_, n)*ds(2)
    flux3 = dot(u_, n)*ds(3)
    Q2 = assemble(flux2) 
    Q3 = assemble(flux3) 
    C =  5e-5
    p0 = 1e-100 #1e-10 
    R2 = (C*Q2 + p0)
    R3 = (C*Q3 + p0)
    p0bc.vector()[:] = R2
    p1bc.vector()[:] = R3
    R = R3/R2
    #l = []
    print R
 
    bcs['p'][0].set_value(p0bc)
    bcs['p'][1].set_value(p1bc)



#----------------------------------------------------------------------------
def pre_solve_hook(mesh, OasisFunction, u_, ds, Q, **NS_namespace):
    Vv = VectorFunctionSpace(mesh, 'CG', 1)
    #bc11 = DirichletBC(Vv, (0,0), ff, 1)
    #bc12 = DirichletBC(Vv, (0,0), ff, 1)
    #bc21 = DirichletBC(Vv, (0,0), ff, 2)
    #bc22 = DirichletBC(Vv, (0,0), ff, 3)
    bc01 = DirichletBC(Vv, (0,0), ff, 4)
    bc02 = DirichletBC(Vv, (0,0), ff, 5)
    bc03 = DirichletBC(Vv, (0,0), ff, 6)
    bc04 = DirichletBC(Vv, (0,0), ff, 7)
    bc05 = DirichletBC(Vv, (0,0), ff, 8)
    bc06 = DirichletBC(Vv, (0,0), ff, 9) 
    bc07 = DirichletBC(Vv, (0,0), ff, 10)
    ds = ds[ff]
    p0bc = Function(Q)
    p1bc = Function(Q)

    bcs_vv = [bc01, bc02, bc03, bc04, bc05, bc06, bc07] # bc21, bc22
    #fu = File("coupled/coupled.pvd")
    fu = File("/media/mo/TOSHIBA EXT/master/2D1e3Re84/2D1e3Re84.pvd")
    #gu = File("Pan2D84.xdmf")
    return dict(Vv=Vv, uv=OasisFunction(u_, Vv), bcs_vv=bcs_vv, fu=fu, p0bc=p0bc, p1bc=p1bc, ds=ds)# gu=gu,
#----------------------------------------------------------------------------
def start_timestep_hook(t, **NS_namespace):
	pass
   
def temporal_hook(tstep, q_, u_, uv, Vv, plot_interval, bcs_vv, fu, **NS_namespace): # gu,

	if tstep % plot_interval == 0:
		uv.assign(project(u_, Vv, bcs=bcs_vv))
        plot(uv, title="Velocity")
        fu << uv




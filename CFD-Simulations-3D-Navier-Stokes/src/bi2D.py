from ..NSfracStep import *

# Mesh
mesh = Mesh("mesh/BiMicro.xml")
#mesh = Mesh("mesh/bi3D.xml")
Re = 84
nu = 1.004 # e+6
#mu = 1.002
NS_parameters.update(dict(
    nu = nu,
    #mu = mu,
    T = 10000.0, #0.0001,
    dt = 100.00, #0.0000001,
    Re = Re,
    folder = "bi2D",
    max_iter = 1,
    plot_interval = 1,
    save_step = 1000000000000,
    checkpoint = 1000000000000,
    velocity_degree = 1,
    use_krylov_solvers = True
  )
)

# inflow and outflow

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
		       between(x[1], (112.5321, 695.)) and on_boundary

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

inlet = inlet()
outlet1 = outlet1()
outlet2 = outlet2()
top = top()
bottom = bottom()
out1L = out1L()
out1R = out1R()
out2L = out2L()
out2R = out2R()

ff = FacetFunction("size_t", mesh, 0)
inlet.mark(ff, 1)
outlet1.mark(ff, 2)
outlet2.mark(ff, 3)
top.mark(ff, 4)
bottom.mark(ff, 5)
out1L.mark(ff, 6)
out1R.mark(ff, 7)
out2L.mark(ff, 8)
out2R.mark(ff, 9)
# Re22=0.008730,Re84=0.333343
u_in = Expression("0.333343*(1 -  pow(x[1],2)/pow(200,2))")

# outflow boundary value for the pressure


def create_bcs(V, Q, sys_comp, **NS_namespace):
	bcs = dict((ui, []) for ui in sys_comp)
    # pressure
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
	
    # collect no slip:
	bcs['u0'] = [bc01, bc02, bc03, bc04, bc05, bc06, bc11]
	bcs['u1'] = [bc01, bc02, bc03, bc04, bc05, bc06, bc12] #bcs['u0']  
    #collect pressure bcs:
	bcs['p'] = [pc21, pc22] 
	return bcs

def pressure_hook(u_, dot, ds, bcs, p0bc, p1bc, **NS_namespace):
    n = FacetNormal(mesh)
    flux2 = dot(u_, n)*ds(2)
    flux3 = dot(u_, n)*ds(3)
    Q2 = assemble(flux2)
    Q3 = assemble(flux3) 
    C = 1e-5 
    p0 = 1.01125e-4 
    R2 = (C*Q2 + p0)
    R3 = (C*Q3 + p0)

    p0bc.vector()[:] = R2
    p1bc.vector()[:] = R3
    #print R2, R3, Q2, Q3
    bcs['p'][0].set_value(p0bc)
    bcs['p'][1].set_value(p1bc)
    return 

def pre_solve_hook(mesh, OasisFunction, u_, ds, Q, **NS_namespace):
    Vv = VectorFunctionSpace(mesh, 'CG', 1)
    bc11 = DirichletBC(Vv, (0,0), ff, 1)
    #bc12 = DirichletBC(Vv, (0,0), ff, 1)
    bc21 = DirichletBC(Vv, (0,0), ff, 2)
    bc22 = DirichletBC(Vv, (0,0), ff, 3)
    bc01 = DirichletBC(Vv, (0,0), ff, 4)
    bc02 = DirichletBC(Vv, (0,0), ff, 5)
    bc03 = DirichletBC(Vv, (0,0), ff, 6)
    bc04 = DirichletBC(Vv, (0,0), ff, 7)
    bc05 = DirichletBC(Vv, (0,0), ff, 8)
    bc06 = DirichletBC(Vv, (0,0), ff, 9)
    ds = ds[ff]
    p0bc = Function(Q)
    p1bc = Function(Q)

    bcs_vv = [bc01, bc02, bc03, bc04, bc05, bc06, bc11, bc21, bc22]
    #fu = File("bi2D84T300000dt100/mcsbi2D84T300000dt100.pvd")
    fu = File("/media/mo/TOSHIBA EXT/master/2D/bi84/bi84.pvd")
    #gu = File("Pbi2D22.xdmf")
    return dict(Vv=Vv, uv=OasisFunction(u_, Vv), bcs_vv=bcs_vv, fu=fu, p0bc=p0bc, p1bc=p1bc, ds=ds)#gu=gu,

def start_timestep_hook(t, **NS_namespace):
    
    #p_in.t = t
	#u_in.t = t
	pass
def temporal_hook(tstep, q_, u_, uv, Vv, plot_interval, bcs_vv, fu, **NS_namespace):#gu,

	if tstep % plot_interval == 0:
		#plot(q_['p'], title="Pressure")
		#gu << q_['p']
		#uv()
		uv.assign(project(u_, Vv, bcs=bcs_vv))
        plot(uv, title="Velocity")
        fu << uv
#u_nodal_values = u_.vector()
"""
u_nodal_values = u.vector()
u_array = u_nodal_values.array()
coor = mesh.coordinates()
if mesh.num_vertices() == len(u_array):
	for i in range(mesh.num_vertices()):
		print "u(%8g,%8g) = %g" % (coor[i][0], coor[i][1], u_array[i])
"""


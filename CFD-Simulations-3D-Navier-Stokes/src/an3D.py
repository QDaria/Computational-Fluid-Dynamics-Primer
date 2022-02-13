from ..NSfracStep import *

# Mesh
mesh = Mesh("mesh/AnMicro3D.xml")
#mesh = Mesh("mesh/an5020.xml")
Re = 170 ; nu = 1.004; mu = Constant(1.002)
#a = 400.; b = 185.
#Dh = (2.*a*b)/(a + b)
#Umax = Re/(nu*Dh) ; Umean = Umax*2./3
NS_parameters.update(dict(
    nu = nu,
    mu = mu,
    T = 20000., 
    dt = 20, 
    Re = Re,
    folder = "an3D",
    max_iter = 1,
    plot_interval = 1,
    save_step = 1000000000000,
    checkpoint = 1000000000000,
    velocity_degree = 1,
    use_krylov_solvers = True
  )
)
# inlet profile
from numpy import cosh, cos
set_log_active(False)
#dpdx = Constant(-1.8e-12) #Re=2.2
#dpdx = Constant(-6.9e-11) #Re84 = 38*Re2_2
#dpdx = Constant(-8.37e-11) # Re=110 = 50*Re2_2
dpdx = Constant(-1.685e-4) #Re170
#dpdx = Constant(-2.782e-10) # Re=340
a = 200. ; b = 92.5
#Umean = Constant(333343)
factor = 16.*a*a/mu(0)/pi**3*(-dpdx(0))
# Re22=8730, Re84=333343, Re170=674624, Re220=873043, Re340=1349249
# Re22=0.00000314,Re84=0.00012,Re110=0.000157, Re170=0.000243, Re220=0.000314, Re340=0.000486
ue_code = """
class U : public Expression
{
  public :

	double a, b, mu , dpdx ;

  void eval(Array <double>& values , const Array <double>& x) const
	{
	  double u = 0.;
      double factor = 16.*a*a/mu/pow(DOLFIN_PI, 3)*(-dpdx);
      for (std::size_t i=1; i< 600 ; i=i+2)
      {
        u += pow(-1, (i-1)/2 % 2)*(1.-cosh(i*DOLFIN_PI*(x[2]-b)/2./a)/
           cosh(i*DOLFIN_PI*b/2./a))*cos(i*DOLFIN_PI*x[1]/2./a)/pow(i, 3);
    
      }
	  values[0] = u*factor;
	}
}; """


u_c = Expression(ue_code)
u_c.a = float(a) 
u_c.b = float(b)
u_c.mu = float(mu(0)) 
u_c.dpdx = float(dpdx(0))
#u_c.Umean = float(Umean(0))
#u_c = Umean(0)*u_c
#u_in = Expression("333343*sin(pi*t)", t=0.)
# plot(u_[0], mesh, interactive=True), plot(u_[1], mesh, interactive=True), plot(u_[2], mesh, interactive=True), plot(u_c, mesh, interactive=True)
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

# No-slip 
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
		   	   between(x[0], (831.659443104, 1232.)) and \
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
		   	   between(x[1], (-721.43, - 296.665527018)) and on_boundary

# lids
class over(SubDomain):
	def inside(self, x, on_boundary):
		return x[2] > 185 - 1e-3 and on_boundary

class under(SubDomain):
	def inside(self, x, on_boundary):
		return x[2] < DOLFIN_EPS and on_boundary

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
over = over()		; over.mark(ff, 10)
under = under()		; under.mark(ff, 11)
#----------------------------------------------------------------------------


#u_in = Expression("Umean*sin(pi*t)",Umean=0.333343, t=0.)
def create_bcs(V, Q, sys_comp, **NS_namespace):
	bcs = dict((ui, []) for ui in sys_comp)
    # pressure
	#pc11 = DirichletBC(Q, 0.0, inlet)
	pc21 = DirichletBC(Q, 0.0 , outlet1)
	pc22 = DirichletBC(Q, 0.0 , outlet2)
	# velocity
	bc11 = DirichletBC(V, u_c, inlet)
	bc12 = DirichletBC(V, 0., inlet)
	bc13 = DirichletBC(V, 0., inlet)
	bc21 = DirichletBC(V, 0., outlet1)
	bc22 = DirichletBC(V, 0., outlet2)
    
	# noslip
	bc01 = DirichletBC(V, 0., top)
	bc02 = DirichletBC(V, 0., bottom)
	bc03 = DirichletBC(V, 0., out1L)
	bc04 = DirichletBC(V, 0., out1R)
	bc05 = DirichletBC(V, 0., out2L)
	bc06 = DirichletBC(V, 0., out2R)
	bc07 = DirichletBC(V, 0., over)
	bc08 = DirichletBC(V, 0., under)
	bc09 = DirichletBC(V, 0., an)
    # collect no slip:
	bcs['u0'] = [bc01, bc02, bc03, bc04, bc05, bc06, bc07, bc08, bc09, bc11]#, bc21, bc22
	bcs['u1'] = [bc01, bc02, bc03, bc04, bc05, bc06, bc07, bc08, bc09, bc12]#, bc21, bc22
	bcs['u2'] = [bc01, bc02, bc03, bc04, bc05, bc06, bc07, bc08, bc09, bc13]#, bc21, bc22
    #collect pressure bcs:
	bcs['p'] = [pc21, pc22] 
	return bcs

def pressure_hook(u_, dot, ds, bcs, p0bc, p1bc, **NS_namespace):
    n = FacetNormal(mesh)
    flux2 = dot(u_, n)*ds(2)
    flux3 = dot(u_, n)*ds(3)
    Q2 = assemble(flux2)
    Q3 = assemble(flux3) 
    C = 1e-5 #
    p0 = 1e-12 #
    R2 = (C*Q2 + p0)
    R3 = (C*Q3 + p0)

    p0bc.vector()[:] = R2
    p1bc.vector()[:] = R3
    R = R3/R2
    print R 
    bcs['p'][0].set_value(p0bc)
    bcs['p'][1].set_value(p1bc)
     

def pre_solve_hook(mesh, OasisFunction, u_, ds, Q, **NS_namespace):
    Vv = VectorFunctionSpace(mesh, 'CG', 1)
    #bc11 = DirichletBC(Vv, (0,0,0), ff, 1)
    #bc12 = DirichletBC(Vv, (0,0), ff, 1)
    #bc21 = DirichletBC(Vv, (0,0,0), ff, 2)
    #bc22 = DirichletBC(Vv, (0,0,0), ff, 3)
    bc01 = DirichletBC(Vv, (0,0,0), ff, 4)
    bc02 = DirichletBC(Vv, (0,0,0), ff, 5)
    bc03 = DirichletBC(Vv, (0,0,0), ff, 6)
    bc04 = DirichletBC(Vv, (0,0,0), ff, 7)
    bc05 = DirichletBC(Vv, (0,0,0), ff, 8)
    bc06 = DirichletBC(Vv, (0,0,0), ff, 9)
    bc07 = DirichletBC(Vv, (0,0,0), ff, 10)
    bc08 = DirichletBC(Vv, (0,0,0), ff, 11)
    bc09 = DirichletBC(Vv, (0,0,0), ff, 12)
    ds = ds[ff]
    p0bc = Function(Q)
    p1bc = Function(Q)

    bcs_vv = [bc01, bc02, bc03, bc04, bc05, bc06, bc07, bc08, bc09]#, bc11, bc21, bc22
    #fu = File("an3D52.xdmf")
    fu = File("/media/mo/TOSHIBA EXT/master/3D/1e5an3DRe170dt20ts1000fine/1e5an3DRe170dt20ts1000fine.pvd")
    #fu = File("/media/mo/TOSHIBA EXT/master/3D/u5020an3D84.xdmf")
    #gu = File("/media/mo/TOSHIBA EXT/master/3D/Uan3D110.xdmf")
    return dict(Vv=Vv, uv=OasisFunction(u_, Vv), bcs_vv=bcs_vv, fu=fu, p0bc=p0bc, p1bc=p1bc, ds=ds)#, ds=ds#gu=gu,

def start_timestep_hook(t, **NS_namespace):
	#u_in.t = t    
	pass
    #p_in.t = t
	#u_in.t = t

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
		print 'u(%8g,%8g) = %g' % (coor[i][0], coor[i][1], u_array[i])
"""


from dolfin import *

# Mesh
mesh = Mesh("mesh/AnMicro3D.xml")

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
#-721.42538074
# lids
class over(SubDomain):
	def inside(self, x, on_boundary):
		return x[2] > 185 - 1e-3 and on_boundary

class under(SubDomain):
	def inside(self, x, on_boundary):
		return x[2] < DOLFIN_EPS and on_boundary

class an(SubDomain):
	def inside(self, x, on_boundary):
		return x[1]**2  > - (x[0] - 767.819926237)**2 + 300.**2 - 40 and \
		   	   between(x[0], (723.215, 1067.82)) and \
		       between(x[1], (-300., 293.13)) and on_boundary
inlet = inlet()
outlet1 = outlet1()
outlet2 = outlet2()
top = top()
bottom = bottom()
out1L = out1L()
out1R = out1R()
out2L = out2L()
out2R = out2R()
over = over()
under = under()
an = an()


ff1 = FacetFunction("size_t", mesh, 0)
ff2 = FacetFunction("size_t", mesh, 0)
ff3 = FacetFunction("size_t", mesh, 0)
ff4 = FacetFunction("size_t", mesh, 0)
ff5 = FacetFunction("size_t", mesh, 0)
ff6 = FacetFunction("size_t", mesh, 0)
ff7 = FacetFunction("size_t", mesh, 0)
ff8 = FacetFunction("size_t", mesh, 0)
ff9 = FacetFunction("size_t", mesh, 0)
ff10 = FacetFunction("size_t", mesh, 0)
ff11= FacetFunction("size_t", mesh, 0)
ff12= FacetFunction("size_t", mesh, 0)

inlet.mark(ff1, 1)
outlet1.mark(ff2, 2)
outlet2.mark(ff3, 3)

# no slip
top.mark(ff4, 4)
bottom.mark(ff5, 5)
out1L.mark(ff6, 6)
out1R.mark(ff7, 7)
out2L.mark(ff8, 8)
out2R.mark(ff9, 9)
over.mark(ff10, 10)
under.mark(ff11, 11)
an.mark(ff12, 12)



plot(ff10, interactive=True, title = "over")
plot(ff11, interactive=True, title = "under")
plot(ff12, interactive=True, title = "aneurysm")

"""
plot(ff1, interactive=True, title = "inlet")
plot(ff2, interactive=True, title = "outlet1")
plot(ff3, interactive=True, title = "outlet2")
plot(ff4, interactive=True, title = "top")
plot(ff5, interactive=True, title = "bottom")
plot(ff6, interactive=True, title = "out1L")
plot(ff7, interactive=True, title = "out1R")
plot(ff8, interactive=True, title = "out2L")
plot(ff9, interactive=True, title = "out2R")

"""

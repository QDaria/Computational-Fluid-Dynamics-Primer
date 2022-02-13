from dolfin import *

# Mesh
mesh = Mesh("mesh/BiMicro.xml")

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
		   	   between(x[0], (597, 1161.)) and \
		       between(x[1], (200., 764.)) and on_boundary

class out1R(SubDomain):
	def inside(self, x, on_boundary):
		return x[1] < -538.53 + x[0] + 1e-3 and \
		   	   between(x[0], (651.062788375, 1232.)) and \
		       between(x[1], (112.5321, 695)) and on_boundary

class out2L(SubDomain):		   
	def inside(self, x, on_boundary):
		return x[1] < 1501.384 - x[0]*tan(4*pi/9) + 1e-3 and \
		   	   between(x[0], (300., 404.2)) and \
		       between(x[1], (-790.9, -200.)) and on_boundary	

class out2R(SubDomain):
	def inside(self, x, on_boundary):
		return x[1] > 3804.89273914 - x[0]*tan(4*pi/9) - 1e-3 and \
		       between(x[0], (651.062788375, 798.113)) and \
		   	   between(x[1], (-721.43, 112.5322 )) and on_boundary #and not between(x[1], (-3., 2.9313)) 

inlet = inlet()
outlet1 = outlet1()
outlet2 = outlet2()
top = top()
bottom = bottom()
out1L = out1L()
out1R = out1R()
out2L = out2L()
out2R = out2R()


ff1 = FacetFunction("size_t", mesh, 0)
ff2 = FacetFunction("size_t", mesh, 0)
ff3 = FacetFunction("size_t", mesh, 0)
ff4 = FacetFunction("size_t", mesh, 0)
ff5 = FacetFunction("size_t", mesh, 0)
ff6 = FacetFunction("size_t", mesh, 0)
ff7 = FacetFunction("size_t", mesh, 0)
ff8 = FacetFunction("size_t", mesh, 0)
ff9 = FacetFunction("size_t", mesh, 0)


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



plot(ff1, interactive=True, title = "inlet")
plot(ff2, interactive=True, title = "outlet1")
plot(ff3, interactive=True, title = "outlet2")
plot(ff4, interactive=True, title = "top")
plot(ff5, interactive=True, title = "bottom")
plot(ff6, interactive=True, title = "out1L")
plot(ff7, interactive=True, title = "out1R")
plot(ff8, interactive=True, title = "out2L")
plot(ff9, interactive=True, title = "out2R")

# BiMN.geo
"""
lc = DefineNumber[ 0.25, Name "Parameters/lc" ];
Point(1) = {0, -2, 0, 0.25};
Point(2) = {3, -2, 0, 0.1};
Point(3) = {4.041889066, -7.90884651807, 0, 0.25};
Point(4) = {7.98112007805, -7.2142538074, 0, 0.25};
Point(5) = {6.51062788375, 1.12532184019, 0, 0.1};
Point(6) = {12.31659443106, 6.9312883875, 0, 0.25};
Point(7) = {11.6094876499, 7.63839516869, 0, 0.25};
Point(8) = {5.97109248121, 2, 0, 0.1};
Point(9) = {0, 2, 0, 0.25};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 1};
Line Loop(11) = {8, 9, 1, 2, 3, 4, 5, 6, 7};


Plane Surface(12) = {11};
"""


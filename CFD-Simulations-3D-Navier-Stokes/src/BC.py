#-------------------------------------BC---------------------------------------
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

class an(SubDomain):
    def inside(self, x, on_boundary):
        return x[1]**2  > - (x[0] - 767.819926237)**2 + 300.**2 - 130 and \
            between(x[0], (723.215, 1067.82)) and \
                between(x[1], (-300., 293.13)) and on_boundary
# lids
class over(SubDomain):
    def inside(self, x, on_boundary):
        return x[2] > 185 - 1e-3 and on_boundary

class under(SubDomain):
    def inside(self, x, on_boundary):
        return x[2] < DOLFIN_EPS and on_boundary



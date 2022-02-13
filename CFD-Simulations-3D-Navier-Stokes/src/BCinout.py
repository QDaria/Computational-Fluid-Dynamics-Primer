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
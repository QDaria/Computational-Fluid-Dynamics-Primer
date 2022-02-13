# lids
class over(SubDomain):
    def inside(self, x, on_boundary):
        return x[2] > 185 - 1e-3 and on_boundary

class under(SubDomain):
    def inside(self, x, on_boundary):
        return x[2] < DOLFIN_EPS and on_boundary


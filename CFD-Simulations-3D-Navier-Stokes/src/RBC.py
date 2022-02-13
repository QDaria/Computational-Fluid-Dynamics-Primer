def pressure_hook(u_, dot, ds, bcs, p0bc, p1bc, **NS_namespace):
    n = FacetNormal(mesh)
    flux2 = dot(u_, n)*ds(2)
    flux3 = dot(u_, n)*ds(3)
    Q2 = assemble(flux2)
    Q3 = assemble(flux3)
    C = 1.2e-5 #
    p0 = 1e-100 #
    R2 = (C*Q2 + p0)
    R3 = (C*Q3 + p0)
    p0bc.vector()[:] = R2
    p1bc.vector()[:] = R3
    R = R3/R2
    print R
    bcs['p'][0].set_value(p0bc)
    bcs['p'][1].set_value(p1bc)
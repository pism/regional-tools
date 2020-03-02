# -*- mode: python -*-

cdef extern from "../src/dbg.hh":
    bint upslope_area(double *x, int Mx, double *y, int My, double *z, int *mask, bint output)
    bint accumulated_flow(double *x, int Mx, double *y, int My, double *z, double *mask, int n_samples)
    bint initialize_mask(int Mx, int My, double *thickness, int *mask)


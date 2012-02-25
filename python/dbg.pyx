# -*- mode: python -*-
#cython: boundscheck=False
#cython: wraparound=False
# Comments above are special. Please do not remove.

cimport numpy as np
import numpy as np

cimport dbg_c

ctypedef np.float64_t double_t
ctypedef np.int32_t int_t

def initialize_mask(np.ndarray[dtype=double_t, ndim=2, mode="c"] thk):
    """
    Use ice thickness to initialize the mask.
    """
    cdef np.ndarray[dtype=int_t, ndim=2, mode="c"] mask
    cdef int Mx, My

    Mx = thk.shape[1]
    My = thk.shape[0]

    mask = np.zeros((thk.shape[0], thk.shape[1]), dtype=np.int32)

    dbg_c.initialize_mask(Mx, My, <double*>thk.data, <int*> mask.data)

    return mask

cdef check_dimensions(x, y, z, mask):
    # z and mask are typed, so z.shape is not a Python object.
    # This means that we have to compare z.shape[0,1] to mask.shape[0,1] 'by hand'.
    if not (z.shape[0] == mask.shape[0] and z.shape[1] == mask.shape[1]):
        raise ValueError("arguments z and mask have to have the same shape: got (%d,%d) and (%d,%d)" %
                         (z.shape[0], z.shape[1], mask.shape[0], mask.shape[1]))

    if y.size != z.shape[0]:
        raise ValueError("the size of y has to match the number of rows in z")

    if x.size != z.shape[1]:
        raise ValueError("the size of x has to match the number of columns in z")

def upslope_area(np.ndarray[dtype=double_t, ndim=1] x,
                 np.ndarray[dtype=double_t, ndim=1] y,
                 np.ndarray[dtype=double_t, ndim=2, mode="c"] z,
                 np.ndarray[dtype=int_t, ndim=2, mode="c"] mask,
                 copy = False, print_output = False):
    """
    Computes the upslope area of points marked in the mask argument.

    Ice-free cells should be marked with -1, icy cells to be processed with -2,
    cells at termini (icy or not) with positive numbers, one per terminus.

    Try initialize_mask(thickness) if you don't know where termini are.

    arguments:
    - x, y: 1D arrays with coordinates
    - z: surface elevation, a 2D NumPy array
    - mask: mask, integers, a 2D NumPy array
    - copy: boolean; False if the mask is to be modified in place
    """
    cdef np.ndarray[dtype=int_t, ndim=2, mode="c"] output

    check_dimensions(x, y, z, mask)

    if copy:
        output = mask.copy()
    else:
        output = mask

    dbg_c.upslope_area(<double*>x.data, x.size, <double*>y.data, y.size,
                       <double*>z.data, <int*>output.data,
                       print_output)

    return output

def accumulated_flow(np.ndarray[dtype=double_t, ndim=1] x,
                     np.ndarray[dtype=double_t, ndim=1] y,
                     np.ndarray[dtype=double_t, ndim=2, mode="c"] z,
                     np.ndarray[dtype=double_t, ndim=2, mode="c"] mask,
                     copy = False, n_samples = 1):
    """
    Computes the accumulated flow map.

    Ice-free cells should be marked with -1, icy cells with 0.

    arguments:
    - x, y: 1D arrays with coordinates
    - z: surface elevation, a 2D NumPy array
    - mask: mask, integers, a 2D NumPy array
    - copy: boolean; False if the mask is to be modified in place
    """
    cdef np.ndarray[dtype=double_t, ndim=2, mode="c"] output

    check_dimensions(x, y, z, mask)

    if copy:
        output = mask.copy()
    else:
        output = mask

    dbg_c.accumulated_flow(<double*>x.data, x.size, <double*>y.data, y.size,
                           <double*>z.data, <double*>output.data, n_samples)

    return output

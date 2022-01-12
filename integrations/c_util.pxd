#cython: language_level=3
"""
definitions file for cython
"""
cimport numpy as cnp
ctypedef cnp.float64_t DTYPE_t

ctypedef struct vec:
    DTYPE_t x
    DTYPE_t y
    DTYPE_t z

cdef vec cvec(DTYPE_t, DTYPE_t, DTYPE_t)
cdef vec vsum(vec, vec)
cdef vec diff(vec, vec)
cdef vec cross(vec, vec)
cdef vec mult(vec, DTYPE_t)
cdef DTYPE_t dot(vec, vec)
cpdef DTYPE_t r(DTYPE_t, DTYPE_t, DTYPE_t)
cpdef DTYPE_t rp(DTYPE_t, DTYPE_t, DTYPE_t)
cpdef DTYPE_t rpspherical(DTYPE_t, DTYPE_t, DTYPE_t)
cpdef DTYPE_t rp2spherical(DTYPE_t, DTYPE_t, DTYPE_t)
cpdef DTYPE_t r2(DTYPE_t, DTYPE_t, DTYPE_t)
cpdef DTYPE_t rp2(DTYPE_t, DTYPE_t, DTYPE_t)
cpdef DTYPE_t magnitude(DTYPE_t, DTYPE_t, DTYPE_t)
cdef DTYPE_t magnitudecv(vec)
cpdef DTYPE_t magnitudev(cnp.ndarray[DTYPE_t, ndim=1])
cpdef DTYPE_t psi_ellipsoid(DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t)
cpdef DTYPE_t magnabla_psi_ellipsoid(DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t)
cpdef DTYPE_t dellipsoiddr(DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t)


#cython: language_level=3
"""
definitions file for cython
"""
cimport c_util as cu
ctypedef cu.vec cv
ctypedef double DTYPE_t

cpdef DTYPE_t roche(DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t)
cdef DTYPE_t rochecv(cv, DTYPE_t)
cdef cv rochegradientc(DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t)
cdef cv rochegradientcv(cv, DTYPE_t)
cpdef DTYPE_t rochespherical(DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t)
cpdef DTYPE_t drochedrspherical(DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t)
cpdef DTYPE_t drochedthetaspherical(DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t,  DTYPE_t)
cpdef DTYPE_t gravityspherical(DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t)
cpdef DTYPE_t invgravityspherical(DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t)
cpdef DTYPE_t rocheprimespherical(DTYPE_t, DTYPE_t, DTYPE_t,  DTYPE_t)
cpdef DTYPE_t drocheprimedrspherical(DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t)
cpdef DTYPE_t drocheprimedthetaspherical(DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t)
cpdef DTYPE_t gravityprimespherical(DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t)
cpdef DTYPE_t invgravityprimespherical(DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t, DTYPE_t)
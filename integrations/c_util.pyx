#cython: language_level=3
"""
various vector functions that are also cythonified for speed, these should be self explanatory
"""
import numpy as np
cimport numpy as cnp  # compile time info from numpy

from libc.math cimport sqrt, atan2

cdef vec cvec(DTYPE_t x, DTYPE_t y, DTYPE_t z):
    cdef vec res
    res.x = x
    res.y = y
    res.z = z
    return res

cdef vec vsum(vec v1, vec v2):
    cdef vec res
    res.x = v1.x + v2.x
    res.y = v1.y + v2.y
    res.z = v1.z + v2.z
    return res

cdef vec diff(vec v1, vec v2):
    cdef vec res
    res.x = v1.x - v2.x
    res.y = v1.y - v2.y
    res.z = v1.z - v2.z
    return res

cdef vec cross(vec v1, vec v2):
    cdef vec res
    res.x = v1.y * v2.z - v1.z * v2.y
    res.y = v1.z * v2.x - v1.x * v2.z
    res.z = v1.x * v2.y - v1.y * v2.x
    return res

cdef vec mult(vec v1, DTYPE_t fac):
    cdef vec res
    res.x = v1.x * fac
    res.y = v1.y * fac
    res.z = v1.z * fac
    return res

cdef DTYPE_t dot(vec v1, vec v2):
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z

cpdef DTYPE_t r(DTYPE_t x, DTYPE_t y, DTYPE_t z):
    return sqrt(r2(x, y, z))

cpdef DTYPE_t rp(DTYPE_t x, DTYPE_t y, DTYPE_t z):
    return r(1 - x, y, z)

cpdef DTYPE_t rpspherical(DTYPE_t ar, DTYPE_t cosp, DTYPE_t sint):
    return sqrt(rp2spherical(ar, cosp, sint))

cpdef DTYPE_t rp2spherical(DTYPE_t ar, DTYPE_t cosp, DTYPE_t sint):
    return 1 - 2 * ar * cosp * sint + ar ** 2.

cpdef DTYPE_t r2(DTYPE_t x, DTYPE_t y, DTYPE_t z):
    return x ** 2. + y ** 2. + z ** 2.

cpdef DTYPE_t rp2(DTYPE_t x, DTYPE_t y, DTYPE_t z):
    return r2(1 - x, y, z)

cpdef DTYPE_t magnitude(DTYPE_t x, DTYPE_t y, DTYPE_t z):
    return r(x, y, z)

cdef DTYPE_t magnitudecv(vec v1):
    return r(v1.x, v1.y, v1.z)

cpdef DTYPE_t magnitudev(cnp.ndarray[DTYPE_t, ndim=1] vec):
    return r(vec[0], vec[1], vec[2])

cpdef cnp.ndarray[DTYPE_t, ndim=2] cart2spher(cnp.ndarray[DTYPE_t, ndim=2] xyz):
    """
    convert list of cartesian points (x, y, z) to points in spherical coordinates (r, theta, phi)
    theta = polar angle
    :param xyz: 2D list of points
    :return: list of points in spherical coordinates (r, theta, phi)
    """
    cdef cnp.ndarray[DTYPE_t, ndim=2] pts = np.empty((xyz.shape[0], 3))
    cdef DTYPE_t[:, :] view = pts, vxyz = xyz
    cdef DTYPE_t XsqPlusYsq
    cdef Py_ssize_t i
    for i in range(xyz.shape[0]):
        XsqPlusYsq = vxyz[i, 0] ** 2 + vxyz[i, 1] ** 2
        view[i, 0] = sqrt(XsqPlusYsq + vxyz[i, 2] ** 2)  # set the memory view (supposedly faster), ...
        view[i, 1] = atan2(sqrt(XsqPlusYsq), vxyz[i, 2])
        view[i, 2] = atan2(vxyz[i, 1], vxyz[i, 0])
    return pts  # ... but return the actual pts array (I think of this like pointers)

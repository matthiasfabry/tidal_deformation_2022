#cython: language_level=3
"""
Functions in the Roche geometry
"""

cimport c_util as cu
from libc.math cimport sqrt

cdef cv rochegradientc(DTYPE_t x, DTYPE_t y, DTYPE_t z, DTYPE_t q):
    """
    return the gradient of Psi
    :param x: cartesian x coordinate
    :param y: cartesian y coordinate
    :param z: cartesian z coordinate
    :param q: mass ratio
    :return: cartesian nabla(Psi)
    """
    cdef DTYPE_t arm3 = cu.r2(x, y, z) ** -1.5
    cdef DTYPE_t arprimem3 = cu.rp2(x, y, z) ** -1.5
    cdef cv ret
    if q == 0:
        ret = cu.cvec(x * arm3 - x, y * (arm3 - 1), z * arm3)
    else:
        ret = cu.cvec(x * arm3 - q * (1 - x) * arprimem3 + q - (q + 1) * x,
                      y * (arm3 + q * arprimem3 - (q + 1)),
                      z * (arm3 + q * arprimem3))
    return ret

cdef cv rochegradientcv(cv point, DTYPE_t q):
    """
    return the gradient of Psi
    :param point: cartesian vector
    :param q: mass ratio
    :return: cartesian nabla(Psi)
    """
    return rochegradientc(point.x, point.y, point.z, q)

cpdef DTYPE_t roche(DTYPE_t x, DTYPE_t y, DTYPE_t z, DTYPE_t q):
    """
    return the value of potential Psi
    :param x: cartesian x coordinate
    :param y: cartesian y coordinate
    :param z: cartesian z coordinate
    :param q: mass ratio
    :return: Psi(x, y, z)
    """
    if q == 0:
        return - cu.r(x, y, z) ** -1. - 0.5 * (x ** 2. + y ** 2.)
    return - cu.r(x, y, z) ** -1. - q * (cu.rp(x, y, z) ** -1. - x) - 0.5 * (1 + q) * (x ** 2. + y ** 2.)

cdef DTYPE_t rochecv(cv point, DTYPE_t q):
    """
    return the value of potential Psi
    :param point: cartesian vector
    :param q: mass ratio
    :return: Psi((x, y, z))
    """
    return roche(point.x, point.y, point.z, q)

cpdef DTYPE_t rochespherical(DTYPE_t ar, DTYPE_t sint, DTYPE_t cosp, DTYPE_t q):
    """
    return the value of potential Psi
    :param ar: spherical coordinate r
    :param sint: sine of spherical polar angle theta
    :param cosp: cosine of spherical azimuthal angle phi
    :param q: mass ratio
    :return: Psi(r, theta, phi)
    """
    return -ar ** -1. - q * (cu.rpspherical(ar, cosp, sint) ** -1. - ar * cosp * sint) - 0.5 * (
            1 + q) * ar ** 2. * sint ** 2.

cpdef DTYPE_t drochedrspherical(DTYPE_t ar, DTYPE_t sint, DTYPE_t cosp, DTYPE_t q):
    """
    return the value of radial component of nabla(Psi)
    :param ar: spherical coordinate r
    :param sint: sine of spherical polar angle theta
    :param cosp: cosine of spherical azimuthal angle phi
    :param q: mass ratio
    :return: dPsi(r, theta, phi)/dr
    """
    if q == 0:
        return ar ** -2. - ar * sint ** 2.
    return ar ** -2. + q * (cu.rp2spherical(ar, cosp, sint) ** -1.5 *
                            (ar - cosp * sint) + cosp * sint) - (1 + q) * ar * sint ** 2.

cpdef DTYPE_t drochedthetaspherical(DTYPE_t ar, DTYPE_t cost, DTYPE_t sint, DTYPE_t cosp, DTYPE_t q):
    """
    return the value of polar component of nabla(Psi)
    :param ar: spherical coordinate r
    :param cost: cosine of spherical polar angle theta
    :param sint: sine of spherical polar angle theta
    :param cosp: cosine of spherical azimuthal angle phi
    :param q: mass ratio
    :return: dPsi(r, theta, phi)/dtheta
    """
    if q == 0:
        return - ar ** 2. * sint * cost
    return -q * ar * cosp * cost * (cu.rp2spherical(ar, cosp, sint) ** -1.5 - 1) - (1 + q) * ar ** 2. * sint * cost

cpdef DTYPE_t gravityspherical(DTYPE_t ar, DTYPE_t cost, DTYPE_t sint, DTYPE_t cosp, DTYPE_t sinp, DTYPE_t q):
    """
    return the value of effective gravity |nabla(Psi)|
    :param ar: spherical coordinate r
    :param cost: cosine of spherical polar angle theta
    :param sint: sine of spherical polar angle theta
    :param cosp: cosine of spherical azimuthal angle phi
    :param sinp: sine of spherical azimuthal angle phi
    :param q: mass ratio
    :return: |nabla(Psi(r, theta, phi))|
    """
    if q == 0:
        return sqrt(drochedrspherical(ar, sint, cosp, q) ** 2. +
                    ar ** -2. * drochedthetaspherical(ar, cost, sint, cosp, q) ** 2.)
    return sqrt(drochedrspherical(ar, sint, cosp, q) ** 2. +
                ar ** -2. * drochedthetaspherical(ar, cost, sint, cosp, q) ** 2. +
                (q * sinp * (1 - cu.rp2spherical(ar, cosp, sint) ** -1.5)) ** 2.)

cpdef DTYPE_t invgravityspherical(DTYPE_t ar, DTYPE_t cost, DTYPE_t sint, DTYPE_t cosp, DTYPE_t sinp, DTYPE_t q):
    """
    return the value of inverse effective gravity 1/|nabla(Psi)|
    :param ar: spherical coordinate r
    :param cost: cosine of spherical polar angle theta
    :param sint: sine of spherical polar angle theta
    :param cosp: cosine of spherical azimuthal angle phi
    :param sinp: sine of spherical azimuthal angle phi
    :param q: mass ratio
    :return: 1/|nabla(Psi(r, theta, phi))|
    """
    return gravityspherical(ar, cost, sint, cosp, sinp, q) ** -1.

# these are functions for the 'displaced' potential x->1-x if you want properties of the secondary, not actively used
cpdef DTYPE_t rocheprimespherical(DTYPE_t ar, DTYPE_t sint, DTYPE_t cosp, DTYPE_t q):
    return -1./sqrt(1 + 2 * ar * cosp * sint + ar ** 2.) - q / ar - ar * cosp * sint - 0.5 * (
            1 + q) * ar ** 2. * sint ** 2. - (1 - q) / 2

cpdef DTYPE_t drocheprimedrspherical(DTYPE_t ar, DTYPE_t sint, DTYPE_t cosp, DTYPE_t q):
    return (ar + cosp * sint) * (1 + 2 * ar * cosp * sint + ar ** 2.) ** -1.5 + q * ar ** -2. - cosp * sint - (
            1 + q) * ar * sint ** 2.

cpdef DTYPE_t drocheprimedthetaspherical(DTYPE_t ar, DTYPE_t cost, DTYPE_t sint, DTYPE_t cosp, DTYPE_t q):
    return ar * cosp * cost * (1 + 2 * ar * cosp * sint * ar ** 2.) ** -1.5 - ar * cosp * cost - (
            1 + q) * ar ** 2. * sint * cost

cpdef DTYPE_t gravityprimespherical(DTYPE_t ar, DTYPE_t cost, DTYPE_t sint, DTYPE_t cosp, DTYPE_t sinp, DTYPE_t q):
    return sqrt(drocheprimedrspherical(ar, sint, cosp, q) ** 2. +
                ar ** -2. * drocheprimedthetaspherical(ar, cost, sint, cosp, q) ** 2. +
                (sinp * (1 - (1 + 2 * ar * cosp * sint + ar ** 2.) ** -1.5)) ** 2.)

cpdef DTYPE_t invgravityprimespherical(DTYPE_t ar, DTYPE_t cost, DTYPE_t sint, DTYPE_t cosp, DTYPE_t sinp, DTYPE_t q):
    return gravityprimespherical(ar, cost, sint, cosp, sinp, q) ** -1.

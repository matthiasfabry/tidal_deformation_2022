"""
pure python version of c_rocheFunctions.pyx. Try the cython version first, this is backup.
"""
import numpy as np
import py_utils as pu


def roche(x, y, z, q):
    if q == 0:
        return - pu.r(x, y, z) ** -1. - 0.5 * (x ** 2. + y ** 2.)
    return - pu.r(x, y, z) ** -1. - q * (pu.rp(x, y, z) ** -1. - x) - 0.5 * (1 + q) * (x ** 2. + y ** 2.)


def rochev(point, q):
    return roche(point[0], point[1], point[2], q)


def rochespherical(ar, sint, cosp, q):
    return -ar ** -1. - q * (pu.rpspherical(ar, cosp, sint) ** -1. - ar * cosp * sint) - 0.5 * (
            1 + q) * ar ** 2. * sint ** 2.


def drochedrspherical(ar, sint, cosp, q):
    if q == 0:
        return ar ** -2. - ar * sint ** 2.
    return ar ** -2. + q * (pu.rp2spherical(ar, cosp, sint) ** -1.5 *
                            (ar - cosp * sint) + cosp * sint) - (1 + q) * ar * sint ** 2.


def eggleton_rl_radius(q):
    return 0.49 / q ** (2. / 3) / (0.6 / q ** (2. / 3) + np.log1p(1. / q ** (1. / 3)))

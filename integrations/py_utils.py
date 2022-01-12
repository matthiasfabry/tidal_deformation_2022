"""
pure python version of c_util.pyx. Try the cython version first. this is backup.
"""

import numpy as np
import c_util as cu


def r(x, y, z):
    return np.sqrt(r2(x, y, z))


def rp(x, y, z):
    return r(1 - x, y, z)


def r2(x, y, z):
    return x ** 2. + y ** 2. + z ** 2.


def rp2(x, y, z):
    return r2(1 - x, y, z)


def rpspherical(ar, cosp, sint):
    return np.sqrt(rp2spherical(ar, cosp, sint))


def rp2spherical(ar, cosp, sint):
    return 1 - 2 * ar * cosp * sint + ar ** 2.


def corr_finder_from_normal(i, j, points, rr, phinum, thetanum):
    if i == 0:
        neighbor1 = np.asarray(points[j * phinum + i])
    else:
        neighbor1 = np.asarray(points[j * phinum + i - 1])
    if i == phinum - 1:
        neighbor2 = np.asarray(points[j * phinum + i])
    else:
        neighbor2 = np.asarray(points[j * phinum + i + 1])
    if j == 0:
        neighbor3 = np.asarray(points[i])
    else:
        neighbor3 = np.asarray(points[(j - 1) * phinum + i])
    if j == thetanum - 1:
        neighbor4 = np.asarray(points[j * phinum + i])
    else:
        neighbor4 = np.asarray(points[(j + 1) * phinum + i])
    v1 = neighbor4 - neighbor3
    v2 = neighbor2 - neighbor1
    cross = np.cross(v1, v2)
    corr = np.dot(points[j * phinum + i] / rr, cross / cu.magnitude(cross))
    return np.abs(corr)

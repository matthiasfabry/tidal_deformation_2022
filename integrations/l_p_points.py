"""
module for computing the location of the lagrangian points and other extrema of a roche shell
these are all root finding algorithms of well known equations (see eg Kopal 1959)
"""

import numpy as np
import scipy.optimize as spopt
import c_util as c


def l1(q):
    def toroot(x):
        return (1 + q) * x ** 5 - (2 + 3 * q) * x ** 4 + (1 + 3 * q) * x ** 3 - x ** 2 + 2 * x - 1

    return spopt.root_scalar(toroot, bracket=[0, 1]).root


def lfar(q):
    def toroot(x):
        return (1 + q) * x ** 5 - (2 + 3 * q) * x ** 4 + (1 + 3 * q) * x ** 3 - (1 + 2 * q) * x ** 2 + 2 * x - 1

    return spopt.root_scalar(toroot, bracket=[1, 2]).root


def lout(q):
    def toroot(x):
        return (1 + q) * x ** 5 - (2 + 3 * q) * x ** 4 + (1 + 3 * q) * x ** 3 + x ** 2 - 2 * x + 1

    return spopt.root_scalar(toroot, bracket=[-1, 0]).root


def p3(q, pot):
    def toroot(x):
        return (1 - pot * x + q * x ** 2 - (1 + q) / 2 * x ** 3) * (1 - x) - q * x

    return spopt.root_scalar(toroot, bracket=[0, lout(q)]).root


def p5(q, pot):
    def toroot(x):
        ar = c.r(x[0], x[1], 0)
        arp = c.r(1 - x[0], x[1], 0)
        return [(pot - q * x[0] + (1 + q) / 2 * (x[0] ** 2 + x[1] ** 2)) * ar * arp + arp + q * ar,
                (q - (q + 1) * x[0]) * ar ** 3 * arp ** 3 + x[0] * arp ** 3 - q * (1 - x[0]) * ar ** 3]

    return spopt.root(toroot, np.array([0, 0])).x


def p7(q, p5x, pot):
    def toroot(z):
        ar = c.r(p5x, 0, z)
        arp = c.r(1 - p5x, 0, z)
        return (pot - q * p5x + (1 + q) / 2 * p5x ** 2) * ar * arp + arp + q * ar

    return spopt.root_scalar(toroot, bracket=[0, 1]).root

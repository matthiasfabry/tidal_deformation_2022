"""
methods relating to the single rotating star potential.
"""
import numpy as np
import scipy.optimize as spopt

pi = np.pi


def omega_single_star(r_eq, q):
    """
    computes the fractional rotation rate of the equivalent single star given the equivalent radius and mass ratio
    See also appendix B
    :param r_eq: equivalent radius of the roche equipotential shell
    :param q: mass ratio M2/M1
    :return: omega, the fractional rotation rate
    """
    if hasattr(r_eq, "__iter__"):  # loop if there are multiple
        ret = np.empty(len(r_eq))
        for i in range(len(r_eq)):
            def to_solve(oom):
                return oom ** (2 / 3) - r_eq[i] * polyequarad(oom) * np.cbrt(1 + q)
            
            try:
                ret[i] = spopt.root_scalar(to_solve, bracket=[0, 1]).root
            except ValueError as e:
                print(e)
                ret[i] = np.nan
        return ret
    else:
        def to_solve(oom):
            try:
                return oom ** (2 / 3) - r_eq * polyequarad(oom) * np.cbrt(1 + q)
            except ValueError as e:
                print(e)
                return np.nan
        
        return spopt.root_scalar(to_solve, bracket=[0, 1]).root


# the following are all the polynomial fits of the single star properties, as function of fractional rotation rate
# see appendix A
def polyequarad(om):  # equatorial radius
    return 1 + om ** 2 / 6 - 0.005124 * om ** 4 + 0.06562 * om ** 6


def polysurftimesgrav(om):  # surface area times average gravity, scales with GM
    return 4 * pi * (1 - 2 / 3 * om ** 2 + 0.2989 * om ** 4 + 0.006978 * om ** 6)


def polyvol(om):  # equipotential volume
    return 4 * pi / 3 * polyequarad(om) ** 3 * (1 - 0.5 * om ** 2 + 0.1180 * om ** 4 - 0.07688 * om ** 6)


def polyvol_pablo(om):  # equipotential volume from Paxton 2019
    return 4 * pi / 3 * polyequarad(om) ** 3 * (1 - 0.5 * om ** 2 + 0.1149 * om ** 4 - 0.07376 * om ** 6)


def polysurf(om):  # surface area
    return 4 * pi * polyequarad(om) ** 2 * (1 - om ** 2 / 3 + 0.08696 * om ** 4 - 0.05079 * om ** 6)


def polysurf_pablo(om):  # surface area from Paxton 2019
    return 4 * pi * polyequarad(om) ** 2 * (1 - om ** 2 / 3 + 0.08525 * om ** 4 - 0.04908 * om ** 6)


def polyequivrad(om):  # equivalent radius, basically inverse from poly_equarad
    return 1 - om ** 2 / 6 + 0.02025 * om ** 4 - 0.03870 * om ** 6


def polysurftimesinvgrav(om):  # area times inverse average gravity, scales with 1/GM
    return 4 * pi * polyequarad(om) ** 4 * A(om)


def A(om):  # auxiliary function
    return 1 + 0.3293 * om ** 4 - 0.4926 * om ** 6 - 0.5600 * np.log1p(-om ** 5.626)


def B(om):  # auxiliary function
    return 1 + 0.2 * om ** 2 + 0.4140 * om ** 4 - 0.8650 * om ** 6 - 1.5 * 0.5600 * np.log1p(-om ** 5.626)


def polyfP(om):  # correction factor f_P
    return (1 - 2 / 3 * om ** 2 + 0.2133 * om ** 4 - 0.1068 * om ** 6) / A(om)


def polyfT(om):  # correction factor f_T
    return (1 - 0.0796 * om ** 4 - 0.2322 * om ** 6) / A(om)


def polyfTovertP(om):  # ratio f_T/f_P (basically to eliminate the pole in the log from A(om)
    return (1 - 0.07955 * om ** 4 - 0.2322 * om ** 6) / \
           (1 - 2 / 3 * om ** 2 + 0.2133 * om ** 4 - 0.1068 * om ** 6)


def polyirot(om):  # moment of inertia
    if om == 1:
        return 2 / 3
    return 2 / 3 * polyequarad(om) ** 2 * B(om) / A(om)


def C(om):  # auxiliary function
    return 1 + 17 / 60 * om ** 2 + 0.3841 * om ** 4 - 0.8202 * om ** 6 - 0.9168 * np.log1p(-om ** 5.558)


def polyjrot(om):  # specific angular momentum, scaled with 1/sqrt(GMr)
    if om == 1:
        return np.sqrt(1 / 0.81488)
    return 2 / 3 * om * C(om) / A(om)

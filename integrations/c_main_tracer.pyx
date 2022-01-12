#cython: language_level=3
"""
integration methods in the Roche geometry
"""

cimport numpy as cnp  # import compile time info for numpy
import numpy as np
cimport c_rocheFunctions as cr
from libc.math cimport sin, cos, cbrt, pi, isnan, NAN, fabs  # noqa

ctypedef double DTYPE_t


def c_roche_tracer(DTYPE_t q, DTYPE_t targetpot, DTYPE_t xl1,
                   cnp.ndarray[DTYPE_t, ndim=1] thetas,
                   cnp.ndarray[DTYPE_t, ndim=1] phis,
                   DTYPE_t dr_stop = 1e-12,
                   cnp.ndarray[DTYPE_t, ndim=2] l1_inter=None,
                   cnp.ndarray[DTYPE_t, ndim=2] lout_inter=None,
                   int save_points=0
                   ):
    """
    finds an equipotential shell and integrates various quantities on them
    :param q: mass ratio m_other/m_this
    :param targetpot: equipotential value Psi you want to find and integrate on
    :param xl1: x coordinate of location of first lagrangian point (used as first step in the root finding,
        can be some other number >0)
    :param phis: samples of phi dividing (0, pi)
    :param thetas: samples of theta dividing (0, pi/2), equally spaced in cos(theta)
    :param dr_stop: smallest dr the bisection algorithm
    :param l1_inter: rs evaluated on the l1 splitting surface on the same grid of (phis, thetas), needed if Psi>Psi_L1
    :param lout_inter: rs evaluated on the Lout splitting surface on the same grid of (phis, thetas)
        needed if Psi>Psi_Lout
    :param save_points: do you want to save the hull points
    :return: list of points on the equipotential (if save_points != 0),
        list of quantities on that equipotential: volume, area,
        equivalent radius, average gravity, average inverse gravity, moment of inertia, f_p, f_t
    """
    # this is C-like, we need to declare stuff to benefit from cython.
    cdef DTYPE_t cosp, sinp, cost, sint, r, \
        vol = 0., area = 0., grav = 0., invgrav = 0., r_eq = 0., f_p = 0., f_t = 0., \
        moi = 0., inv_grav_here, g_here, corr, element, surf_here, r1cross = NAN, r3cross = NAN, dphi, dcostheta, dpsidr
    cdef Py_ssize_t i, j, phinum = len(phis), thetanum = len(thetas)  # these are array indices/lengths
    cdef int crosses = 0
    cdef int *crosses_ref = &crosses  # pointer to the location of crosses
    cdef list points1 = list()
    
    dphi = pi / phinum
    dcostheta = 1.0 / thetanum
    
    for j in range(thetanum): # find hull first
        cost = cos(thetas[j])
        sint = sin(thetas[j])
        for i in range(phinum):
            cosp = cos(phis[i])
            sinp = sin(phis[i])
            
            r1cross = l1_inter[i, j]
            r3cross = lout_inter[i, j]
            
            r = r_finder(targetpot, xl1, dr_stop, sint, cosp,
                         cost, sinp, q, r3cross, r1cross, crosses_ref)  # find r(theta, psi)
            
            if save_points:
                points1.append(np.array([r*sint*cosp, r*sint*sinp, r*cost]))
                
            element = r ** 2 * dphi * dcostheta  # base integrand
            vol += r * element  # volume is r cubed
            
            if not crosses:  # if on a splitting surface, it doesn't contribute to either of the following:
                g_here = cr.gravityspherical(r, cost, sint, cosp, sinp, q)
                dpsidr = cr.drochedrspherical(r, sint, cosp, q)
                corr = dpsidr / g_here  # == nhat.rhat
                surf_here = element / corr  # projection correction for surface
                inv_grav_here = surf_here / g_here
                invgrav += inv_grav_here
                grav += surf_here * g_here
                area += surf_here
                moi += inv_grav_here * r ** 2 * sint ** 2
                
    area *= 4  # normalize with correct factors
    grav *= 4 / area
    invgrav *= 4 / area
    vol *= 4 / 3
    r_eq = cbrt(3 * vol / (4 * pi))
    moi *= 4 / (area * invgrav)

    f_p = 4 * pi * r_eq ** 4 / (area * invgrav)  # calculate fp, ft
    f_t = (4 * pi * r_eq ** 2 / area) ** 2 / (grav * invgrav)
    return points1, [vol, area, grav, invgrav, r_eq, moi, f_p, f_t]


cdef DTYPE_t r_finder(DTYPE_t targetpot, DTYPE_t dr_start, DTYPE_t dr_stop,
                       DTYPE_t sint, DTYPE_t cosp, DTYPE_t cost, DTYPE_t sinp, DTYPE_t q,
                       DTYPE_t loutcross, DTYPE_t l1cross, int *crossesref, DTYPE_t r_start=0.):
    """
    root finder that solves Psi = Psi(r, theta, phi) to r for a given Psi, theta and phi.
    :param targetpot: equipotential Psi considered
    :param dr_start: first step that will be taken from r=0
    :param dr_stop: smallest step the bisection is allowed to take
    :param sint: sin(theta)
    :param cosp: cos(phi)
    :param cost: cos(theta)
    :param sinp: sin(phi)
    :param q: mass ratio m_other/m_this
    :param loutcross: r(theta, phi) at which the Lout splitting surface is crossed
        (can be NAN if it doesn't for this (theta, phi))
    :param l1cross: r(theta, phi) at which the L1 splitting surface is crossed
        (can be NAN if it doesn't for this (theta, phi))
    :param crossesref: the address where the int crosses is stored, the integrator above needs to know whether we
        crossed a splitting surface
    :param r_start: start the bisection on another r (optional, default=0)
    :return: the radius r at which Psi=Psi(r, theta, phi) with at most dr = dr_stop
        OR the radius at the splitting surface if we crossed one
    """
    cdef DTYPE_t dr = dr_start, r = r_start, newpot, newr
    cdef int reduce
    
    while True:
       
        newr = r + dr # step radius
        reduce = 0
        crossesref[0] = 0  # we did not cross a splitting surface yet
        
        newpot = cr.rochespherical(newr, sint, cosp, q)
            
        if not isnan(l1cross):  # check if you crossed L1 or Lout surface, if applicable
            if newr > l1cross:
                reduce = 1
                crossesref[0] = 1
        elif not isnan(loutcross):
            if newr > loutcross:
                reduce = 1
                crossesref[0] = 1

        if reduce == 0 and newpot > targetpot:  # check if you overshot potential
            reduce = 1

        if reduce == 0 and cr.drochedrspherical(newr, sint, cosp, q) < 0:  # check if derivative is still positive
            reduce = 1
            if not (isnan(l1cross) and isnan(loutcross)):
                crossesref[0] = 1
        
        if reduce == 1:  # do reducing procedure
            if dr <= dr_stop:
                break  # stop the loop if dr is too small, ie we have the r we need;
            dr = tighten_dr(dr)  # tighten stepping
            continue  # go to the top; don't update r if reducing was triggered
        
        dr = caution_dr(dr, dr_stop, newpot, targetpot) # check cautioning condition
        r = newr # update r
    return r

cdef DTYPE_t caution_dr(DTYPE_t dr, DTYPE_t dr_stop, DTYPE_t newpot, DTYPE_t targetpot):
    # tighten stepping if you're getting close to targetpot (if dr is not already small)
    if dr > 4./3 * dr_stop and fabs(targetpot - newpot) < dr:
        dr = max(0.75 * dr, fabs(targetpot - newpot))
    return dr

cdef DTYPE_t tighten_dr(DTYPE_t dr):
    dr *= 0.5  # tighten stepping size if we overshot
    return dr

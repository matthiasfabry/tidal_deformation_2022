#cython: language_level=3
"""
Tracer that finds the splitting surface between stars in the Roche geometry
"""

cimport numpy as cnp  # import compile type info for numpy
import numpy as np
cimport c_rocheFunctions as cr
cimport c_util as cu
from libc.math cimport sin, cos, pi, log10

ctypedef cu.vec cv
ctypedef cnp.float64_t DTYPE_t

def c_gradient_tracer(DTYPE_t x_start, DTYPE_t stop_pot, DTYPE_t q):
    """
    follows the gradient out from x_start to construct a splitting surface s.t. DelPsi is perpendicular to dS
    :param x_start: starting x coordinate (usually x_l1 or x_lout)
    :param stop_pot: potential Psi at which to stop tracing
    :param q: mass ratio m_other/m_this
    :return: dict of points with key 'bulk' specifying points that fall on the splitting surface and
        key 'edge' specifying points on the outermost (theta, phi) line.
    """
    cdef DTYPE_t start_pot = cr.roche(x_start, 0, 0, q)
    if start_pot >= stop_pot: # quick check if we have anything to do at all
        return np.array([[x_start, 0, 0]])
    
    cdef cv p_start
    cdef dict poss
    cdef int raynum = 1000
    
    p_start = cu.cvec(x_start, 0, 0)
    if x_start < 0:  # Lout
        if q < 0.05:  # values that determine which algorithm is used to trace a ray; these are somewhat arbitrary, YMMV
            poss = split_tracer_from_equator(q, p_start, stop_pot, raynum)
        elif q < 10:
            poss = split_tracer_decreasing_angle(q, p_start, stop_pot, raynum)
        else:
            poss = split_tracer_grid(q, p_start, stop_pot, raynum)
    else:  # L1
        if q < 0.05:  # same here
            poss = split_tracer_decreasing_angle(q, p_start, stop_pot, raynum)
        else:
            poss = split_tracer_grid(q, p_start, stop_pot, raynum)
    return poss

cdef dict split_tracer_grid(DTYPE_t q, cv p_start, DTYPE_t stop_pot, int raynum):
    """
    traces the points on the splitting surface in a (pseudo) regular grid of polar angles
    :param q: mass ratio m_other/m_this
    :param p_start: point at which to start
    :param stop_pot: potential Psi at which to stop
    :param raynum: how many rays (=angles thetas) need to be sampled
    :return: dict of points on the splitting surface
    """
    cdef dict poss = {'bulk':[], 'edge':[]}
    cdef list rayhere
    cdef cnp.ndarray[DTYPE_t, ndim=1] thetas = pi / 2 * np.linspace(0, 1, raynum) ** 4  # you can adapt this grid
    cdef Py_ssize_t i, raynumsize = raynum
    
    poss['bulk'].append(np.array([p_start.x, p_start.y, p_start.z]))  # starting point is part of the surface
    for i in range(raynum):  # compute single rays for all thetas
        rayhere = single_ray_gradient(q, p_start, thetas[i], stop_pot)
        poss['bulk'].extend(rayhere)
        poss['edge'].append(rayhere[-1])
    if pi/2 not in thetas:
        rayhere = single_ray_gradient(q, p_start, pi/2, stop_pot)  # make sure you have the meridian,
                                                                    # linspace can sometimes miss
        poss['bulk'].extend(rayhere)
        poss['edge'].append(rayhere[-1])
    return poss

cdef dict split_tracer_decreasing_angle(DTYPE_t q, cv p_start, DTYPE_t stop_pot, int raynum):
    """
    traces the points on the splitting surface by progressively decreasing the starting polar angle
    :param q: mass ratio m_other/m_this
    :param p_start: point at which to start
    :param stop_pot: potential Psi at which to stop
    :param raynum: how many rays (=angles thetas) need to be sampled
    :return: list of points on the splitting surface
    """
    cdef dict poss = {'bulk':[], 'edge':[]}, eqrun
    cdef list posshere
    cdef DTYPE_t angle = 0.5 * pi, dist, dist0, disthere  # start at theta=pi/2
    cdef cv equatorpos, meridianposs, lastpos
    cdef int i = 0
    
    poss['bulk'].append(np.array([p_start.x, p_start.y, p_start.z]))
    posshere = single_ray_gradient(q, p_start, 0, stop_pot)  # trace equator (theta=0) first
    poss['bulk'].extend(posshere)
    poss['edge'].append(posshere[-1])
    equatorpos = cu.cvec(posshere[-1][0], posshere[-1][1], posshere[-1][2])  # record furthest point on the traced equator ray
    posshere = single_ray_gradient(q, p_start, angle, stop_pot)  # trace meridian (theta=pi/2)
    poss['bulk'].extend(posshere)
    poss['edge'].append(posshere[-1])
    meridianposs = cu.cvec(posshere[-1][0], posshere[-1][1], posshere[-1][2])
    # distance from furthest point on meridian to extreme of equator
    dist0 = cu.magnitudecv(cu.diff(meridianposs, equatorpos))
    
    dist = dist0
    # print('enter angle loop')
    while dist > dist0 / raynum:
        if angle < 1e-14:  # stop criterion
            break
        else:
            angle *= 0.95  # step down angle
        posshere = single_ray_gradient(q, p_start, angle, stop_pot)  # trace a single ray
        lastpos = cu.cvec(posshere[-1][0], posshere[-1][1], posshere[-1][2])
        disthere = cu.magnitudecv(cu.diff(lastpos, equatorpos))
        # print(disthere)
        if dist - disthere > dist0 / raynum:  # only add if the last ray is significantly different
            # print('adding ray!')
            poss['bulk'].extend(posshere)
            poss['edge'].append(posshere[-1])
            dist = disthere
    # print('exit angle loop')
    if dist > dist0 / raynum:  # if theta=1e-14 still leaves too much space to equatorpos, do an equator run
        # print('finishing with equator run')
        eqrun = split_tracer_from_equator(q, p_start, stop_pot, raynum)
        poss['bulk'].extend(eqrun['bulk'])
        poss['edge'].extend(eqrun['edge'])
    return poss

cdef dict split_tracer_from_equator(DTYPE_t q, cv p_start, DTYPE_t stop_pot, int raynum):
    """
    traces the points on the splitting surface by starting from the equator, taking a tiny step up, and then tracing
    normally.
    :param q: mass ratio m_other/m_this
    :param p_start: point at which to start
    :param stop_pot: potential Psi at which to stop
    :param raynum: how many rays (=angles thetas) need to be sampled
    :return: list of points on the splitting surface
    """
    cdef dict poss = {'bulk':[], 'edge':[]}
    cdef list posshere, equator
    cdef Py_ssize_t i
    cdef cv start
    
    poss['bulk'].append(np.array([p_start.x, p_start.y, p_start.z]))
    equator = single_ray_gradient(q, p_start, 0, stop_pot)
    poss['bulk'].extend(equator)
    poss['edge'].append(equator[-1])
    posshere = single_ray_gradient(q, p_start, pi / 2, stop_pot)
    poss['bulk'].extend(posshere)
    poss['edge'].append(posshere[-1])
    for i in range(len(equator[:-1])):
        start = cu.cvec(equator[i][0], equator[i][1], equator[i][2])
        posshere = single_ray_gradient(q, start, pi / 2, stop_pot)
        poss['bulk'].extend(posshere)
        poss['edge'].append(posshere[-1])
    return poss

cdef list single_ray_gradient(DTYPE_t q, cv p_start, DTYPE_t theta, DTYPE_t stop_pot,
                              DTYPE_t start_step = 1e-14, DTYPE_t dr_save = 1e-3):
    """
    follows the gradient from a starting point until stop_pot is reached
    :param q: mass ratio m_other/m_this
    :param p_start: point where you start
    :param theta: polar angle at which you make the first tiny step
    :param stop_pot: when to stop tracing
    :param start_step: size of the first tiny step away from p_start
    :param dr_save: how far spaced must points be for them to be saved
    :return: list of points on the traced ray
    """
    cdef cv point, lastpoint, gradhere, lastgrad, ghat
    cdef DTYPE_t pot, dr_here=0., step_here, change=1e3
    cdef list poss = []
    cdef int i = 0
    
    dr_save = min(dr_save, 10**(log10(dr_save) - 1/7*log10(q)))  # extra care if q >> 1
    # first point: step away from the symmetry axis by start_step in direction theta
    point = cu.cvec(p_start.x, p_start.y + start_step * cos(theta), p_start.z + start_step * sin(theta))
    lastpoint = point
    gradhere = cr.rochegradientcv(point, q)
    ghat = cu.mult(gradhere, 1./cu.magnitudecv(gradhere))  # next direction to step in
    while True:
        pot = cr.rochecv(point, q)
        if pot > stop_pot:  # check to stop loop, and add last point
            poss.append(np.array([point.x, point.y, point.z]))
            break
        dr_here = cu.magnitudecv(cu.diff(point, lastpoint))
        if dr_here > dr_save:  # save a point to the list
            i += 1
            lastpoint = point
            poss.append(np.array([point.x, point.y, point.z]))
        step_here = 1e-6 / max(change, 1e-2)  # if gradient changes a lot, make step smaller than 1e-4
        point = cu.vsum(point, cu.mult(ghat, step_here))  # newpoint
        lastgrad = gradhere
        gradhere = cr.rochegradientcv(point, q)
        ghat = cu.mult(gradhere, 1./cu.magnitudecv(gradhere))  # next direction
        change = cu.magnitudecv(cu.diff(gradhere, lastgrad))
    return poss

"""
methods for interpolating the splitting surfaces on the (theta, phi) rays used in by the main roche integrations module
"""

import numpy as np
import scipy.interpolate as spinp
import scipy.spatial as scpa


def interpolate_split(points, edgepoints, thetas, phis, q, THETAS=None, PHIS=None, try_cubic=True):
    """
    given data from the splitting tracer, and chosen thetas and phis, compute interpolated values of the splitting
    surfaces on the (theta, phi) grid. First a cubic interpolation is tried. However for some q's cubic produces
    nan values on the splitting surface, probably because it's a non-convex surface although not sure. Then a linear
    interpolation is attempted, this is much more robust.
    :param points: points on the splitting surface (as computed by split_tracer), in spherical coordinates
    :param edgepoints: points that make up the outer edge of the splitting surface (also delivered by split_tracer),
        in spherical coordinates
    :param thetas: list of theta values to sample in the roche integration
    :param phis: list of phi values to sample in the roche integration
    :param q: mass ratio M2/M1
    :param THETAS: meshgrid of thetas, if None it's computed here
    :param PHIS: meshgrid of phis, if None it's computed here
    :param try_cubic: bool whether to try a cubic interpolation, default True.
    :return:
    """
    if THETAS is None or PHIS is None:
        THETAS, PHIS = np.meshgrid(thetas, phis)
    bulk = None
    edge = None
    if try_cubic:
        edge = interpolate_edge(edgepoints)
        bulk = interpolate_tracks(points, THETAS, PHIS)
        bulk = check_inter_and_clean(bulk, edge, thetas, phis)
        if bulk is None:
            print('{}: error detected in cubic interpolation, doing linear'.format(q))
            bulk = interpolate_tracks(points, THETAS, PHIS, method='linear')
            bulk = check_inter_and_clean(bulk, edge, thetas, phis)
    if not try_cubic or bulk is None:
        edge = interpolate_edge(edgepoints, kind='linear')
        bulk = interpolate_tracks(points, THETAS, PHIS, method='linear')
        bulk = check_inter_and_clean(bulk, edge, thetas, phis)
        if bulk is None:
            raise ValueError('{}: after linear edge + bulk, still no good splitting interpolation, best to skip this q'.format(q))
    return bulk


def interpolate_edge(edge, kind='cubic'):
    """
    builds a 1D interpolant of the edge of the splitting surface
    :param edge: points on the edge of the splitting surface
    :param kind: either 'cubic' or 'linear'
    :return: an interpolant function f that can be evaluated as f(phi) giving the theta value of that phi that lies on
        the edge
    """
    try:
        uphis = np.unique(edge[:, 2], return_index=True)  # indices where phi values are unique
        uniquel1edge = edge[uphis[1]]  # unique edge points
        uniquel1edge = uniquel1edge[uniquel1edge[:, 2].argsort()]
        # make interpolant for theta(phi)
        f = spinp.interp1d(uniquel1edge[:, 2], uniquel1edge[:, 1], kind=kind, bounds_error=False)
    except scpa.qhull.QhullError as e:
        print(' QhullError: probably cannot tesselate tracks')
        print(e, edge)
        f = None
    except IndexError as e:
        print(' Cannot index tracks, probably no tracks at all')
        print(e, edge)
        f = None
    return f


def interpolate_tracks(points, TTHETAS, PPHIS, method='cubic'):
    """
    interpolate points on a grid of (THETA, PHIS)
    :param points: points on the splitting surface
    :param TTHETAS: meshgrid of thetas
    :param PPHIS: meshgrid of phis
    :param method: either 'cubic', 'linear' or 'nearest'
    :return: 2d list of interpolated values of the splitting surface on the theta, phi grid
    """
    try:
        inter = spinp.griddata((points[:, 1], points[:, 2]), points[:, 0], (TTHETAS, PPHIS), method=method)
    except scpa.qhull.QhullError as e:
        print(' QhullError: probably cannot tesselate tracks')
        print(e, points)
        return None
    except IndexError as e:
        print(' Cannot index tracks, probably no tracks at all')
        print(e, points)
        return None
    return inter


def check_inter_and_clean(bulk_inter, edge_inter, tthetas, pphis):
    """
    check whether the interpolation is good given the top edge of the splitting surface and put all interpolated
    values outside of the relevant domain to NaN.
    :param bulk_inter: bulk of the interpolation (from interpolate_tracks)
    :param edge_inter: edge interpolant (from interpolate_edge)
    :param tthetas: thetas that are sampled
    :param pphis: phis that are sampled
    :return: 2d list of interpolated values or None if the original interpolation was no good.
    """
    ok = True
    i = 0
    while ok and i < len(bulk_inter[:, 1]):  # go over all phis
        thetamin = edge_inter(pphis[i])  # get the minimal theta for this phi from interpolant
        if np.isnan(thetamin):
            break  # if here, we have checked all phis as we are outside domain on which the edge was interpolated
        j = 0
        while ok and j < len(bulk_inter[i]) and tthetas[j] >= thetamin:  # now go over thetas for this phi
            if np.isnan(bulk_inter[i, j]):
                ok = False
                break  # if here, the interpolation of bulk has produced a nan within the relevant domain, not ok!
            j += 1
        bulk_inter[i, j:] = np.nan  # set all values in the interpolation outside of the edge to nan
        i += 1
    if not ok:
        return None
    return bulk_inter

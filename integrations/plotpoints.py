"""
convenience methods for plotting point clouds of calculated equipotential shells, uses mayavi
"""
import _tkinter  # noqa
import mayavi.mlab as mlab
import numpy as np
import scipy.spatial as scpa


@mlab.animate(delay=50)
def anim():
    f = mlab.gcf()
    while 1:
        f.scene.camera.azimuth(1)
        f.scene.render()
        yield


def plot_points(points, *mlabargs, x_offset=0, **mlabkwargs):
    mlab.points3d(points[:, 0] + x_offset, points[:, 1], points[:, 2], *mlabargs, **mlabkwargs)


def symmetrize(endpoints):
    endpoints = np.vstack(
        (endpoints, np.array([endpoints[:, 0], -endpoints[:, 1], endpoints[:, 2]]).T))
    endpoints = np.vstack(
        (endpoints, np.array([endpoints[:, 0], endpoints[:, 1], -endpoints[:, 2]]).T))
    return endpoints


def sym_hullify_and_plot(endpoints, **mlabkwargs):
    endpoints = symmetrize(endpoints)
    hull = scpa.ConvexHull(endpoints)
    sim = hull.simplices
    mask = [True] * sim.shape[0]
    for i in range(len(sim)):
        for pt in sim[i]:
            if 0. in hull.points[pt]:
                mask[i] = False
    sim = sim[mask]
    mlab.triangular_mesh(endpoints[:, 0], endpoints[:, 1], endpoints[:, 2], sim,
                         **mlabkwargs)


def sym_and_plot(endpoints, *mlabargs, **mlabkwargs):
    endpoints = symmetrize(endpoints)
    plot_points(endpoints, *mlabargs, **mlabkwargs)

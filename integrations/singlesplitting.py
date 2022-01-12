"""
script to make a splitting surfaces for a single mass ratio q. Useful to see what is going on.
"""
# %%
import matplotlib

matplotlib.use('Qt5Agg')
import _tkinter  # noqa
import numpy as np
import mayavi.mlab as mlab
import l_p_points as lp
import c_rocheFunctions as cf
import c_split_tracer as ct
import c_util as cu
import interpolatesplits as inter
import py_rocheFunctions as rf
import plotpoints as pp
import matplotlib.pyplot as plt

# %%  # info about the selected q
q = 1.0
x_l1 = lp.l1(q)
x_lfar = lp.lfar(q)
x_lout = lp.lout(q)
lout_pot = cf.roche(x_lout, 0, 0, q)
l4_pot = cf.roche(0.5, 0.5 * np.sqrt(3), 0, q)

# %%  # calculate the splitting surface around L1
l1_pos = ct.c_gradient_tracer(x_l1, 0.5 * lout_pot + 0.5 * l4_pot, q)
print('l1 found')
l1_bulk = np.array(l1_pos['bulk'])
l1_edge = np.array(l1_pos['edge'])

# %%  # check out the point cloud of the splitting surface
trim = 1
mlab.figure()
pp.plot_points(l1_bulk[::trim], lambda x, y, z: rf.roche(x, y, z, q))

mlab.axes()
mlab.colorbar()
mlab.show()

# %%  # do a mock interpolation on a theta, phi grid
phinum = 660
dphi = np.pi / phinum
phis = np.linspace(dphi / 2, np.pi - dphi / 2, phinum)
thetanum = int(phinum / 2)
dcostheta = 1. / thetanum
thetas = np.arccos(np.linspace(dcostheta / 2, 1 - dcostheta / 2, thetanum))
THETAS, PHIS = np.meshgrid((thetas, phis))
l1_bulk = cu.cart2spher(l1_bulk)
l1_edge = cu.cart2spher(l1_edge)
l1inter = inter.interpolate_split(l1_bulk, l1_edge, thetas, phis, q)

plt.contourf(THETAS, PHIS, l1inter)
plt.colorbar()
plt.show()

# %%  # compute lout splitting surface
lout_pos = ct.c_gradient_tracer(x_lout, 0.5 * lout_pot + 0.5 * l4_pot, q)
lout_pos = np.array(lout_pos)
print('lout found')
lout_bulk = np.array(lout_pos['bulk'])
lout_edge = np.array(lout_pos['edge'])

# %%  # check out point cloud
trim = 1
pp.plot_points(lout_bulk[::trim], lambda x, y, z: rf.roche(x, y, z, q))

mlab.axes()
mlab.colorbar()
mlab.show()

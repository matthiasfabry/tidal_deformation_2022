"""
script to make a Roche integration for a single mass ratio q. Useful to see what is going on.
"""
# %%
import matplotlib
matplotlib.use('Qt5Agg')
import _tkinter  # noqa
import time as t
import mayavi.mlab as mlab
import numpy as np
import loadsplitting as split
import c_main_tracer as ct
import c_util as cu
import c_rocheFunctions as cf
import interpolatesplits as inter
import l_p_points as lp
import matplotlib.pyplot as plt
import plotpoints as pp
import py_rocheFunctions as pr


# %%  # select q and compute Lagrangian points and stuff
q = 10**1
q = split.closest_q(q)  # make sure q has computed splitting surfaces
print(q)
# determine Lagrange points and set potential to determine volume/surface/etc...
x_l1 = lp.l1(q)
x_lfar = lp.lfar(q)
x_lout = lp.lout(q)
print('x-lagrangians: ', x_lout, x_l1, x_lfar)
L1pot = cf.roche(x_l1, 0, 0, q)
Loutpot = cf.roche(x_lout, 0, 0, q)
Lfarpot = cf.roche(x_lfar, 0, 0, q)
L4pot = cf.roche(0.5, 0.5 * np.sqrt(3), 0, q)
print('potential at Lout =', Loutpot, '=', Loutpot / L1pot, 'times pot @ L1')
print('potential at L1 =', L1pot)
print('potential at Lfar =', Lfarpot, '=', Lfarpot / L1pot, 'times pot @ L1')
print('potential at L4 =', L4pot, '=', L4pot / L1pot, 'times pot @ L1')
targetpot = 0.999999 * Loutpot
if targetpot > Lfarpot:
    print('system will lose mass through Lfar!')
if targetpot > Loutpot:
    print('system will lose mass through Lout!')
print('target potential:', targetpot, "=", targetpot / L1pot, 'times pot @ L1')
if targetpot < Loutpot:
    # determine rest of extrema
    x_p3 = lp.p3(q, targetpot)
    x_p5, y_p5 = lp.p5(q, targetpot)
    z_p7 = lp.p7(q, x_p5, targetpot)
    print('P3 =', [x_p3, 0, 0])
    print('radius to P5_x =', x_p5 - x_p3)
    print('P5 =', [x_p5, y_p5, 0])
    print('P7 =', [x_p5, 0, z_p7])
    
# %%  # load the respective splitting surfaces
print('loading splitting surfaces')
splits = split.load_split(q)
sl1 = splits['l1_bulk']
sl3 = splits['l3_bulk']
l1edge = splits['l1_edge']
l3edge = splits['l3_edge']
l1bulkspher = cu.cart2spher(sl1)
loutbulkspher = cu.cart2spher(sl3)
l1edgespher = cu.cart2spher(l1edge)
loutedgespher = cu.cart2spher(l3edge)

#%%  # check out the splitting surface
mlab.figure()
pp.plot_points(sl1, lambda x, y, z: pr.roche(x, y, z, q), scale_factor=0.001, colormap='inferno')
pp.plot_points(l1edge, scale_factor=0.001)
mlab.axes()
mlab.show()

# %%  # make grid of phi, theta we're gonna integrate over
phinum = 4*660
dphi = np.pi / phinum
phis = np.linspace(dphi / 2, np.pi - dphi / 2, phinum)
thetanum = int(phinum / 2)
dcostheta = 1. / thetanum
thetas = np.arccos(np.linspace(dcostheta / 2, 1 - dcostheta / 2, thetanum))
THETAS, PHIS = np.meshgrid(thetas, phis)

#%%  # interpolate the splitting surfaces on that grid
l1_inter = inter.interpolate_split(l1bulkspher, l1edgespher, thetas, phis, q)
lout_inter = inter.interpolate_split(loutbulkspher, loutedgespher, thetas, phis, q)

#%%  # check out interpolation result
plt.figure()
plt.contourf(THETAS, PHIS, l1_inter, levels=50)
plt.colorbar()
plt.show()

# %%  # do the roche integration
print('rays to integrate:', phinum * thetanum)
time = t.time()
# save_points must be True in order to visualize the hull
endpoints1, vals = ct.c_roche_tracer(q, targetpot, x_l1, thetas, phis,
                                     l1_inter=l1_inter, lout_inter=lout_inter, save_points=True)
endpoints1 = np.array(endpoints1)
print('main tracer done in', t.time() - time, 's\n')
print('Volume =', vals[0])
print('Surface area =', vals[1])
print('Avg effective grav =', vals[2])
print('Avg inv effective grav =', vals[3])
print('Volume equivalant radius =', vals[4])
print('Moment of inertia =', vals[5])
print('f_P =', vals[6])
print('f_T =', vals[7])

# %%  # plot the resulting point cloud on top of the splitting surfaces
mlab.figure()
trim = 10
pp.plot_points(endpoints1[::trim], color=(1, 0.1, 0), scale_factor=0.001)
pp.plot_points(splits['l1_bulk'], scale_factor=0.01, color=(0, 0, 0))
pp.plot_points(splits['l3_bulk'], scale_factor=0.01, color=(0, 0, 0))

pp.sym_and_plot(endpoints1, scale_factor=0.01, color=(0, 0.1, 1))
mlab.axes()
mlab.orientation_axes()
mlab.show()

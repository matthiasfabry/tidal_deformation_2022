"""
module to pre-interpolate the results on a square grid that the MESA interpolators can work with
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as spinp

import resultloader as rl

# %%  # get relevant integration results
vals = rl.get_results(rebuild=True)

qs = np.unique(vals[:, 0])
toadd = np.empty((qs.size, vals[0].size))
for i in range(len(toadd)):  # make extra entries for the r_psi = 0 case for all qs, we know f_P=f_T=1 there and i_rot=0
    toadd[i] = np.array([qs[i], 0, np.NINF, 0, 0, np.PINF, 0, 0, 0, 1, 1, 0, np.NINF, 0, 0, np.NINF])
vals = np.vstack((vals, toadd))

max_requiv = np.max(vals[:, 7])
max_rrRL = np.max(vals[:, 14])

# make one extra entry in the top left of the q, r_psi space with standard spherical data, to make the interpolation
# not give NaNs there, we shouldn't need this part of the parameter space, but this is just to make MESA behave
toadd = np.array([min(qs), 0, 0, 4 * np.pi / 3 * max_requiv ** 3, 4 * np.pi * max_requiv ** 2, 0, 0,
                  max_requiv, max_requiv ** 2, 0, 0, 0, 0, 0, max_rrRL, 0])
vals = np.vstack((vals, toadd))

# %%  # grid of log q and r_psi/r_RL to interpolate on
lqs = np.linspace(-7, 7, 401)
roverrRL = np.linspace(0, max_rrRL, 5000)
LQS, RS = np.meshgrid(lqs, roverrRL)

logqs = np.log10(vals[:, 0])
# interpolate!
fp_table = spinp.griddata((logqs, vals[:, 14]), vals[:, 9], (LQS, RS))
ft_table = spinp.griddata((logqs, vals[:, 14]), vals[:, 10], (LQS, RS))
irot_table = spinp.griddata((logqs, vals[:, 14]), vals[:, 8], (LQS, RS))
area_table = spinp.griddata((logqs, vals[:, 14]), vals[:, 4], (LQS, RS))

# %%  # check out the result
plt.contourf(LQS, RS, fp_table)


# %%
# fill out any remaining nans in the grid with extrapolated values
def fill_out(grid):
    for i in range(grid[0, :].size):
        size = grid[:, 0].size
        j = size - 1
        while j >= 0 and np.isnan(grid[j, i]):  # nans should only appear at the top, we bounded with the r=0 case above
            j -= 1
        if j != size - 1:
            for k in range(j, size - 1):  # extrapolate outward again
                grid[k, i] = grid[j, i] + (grid[j, i] - grid[j - 20, i]) / (roverrRL[j] - roverrRL[j - 20]) \
                             * (roverrRL[k] - roverrRL[j])


fill_out(fp_table)
fill_out(ft_table)
fill_out(irot_table)
fill_out(area_table)

# %%  # check out fp(r) for a selected log q, mainly to verify the extrapolation behaves well
plt.figure()
num = 110
print(lqs[num])
plt.plot(roverrRL, fp_table[:, num])


# %% # write interpolated tables to a file readable for MESA, run_star_extras.f90 setup_interpolator() reads in these
# files
def write_table_to_file(table, filename):
    with open(filename, 'w') as f:
        f.write(str(len(lqs)))
        f.write('\n')
        for lq in lqs:
            f.write(str(np.round(lq, 5)))
            f.write('\n')
        f.write(str(len(roverrRL)))
        f.write('\n')
        for r in roverrRL:
            f.write(str(np.round(r, 7)))
            f.write('\n')
        for i in range(len(lqs)):
            line = ""
            for j in range(len(roverrRL)):
                line += " " + str(np.round(table[j, i], 10))
            f.write(line)
            f.write('\n')


write_table_to_file(fp_table, 'fp_data.txt')
write_table_to_file(ft_table, 'ft_data.txt')
write_table_to_file(irot_table, 'irot_data.txt')
write_table_to_file(area_table, 'area_data.txt')

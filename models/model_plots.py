import glob

import matplotlib.legend_handler as handler
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1.inset_locator as inset
import numpy as np

import mkipp.mesadata as md

plt.style.use('mesa')

# %%  # overcontact models
newbctidalh = md.MesaData(glob.glob('methodspaper/overcontact/newBCtidal/LOGS1/history.data')[0])
oldbctidalh = md.MesaData(glob.glob('methodspaper/overcontact/oldBCtidal/LOGS1/history.data')[0])
newbcsingleh = md.MesaData(glob.glob('methodspaper/overcontact/newBCsingle/LOGS1/history.data')[0])
oldbcsingleh = md.MesaData(glob.glob('methodspaper/overcontact/oldBCsingle/LOGS1/history.data')[0])
oldbcnoroth = md.MesaData(glob.glob('methodspaper/overcontact/oldBCnorot/LOGS1/history.data')[0])
newbcsingleh.trim_PMS()
newbctidalh.trim_PMS()
oldbctidalh.trim_PMS()
oldbcsingleh.trim_PMS()
oldbcnoroth.trim_PMS()
# %%  ## contact plots
oldbcnorotc = np.argmin(abs(oldbcnoroth.get('rl_relative_overflow_1')))
oldbcsinglec = np.argmin(abs(oldbcsingleh.get('rl_relative_overflow_1')))
newbcsinglec = np.argmin(abs(newbcsingleh.get('rl_relative_overflow_1')))
oldbctidalc = np.argmin(abs(oldbctidalh.get('rl_relative_overflow_1')))
newbctidalc = np.argmin(abs(newbctidalh.get('rl_relative_overflow_1')))
md.plot(newbctidalh, 'star_age', 'log_Teff', xscale=1e6, color='r', label='new BC, tidal')
plt.plot(newbctidalh.get('star_age')[newbctidalc] / 1e6, newbctidalh.get('log_Teff')[newbctidalc], 'r.')
md.plot(oldbctidalh, 'star_age', 'log_Teff', xscale=1e6, color='b', label='def. BC, tidal')
plt.plot(oldbctidalh.get('star_age')[oldbctidalc] / 1e6, oldbctidalh.get('log_Teff')[oldbctidalc], 'b.')
md.plot(newbcsingleh, 'star_age', 'log_Teff', 'r--', dashes=[2, 2], xscale=1e6, label='new BC, single')
plt.plot(newbcsingleh.get('star_age')[newbcsinglec] / 1e6, newbcsingleh.get('log_Teff')[newbcsinglec], 'r.')
md.plot(oldbcsingleh, 'star_age', 'log_Teff', 'b--', dashes=[2, 2], xscale=1e6, label='def. BC, single')
plt.plot(oldbcsingleh.get('star_age')[oldbcsinglec] / 1e6, oldbcsingleh.get('log_Teff')[oldbcsinglec], 'b.')
md.plot(oldbcnoroth, 'star_age', 'log_Teff', xscale=1e6, color='k', label='no rotation')
plt.plot(oldbcnoroth.get('star_age')[oldbcnorotc] / 1e6, oldbcnoroth.get('log_Teff')[oldbcnorotc], 'k.')

plt.xlabel('age (Myr)')
plt.ylabel('$\log T_{\\rm eff}$')
handles, _ = plt.gca().get_legend_handles_labels()


class HandlerCircle(handler.HandlerPatch):
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
        p = mpatch.Circle(xy=center, radius=1.5)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]


toadd = mpatch.Circle((0., 0.), radius=0.1, color='gray', label='$L_1$ contact')
handles.append(toadd)
plt.legend(handles=handles, handler_map={mpatch.Circle: HandlerCircle()})
plt.show()
# %%  # single star
newbcsingleh = md.MesaData(glob.glob('methodspaper/singlerot/newBC/LOGS/history.data')[0])
oldbcsingleh = md.MesaData(glob.glob('methodspaper/singlerot/oldBC/LOGS/history.data')[0])
oldbcsingleh.trim_PMS()
newbcsingleh.trim_PMS()
# %%  ##  single star plots
plt.figure()
oldbcsingleh.trim('center_h1', 0.001, 1)
newbcsingleh.trim('center_h1', 0.001, 1)

md.plot_HRD(oldbcsingleh, 'b', label='default BC')
md.plot_HRD(newbcsingleh, 'r', label='new BC')
plt.legend()
axins2 = inset.inset_axes(plt.gca(), width="1000%", height="1000%",
                          bbox_to_anchor=(4.055, 2.75, 0.0055, 0.0175), loc=4,
                          bbox_transform=plt.gca().transData, borderpad=0)
md.plot_HRD(oldbcsingleh, 'b', label='default BC')
md.plot_HRD(newbcsingleh, 'r', label='new BC')
plt.xlabel('')
plt.ylabel('')
plt.xlim((4.065, 4.053))
plt.ylim((3.085, 3.14))
plt.show()
# %%  # single star area ratio
model = newbcsingleh
fig = plt.figure(constrained_layout=True)
gs = fig.add_gridspec(2, 1, height_ratios=(2, 7),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      hspace=0.05)
fig.add_subplot(gs[0, 0])
plt.plot(model.get('star_age') / 1e6, model.get('omega_surf'), 'r')
plt.ylabel('$\Omega/\Omega_{\\rm crit}$')
plt.gca().axes.xaxis.set_ticklabels([])
ax = fig.add_subplot(gs[1, 0])
plt.plot(model.get('star_age') / 1e6, 4 * np.pi * model.get('radius') ** 2 * Rsun ** 2 / model.get('surf_area'), 'r')
plt.xlabel('age (Myr)')
plt.ylabel(r'$4\pi r_\Psi^2 / S_\Psi$')
plt.tight_layout()
plt.show()

# %%  # semi detached EB
newbctidalh = md.MesaData(glob.glob('methodspaper/semidetached/newBCtidal/LOGS1/history.data')[0])
oldbctidalh = md.MesaData(glob.glob('methodspaper/semidetached/oldBCtidal/LOGS1/history.data')[0])
newbcsingleh = md.MesaData(glob.glob('methodspaper/semidetached/newBCsingle/LOGS1/history.data')[0])
oldbcsingleh = md.MesaData(glob.glob('methodspaper/semidetached/oldBCsingle/LOGS1/history.data')[0])
oldbcnoroth = md.MesaData(glob.glob('methodspaper/semidetached/oldBCnorot/LOGS1/history.data')[0])
newbcsingleh.trim_PMS()
newbctidalh.trim_PMS()
oldbctidalh.trim_PMS()
oldbcsingleh.trim_PMS()
oldbcnoroth.trim_PMS()
# %%  ## detached plots
plt.figure()
md.plot(oldbcnoroth, 'center_h1', 'radius', 'k', label='no rotation')
md.plot(oldbcsingleh, 'center_h1', 'radius', 'b--', dashes=[1, 2], label='def. BC, single')
md.plot(newbcsingleh, 'center_h1', 'radius', 'r--', dashes=[1, 2], label='new BC, single')
md.plot(oldbctidalh, 'center_h1', 'radius', 'b', label='def. BC, tidal')
md.plot(newbctidalh, 'center_h1', 'radius', 'r', label='new BC, tidal')
plt.xlabel('$X_c$')
plt.ylabel('$R/R_\\odot$')
plt.gca().invert_xaxis()
plt.legend(loc=4, fontsize=7)
axins = inset.inset_axes(plt.gca(), width="1000%", height="1000%",
                         bbox_to_anchor=(0.625, 6.8, 0.025, 0.2), loc=2,
                         bbox_transform=plt.gca().transData, borderpad=0)
md.plot(oldbcnoroth, 'center_h1', 'radius', 'k', label='no rotation')
md.plot(oldbctidalh, 'center_h1', 'radius', 'b', label='def. BC, tidal')
md.plot(newbctidalh, 'center_h1', 'radius', 'r', label='new BC, tidal')
md.plot(oldbcsingleh, 'center_h1', 'radius', 'b--', dashes=[1, 3], label='def. BC, single')
md.plot(newbcsingleh, 'center_h1', 'radius', 'r--', dashes=[0, 2, 1, 1], label='new BC, single')
plt.xlabel('')
plt.ylabel('')
plt.xlim((0.12, 0.06))
plt.ylim((6.7, 7.4))
plt.show()

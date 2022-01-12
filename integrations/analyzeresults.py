import numpy as np

import resultloader as rl

cols = rl.get_cols()

import _tkinter  # noqa
import matplotlib.pyplot as plt
import singlestarcheckingfuncs as single
import loadsplitting as split
import matplotlib.tri as tri


# custom plotting routines
def make_tricontour(vals, x_ind, y_ind, z_ind, **kwargs):
    plt.tricontourf(vals[:, x_ind], vals[:, y_ind], vals[:, z_ind], cmap='inferno', **kwargs)
    plt.xlabel(cols[x_ind])
    # plt.ylabel(cols[y_ind])
    plt.colorbar()
    plt.title(cols[z_ind])


def make_tricontour_semilogx(vals, x_ind, y_ind, z_ind, **kwargs):
    plt.tricontourf(np.log10(vals[:, x_ind]), vals[:, y_ind], vals[:, z_ind], cmap='inferno',
                    **kwargs)
    plt.xlabel(r'$\log\,$' + cols[x_ind])
    # plt.ylabel(cols[y_ind])
    plt.colorbar()


def scatter_data(vals, x_ind, y_ind, **kwargs):
    plt.scatter(vals[:, x_ind], vals[:, y_ind], **kwargs)
    plt.xlabel(cols[x_ind])
    # plt.ylabel(cols[y_ind])


def scatter_data_semilogx(vals, x_ind, y_ind, **kwargs):
    plt.scatter(np.log10(vals[:, x_ind]), vals[:, y_ind], **kwargs)
    plt.xlabel(r'$\log\,$' + cols[x_ind])
    # plt.ylabel(cols[y_ind])


def scatter_marginalized(vals, x_ind, y_ind, constraint, label=None, **kwargs):
    mask = vals[:, constraint[0]] == constraint[1]
    x = vals[:, x_ind][mask]
    y = vals[:, y_ind][mask]
    if label is not None:
        label = label
    else:
        label = '{}({}); {} = {}'.format(cols[y_ind], cols[x_ind], cols[constraint[0]],
                                         constraint[1])
    plt.scatter(x, y, label=label, **kwargs)
    plt.xlabel(cols[x_ind])
    # plt.ylabel(cols[y_ind])


def scatter_marginalized_semilogx(vals, x_ind, y_ind, constraint, **kwargs):
    mask = vals[:, constraint[0]] == constraint[1]
    x = vals[:, x_ind][mask]
    y = vals[:, y_ind][mask]
    plt.scatter(np.log10(x), y,
                label='{}({}); {} = {}'.format(cols[y_ind], r'$\log\,$' + cols[x_ind],
                                               cols[constraint[0]], constraint[1]), **kwargs)
    plt.xlabel(r'$\log\,$' + cols[x_ind])
    # plt.ylabel(cols[y_ind])


# %%
if __name__ == '__main__':
    # %%  load in the set of integrations
    vals = rl.get_results(rebuild=True)
    # %%  # plot fp curves + single comparison (Fig. 2)
    plt.figure()
    vals = vals[vals[:, 14].argsort()]
    cs = ['k', 'b', 'r']
    qs = [0.1, 1, 10]
    for i in range(len(qs)):
        q = split.closest_q(qs[i])
        plt.plot(vals[vals[:, 0] == q, 14], vals[vals[:, 0] == q, 9], cs[i],
                 label='$\\log q={}$'.format(round(np.log10(q))))
        oms = vals[vals[:, 0] == q, 11]
        plt.plot(vals[vals[:, 0] == q, 14], single.polyfP(oms), cs[i], ls='--')
    lgd = plt.legend()
    plt.xlabel('$r_\\Psi/r_{\\rm RL}$')
    plt.ylabel('$f_P$')
    plt.xlim((0, 1.45))
    plt.ylim((0.3, 1.05))

    # single equivalent legend entry
    def add_patch(legend):
        from matplotlib.lines import Line2D
        ax = legend.axes
    
        handles, labels = ax.get_legend_handles_labels()
        handles.append(Line2D([0], [0], linestyle='--', color='gray'))
        labels.append("Single equivalent")
    
        legend._legend_box = None
        legend._init_legend_box(handles, labels)
        legend._set_loc(legend._loc)
        legend.set_title(legend.get_title().get_text())


    add_patch(lgd)
    plt.show()
    
    # %%  # make Eddington factor comparison in the single rotating star case (Fig. 4)
    plt.figure()
    oms = np.linspace(0, 1, 200)
    plt.plot(oms, 1 - oms ** 2, 'b', label='Langer')
    plt.plot(oms, 1. / single.polyfTovertP(oms), 'r', label='$f_P/f_T$')
    plt.ylabel('$\Gamma_{\\rm max}$')
    plt.xlabel('$\omega$')
    plt.legend()
    plt.show()
    
    # %%  # make eddington factor reduction plot in tidal case (Fig. 5)
    plt.figure()
    
    # all this is to mask the bad part of the triangulation in tricontourf
    maxrs = list()
    for q in np.unique(vals[:, 0]):
        maxr = np.max(vals[vals[:, 0] == q, 14])
        maxrs.append((q, maxr))
    maxrs = np.array(maxrs)
    isbad = np.empty(len(vals), dtype=bool)
    for i in range(len(vals)):
        qhere = vals[i, 0]
        isbad[i] = vals[i, 14] == maxrs[maxrs[:, 0] == qhere, 1]
    t = tri.Triangulation(np.log10(vals[:, 0]), vals[:, 14])
    mask = np.any(np.where(isbad[t.triangles], True, False), axis=1)
    t.set_mask(mask)
    z = np.ma.asarray(vals[:, 9] / vals[:, 10])
    
    # finally plot
    fig = plt.figure()
    gs = fig.add_gridspec(2, 2, width_ratios=(7, 2), height_ratios=(2, 7),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)
    fig.add_subplot(gs[0, 0])
    plt.plot(np.log10(vals[vals[:, 1] == 1, 0]), vals[vals[:, 1] == 1, 9] / vals[vals[:, 1] == 1, 10], label='$L_1$')
    plt.plot(np.log10(vals[vals[:, 12] == 1, 0]), vals[vals[:, 12] == 1, 9] / vals[vals[:, 12] == 1, 10],
             label='$L_{\\rm out}$')
    plt.ylabel('$f_P/f_T$')
    plt.xlim((-2, 2))
    
    ax = fig.add_subplot(gs[1, 0], autoscalex_on=False)
    plt.gca().sharex(ax)
    plt.tricontourf(t, z, cmap='inferno', levels=10)
    cbar = plt.colorbar()
    
    plt.tricontour(t, z, colors='gray')
    l3s = list()
    for q in np.unique(vals[:, 0]):
        l3r = vals[np.logical_and(vals[:, 0] == q, vals[:, 12] == 1), 14]
        l3s.append((np.log10(q), l3r))
    maxs = np.array(l3s)
    plt.plot(maxs[:, 0], maxs[:, 1], 'g', label='$r_{L_{\\rm out}}$')
    
    plt.xlabel('$\\log q$')
    plt.ylabel('$r_\\Psi/r_{\\rm RL}$')
    plt.xlim((-2, 2))
    plt.ylim((0.7, 1.5))
    plt.legend(loc=4)
    cbar.set_label('$f_P/f_T$', rotation=270, labelpad=15)
    plt.show()
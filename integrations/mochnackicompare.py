"""
script to compare integration results against Mochnacki 1984
"""

import resultloader as rl
import matplotlib.pyplot as plt
import mochnackicheckingfuncs as moch
import numpy as np
plt.style.use('mesa')

#%%  # load in results for various resolutions
vals = rl.get_results('resultsmoch', rebuild=True)
vals2 = rl.get_results('resultsmoch2', rebuild=True)
vals4 = rl.get_results('resultsmoch4', rebuild=True)
vals8 = rl.get_results('resultsmoch8', rebuild=True)
q = 1.0
kwargs = dict(nrows=2, gridspec_kw={'height_ratios': [4, 1]}, sharex=True, constrained_layout=True)


#%%  # make Mochnacki comparison figure (Fig. 3)
def plotmoch(mochs, column):
    fig, axs = plt.subplots(**kwargs)
    vvals = vals[vals[:, 14].argsort()]
    As = vvals[vvals[:, 0] == q, column]
    Fs = vvals[vvals[:, 0] == q, 14]
    vvals2 = vals2[vals2[:, 14].argsort()]
    As2 = vvals2[vvals2[:, 0] == q, column]
    Fs2 = vvals2[vvals2[:, 0] == q, 14]
    vvals4 = vals4[vals4[:, 14].argsort()]
    As4 = vvals4[vvals4[:, 0] == q, column]
    Fs4 = vvals4[vvals4[:, 0] == q, 14]
    vvals8 = vals8[vals8[:, 14].argsort()]
    As8 = vvals8[vvals8[:, 0] == q, column]
    Fs8 = vvals8[vvals8[:, 0] == q, 14]

    axs[0].plot(moch.rs[15]/moch.rs[15][9], mochs[15], '--', label=r'Mochnacki (1984)', color='orange')
    axs[0].plot(Fs, As, 'bd', label=r'default resolution')
    axs[0].plot(Fs2, As2, 'rx', label=r'double resolution')
    axs[0].plot(Fs4, As4, 'g.', label=r'quadruple resolution')
    axs[0].plot(Fs8, As8, 'm+', label=r'octuple resolution')

    axs[1].plot(Fs, (mochs[15] - As) / As, 'bd')
    axs[1].plot(Fs2, (mochs[15] - As2) / As2, 'rx')
    axs[1].plot(Fs4, (mochs[15] - As4) / As4, 'g.')
    axs[1].plot(Fs4, (mochs[15] - As8) / As8, 'm+')
    
    axs[1].set_xlabel('$r_\\Psi/r_{\\rm RL}$')
    axs[0].set_ylabel(r'${}$'.format(rl.get_cols()[column]))
    axs[0].legend()
    axs[1].set_ylabel(r'rel. error')

#%%
plotmoch(4*np.pi*moch.rs**4/(moch.areas*moch.gminus1s), 9)


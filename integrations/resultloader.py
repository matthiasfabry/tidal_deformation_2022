"""
methods for loading/manipulating results saved by tracertreads.py
"""
import numpy as np


def get_cols():
    """
    convenience list containing names for all columns in a result. useful for plotting.
    i.e.: vals[i, 4] will contain an Area and vals[j, 10] will contain the value of f_T
    :return: list of column names
    """
    return [r'$q$',  # 0  mass ratio
            r'$\Psi/\Psi_{L1}$',  # 1  fractional potential wrt potential at L1
            r'$\Psi$',  # 2  equipotential value
            r'$V_\Psi$',  # 3   equipotential volume
            r'$S_\Psi$',  # 4   equipotential surface area
            r'$\langle g\rangle$',  # 5   average effective gravity
            r'$\langle g^{-1}\rangle$',  # 6   average inverse effective gravity
            r'$r_\Psi$',  # 7   equivalent radius
            r'$i_{\rm rot}$',  # 8   moment of inertia
            r'$f_{P}$',  # 9   f_P correction factor
            r'$f_{T}$',  # 10   f_T correction factor
            r'$\omega_{\rm single}$',  # 11   single equivalent rotation rate
            r'$\Psi/\Psi_{L_{\rm out}}$',  # 12   fractional potential wrt potential at outer L point
            r'Filling',  # 13   Mochnacki's filling factor
            r'$r_\Psi/r_{\rm RL}$',  # 14   fractional roche lobe radius
            r'$\Psi/\Psi_{L_{\rm far}}$',  # 15   fractional potential wrt potential at far L point
            ]


def get_results(folder='results', rebuild=False):
    """
    loads all available results from saved data
    :param folder: where to look for data
    :param rebuild: If true, go over all txt files to rebuild array to save in a numpy array
        If false, simply load in the already the numpified array
    :return: ndarray containing all integration data
    """
    import glob

    if rebuild:
        files = glob.glob(folder + '/q*/super*.txt')
        vals = np.empty((len(files), len(np.loadtxt(files[0]))))
        for i in range(len(files)):
            vals[i, :] = np.loadtxt(files[i])
        vals = add_post_calcs(vals)
        np.save(folder + '/consolidated_results', vals)
    vals = np.load(folder + '/consolidated_results.npy')
    return vals[vals[:, 0].argsort()]


def add_post_calcs(vals):
    """
    does post processing to add useful columns to the result array, they are listed in get_cols()
    :param vals: current ndarray of all integration results
    :return: new ndarray containing more columns of data
    """
    import singlestarcheckingfuncs as single
    import c_rocheFunctions as cf
    import l_p_points as lp
    import mochnackicheckingfuncs as moch
    
    n = len(vals)
    omegas_single = np.empty((n, 1))  # have to be 2d in order to hstack them later
    fillings = np.empty((n, 1))
    l2supers = np.empty((n, 1))
    l3supers = np.empty((n, 1))
    rover_rl = np.empty((n, 1))
    for i in range(n):
        # print(vals[i, 0], vals[i, 1])
        try:  # what would the fractional rotation rate be if the binary lobe is interpreted as single star rotation
            omegas_single[i] = single.omega_single_star(vals[i, 7], vals[i, 0])
        except (ValueError, SystemError) as e:  # can be undetermined
            # print(vals[i, 0], vals[i, 1], vals[i, 7], e)
            omegas_single[i] = np.nan

        fillings[i] = moch.filling(vals[i, 2], vals[i, 0])  # Moch's Filling factor
        l2supers[i] = vals[i, 2] / cf.roche(lp.lfar(vals[i, 0]), 0, 0, vals[i, 0])  # fractional potential vs that of L2
        l3supers[i] = vals[i, 2] / cf.roche(lp.lout(vals[i, 0]), 0, 0, vals[i, 0])  # fractional potential vs that of L3
        
        temp = vals[np.logical_and(vals[:, 1] == 1., vals[:, 0] == vals[i, 0]), 7]
        # print(i, temp, vals[i, 0], vals[i, 1])
        rover_rl[i] = vals[i, 7] / temp  # fractional roche lobe radius
    return np.hstack((vals, omegas_single, l3supers, fillings, rover_rl, l2supers))  # add them all to the vals array


def roverrl(vals):
    """
    convenience method to recompute the fractional rl sizes
    :param vals: ndarray with integration data
    :return: new ndarray with modified fractional rl sizes
    """
    for i in range(len(vals)):
        temp = vals[np.logical_and(vals[:, 1] == 1., vals[:, 0] == vals[i, 0]), 7]
        # print(vals[i, 0], temp)
        vals[i, 14] = vals[i, 7] / temp  # overwrite the 14th column with new computed rlsize
    return vals


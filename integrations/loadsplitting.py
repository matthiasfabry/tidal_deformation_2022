"""
methods to load splitting surfaces into memory
"""
import glob

import numpy as np


def closest_q(qq):
    """
    selects the closest mass ratio available from all constructed splitting surfaces
    :param qq: input mass ratio M2/M1
    :return: output mass ratio M2/M1 that is available
    """
    split_files = glob.glob('splits/**')
    n = len(split_files)
    qs = np.empty(n)
    for i in range(n):
        qs[i] = split_files[i].rsplit('.', maxsplit=1)[0].split('_', maxsplit=1)[-1]
    index = 0
    for ii in range(n):
        if abs(qq - qs[ii]) < abs(qq - qs[index]):
            index = ii
    return qs[index]


def load_split(qq):
    """
    loads in the splitting surfaces of the given mass ratio q
    :param qq: mass ratio M2/M1
    :return: dict of splitting surface data, as computed in splittingthreads.py
    """
    return np.load('splits/splittracks_' + str(qq) + '.npz')

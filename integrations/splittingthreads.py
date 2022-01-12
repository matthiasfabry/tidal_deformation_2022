"""
Main script for setting up the computations of the splitting surfaces of Sect 2.4
This uses multiprocess classes to parallelize the work.
"""

import multiprocessing
import pathlib
import time as t

import numpy as np

import c_rocheFunctions as cf
import c_split_tracer as ct
import l_p_points as rf


# %%
class SplittingThread(multiprocessing.Process):
    def __init__(self, thread_no, qqs, wd, overwrite=True):
        super().__init__()
        self.thread_no = thread_no
        self.qs = qqs
        self.wd = wd
        self.overwrite = overwrite
        self.starttime = t.time()
    
    def run(self):
        count = 0
        for q in self.qs:
            count += 1
            if self.overwrite and pathlib.Path('{}/splittracks_{}'.format(self.wd, q)).exists():
                print('q = {} already exists, skipping...')
                continue
            if count == 1:
                print('thread no {} starting with first iteration!'.format(self.thread_no))
            else:
                print('thread no {} doing iter {} of {}; estimated time remaining: {} h'
                      .format(self.thread_no, count, self.qs.size, np.round((t.time() - self.starttime) *
                                                                            (self.qs.size - count + 1)
                                                                            / ((count - 1) * 3600), 2)))
            
            x_l1 = rf.l1(q)  # lagrangian points
            x_lout = rf.lout(q)
            lout_pot = cf.roche(x_lout, 0, 0, q)
            l4_pot = cf.roche(0.5, 0.5 * np.sqrt(3), 0, q)
            stop = 0.2 * lout_pot + 0.8 * l4_pot  # furthest potential
            l1_pos = ct.c_gradient_tracer(x_l1, stop, q)  # heavy lifting
            lout_pos = ct.c_gradient_tracer(x_lout, stop, q)  # heavy lifting
            l1_bulk = np.array(l1_pos['bulk'])
            l1_edge = np.array(l1_pos['edge'])
            lout_bulk = np.array(lout_pos['bulk'])
            lout_edge = np.array(lout_pos['edge'])
            np.savez_compressed('{}/splittracks_{}'.format(self.wd, q), l1_bulk=l1_bulk, l1_edge=l1_edge,
                                l3_bulk=lout_bulk, l3_edge=lout_edge)  # save the surface points in a binary file


# we need to protect this module from being imported recursively by subprocesses
if __name__ == '__main__':
    
    folder = 'splits'
    # q array
    N_q = 20 * 14
    qs = np.logspace(-7, 7, N_q)
    qs = np.round(qs, 8)
    qs = np.append(qs, 1.0)
    pathlib.Path(folder).mkdir(exist_ok=True, parents=True)
    
    # multiprocessing params
    threadno = min(16, qs.size)
    threads = list()
    for i in range(threadno):  # build process objects
        threads.append(SplittingThread(i + 1, qs[i::threadno], folder))
    start = t.time()
    for thread in threads:
        thread.start()  # runs the process' computations
    
    for thread in threads:
        thread.join()
    print('tracing finished in', t.time() - start)

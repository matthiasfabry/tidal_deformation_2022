"""
module for setting up Roche integration computation of Sect 2.6
this uses multiprocess to parallelize the work
"""
import multiprocessing
import os
import pathlib
import time as t

import numpy as np

import c_main_tracer as ct
import c_rocheFunctions as cf
import c_util as cu
import interpolatesplits as inter
import l_p_points as lp
import loadsplitting as split
import mochnackicheckingfuncs as moch  # noqa

# global variables
folder = 'resultsmoch8'
overwrite = True
pointsfolder = folder + '/points'
drstop = 1e-12
phinum = 8 * 660
trim = 50

# setup qs to be sampled
Nqs = 14 * 20
qs = np.logspace(-7, 7, Nqs)
qs = np.round(qs, 8)
qs = np.append(qs, [1.0])
# qs = np.array([1.0])

# grid of theta, phi
dphi = np.pi / phinum
phis = np.linspace(dphi / 2, np.pi - dphi / 2, phinum)
thetanum = int(phinum / 2)
dcostheta = 1. / thetanum
thetas = np.arccos(np.linspace(dcostheta / 2, 1 - dcostheta / 2, thetanum))
THETAS, PHIS = np.meshgrid(thetas, phis)


class TracerProcess(multiprocessing.Process):
    """
    class for parallelizing tracing integrations
    """
    
    def __init__(self, thread_no, qqs, wd, ovrwrt=False):
        """
        constructs a tracing object.
        :param thread_no: name for this object
        :param qqs: list of mass ratios to integrate
        :param wd: working directory to save the results to
        :param ovrwrt: can
        """
        super().__init__()
        self.thread_no = thread_no
        self.qs = qqs
        self.wd = wd
        self.iterations = np.nan
        self.tick = t.time()
        self.overwrite = ovrwrt
    
    def run(self):
        count = 0
        for q in self.qs:
            # set up common items for this mass ratio q
            q = split.closest_q(q)
            xl1 = lp.l1(q)
            l1pot = cf.roche(xl1, 0, 0, q)
            l3pot = cf.roche(lp.lout(q), 0, 0, q)
            l4pot = cf.roche(0.5, 0.5 * np.sqrt(3), 0, q)
            bound = 0.5 * (l3pot + l4pot)
            
            # set up sampled potential values
            p1 = 0.9
            p2 = 0.99
            p3 = 0.999
            p4 = 0.9999
            p5 = 0.99999
            
            # pots = moch.psi_from_F(moch.Fs, q)  # if you want to compare to Mochnacki
            
            # all potential values we sample
            pots = np.concatenate((
                np.geomspace(50 * l1pot, 2 * l1pot, 30, endpoint=False),
                np.geomspace(2 * l1pot, (2 - p1) * l1pot, 10, endpoint=False),
                np.geomspace((2 - p1) * l1pot, (2 - p2) * l1pot, 10, endpoint=False),
                np.geomspace((2 - p2) * l1pot, (2 - p3) * l1pot, 10, endpoint=False),
                np.geomspace((2 - p3) * l1pot, (2 - p4) * l1pot, 10, endpoint=False),
                np.geomspace((2 - p4) * l1pot, (2 - p5) * l1pot, 10, endpoint=False),
                np.geomspace((2 - p5) * l1pot, l1pot, 10, endpoint=False),
                np.geomspace(l1pot, p4 * l1pot + (1 - p4) * l3pot, 10, endpoint=False),
                np.geomspace(p4 * l1pot + (1 - p4) * l3pot, p3 * l1pot + (1 - p3) * l3pot, 10, endpoint=False),
                np.geomspace(p3 * l1pot + (1 - p3) * l3pot, p2 * l1pot + (1 - p2) * l3pot, 10, endpoint=False),
                np.geomspace(p2 * l1pot + (1 - p2) * l3pot, p1 * l1pot + (1 - p1) * l3pot, 10, endpoint=False),
                np.geomspace(p1 * l1pot + (1 - p1) * l3pot, (1 - p1) * l1pot + p1 * l3pot, 20, endpoint=False),
                np.geomspace((1 - p1) * l1pot + p1 * l3pot, (1 - p2) * l1pot + p2 * l3pot, 10, endpoint=False),
                np.geomspace((1 - p2) * l1pot + p2 * l3pot, (1 - p3) * l1pot + p3 * l3pot, 10, endpoint=False),
                np.geomspace((1 - p3) * l1pot + p3 * l3pot, (1 - p4) * l1pot + p4 * l3pot, 10, endpoint=False),
                np.geomspace((1 - p4) * l1pot + p4 * l3pot, (1 - p5) * l1pot + p5 * l3pot, 10, endpoint=False),
                np.geomspace((1 - p5) * l1pot + p5 * l3pot, l3pot, 10, endpoint=False),
                np.geomspace(l3pot, p4 * l3pot + (1 - p4) * l4pot, 10, endpoint=False),
                np.geomspace(p4 * l3pot + (1 - p4) * l4pot, p3 * l3pot + (1 - p3) * l4pot, 10, endpoint=False),
                np.geomspace(p3 * l3pot + (1 - p3) * l4pot, p2 * l3pot + (1 - p2) * l4pot, 10, endpoint=False),
                np.geomspace(p2 * l3pot + (1 - p2) * l4pot, p1 * l3pot + (1 - p1) * l4pot, 10, endpoint=False),
                np.geomspace(p1 * l3pot + (1 - p1) * l4pot, bound, 10, endpoint=False),
            ))
            if l1pot not in pots:  # make sure you have l1 (geomspace can miss due to floating point errors)
                pots = np.append(pots, l1pot)
            
            # pots = np.array([l1pot])  # if you only want to do the RL case

            if not self.overwrite:
                skip = True
                ii = 0
                n = len(pots)
                while ii < n and skip:
                    skip = pathlib.Path('{}/q{}/super{}.txt'.format(self.wd, q, pots[ii] / l1pot)).exists()
                    ii += 1
                if skip:
                    print('q {}, all pots already exist, skipping...'.format(q))
                    continue
                
            # load the splitting surfaces for this q (must have exactly the same q value)
            splits = np.load('splits/splittracks_' + str(q) + '.npz')
            pos_tracks_l1 = cu.cart2spher(splits['l1_bulk'])
            l1_edge = cu.cart2spher(splits['l1_edge'])
            pos_tracks_l3 = cu.cart2spher(splits['l3_bulk'])
            l3_edge = cu.cart2spher(splits['l3_edge'])

            # interpolate the crossing rs on the used grid of thetas and phis
            try:
                l1inter, _ = inter.interpolate_split(pos_tracks_l1, l1_edge, thetas, phis, q,
                                                     THETAS=THETAS, PHIS=PHIS, try_cubic=False)
                l3inter, _ = inter.interpolate_split(pos_tracks_l3, l3_edge, thetas, phis, q,
                                                     THETAS=THETAS, PHIS=PHIS, try_cubic=False)
            except ValueError as e:
                print(e)
                count += pots.size
                continue
            
            self.iterations = self.qs.size * pots.size
            pathlib.Path(self.wd + '/q{}'.format(q)).mkdir(exist_ok=True, parents=True)
            # now go over all selected target equipotentials
            for targetpot in pots:
                count += 1
                ssuper = targetpot / l1pot
                if targetpot >= bound:
                    print('thread {}, iter {}, q = {}; potential above 1/2(L3+L4)! skipping'.format(
                        self.thread_no, count, q))
                    continue
                if not self.overwrite and pathlib.Path(
                        '{}/q{}/super{}.txt'.format(self.wd, q, ssuper)).exists():
                    print('q {}, super {}, already exists, skipping...'.format(q, ssuper))
                    continue
                if count == 1:
                    print('thread no {} starting with first iteration!'.format(self.thread_no))
                else:
                    print('thread no {} doing iter {} of {}; estimated time remaining: {} h'
                          .format(self.thread_no, count, self.iterations,
                                  np.round((t.time() - self.tick) * (self.iterations - count + 1)
                                           / 3600, 3)))
                self.tick = t.time()
                # heavy lifting, integrate in the roche geometry
                try:
                    _, vals1 = ct.c_roche_tracer(q, targetpot, xl1, thetas, phis,
                                                 l1_inter=l1inter, lout_inter=l3inter,
                                                 dr_stop=drstop)
                except ZeroDivisionError as e:
                    print(e, q, ssuper)
                    break
                vals1.insert(0, targetpot)
                vals1.insert(0, ssuper)
                vals1.insert(0, q)
                vals1 = np.array(vals1)
                # save the results
                np.savetxt('{}/q{}/super{}.txt'.format(self.wd, q, ssuper), vals1)
                # save a trimmed version of the points on the equipotential for future reference/plotting
                # np.savez_compressed('{}/points/primarypoints{}_{}'.format(self.wd, q, ssuper),
                #                     np.array(points1)[::trim])
            try:  # if empty, remove the q dir to avoid clutter
                pathlib.Path(self.wd + '/q{}'.format(q)).rmdir()
            except OSError:
                pass


# we need to protect this from being imported recursively by subprocesses
if __name__ == '__main__':
    pathlib.Path(pointsfolder).mkdir(exist_ok=True, parents=True)
    print(len(qs))
    
    # multiprocessing params
    procno = min(len(qs), os.cpu_count())
    procs = list()
    for i in range(procno):
        # make the processes
        procs.append(TracerProcess(i + 1, qs[i::procno], folder, ovrwrt=overwrite))
    
    # save some params for reference
    with open(folder + '/paramsextra.txt', 'a') as file:
        file.write('qs \t {}\n'.format(qs))
        file.write('phinum \t {}\n'.format(phinum))
        file.write('drstop \t {}\n'.format(drstop))
    
    # start the run!
    start = t.time()
    for proc in procs:
        proc.start()
    
    # gather processes
    for proc in procs:
        proc.join()
    print('tracing finished in', t.time() - start)

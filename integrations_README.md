Roche integrations
==================

This folder contains the scripts and code to reproduce the Roche geometry integrations results, 
such as those shown in Figs. 2, 3, 5.
Additionally the fitting procedures for the single rotating star done in appendix A are included 
here in `singlestarfitting.py`

Typical usage is as follows:
0. compile the `cython` code for your machine with `python setup.py build_ext --inplace`
1. play around in the `singlesplitting.py` and `singlevolume.py` scripts to get a feel for the 
   geometry, edge cases, maybe make cool peanut plots, etc...
2. compute all splitting surfaces for your choice of mass ratio sampling in `spittingthreads.py`
3. compute all integrations for your choice of equipotential values in `tracerthreads.py`
4. analyze and plot the results with `analyzeresults.py`
5. interpolate on a square grid to use results computed here in MESA with `squareinterpolate.py`

To compare with Mochnacki 1984:
1. use `tracerthreads.py` repeatedly with different resolutions, saving in different folders of 
   course
2. use `mochnackicompare.py` to load in these results and compare to Mochnacki's 1984 own 
   results, which are codified in `mochnackicheckingfuncs.py`
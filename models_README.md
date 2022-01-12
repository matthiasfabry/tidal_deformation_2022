Evolutionary models
===================

This folder contains all the files necessary to reproduce the results presented in Sect. 5.
It contains MESA working directories that can be used in conjunction with a MESA installation 
for the three situation explored:
1. `singlerot` for the single star rotation case (Sect. 5.2)
2. `semidetached` for the semidetached eclipsing binary-like case (Sect. 5.3)
3. `overcontact` for the twin binary that evolves into contact (VFTS 352-like) (Sect. 5.4)

In each of these folders, subfolders indicate what boundary conditions and what deformation 
model is chosen, ie `newBC` are the new boundary conditions developed in Sect. 4, wheras `oldBC` 
are the default MESA boundary conditions.
The distortion model is either `norot` (star has no angular momentum and is thus perfectly 
spherical), `single` which uses the single rotating star corrections of Appendix A and `tidal`, 
which uses the corrections calculated with the integrations in the Roche geometry

the files `area_data.txt, fp_data.txt, ft_data.txt` and `irot_data.txt` contain grids of 
equipotential surface area, fp, ft correction and moment of inertia respectively as function of 
mass ratio and fractional roche lobe size. They are the direct output of `squareinterpolate.py` 
of the `integrations` section of this repository.
These files are then read in here in `tidal`ly distorted models in `run_star_extras.f90` to 
compute the required corrections to the stellar structure equations and boundary conditions (if 
`newBC`)

`model_plots.py` allows to recreate Figs. 6-9 of Sect. 5 (`mesadata.py` is a helper file to read 
in the logs from MESA)

f-kriging
=============================

This is a fortran code to perform the Cokriging method, Kriging method, for a predetermined
grid or for a grid constructed via adaptive sequential sampling. 
 
Installation Prerequesites
-----------------------------

This setup requires the following libraries to be installed:
* `BLAS`
* `LAPACK`
* `LAPACK95` 
* `FortArrange` (can be found on github @ https://github.com/arinrb/FortArrange)

The makefile is setup for ubuntu with the libraries in the normal path 
plus the following modules folders : 
* `/usr/lib/lapack95_modules/`
* `/usr/lib/fortarrange_header/`


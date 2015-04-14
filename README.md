f-kriging
=============================
### NOT MAINTAINED. OLD CODE. I'M NOT PROUD OF THIS.

This is a fortran code to perform the Cokriging method, Kriging method, for a predetermined
grid or for a grid constructed via adaptive sequential sampling. I am also currently working 
on coupling it with an optimization scheme. 
 
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

Current Goals
------------------------------
* Investigate optimization of objective functions with this indirect cokriging solver. 
* Provide an easy interface to the code
* Turn this whole thing into a library
* Document it

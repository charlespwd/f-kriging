FortLinesearch
=====================
This is a simple module that performs a linesearch on a multivariable double precision 
function with a sufficient decrease and backtracking step algorithm and a gradient 
calculated with a central finite difference scheme.

Requirements
---------------------
The BLAS library must be linked to the project

Module Globals
---------------------
* `iimax`, maximum number of iteration to converge to a solution. defaults to 100
* `alpha_0`, the initial alpha chosen by the backtracking algorithm. defaults to 10
* `c_1`, the wolve condition constant. defaults to 10^-4
* `tol`, the linesearch tolerance. A minimum is assumed to be found when sum(abs(grad)) is
smaller than this one. defaults to 10^-8.
* `h`, the stepsize used in the central finite difference. defaults to 10^-10

Usage
---------------------
`linesearch(x_star, x_initial, f, D)`
* `x_star`, `(OUT)`, is the minimum
* `x_initial`, `(IN)`, is the initial guess
* `f`, description below
* `D`, `(IN)` the dimension of the domain

`f` must be a pure double precision function with the following interface:

		interface func
			double precision pure function f(x,D)
				integer, intent(in) :: D
				double precision, intent(in) :: x(D,1)	
			end function
		end interface func

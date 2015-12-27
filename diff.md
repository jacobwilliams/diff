project: diff
project_dir: ./src
output_dir: ./doc
project_github: https://github.com/jacobwilliams/diff
summary: Numerical Differentiation of a User Defined Function
author: Jacob Williams
github: https://github.com/jacobwilliams
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
display: private
display: protected
source: true
graph: true
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html

Brief description
---------------

This code is a modern Fortran update of the DIFF subroutine found [here](ftp://math.nist.gov/pub/repository/diff/).

The DIFF subroutine computes the first, second or third derivative of a real function of a single real variable.  The user provides the function, a real interval `[xmin,xmax]` on which the function is continuous, and a point `x0` lying in `[xmin,xmax]`. Optionally, the user may provide an estimate of the absolute or relative accuracy of function evaluation in `[xmin,xmax]` and the absolute or relative error tolerance that is acceptable in the computed derivative. The subroutine returns the computed derivative and an estimated upper bound on the absolute error of the computed derivative. The method used is Neville's process of extrapolating from a sequence of interpolating polynomials with interpolating points distributed symmetrically about `x0` or, if this is not possible, to one side of `x0`.

The original sourcecode was produced by the National Bureau of Standards, and is presumed to be in the public domain.

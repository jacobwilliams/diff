diff
============

Numerical Differentiation of a User Defined Function

## Overview

This code is a modern Fortran update of the DIFF subroutine found [here](ftp://math.nist.gov/pub/repository/diff/).

The DIFF subroutine computes the first, second or third derivative of a real function of a single real variable.  The user provides the function, a real interval `[xmin,xmax]` on which the function is continuous, and a point `x0` lying in `[xmin,xmax]`. Optionally, the user may provide an estimate of the absolute or relative accuracy of function evaluation in `[xmin,xmax]` and the absolute or relative error tolerance that is acceptable in the computed derivative. The subroutine returns the computed derivative and an estimated upper bound on the absolute error of the computed derivative. The method used is Neville's process of extrapolating from a sequence of interpolating polynomials with interpolating points distributed symmetrically about `x0` or, if this is not possible, to one side of `x0`.

The original sourcecode was produced by the National Bureau of Standards, and is presumed to be in the public domain.

## Example Useage

```fortran
    program example

    use diff_module
    use iso_fortran_env, only: wp => real64 !use double precision

    implicit none

    integer,parameter  :: iord  = 1
    real(wp),parameter :: x0    = 0.12345_wp
    real(wp),parameter :: xmin  = 0.0_wp
    real(wp),parameter :: xmax  = 1.0_wp
    real(wp),parameter :: eps   = 1.0e-9_wp
    real(wp),parameter :: acc   = 0.0_wp

    real(wp) :: deriv, error
    integer  :: ifail
    type(diff_func) :: d

    call d%set_function(sin_func) !set function
    call d%compute_derivative(iord,x0,xmin,xmax,eps,acc,deriv,error,ifail)

    write(*,'(A)') ''
    write(*,'(A,I5)')     'ifail                :', ifail
    write(*,'(A,E25.16)') 'estimated derivative :', deriv
    write(*,'(A,E25.16)') 'actual derivative    :', cos(x0)
    write(*,'(A,E25.16)') 'estimated error      :', error
    write(*,'(A,E25.16)') 'actual error         :', cos(x0) - deriv
    write(*,'(A)') ''

    contains

        function sin_func(me,x) result(fx)

        implicit none

        class(diff_func),intent(inout) :: me
        real(wp),intent(in) :: x
        real(wp) :: fx

        fx = sin(x)

        end function sin_func

    end program example
```

Which produces:

```
ifail                :   0
estimated derivative :   0.9923897210998529E+00
actual derivative    :   0.9923897211114882E+00
estimated error      :   0.5805992701806732E-09
actual error         :   0.1163524832037410E-10
```

## Notes

Alan Miller's [to_f90](http://jblevins.org/mirror/amiller/to_f90.f90) program was used to assist in the conversion to modern Fortran.

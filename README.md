diff
============

Numerical Differentiation of a User Defined Function

Reference
---------------

This code is a modern Fortran update of the diff subroutine found here:

<ftp://math.nist.gov/pub/repository/diff/src/DIFF>

Example
---------------

	program example

	use diff_module
	
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
	
	d%f => sin_func  !set function
	
	call d%compute_derivative(iord,x0,xmin,xmax,eps,acc,deriv,error,ifail)
	
	write(*,'(A)') ''
	write(*,'(A,I5)')     'ifail                :', ifail
	write(*,'(A,E25.16)') 'estimated derivative :', deriv
	write(*,'(A,E25.16)') 'actual derivative    :', cos(x0)
	write(*,'(A,E25.16)') 'estimated error      :', error
	write(*,'(A,E25.16)') 'actual error         :', cos(x0) - deriv
	write(*,'(A)') ''
	
	contains
			
	!***********************************************************
		function sin_func(me,x) result(fx)
	!***********************************************************
	
		implicit none
	
		class(diff_func),intent(inout) :: me
		real(wp),intent(in) :: x
		real(wp) :: fx
				
		fx = sin(x)
				
	!***********************************************************
		end function sin_func
	!***********************************************************
				
	end program example

Notes
---------------

Alan Miller's to_f90 program was used to assist in the conversion to modern Fortran: 
<http://jblevins.org/mirror/amiller/to_f90.f90>



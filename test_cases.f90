!*****************************************************************************************
	program test_cases
!*****************************************************************************************
!
!	Test cases for the diff_module
!
!*****************************************************************************************
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
	
	type,extends(diff_func) :: my_func
		integer :: ifunc = 0
	end type my_func
	type(my_func) :: d
	
	d%ifunc = 0
	d%f => sin_func
	call d%compute_derivative(iord,x0,xmin,xmax,eps,acc,deriv,error,ifail)
	call write_results(cos(x0))

	d%ifunc = 0
	d%f => test_func_1
	call d%compute_derivative(iord,x0,xmin,xmax,eps,acc,deriv,error,ifail)
	call write_results(cos(x0) - sin(x0) + 2.0_wp*x0)

	d%ifunc = 0
	d%f => test_func_2
	call d%compute_derivative(iord,x0,xmin,xmax,eps,acc,deriv,error,ifail)
	call write_results(6.0_wp*x0**5 + 5.0_wp*x0**4 + 4.0_wp*x0**3 + 2.0_wp*x0 + 1.0_wp)
	
	contains
!*****************************************************************************************
	
	!***********************************************************
		subroutine write_results(truth)
	!***********************************************************
		
		implicit none
		
		real(wp),intent(in) :: truth
		
		write(*,'(A)') ''
		write(*,'(A,E25.16)') 'solution          :', deriv
		write(*,'(A,E25.16)') 'actual derivative :', truth
		write(*,'(A,E25.16)') 'estimated error   :', error
		write(*,'(A,E25.16)') 'actual error      :', truth - deriv
		write(*,'(A,I5)')     'ifail             :', ifail
		write(*,'(A,I5)')     'func evaluations  :', d%ifunc
		write(*,'(A)') ''
		
	!***********************************************************
		end subroutine write_results
	!***********************************************************
		
	!***********************************************************
		function sin_func(me,x) result(fx)
	!***********************************************************
	
		implicit none
	
		class(diff_func),intent(inout) :: me
		real(wp),intent(in) :: x
		real(wp) :: fx
		
		select type (me)
		type is (my_func)
			me%ifunc = me%ifunc + 1
		end select
		
		fx = sin(x)
				
	!***********************************************************
		end function sin_func
	!***********************************************************
		
	!***********************************************************
		function test_func_1(me,x) result(fx)
	!***********************************************************
	
		implicit none
	
		class(diff_func),intent(inout) :: me
		real(wp),intent(in) :: x
		real(wp) :: fx
		
		select type (me)
		type is (my_func)
			me%ifunc = me%ifunc + 1
		end select

		fx = sin(x) + cos(x) + x**2
				
	!***********************************************************
		end function test_func_1		
	!***********************************************************

	!***********************************************************
		function test_func_2(me,x) result(fx)
	!***********************************************************
	
		implicit none
	
		class(diff_func),intent(inout) :: me
		real(wp),intent(in) :: x
		real(wp) :: fx
		
		select type (me)
		type is (my_func)
			me%ifunc = me%ifunc + 1
		end select
		
		fx = x**6 + x**5 + x**4 + x**2 + x
				
	!***********************************************************
		end function test_func_2
	!***********************************************************
				
!*****************************************************************************************
	end program test_cases
!*****************************************************************************************
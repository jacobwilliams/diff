!*****************************************************************************************
!> author: Jacob Williams
!  license: BSD
!
!  The example in the readme file.

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
!*****************************************************************************************
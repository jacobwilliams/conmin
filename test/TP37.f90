program tp37

!! this program executes the example problem from
!! https://mdolab-pyoptsparse.readthedocs-hosted.com/en/latest/quickstart.html
!!
!! This is Schittkowski TP37 test problem.
!! see also: https://klaus-schittkowski.de/test_problems.pdf

use conmin_module, only: conmin_class, wp => conmin_wp
use pyplot_module

implicit none

integer,parameter :: n = 3  !! 3 optimization variables (ndv)
integer,parameter :: m = 2 !! 2 constraint functions (ncom)

integer,parameter :: n1 = n + 2
integer,parameter :: n2 = m + 2*n
integer,parameter :: n3 = n + 2
integer,parameter :: n4 = max(n3,n)
integer,parameter :: n5 = 2*n4

real (wp) :: s(n1), g1(n2), g2(n2), b(n3,n3), c(n4), vlb(n1), vub(n1), &
             scal(n1), df(n1), a(n1,n3), x(n1), g(n2)
integer   :: ms1(n5), isc(n2), ic(n3)
integer   :: i, nlim
type(conmin_class) :: solver
type(pyplot) :: plt

real(wp),dimension(:),allocatable :: iter_hist
real(wp),dimension(:),allocatable :: g1_hist
real(wp),dimension(:),allocatable :: g2_hist
real(wp),dimension(:),allocatable :: obj_hist

open (newunit=solver%iunit,file='tp37.txt',status='REPLACE')

! defaults:
solver%infog  = 0
solver%info   = 0
solver%icndir = 0
solver%itrm   = 0
solver%nscal  = 0
solver%ctl    = 0.0_wp
solver%ctlmin = 0.0_wp
solver%theta  = 0.0_wp
solver%phi    = 0.0_wp
solver%linobj = 0.0_wp
solver%alphax = 0.0_wp
solver%abobj1 = 0.0_wp
solver%ctl    = 0.0_wp

! modified values from default:
! solver%nfdg   = 0           ! use finite diff gradients
solver%nfdg   = 1              ! user-computed gradients
solver%fdch   = 1.0e-4_wp
solver%fdchm  = 1.0e-5_wp
solver%iprint = 4
solver%itmax  = 100          ! maximum number of iterations
solver%ct     = 1.0e-5_wp    ! Constraint thickness parameter
solver%ctmin  = 1.0e-6_wp    ! constraint acceptable tolerance
solver%delfun = 1.0e-12_wp   ! objective relative tolerance
solver%dabfun = 1.0e-13_wp   ! objective absolute tolerance
solver%ncon   = m            ! number of constraints
solver%ndv    = n            ! number of optimization variables
solver%nside  = 1            ! upper and lower bounds
do i = 1, n
    x(i) = 10.0_wp   ! initial values
    vlb(i) = 0.0_wp  ! lower bounds
    vub(i) = 42.0_wp ! upper bounds
    scal(i) = 1.0_wp ! scale  - not used if nscal = 0
end do
isc(1:m) = 0
nlim = solver%itmax * (n+5)

! non-iterative part of analysis
solver%igoto = 0

solver%report => report  ! to report iterations

! iterative part of analysis
do  i = 1, nlim
    ! call the optimization routine conmin
    call solver%solve(x,vlb,vub,g,scal,df,a,s,g1,g2,b,c,isc,ic,ms1,n1,n2,n3,n4,n5)

    ! info = 1:  calculate obj and g(i), i = 1, ncon
    !
    ! info = 2:  calculate nac, ic(i), i = 1,nac, the gradient of obj, and
    !            the gradient of g(j), where j = ic(i), i = 1,nac.
    !            store the gradients of g in columns of a.

    !...JW : need to clean up the docstrings to make sure
    !        they are right. i think some of them are for
    !        the old code.

    select case (solver%info)
    case (1)
        ! objective function
        solver%obj = -x(1) * x(2) * x(3)
        ! constraint values
        g(1) =  x(1) + 2.0_wp * x(2) + 2.0_wp * x(3) - 72.0_wp
        g(2) = -x(1) - 2.0_wp * x(2) - 2.0_wp * x(3)

    case (2)
        ! objective gradient:
        df(1) = -x(2) * x(3)
        df(2) = -x(1) * x(3)
        df(3) = -x(1) * x(2)

        ! constraint gradient:
        ! [note that they are stored by columns]
        solver%nac = 0
        if (g(1) >= solver%ct) then
            solver%nac = 1
            ic(1) = 1
            a(1,1) = 1.0_wp
            a(2,1) = 2.0_wp
            a(3,1) = 2.0_wp
        end if
        if (g(2) >= solver%ct) then
            solver%nac = solver%nac + 1
            ic(solver%nac) = 2
            a(1,2) = -1.0_wp
            a(2,2) = -2.0_wp
            a(3,2) = -2.0_wp
        end if

    end select

    if (solver%igoto == 0) exit
end do

close (solver%iunit)

! make a plot of the constraint iteration history
call plt%initialize(grid=.true.,xlabel='Iteration',ylabel='Constraint Violation',&
                    title='TP37 Constraints',legend=.true.)
call plt%add_plot(iter_hist,g1_hist,label='g(1)',linestyle='b.-',markersize=5,linewidth=2)
call plt%add_plot(iter_hist,g2_hist,label='g(2)',linestyle='r.-',markersize=5,linewidth=2)
call plt%savefig('iterations.png')

call plt%initialize(grid=.true.,xlabel='Iteration',ylabel='Objective Function Value',&
                    title='TP37 Objective Function',legend=.true.)
call plt%add_plot(iter_hist,obj_hist,label='obj',linestyle='k.-',markersize=5,linewidth=2)
call plt%savefig('obj.png')

contains

    subroutine report(me, iter, x, obj, g)
        !! example iteration reporting function
        class(conmin_class), intent(inout) :: me
        integer, intent(in)                :: iter !! Iteration number
        real(wp),dimension(:), intent(in)  :: x    !! Optimization variables
        real(wp),intent(in)                :: obj  !! Objective function value
        real(wp),dimension(:),intent(in)   :: g    !! Constraint functions

        write(*,'(a,i2,a,i1,a,i1,a,f26.12,1x,a,3(f20.6),1x,a,2(f20.6))') &
            '(', iter, ',', solver%info, ',', solver%igoto, ') '// &
            'j:', obj, &
            'x:', x(1), x(2), x(3), &
            'g:', g(1), g(2)

        ! for the plot:
        if (.not. allocated(iter_hist)) then
            allocate(iter_hist(0))
            allocate(g1_hist(0))
            allocate(g2_hist(0))
            allocate(obj_hist(0))
        end if
        iter_hist = [iter_hist, real(iter,wp)]
        g1_hist = [g1_hist, g(1)]
        g2_hist = [g2_hist, g(2)]
        obj_hist = [obj_hist, obj]
    end subroutine report

end program tp37
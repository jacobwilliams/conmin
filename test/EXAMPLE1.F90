PROGRAM exampl1

! This program executes the example problem one of the conmin manual.

USE conmin_module, only: conmin_class, wp => conmin_wp

IMPLICIT NONE

REAL (wp) :: s(6), g1(11), g2(11), b(11,11), c(11), vlb(6), vub(6), &
             scal(6), df(6), a(6,11), aobj, x(6), g(11)
INTEGER   :: ms1(22), isc(11), ic(11)
INTEGER   :: i, n1, n2, n3, n4, n5, nlim
type(conmin_class) :: solver

OPEN (newunit=solver%iunit,FILE='EXOUT1.TXT',STATUS='REPLACE')

!  INITIALIZE
solver%infog = 0
solver%info = 0
solver%nfdg = 0
solver%iprint = 2
solver%ndv = 4
solver%itmax = 40
solver%ncon = 3
solver%nside = 0
solver%icndir = 0
solver%nscal = 0
solver%fdch = 0.0_wp
solver%fdchm = 0.0_wp
solver%ct = 0.0_wp
solver%ctmin = 0.0_wp
solver%ctl = 0.0_wp
solver%ctlmin = 0.0_wp
solver%theta = 0.0_wp
solver%phi = 0.0_wp
solver%delfun = 0.0_wp
solver%dabfun = 0.0_wp
solver%linobj = 0.0_wp
solver%itrm = 0
n1 = 6
n2 = 11
n3 = 11
n4 = 11
n5 = 22
solver%alphax = 0.0_wp
solver%abobj1 = 0.0_wp
solver%ctl = 0.0_wp
DO  i = 1, solver%ndv
  x(i) = 1.0_wp
  vlb(i) = -99999.0_wp
  vub(i) = 99999.0_wp
END DO

isc(1:solver%ncon) = 0
nlim = solver%itmax * (solver%ndv+5)

! NON-ITERATIVE PART OF ANALYSIS
solver%igoto = 0

! ITERATIVE PART OF ANALYSIS

DO  i = 1, nlim
  ! CALL THE OPTIMIZATION ROUTINE CONMIN
  CALL solver%solve(x,vlb,vub,g,scal,df,a,s,g1,g2,b,c,isc,ic,ms1,n1,n2,n3,n4,n5)
  ! ANALYSIS MODULE
  CALL analys()
  solver%obj = aobj
  IF (solver%igoto == 0) EXIT
END DO

CLOSE (solver%iunit)

CONTAINS

  SUBROUTINE analys()

    !! Routine to calculate objective function and constraints

    ! OBJECTIVE FUNCTION
    aobj = x(1) ** 2 - 5.0_wp * x(1) + x(2) ** 2 - 5.0_wp * x(2) + 2.0_wp * x(3) ** 2  &
          - 21.0_wp * x(3) + x(4) ** 2 + 7.0_wp * x(4) + 50.0_wp

    ! CONSTRAINT VALUES
    g(1) = x(1) ** 2 + x(1) + x(2) ** 2 - x(2) + x(3) ** 2 + x(3) +  &
          x(4) ** 2 - x(4) - 8.0_wp

    g(2) = x(1) ** 2 - x(1) + 2.0_wp * x(2) ** 2 + x(3) ** 2 + 2.0_wp * x(4) ** 2 &
          - x(4) - 10.0_wp

    g(3) = 2.0_wp * x(1) ** 2 + 2.0_wp * x(1) + x(2) ** 2 - x(2) + x(3) ** 2 - x(4) - 5.0_wp

  END SUBROUTINE analys

END PROGRAM exampl1
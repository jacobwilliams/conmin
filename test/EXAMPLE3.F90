PROGRAM exampl3

  !! THIS PROGRAM EXECUTES THE EXAMPLE PROBLEM THREE OF THE CONMIN MANUAL

  USE conmin_module, only: conmin_class, wp => conmin_wp

  IMPLICIT NONE

  REAL (wp)  :: s(4), g1(10), g2(10), b(10,10), c(10), vlb(4), vub(4),  &
                scal(4), df(4), a(4,10)
  INTEGER    :: ms1(20), isc(10), ic(10)
  REAL (wp)  :: aobj, x(6), g(11)
  INTEGER    :: i, n1, n2, n3, n4, n5, nlim
  type(conmin_class) :: solver

OPEN (UNIT=6,FILE='EXOUT3.TXT',STATUS='REPLACE')

!  INITIALIZE
solver%infog = 0
solver%info = 0
solver%nfdg = 2
solver%iprint = 1
solver%ndv = 2
solver%itmax = 40
solver%ncon = 6
solver%nside = 1
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
solver%linobj = 1
solver%itrm = 0
n1 = 4
n2 = 10
n3 = 10
n4 = 10
n5 = 20
solver%alphax = 0.0_wp
solver%abobj1 = 0.0_wp
solver%ctl = 0.0_wp
DO  i = 1, solver%ndv
  x(i) = 1.0_wp
  vlb(i) = 0.001_wp
  vub(i) = 1.0e+10_wp
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

CLOSE (6)

CONTAINS

SUBROUTINE analys()

  !! ROUTINE TO CALCULATE OBJECTIVE FUNCTION AND CONSTRAINT VALUES
  !! FOR OPTIMIZATION OF CONSTRAINED ROSEN-SUZUKI FUNCTION.

  REAL (wp) :: a1, a2, denom, rho, sig11, sig21, sig31

  rho = 0.1_wp
  a1 = x(1)
  a2 = x(2)

  IF (solver%info < 2) THEN

    ! OBJECTIVE FUNCTION
    aobj = 10.0_wp * rho * (2.0_wp*SQRT(2.0_wp)*a1+a2)

    ! CONSTRAINT VALUES
    denom = 2.0_wp * a1 * a2 + SQRT(2.0_wp) * a1 * a1
    sig11 = 20.0_wp * (SQRT(2.0_wp)*a1+a2) / denom
    sig21 = 20.0_wp * SQRT(2.0_wp) * a1 / denom
    sig31 = -20.0_wp * a2 / denom

    g(1) = -sig11 / 15.0_wp - 1.0_wp
    g(2) =  sig11 / 20.0_wp - 1.0_wp
    g(3) = -sig21 / 15.0_wp - 1.0_wp
    g(4) =  sig21 / 20.0_wp - 1.0_wp
    g(5) = -sig31 / 15.0_wp - 1.0_wp
    g(6) =  sig31 / 20.0_wp - 1.0_wp

  ELSE
    ! GRADIENT INFORMATION
    df(1) = 20.0_wp * SQRT(2.0_wp) * rho
    df(2) = 10.0_wp * rho
  END IF

  !  GRADIENTS OF ACTIVE AND VIOLATED CONSTRAINTS
  !  WILL BE CALCULATED BY FINITE DIFFERENCE WITHIN
  !  CONMIN

END SUBROUTINE analys

END PROGRAM exampl3


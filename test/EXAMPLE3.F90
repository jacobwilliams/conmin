PROGRAM exampl3
!
!     THIS PROGRAM EXECUTES THE EXAMPLE PROBLEM THREE OF THE CONMIN MANUAL

  USE conmin_module, only: conmin_class

  IMPLICIT NONE

INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp)  :: s(4), g1(10), g2(10), b(10,10), c(10), vlb(4), vub(4),  &
              scal(4), df(4), a(4,10)
INTEGER    :: ms1(20), isc(10), ic(10)
type(conmin_class) :: solver

! COMMON /varable/ aobj, x(6), g(11)
REAL (dp)  :: aobj, x(6), g(11)

! COMMON /andata/ loopcnt

! COMMON /grad/ isc(11), ic(11), df(6), a(6,11)

! COMMON /consav/ rnum(50), inum(25)

INTEGER    :: i, n1, n2, n3, n4, n5, nlim
! REAL (dp)  :: phi

! NAMELIST /conpar/ infog, info, nfdg, iprint, ndv, itmax, ncon,  &
!     nside, icndir, nscal, fdch, fdchm, ct, ctmin, ctlmin, theta,  &
!     phi, delfun, dabfun, linobj, itrm, x, vlb, vub, n1, n2, n3,  &
!     n4, n5, alphax, abobj1, ctl, isc, scal
! !
OPEN (UNIT=6,FILE='EXOUT3.TXT',STATUS='REPLACE')
!
!
!  INITIALIZE
!
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
solver%fdch = 0.0
solver%fdchm = 0.0
solver%ct = 0.0
solver%ctmin = 0.0
solver%ctl = 0.0
solver%ctlmin = 0.0
solver%theta = 0.0
solver%phi = 0.0
solver%delfun = 0.0
solver%dabfun = 0.0
solver%linobj = 1
solver%itrm = 0
n1 = 4
n2 = 10
n3 = 10
n4 = 10
n5 = 20
solver%alphax = 0.0
solver%abobj1 = 0.0
solver%ctl = 0.0
DO  i = 1, solver%ndv
  x(i) = 1.0
  vlb(i) = 0.001
  vub(i) = 1.0E+10
END DO
!
isc(1:solver%ncon) = 0
!
!     READ THE PARAMETERS FOR CONMIN
!
!CC   READ(5,CONPAR)      USE DEFAULT VALUES
! WRITE (6,conpar)
nlim = solver%itmax * (solver%ndv+5)
!
!     NON-ITERATIVE PART OF ANALYSIS
!
solver%igoto = 0
!
!     ITERATIVE PART OF ANALYSIS
!
DO  i = 1, nlim
!
!       CALL THE OPTIMIZATION ROUTINE CONMIN
!
  CALL solver%solve(x,vlb,vub,g,scal,df,a,s,g1,g2,b,c,isc,ic,ms1,n1,n2,n3,n4,n5)
!
!
!       ANALYSIS MODULE
!
  CALL analys()
  solver%obj = aobj
  IF (solver%igoto == 0) EXIT
END DO
!
!
CLOSE (6)


CONTAINS


SUBROUTINE analys()
!
!   ROUTINE TO CALCULATE OBJECTIVE FUNCTION AND CONSTRAINT VALUES
!   FOR OPTIMIZATION OF CONSTRAINED ROSEN-SUZUKI FUNCTION.
!
!
REAL (dp)  :: a1, a2, denom, rho, sig11, sig21, sig31

rho = 0.1
a1 = x(1)
a2 = x(2)
!
IF (solver%info < 2) THEN
!
!  OBJECTIVE FUNCTION
!
  aobj = 10. * rho * (2.*SQRT(2.)*a1+a2)
!
!
!   CONSTRAINT VALUES
!
  denom = 2. * a1 * a2 + SQRT(2.) * a1 * a1
  sig11 = 20. * (SQRT(2.)*a1+a2) / denom
  sig21 = 20. * SQRT(2.) * a1 / denom
  sig31 = -20. * a2 / denom
!
  g(1) = -sig11 / 15. - 1.
  g(2) = sig11 / 20. - 1.
  g(3) = -sig21 / 15. - 1.
  g(4) = sig21 / 20. - 1.
  g(5) = -sig31 / 15. - 1.
  g(6) = sig31 / 20. - 1.
!
ELSE
!
!
!    GRADIENT INFORMATION
!
  df(1) = 20. * SQRT(2.) * rho
  df(2) = 10. * rho
END IF
!
!  GRADIENTS OF ACTIVE AND VIOLATED CONSTRAINTS
!  WILL BE CALCULATED BY FINITE DIFFERENCE WITHIN
!  CONMIN
!
!
END SUBROUTINE analys

END PROGRAM exampl3


PROGRAM exampl1
!
!     THIS PROGRAM EXECUTES THE EXAMPLE PROBLEM ONE OF THE CONMIN MANUAL.

USE conmin_module, only: conmin_class
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp)  :: s(6), g1(11), g2(11), b(11,11), c(11), vlb(6), vub(6),  &
              scal(6), df(6), a(6,11)
INTEGER    :: ms1(22), isc(11), ic(11)
type(conmin_class) :: solver

! COMMON /varable/ aobj, x(6), g(11)
REAL (dp)  :: aobj, x(6), g(11)

! COMMON /andata/ loopcnt

INTEGER    :: i, n1, n2, n3, n4, n5, nlim
! REAL (dp)  :: phi

! NAMELIST /conpar/ infog, info, nfdg, iprint, ndv, itmax, ncon,  &
!     nside, icndir, nscal, fdch, fdchm, ct, ctmin, ctlmin, theta,  &
!     phi, delfun, dabfun, linobj, itrm, x, vlb, vub, n1, n2, n3,  &
!     n4, n5, alphax, abobj1, ctl, isc, scal

OPEN (UNIT=6,FILE='EXOUT1.TXT',STATUS='REPLACE')

!
!
!  INITIALIZE
!
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
solver%linobj = 0.0
solver%itrm = 0
n1 = 6
n2 = 11
n3 = 11
n4 = 11
n5 = 22
solver%alphax = 0.0
solver%abobj1 = 0.0
solver%ctl = 0.0
DO  i = 1, solver%ndv
  x(i) = 1.0_dp
  vlb(i) = -99999.0_dp
  vub(i) = 99999.0_dp
END DO
!
isc(1:solver%ncon) = 0
!
!     READ THE PARAMETERS FOR CONMIN
!
!CC   READ(5,CONPAR)                  USE DEFAULT VALUES
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
!       ANALYSIS MODULE
!
  CALL analys()
  solver%obj = aobj
  IF (solver%igoto == 0) EXIT
END DO

CLOSE (6)

CONTAINS



SUBROUTINE analys()
!
!   ROUTINE TO CALCULATE OBJECTIVE FUNCTION AND
!   CONSTRAINTS
!
!
!  OBJECTIVE FUNCTION
!
aobj = x(1) ** 2 - 5. * x(1) + x(2) ** 2 - 5. * x(2) + 2. * x(3) ** 2  &
       - 21. * x(3) + x(4) ** 2 + 7.0 * x(4) + 50.
!
!
!   CONSTRAINT VALUES
!
g(1) = x(1) ** 2 + x(1) + x(2) ** 2 - x(2) + x(3) ** 2 + x(3) +  &
       x(4) ** 2 - x(4) - 8.0
!
g(2) = x(1) ** 2 - x(1) + 2. * x(2) ** 2 + x(3) ** 2 + 2. * x(4) ** 2 &
     - x(4) - 10.0
!
g(3) = 2. * x(1) ** 2 + 2. * x(1) + x(2) ** 2 - x(2) + x(3) ** 2 - x(4) - 5.0
!
END SUBROUTINE analys

END PROGRAM exampl1


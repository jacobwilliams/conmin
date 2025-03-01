PROGRAM exampl2
!
!     THIS PROGRAM EXECUTES THE EXAMPLE PROBLEM TWO OF THE CONMIN MANUAL

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

! COMMON /grad/ isc(11), ic(11), df(6), a(6,11)

INTEGER    :: i, n1, n2, n3, n4, n5, nlim
! REAL (dp)  :: phi

! NAMELIST /conpar/ infog, info, nfdg, iprint, ndv, itmax, ncon,  &
!     nside, icndir, nscal, fdch, fdchm, ct, ctmin, ctlmin, theta,  &
!     phi, delfun, dabfun, linobj, itrm, x, vlb, vub, n1, n2, n3,  &
!     n4, n5, alphax, abobj1, ctl, isc, scal

OPEN (UNIT=6,FILE='EXOUT2.TXT',STATUS='REPLACE')
!
!
!  INITIALIZE
!
solver%infog = 0
solver%info = 0
solver%nfdg = 1
solver%iprint = 1
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
  x(i) = 1.0
  vlb(i) = -99999.
  vub(i) = 99999.
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

CONTAINS



SUBROUTINE analys()
!
!   ROUTINE TO CALCULATE OBJECTIVE FUNCTION AND CONSTRAINT VALUES
!   FOR OPTIMIZATION OF CONSTRAINED ROSEN-SUZUKI FUNCTION.
!
!
IF (solver%info < 2) THEN
!
!  OBJECTIVE FUNCTION
!
  aobj = x(1) ** 2 - 5. * x(1) + x(2) ** 2 - 5. * x(2) + 2. *  &
         x(3) ** 2 - 21. * x(3) + x(4) ** 2 + 7.0 * x(4) + 50.
!
!
!   CONSTRAINT VALUES
!
  g(1) = x(1) ** 2 + x(1) + x(2) ** 2 - x(2) + x(3) ** 2 + x(3) +  &
         x(4) ** 2 - x(4) - 8.0
!
  g(2) = x(1) ** 2 - x(1) + 2. * x(2) ** 2 + x(3) ** 2 + 2. *  &
         x(4) ** 2 - x(4) - 10.0
!
  g(3) = 2. * x(1) ** 2 + 2. * x(1) + x(2) ** 2 - x(2) + x(3) ** 2  &
         - x(4) - 5.0
!
ELSE
!
!
!    GRADIENT INFORMATION
!
  df(1) = 2.0 * x(1) - 5.0
  df(2) = 2.0 * x(2) - 5.0
  df(3) = 4.0 * x(3) - 21.
  df(4) = 2.0 * x(4) + 7.
!
!  GRADIENTS OF ACTIVE AND VIOLATED CONSTRAINTS
!
  solver%nac = 0
  IF (g(1) >= solver%ct) THEN
    solver%nac = 1
    ic(1) = 1
    a(1,1) = 2. * x(1) + 1.
    a(2,1) = 2. * x(2) - 1.
    a(3,1) = 2. * x(3) + 1.
    a(4,1) = 2. * x(4) - 1.
  END IF
!
  IF (g(2) >= solver%ct) THEN
    solver%nac = solver%nac + 1
    ic(solver%nac) = 2
    a(1,solver%nac) = 2. * x(1) - 1.0
    a(2,solver%nac) = 4. * x(2)
    a(3,solver%nac) = 2. * x(3)
    a(4,solver%nac) = 4. * x(4) - 1.0
  END IF
!
  IF (g(3) >= solver%ct) THEN
    solver%nac = solver%nac + 1
    ic(solver%nac) = 3
    a(1,solver%nac) = 4. * x(1) + 2.
    a(2,solver%nac) = 2. * x(2) - 1.
    a(3,solver%nac) = 2. * x(3)
    a(4,solver%nac) = -1.
  END IF
END IF

END SUBROUTINE analys

END PROGRAM exampl2


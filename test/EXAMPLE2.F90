PROGRAM exampl2

  !! THIS PROGRAM EXECUTES THE EXAMPLE PROBLEM TWO OF THE CONMIN MANUAL

  USE conmin_module, only: conmin_class, wp => conmin_wp

  IMPLICIT NONE

REAL (wp)  :: s(6), g1(11), g2(11), b(11,11), c(11), vlb(6), vub(6),  &
              scal(6), df(6), a(6,11)
INTEGER    :: ms1(22), isc(11), ic(11)
REAL (wp)  :: aobj, x(6), g(11)
INTEGER    :: i, n1, n2, n3, n4, n5, nlim
type(conmin_class) :: solver

OPEN (UNIT=6,FILE='EXOUT2.TXT',STATUS='REPLACE')

!  INITIALIZE
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

! iterative part of analysis
do  i = 1, nlim
  ! call the optimization routine conmin
  call solver%solve(x,vlb,vub,g,scal,df,a,s,g1,g2,b,c,isc,ic,ms1,n1,n2,n3,n4,n5)
  ! analysis module
  call analys()
  solver%obj = aobj
  if (solver%igoto == 0) exit
end do

contains

SUBROUTINE analys()

  !! routine to calculate objective function and constraint values
  !! for optimization of constrained rosen-suzuki function.

  if (solver%info < 2) then

    ! objective function

    aobj = x(1) ** 2 - 5.0_wp * x(1) + x(2) ** 2 - 5.0_wp * x(2) + 2.0_wp *  &
      x(3) ** 2 - 21.0_wp * x(3) + x(4) ** 2 + 7.0_wp * x(4) + 50.0_wp

    ! constraint values

    g(1) = x(1) ** 2 + x(1) + x(2) ** 2 - x(2) + x(3) ** 2 + x(3) +  &
      x(4) ** 2 - x(4) - 8.0_wp

    g(2) = x(1) ** 2 - x(1) + 2.0_wp * x(2) ** 2 + x(3) ** 2 + 2.0_wp *  &
      x(4) ** 2 - x(4) - 10.0_wp

    g(3) = 2. * x(1) ** 2 + 2.0_wp * x(1) + x(2) ** 2 - x(2) + x(3) ** 2  &
      - x(4) - 5.0_wp

  else

    ! gradient information

    df(1) = 2.0_wp * x(1) - 5.0_wp
    df(2) = 2.0_wp * x(2) - 5.0_wp
    df(3) = 4.0_wp * x(3) - 21.0_wp
    df(4) = 2.0_wp * x(4) + 7.0_wp

    ! gradients of active and violated constraints

    solver%nac = 0
    if (g(1) >= solver%ct) then
      solver%nac = 1
      ic(1) = 1
      a(1,1) = 2.0_wp * x(1) + 1.0_wp
      a(2,1) = 2.0_wp * x(2) - 1.0_wp
      a(3,1) = 2.0_wp * x(3) + 1.0_wp
      a(4,1) = 2.0_wp * x(4) - 1.0_wp
    end if

    if (g(2) >= solver%ct) then
      solver%nac = solver%nac + 1
      ic(solver%nac) = 2
      a(1,solver%nac) = 2.0_wp * x(1) - 1.0_wp
      a(2,solver%nac) = 4.0_wp * x(2)
      a(3,solver%nac) = 2.0_wp * x(3)
      a(4,solver%nac) = 4.0_wp * x(4) - 1.0_wp
    end if

    if (g(3) >= solver%ct) then
      solver%nac = solver%nac + 1
      ic(solver%nac) = 3
      a(1,solver%nac) = 4.0_wp * x(1) + 2.0_wp
      a(2,solver%nac) = 2.0_wp * x(2) - 1.0_wp
      a(3,solver%nac) = 2.0_wp * x(3)
      a(4,solver%nac) = -1.0_wp
    end if
  end if

END SUBROUTINE analys

END PROGRAM exampl2


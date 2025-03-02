!***********************************************************************
!>
!  CONMIN is a Fortran package for the solution
!  of linear or nonlinear constrained optimization problems.  The basic
!  optimization algorithm is the Method of Feasible Directions.  The user
!  must provide a main calling program and an external routine to evaluate
!  the objective and constraint functions and to provide gradient
!  information.  If analytic gradients of the objective or constraint
!  functions are not available, this information is calculated by finite
!  difference.  While the program is intended primarily for efficient
!  solution of constrained problems, unconstrained function minimization
!  problems may also be solved, and the conjugate direction method of
!  Fletcher and Reeves is used for this purpose.
!
!### History
!   * BY G. N. VANDERPLAATS
!     NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF.
!
!### See also
!  * CONMIN - A FORTRAN PROGRAM FOR CONSTRAINED FUNCTION
!    MINIMIZATION:  USER'S MANUAL,  BY G. N. VANDERPLAATS,
!    NASA TM X-62,282, AUGUST, 1973.
!  * [conmin.zip](https://jblevins.org/mirror/amiller/conmin.zip) from Alan Miller.

module conmin_module

    use iso_fortran_env

    implicit none

    private

#ifdef REAL32
    integer,parameter :: wp = real32   !! Real working precision [4 bytes]
#elif REAL64
    integer,parameter :: wp = real64   !! Real working precision [8 bytes]
#elif REAL128
    integer,parameter :: wp = real128  !! Real working precision [16 bytes]
#else
    integer,parameter :: wp = real64   !! Real working precision if not specified [8 bytes]
#endif

    integer,parameter,public :: conmin_wp = wp   !! Working real precision

    type, public :: conmin_class

    !! main class.

        ! private

        ! cnmn1 variables
        real(wp) :: delfun = 0.001_wp   !! Minimum relative change in the objective
                                        !! function to indicate convergence.  If in `ITRM` consecutive
                                        !! iterations, `ABS(1.0-OBJ(J-1)/OBJ(J))<DELFUN` and the current
                                        !! design is feasible (all `G(J)<=ABS(CT)`), the minimization
                                        !! process is terminated.  If the current design is infeasible
                                        !! (some `G(J)>ABS(CT)`), five iterations are required to
                                        !! terminate and this situation indicates that a feasible design
                                        !! may not exist.
        real(wp) :: dabfun = 0.0_wp !! Default value = 0.001 times the initial function value.  Same
                                    !! as DELFUN except comparison is on absolute change in the
                                    !! objective function, `ABS(OBJ(J)-OBJ(J-1))`, instead of relative
                                    !! change.
        real(wp) :: fdch = 0.01_wp   !! Not used if `NFDG = 0`.  Relative change in
                                     !! decision variable `X(I)` in calculating finite difference
                                     !! gradients.  For example, `FDCH = 0.01` corresponds to a finite
                                     !! difference step of one percent of the value of the decision
                                     !! variable.
        real(wp) :: fdchm = 0.01_wp !! Not used if `NFDG = 0`.  Minimum absolute
                                    !! step in finite difference gradient calculations.  FDCHM applies
                                    !! to the unscaled variable values.
        real(wp) :: ct = -0.1_wp    !! Not used if `NCON = NSIDE = 0`.
                                    !! Constraint thickness parameter.  If `CT<=G(J)<=ABS(CT)`,
                                    !! `G(J)` is defined as active.  If `G(J)>ABS(CT)`, `G(J)` is said to
                                    !! be violated.  If `G(J)<CT`, `G(J)` is not active.  CT is
                                    !! sequentially reduced in magnitude during the optimization
                                    !! process.  If `ABS(CT)` is very small, one or more constraints
                                    !! may be active on one iteration and inactive on the next,
                                    !! only to become active again on a subsequent iteration.
                                    !! This is often referred to as "zigzagging" between constraints.
                                    !! A wide initial value of the constraint thickness is desirable
                                    !! for highly nonlinear problems so that when a constraint
                                    !! becomes active it tends to remain active, thus reducing the
                                    !! zigzagging problem.  The default value is usually adequate.
        real(wp) :: ctmin = 0.004_wp    !! Not used if `NCON = NSIDE = 0`.  Minimum
                                        !! absolute value of CT considered in the optimization process.
                                        !! CTMIN may be considered as "numerical zero" since it may not be
                                        !! meaningful to compare numbers smaller than CTMIN.  The value of
                                        !! CTMIN is chosen to indicate that satisfaction of a constraint
                                        !! within this tolerance is acceptable.  The default value is usually
                                        !! adequate.
        real(wp) :: ctl = -0.01_wp  !! Not used if `NCON = NSIDE = 0`.
                                    !! Constraint thickness parameter for linear and side constraints.
                                    !! CTL is smaller in magnitude than CT because the zigzagging
                                    !! problem is avoided with linear and side constraints.  The default
                                    !! value is usually adequate.
        real(wp) :: ctlmin = 0.001_wp   !! Not used if `NCON = NSIDE = 0`.  Minimum
                                        !! absolute value of CTL considered in the optimization process.
                                        !! The default value is usually adequate.
        real(wp) :: alphax = 0.1_wp  !! the maximum fractional change in any
                                     !! component of X as an initial estimate for ALPHA in the one-dimensional
                                     !! search. That is, the initial ALPHA will be such that no component of X is changed by
                                     !! more than this amount. This only applies to those X(i) of magnitude greater than 0.1.
                                     !! If an optimization run shows numerous ALPHA = 0 results for the one-dimensional search,
                                     !! it may help to try ALPHAX less than the default. ALPHAX is changed by CONMIN depending
                                     !! on the progress of the optimization.
        real(wp) :: abobj1 = 0.1_wp  !! the fractional change attempted as a first step in
                                     !! the one-dimensional search and is based on a linear approximation.
                                     !! ABOBJ1 is updated during the optimization, depending on progress.
                                     !! The initial step in the one-dimensional search is taken as the amount
                                     !! necessary to change `OBJ` by `ABOBJ1*ABS(OBJ)` or to change some `X(i)` by `ALPHAX*ABS( X(i) )`,
                                     !! whichever is less.
        real(wp) :: theta = 1.0_wp !! Not used if `NCON = NSIDE = 0`.  Mean value
                                    !! of the push-off factor in the method of feasible directions.
                                    !! A larger value of THETA is desirable if the constraints, `G(J)`,
                                    !! are known to be highly nonlinear, and a smaller value may be
                                    !! used if all `G(J)` are known to be nearly linear.  The actual
                                    !! value of the push-off factor used in the program is a quadratic
                                    !! function of each `G(J)`, varying from 0.0 for `G(J) = CT` to `4.0*THETA`
                                    !! for `G(J) = ABS(CT)`.  A value of `THETA = 0.0` is used in the
                                    !! program for constraints which are identified by the user to be
                                    !! strictly linear.  THETA is called a "push-off" factor because
                                    !! it pushes the design away from the active constraints into the
                                    !! feasible region.  The default value is usually adequate.
        real(wp) :: obj !! Value of objective function for the current decision variables,
                        !! `X(I), I = 1, NDV` contained in vector `X`.  Calculate OBJ if
                        !! INFO = 1 or INFO = 2.
        integer :: ndv !! Number of decision variables, `X(I)`, contained in vector `X`.
        integer :: ncon !! Number of constraint functions, `G(J)`.  NCON may be zero.
        integer :: nside    !! Side constraint parameter.  `NSIDE = 0` signifies that the
                            !! variables `X(I)` do not have lower or upper bounds.  `NSIDE>0`
                            !! signifies that all variables `X(I)` have lower and upper bounds
                            !! defined by `VLB(I)` and `VUB(I)` respectively.  If one or more
                            !! variables are not bounded while others are, the values of the
                            !! lower and upper bounds on the unbounded variables must be taken
                            !! as very large negative and positive values respectively
                            !! (i.e., `VLB(I) = -1.0E+10`, `VUB(I) = 1.0E+10`).
        integer :: iprint = 0   !! Print control.  All printing is done on unit number 6.
                                !!
                                !!  * `iprint=0`:  Print nothing.
                                !!  * `iprint=1`:  Print initial and final function information.
                                !!  * `iprint=2`:  1st debug level.  Print all of above plus control
                                !!      parameters.  Print function value and X-vector at each
                                !!      iteration.
                                !!  * `iprint=3`:  2nd. debug level.  Print all of above plus all constraint
                                !!      values, numbers of active or violated constraints, direction
                                !!      vectors, move parameters and miscellaneous information.  The
                                !!      constraint parameter, BETA, printed under this option
                                !!      approaches zero as the optimum objective is achieved.
                                !!  * `iprint=4`:  Complete debug.  Print all of above plus gradients of
                                !!      objective function, active or violated constraint functions
                                !!      and miscellaneous information.
                                !!  * `iprint=5`:  all of above plus each proposed design vector, objective
                                !!      and constraints during the one-dimensional search.
        integer :: nfdg !! the finite difference gradient parameter:
                        !!
                        !!  * `NFDG = 0`:  all gradient information is calculated by finite difference
                        !!    within [[conmin]].
                        !!  * `NFDG = 1`:  all gradient information is supplied by the user.
                        !!  * `NFDG = 2`:  the gradient of `OBJ` is supplied by the user and the
                        !!    gradients of constraints are calculated by finite
                        !!    difference within [[conmin]].
        integer :: nscal    !! Scaling control parameter.  The decision variables will be
                            !! scaled linearly.
                            !!
                            !!  * `NSCAL<0`:  Scale variables `X(I) `by dividing by `SCAL(I)`, where
                            !!                vector SCAL is defined by the user.
                            !!  * `NSCAL==0`: Do not scale the variables.
                            !!  * `NSCAL>0`:  Scale the variables every NSCAL iterations.
                            !!                Variables are normalized so that scaled
                            !!                `X(I) = X(I)/ABS(X(I))`.  When using this option, it
                            !!                is desirable that `NSCAL = ICNDIR` if ICNDIR is input
                            !!                as nonzero, and `NSCAL = NDV + 1` in ICNDIR is input
                            !!                as zero.
        integer :: linobj = 0   !! Not used if NCON = NSIDE = 0.  Linear objective function
                                !! identifier.  If the objective, OBJ, is specifically known to
                                !! be a strictly linear function of the decision variables, `X(I)`,
                                !! set `LINOBJ = 1`.  If OBJ is a general nonlinear function, set
                                !! `LINOBJ = 0`.
        integer :: itmax = 10 !! Maximum number of iterations in the
                              !! minimization process.  If `NFDG==0` each iteration requires one
                              !! set of gradient computations (INFO = 3 or 4) and approximately
                              !! three function evaluations (INFO = 1 or 2).  If `NFDG>0`
                              !! each iteration requires approximately NDV + 3 function
                              !! evaluations (INFO = 1 or 2).
        integer :: itrm = 3 !! Number of consecutive iterations to indicate
                            !! convergence by relative or absolute changes, `DELFUN` or `DABFUN`.
        integer :: icndir !! Default value = `NDV + 1`.  Conjugate direction restart parameter.
                          !! If the function is currently unconstrained, (all `G(J)<CT` or
                          !! `NCON = NSIDE = 0`), Fletcher-Reeves conjugate direction method will
                          !! be restarted with a steepest descent direction every ICNDIR
                          !! iterations.  If `ICNDIR = 1` only steepest descent will be used.
        integer :: igoto = 0 !! Reverse communication flag.
                             !! CONMIN executes according to the parameter IGOTO which must be initialized to zero.
        integer :: nac  !! Number of active and violated constraints (`G(J)>=CT`).
                        !! Calculate `NAC` if `INFO = 4` and `NFDG = 0`.
        integer :: info !! * `INFO = 1`:  calculate OBJ and G(I), I = 1, NCON
                        !! * `INFO = 2`:  calculate NAC, IC(I), I = 1,NAC, the gradient of OBJ, and
                        !!   the gradient of G(J), where J = IC(I), I = 1,NAC.
                        !!   Store the gradients of G in columns of A.
        integer :: infog !! * `INFOG = 0`:   same as when INFOG was not used.
                         !! * `INFOG = 1`:   only those constraints identified as active or violated
                         !!                in array IC(I), I = 1, NAC need be evaluated.  This is
                         !!                only meaningful if finite difference gradients are
                         !!                calculated, and allows the user to avoid
                         !!                calculating non-essential information.  If it is
                         !!                convenient to evaluate all constraints each time,
                         !!                variable INFOG may be ignored.
        integer :: iter !! Iteration number.  The optimization process is iterative so
                        !! that the vector of decision variables at the Kth iteration
                        !! is defined by `X(K) = X(K - 1) + ALPHA*S(K)`, where in this case
                        !! `K` refers to the iteration number and the components `X(I)` are
                        !! all changed simultaneously.  `ALPHA` is defined as the move
                        !! parameter and is printed if the print control `IPRINT>=3`.
                        !! `S` is the move direction.

        ! consav variables
        real(wp) :: dm1, dm2, dm3, dm4, dm5, dm6, dm7, dm8, dm9, dm10, dm11, &
                    dm12, dct, dctl, abobj, cta, ctam, ctbm, obj1, &
                    slope, dx, dx1, fi, xi, dftdf1, alp, fff, a1, a2, a3, &
                    a4, f1, f2, f3, f4, cv1, cv2, cv3, cv4, app, alpca, &
                    alpfes, alpln, alpmin, alpnc, alpsav, alpsid, alptot, &
                    rspace
        real(wp) :: phi = 5.0_wp    !! Not used if `NCON = NSIDE = 0`.
                                    !! Participation coefficient, used if a design is infeasible
                                    !! (one or more `G(J)>ABS(CT)`).  PHI is a measure of how hard
                                    !! the design will be "pushed" towards the feasible region and
                                    !! is, in effect, a penalty parameter.  If in a given problem, a
                                    !! feasible solution cannot be obtained with the default value,
                                    !! PHI should be increased, and the problem run again.  If a
                                    !! feasible solution cannot be obtained with `PHI = 100`, it is
                                    !! probable that no feasible solution exists.  The default value
                                    !! is usually adequate.
        integer :: idm1, idm2, idm3, jdir, iobj, kobj, kcount, &
                   nfeas, mscal, ncobj, nvc, kount, icount, igood1, igood2, &
                   igood3, igood4, ibest, iii, nlnc, jgoto, ispace(2)
        integer :: ncal(2)  !! Bookkeeping information.  NCAL(1) gives the number of times
                            !! external routine SUB1 was called with INFO = 1.  NCAL(2) gives
                            !! the number of times INFO = 2.  NCAL(3) gives the number of times
                            !! INFO = 3 and NCAL(4) gives the number of times INFO = 4.

    contains
        private
        procedure, public :: solve => conmin
        procedure :: cnmn01
        procedure :: cnmn02
        procedure :: cnmn03
        procedure :: cnmn05
        procedure :: cnmn06
    end type conmin_class

contains

    subroutine conmin(me, x, vlb, vub, g, scal, df, a, s, g1, g2, b, c, isc, ic, ms1, n1, n2, n3, n4, n5)

        !!  Routine to solve constrained or unconstrained function minimization.
        !!
        !!  JUNE, 1979 VERSION
        !!
        !!  STORAGE REQUIREMENTS:
        !!      PROGRAM - 7000 DECIMAL WORDS (CDC COMPUTER)
        !!      ARRAYS  - APPROX. 2*(NDV**2)+26*NDV+4*NCON, WHERE N3 = NDV+2.
        !!  RE-SCALE VARIABLES IF REQUIRED.

        class(conmin_class), intent(inout) :: me
        integer, intent(in)       :: n1 !! `N1 = NDV + 2`
        integer, intent(in)       :: n2 !! `N2 = NCON + 2*NDV`
        integer, intent(in)       :: n3 !! `N3 = NACMX1`
        integer, intent(in)       :: n4 !! `N4 = MAX (N3,NDV)`
        integer, intent(in)       :: n5 !! `N5 = 2*N4`
        real(wp), intent(inout)   :: x(n1) !! Vector of decision variables, `X(I), I = 1, NDV`.  The initial
                                           !! X-vector contains the user's best estimate of the set of optimum
                                           !! design variables.
        real(wp), intent(inout)   :: vlb(n1) !! Used only if `NSIDE/=0`.  `VLB(I)` is the lower allowable value
                                             !! (lower bound) of variable `X(I)`.  If one or more variables, `X(I)`,
                                             !! do not have lower bounds, the corresponding `VLB(I)` must be
                                             !! initialized to a very large negative number (say `-1.0E+10`).
        real(wp), intent(inout)   :: vub(n1) !! Used only if `NSIDE/=0`.  `VUB(I)` is the maximum allowable value
                                             !! (upper bound) of `X(I)`.  If one or more variables, `X(I)`, do not
                                             !! have upper bounds, the corresponding `VUB(I)` must be initialized
                                             !! to a very large positive number (say `1.0E+10`).
        real(wp), intent(inout)   :: g(n2) !! Not used if `NCON = NSIDE = 0`.  Vector containing all constraint
                                           !! functions, `G(J), J = 1, NCON` for current decision variables, `X`.
                                           !! Calculate `G(J), J = 1, NCON` if `INFO = 2`.
        real(wp), intent(inout)   :: scal(n1) !! Not used if `NSCAL = 0`.  Vector of scaling parameters.  If
                                              !! `NSCAL>0` vector SCAL need not be initialized since SCAL will
                                              !! be defined in CONMIN and its associated routines.  If `NSCAL<0`,
                                              !! vector SCAL is initialized in the main program, and the scaled
                                              !! variables `X(I) = X(I)/SCAL(I)`.  Efficiency of the optimization
                                              !! process can sometimes be improved if the variables are either
                                              !! normalized or are scaled in such a way that the partial deri-
                                              !! vative of the objective function, OBJ, with respect to variable
                                              !! `X(I)` is of the same order of magnitude for all `X(I)`.  `SCAL(I)`
                                              !! must be greater than zero because a negative value of `SCAL(I)`
                                              !! will result in a change of sign of `X(I)` and possibly yield
                                              !! erroneous optimization results.  The decision of if, and how, the
                                              !! variables should be scaled is highly problem dependent, and some
                                              !! experimentation is desirable for any given class of problems.
        real(wp), intent(inout)   :: df(n1) !! Analytic gradient of the objective function for the current
                                            !! decision variables, `X(I)`.  `DF(I)` contains the partial derivative
                                            !! of `OBJ` with respect to `X(I)`.  Calculate `DF(I), I = 1,
                                            !! NDV` if `INFO = 3` or `INFO = 4` and if `NFDG = 0` or `NFDG = 2`.
        real(wp), intent(inout)   :: a(n1, n3) !! Not used if `NCON = NSIDE = 0`.  Gradients of active or violated
                                               !! constraints, for current decision variables, `X(I)`.
                                               !! `A(J,I)` contains the gradient of the Jth active or violated
                                               !! constraint, `G(J)`, with respect to the Ith decision variable,
                                               !! `X(I)` for `J = 1, NAC` and `I = 1, NDV`.  Calculate if `INFO = 4`
                                               !! and `NFDG = 0`.
        real(wp), intent(inout)   :: s(n1) !! Move direction in the NDV-dimensional optimization space.  `S(I)`
                                           !! gives the rate at which variable `X(I)` changes with respect to
                                           !! `ALPHA`.
        real(wp), intent(inout)   :: g1(n2) !! Not used if `NCON = NSIDE = NSCAL = 0`.  Used for temporary
                                            !! storage of constraint values `G(J), J = 1, NCON` and decision
                                            !! variables `X(I), I = 1, NDV`.
        real(wp), intent(inout)   :: g2(n2) !! Not used if `NCON = NSIDE = 0`.  Used for temporary storage of
                                            !! constraint values `G(J), J = 1, NCON`.
        real(wp), intent(inout)   :: b(n3, n3) !! Not used if `NCON = NSIDE = 0`.  Used in determining direction
                                               !! vector S for constrained minimization problems.  Array `B` may
                                               !! be used for temporary storage in external routine SUB1.
        real(wp), intent(inout)   :: c(n4) !! Not used in `NCON = NSIDE = 0`.  Used with array B in determining
                                           !! direction vector S for constrained minimization problems.  Used
                                           !! for temporary storage of vector X if `NSCAL/=0`. routine SUB1.
        integer, intent(inout)    :: isc(n2) !! Not used if `NCON = 0`.  Linear constraint identification vector.
                                             !! If constraint `G(J)` is known to be a linear function of the
                                             !! decision variables, `X(I)`, `ISC(I)` should be initialized to
                                             !! `ISC(I) = 1`.  If constraint `G(J)` is nonlinear `ISC(I)` is initialized
                                             !! to `ISC(I) = 0`.  Identification of linear constraints may improve
                                             !! efficiency of the optimization process and is therefore desirable,
                                             !! but is not essential.  If `G(J)` is not specifically known to be
                                             !! linear, set `ISC(I) = 0`.
        integer, intent(inout)    :: ic(n3) !! Identifies which constraints are active or violated.  `IC(J)`
                                            !! contains the number of the Jth active or violated constraint
                                            !! for `J = 1, NAC`.  For example, if `G(10)` is the first active
                                            !! or violated constraint (`G(J)<CT, J = 1,9`), set `IC(1) = 10`.
                                            !! Calculate if `INFO = 4` and `NFDG = 0`.
        integer, intent(inout)    :: ms1(n5) !! Not used if `NCON = NSIDE = 0`.  Used with array `B` in determining
                                             !! direction vector `S` for constrained minimization problems.  Array
                                             !! MS1 may be used for temporary storage in external routine SUB1.

        real(wp) :: alp1, alp11, alp12, c1, ct1, ctc, ff1, gi, objb, objd, &
                    scj, si, sib, x1, x12, xid, xx
        integer :: i, ii, j, k, m1, m2, m3, mcn1, nci, ndv1, ndv2, nfeasct, nic, nnac

        if (me%nscal /= 0 .and. me%igoto /= 0) then
            do i = 1, me%ndv
                x(i) = c(i)
            end do
        end if
        ! CONSTANTS.
        ndv1 = me%ndv + 1
        ndv2 = me%ndv + 2
        if (me%igoto /= 0) then
            ! ------------------------------------------------------------------
            !                 CHECK FOR UNBOUNDED SOLUTION
            ! ------------------------------------------------------------------
            ! STOP IF OBJ IS LESS THAN -1.0D+40
            if (me%obj <= -1.0e+40_wp) then
                write (6, 5100)
                go to 520
            end if
            select case (me%igoto)
            case (1); go to 60
            case (2); go to 210
            case (3); go to 200
            case (4); go to 410
            case (5); go to 430
            end select
        end if
        ! ------------------------------------------------------------------
        !                  SAVE INPUT CONTROL PARAMETERS
        ! ------------------------------------------------------------------
        if (me%iprint > 0) write (6, 7500)
        if (.not. (me%linobj == 0 .or. (me%ncon > 0 .or. me%nside > 0))) then
            ! TOTALLY UNCONSTRAINED FUNCTION WITH LINEAR OBJECTIVE.
            ! SOLUTION IS UNBOUNDED.
            write (6, 5000) me%linobj, me%ncon, me%nside
            return
        end if
        me%idm1 = me%itrm
        me%idm2 = me%itmax
        me%idm3 = me%icndir
        me%dm1 = me%delfun
        me%dm2 = me%dabfun
        me%dm3 = me%ct
        me%dm4 = me%ctmin
        me%dm5 = me%ctl
        me%dm6 = me%ctlmin
        me%dm7 = me%theta
        me%dm8 = me%phi
        me%dm9 = me%fdch
        me%dm10 = me%fdchm
        me%dm11 = me%abobj1
        me%dm12 = me%alphax
        ! ------------------------------------------------------------------
        !                            DEFAULTS
        ! ------------------------------------------------------------------
        if (me%itrm <= 0) me%itrm = 3
        if (me%itmax <= 0) me%itmax = 20
        ndv1 = me%ndv + 1
        if (me%icndir == 0) me%icndir = ndv1
        if (me%delfun <= 0.0_wp) me%delfun = 0.0001_wp
        me%ct = -abs(me%ct)
        if (me%ct >= 0.0_wp) me%ct = -0.1_wp
        me%ctmin = abs(me%ctmin)
        if (me%ctmin <= 0.0_wp) me%ctmin = 0.004_wp
        me%ctl = -abs(me%ctl)
        if (me%ctl >= 0.0_wp) me%ctl = -0.01_wp
        me%ctlmin = abs(me%ctlmin)
        if (me%ctlmin <= 0.0_wp) me%ctlmin = 0.001_wp
        if (me%theta  <= 0.0_wp) me%theta  = 1.0_wp
        if (me%abobj1 <= 0.0_wp) me%abobj1 = 0.1_wp
        if (me%alphax <= 0.0_wp) me%alphax = 0.1_wp
        if (me%fdch   <= 0.0_wp) me%fdch   = 0.01_wp
        if (me%fdchm  <= 0.0_wp) me%fdchm  = 0.01_wp
        ! ------------------------------------------------------------------
        !                 INITIALIZE INTERNAL PARAMETERS
        ! ------------------------------------------------------------------
        me%infog = 0
        me%iter = 0
        me%jdir = 0
        me%iobj = 0
        me%kobj = 0
        ndv2 = me%ndv + 2
        me%kcount = 0
        me%ncal(1) = 0
        me%ncal(2) = 0
        me%nac = 0
        me%nfeas = 0
        me%mscal = me%nscal
        ct1 = me%itrm
        ct1 = 1.0_wp/ct1
        me%dct = (me%ctmin/abs(me%ct))**ct1
        me%dctl = (me%ctlmin/abs(me%ctl))**ct1
        me%phi = 5.0_wp
        me%abobj = me%abobj1
        me%ncobj = 0
        me%ctam = abs(me%ctmin)
        me%ctbm = abs(me%ctlmin)
        ! CALCULATE NUMBER OF LINEAR CONSTRAINTS, NLNC.
        me%nlnc = 0
        if (me%ncon /= 0) then
            do i = 1, me%ncon
                if (isc(i) > 0) me%nlnc = me%nlnc + 1
            end do
        end if
        ! ------------------------------------------------------------------
        !      CHECK TO BE SURE THAT SIDE CONSTRAINTS ARE SATISFIED
        ! ------------------------------------------------------------------
        if (me%nside /= 0) then
            do i = 1, me%ndv
                if (vlb(i) > vub(i)) then
                    xx = .5*(vlb(i) + vub(i))
                    x(i) = xx
                    vlb(i) = xx
                    vub(i) = xx
                    write (6, 6500) i
                end if
                xx = x(i) - vlb(i)
                if (xx < 0.0_wp) then
                    ! LOWER BOUND VIOLATED.
                    write (6, 6600) x(i), vlb(i), i
                    x(i) = vlb(i)
                else
                    xx = vub(i) - x(i)
                    if (xx < 0.0_wp) then
                        write (6, 6700) x(i), vub(i), i
                        x(i) = vub(i)
                    end if
                end if
            end do
        end if
        ! ------------------------------------------------------------------
        !                    INITIALIZE SCALING VECTOR, SCAL
        ! ------------------------------------------------------------------
        if (me%nscal /= 0) then
            if (me%nscal >= 0) then
                scal(1:me%ndv) = 1.0_wp
            else
                do i = 1, me%ndv
                    si = abs(scal(i))
                    if (si < 1.0e-20_wp) si = 1.0e-5_wp
                    scal(i) = si
                    si = 1.0_wp/si
                    x(i) = x(i)*si
                    if (me%nside /= 0) then
                        vlb(i) = vlb(i)*si
                        vub(i) = vub(i)*si
                    end if
                end do
            end if
        end if
        ! ------------------------------------------------------------------
        ! ***** CALCULATE INITIAL FUNCTION AND CONSTRAINT VALUES  *****
        ! ------------------------------------------------------------------
        me%info = 1
        me%ncal(1) = 1
        me%igoto = 1
        go to 580

60      me%obj1 = me%obj
        if (me%dabfun <= 0.0_wp) me%dabfun = 0.001_wp*abs(me%obj)
        if (me%dabfun < 1.0e-10_wp) me%dabfun = 1.0e-10_wp
        if (me%iprint > 0) then
            ! ------------------------------------------------------------------
            !                PRINT INITIAL DESIGN INFORMATION
            ! ------------------------------------------------------------------
            if (me%iprint > 1) then
                if (me%nside == 0 .and. me%ncon == 0) write (6, 8200)
                if (me%nside /= 0 .or. me%ncon > 0) write (6, 7600)
                write (6, 7700) me%iprint, me%ndv, me%itmax, me%ncon, me%nside, me%icndir, me%nscal, &
                                me%nfdg, me%linobj, me%itrm, n1, n2, n3, n4, n5
                write (6, 7900) me%ct, me%ctmin, me%ctl, me%ctlmin, me%theta, me%phi, me%delfun, me%dabfun
                write (6, 7800) me%fdch, me%fdchm, me%alphax, me%abobj1
                if (me%nside /= 0) then
                    write (6, 8000)
                    do i = 1, me%ndv, 6
                        m1 = min(me%ndv, i + 5)
                        write (6, 5400) i, vlb(i:m1)
                    end do
                    write (6, 8100)
                    do i = 1, me%ndv, 6
                        m1 = min(me%ndv, i + 5)
                        write (6, 5400) i, vub(i:m1)
                    end do
                end if
                if (me%nscal < 0) then
                    write (6, 8300)
                    write (6, 9900) scal(1:me%ndv)
                end if
                if (me%ncon /= 0) then
                    if (me%nlnc /= 0 .and. me%nlnc /= me%ncon) then
                        write (6, 5500)
                        do i = 1, me%ncon, 15
                            m1 = min(me%ncon, i + 14)
                            write (6, 5600) i, isc(i:m1)
                        end do
                    else
                        if (me%nlnc == me%ncon) write (6, 5700)
                        if (me%nlnc == 0) write (6, 5800)
                    end if
                end if
            end if
            write (6, 9700) me%obj
            write (6, 9800)
            do i = 1, me%ndv
                x1 = 1.
                if (me%nscal /= 0) x1 = scal(i)
                g1(i) = x(i)*x1
            end do
            do i = 1, me%ndv, 6
                m1 = min(me%ndv, i + 5)
                write (6, 5400) i, g1(i:m1)
            end do
            if (me%ncon /= 0) then
                write (6, 10000)
                do i = 1, me%ncon, 6
                    m1 = min(me%ncon, i + 5)
                    write (6, 5400) i, g(i:m1)
                end do
            end if
        end if
        if (me%iprint > 1) write (6, 8900)
        ! ------------------------------------------------------------------
        ! ********************  BEGIN MINIMIZATION  ************************
        ! ------------------------------------------------------------------
130     me%iter = me%iter + 1
        if (me%abobj1 < 0.0001_wp) me%abobj1 = 0.0001_wp
        if (me%abobj1 > 0.2_wp) me%abobj1 = 0.2_wp
        if (me%alphax > 1.0_wp) me%alphax = 1.0_wp
        if (me%alphax < 0.001_wp) me%alphax = 0.001_wp

        ! THE FOLLOWING TWO LINES OF CODE WERE COMMENTED OUT ON 3/5/81

!       NFEAS=NFEAS+1
!       IF (NFEAS>10) GO TO 810
        if (me%iprint > 2) write (6, 8400) me%iter
        if (me%iprint > 3 .and. me%ncon > 0) write (6, 8500) me%ct, me%ctl, me%phi
        me%cta = abs(me%ct)
        if (me%ncobj /= 0) then
            ! ------------------------------------------------------------------
            ! NO MOVE ON LAST ITERATION.  DELETE CONSTRAINTS THAT ARE NO LONGER ACTIVE.
            ! ------------------------------------------------------------------
            nnac = me%nac
            do i = 1, nnac
                if (ic(i) > me%ncon) me%nac = me%nac - 1
            end do
            if (me%nac <= 0) go to 250
            nnac = me%nac
            do i = 1, nnac
150             nic = ic(i)
                ct1 = me%ct
                if (isc(nic) > 0) ct1 = me%ctl
                if (g(nic) <= ct1) then
                    me%nac = me%nac - 1
                    if (i > me%nac) go to 250
                    do k = i, me%nac
                        ii = k + 1
                        do j = 1, ndv2
                            a(j, k) = a(j, ii)
                        end do
                        ic(k) = ic(ii)
                    end do
                    go to 150
                end if
            end do
            go to 250
        end if
        if (me%mscal >= me%nscal .and. me%nscal /= 0) then
            if (me%nscal >= 0 .or. me%kcount >= me%icndir) then
                me%mscal = 0
                me%kcount = 0
                ! ------------------------------------------------------------------
                !                      SCALE VARIABLES
                ! ------------------------------------------------------------------
                do i = 1, me%ndv
                    si = scal(i)
                    me%xi = si*x(i)
                    sib = si
                    if (me%nscal > 0) si = abs(me%xi)
                    if (si >= 1.0e-10_wp) then
                        scal(i) = si
                        si = 1.0_wp/si
                        x(i) = me%xi*si
                        if (me%nside /= 0) then
                            vlb(i) = sib*si*vlb(i)
                            vub(i) = sib*si*vub(i)
                        end if
                    end if
                end do
                if (.not. (me%iprint < 4 .or. (me%nscal < 0 .and. me%iter > 1))) then
                    write (6, 8600)
                    write (6, 9900) scal(1:me%ndv)
                end if
            end if
        end if
        me%mscal = me%mscal + 1
        me%nac = 0
        ! ------------------------------------------------------------------
        !      OBTAIN GRADIENTS OF OBJECTIVE AND ACTIVE CONSTRAINTS
        ! ------------------------------------------------------------------
        me%info = 2
        me%ncal(2) = me%ncal(2) + 1
        if (me%nfdg == 1) then
            me%igoto = 2
            go to 580
        end if
        me%jgoto = 0

200     call me%cnmn01(me%jgoto, x, df, g, isc, ic, a, g1, vub, scal, me%ncal, me%dx, me%dx1, &
                       me%fi, me%xi, me%iii, n1, n2, n3)
        me%igoto = 3
        if (me%jgoto > 0) go to 580

210     me%info = 1
        if (me%nac >= n3) go to 520
        if (me%nscal /= 0 .and. me%nfdg /= 0) then
            ! ------------------------------------------------------------------
            !                          SCALE GRADIENTS
            ! ------------------------------------------------------------------
            ! SCALE GRADIENT OF OBJECTIVE FUNCTION.
            df(1:me%ndv) = df(1:me%ndv)*scal(1:me%ndv)
            if (me%nfdg /= 2 .and. me%nac /= 0) then
                ! SCALE GRADIENTS OF ACTIVE CONSTRAINTS.
                do j = 1, me%ndv
                    scj = scal(j)
                    a(j, 1:me%nac) = a(j, 1:me%nac)*scj
                end do
            end if
        end if

250     if (me%iprint >= 3 .and. me%ncon /= 0) then
            ! ------------------------------------------------------------------
            !                               PRINT
            ! ------------------------------------------------------------------
            ! PRINT ACTIVE AND VIOLATED CONSTRAINT NUMBERS.
            m1 = 0
            m2 = n3
            if (me%nac /= 0) then
                do i = 1, me%nac
                    j = ic(i)
                    if (j <= me%ncon) then
                        gi = g(j)
                        c1 = me%ctam
                        if (isc(j) > 0) c1 = me%ctbm
                        gi = gi - c1
                        if (gi <= 0.0_wp) then
                            ! ACTIVE CONSTRAINT.
                            m1 = m1 + 1
                            ms1(m1) = j
                        else
                            m2 = m2 + 1
                            ! VIOLATED CONSTRAINT.
                            ms1(m2) = j
                        end if
                    end if
                end do
            end if
            m3 = m2 - n3
            write (6, 5900) m1
            if (m1 /= 0) then
                write (6, 6000)
                write (6, 10100) ms1(1:m1)
            end if
            write (6, 6100) m3
            if (m3 /= 0) then
                write (6, 6000)
                m3 = n3 + 1
                write (6, 10100) ms1(m3:m2)
            end if
        end if
        ! ------------------------------------------------------------------
        !        CALCULATE GRADIENTS OF ACTIVE SIDE CONSTRAINTS
        ! ------------------------------------------------------------------
        if (me%nside /= 0) then
            mcn1 = me%ncon
            m1 = 0
            do i = 1, me%ndv
                ! LOWER BOUND.
                me%xi = x(i)
                xid = vlb(i)
                x12 = abs(xid)
                if (x12 < 1.0_wp) x12 = 1.0_wp
                gi = (xid - me%xi)/x12
                if (gi >= -1.0e-6_wp) then
                    m1 = m1 + 1
                    ms1(m1) = -i
                    me%nac = me%nac + 1
                    if (me%nac >= n3) go to 520
                    mcn1 = mcn1 + 1
                    do j = 1, me%ndv
                        a(j, me%nac) = 0.0_wp
                    end do
                    a(i, me%nac) = -1.0_wp
                    ic(me%nac) = mcn1
                    g(mcn1) = gi
                    isc(mcn1) = 1
                end if
                ! UPPER BOUND.
                xid = vub(i)
                x12 = abs(xid)
                if (x12 < 1.0_wp) x12 = 1.0_wp
                gi = (me%xi - xid)/x12
                if (gi >= -1.0e-6_wp) then
                    m1 = m1 + 1
                    ms1(m1) = i
                    me%nac = me%nac + 1
                    if (me%nac >= n3) go to 520
                    mcn1 = mcn1 + 1
                    do j = 1, me%ndv
                        a(j, me%nac) = 0.0_wp
                    end do
                    a(i, me%nac) = 1.0_wp
                    ic(me%nac) = mcn1
                    g(mcn1) = gi
                    isc(mcn1) = 1
                end if
            end do
            ! ------------------------------------------------------------------
            !                              PRINT
            ! ------------------------------------------------------------------
            ! PRINT ACTIVE SIDE CONSTRAINT NUMBERS.
            if (me%iprint >= 3) then
                write (6, 6200) m1
                if (m1 /= 0) then
                    write (6, 6300)
                    write (6, 10100) (ms1(j), j=1, m1)
                end if
            end if
        end if
        ! PRINT GRADIENTS OF ACTIVE AND VIOLATED CONSTRAINTS.
        if (me%iprint >= 4) then
            write (6, 8700)
            do i = 1, me%ndv, 6
                m1 = min(me%ndv, i + 5)
                write (6, 5400) i, df(i:m1)
            end do
            if (me%nac /= 0) then
                write (6, 8800)
                do i = 1, me%nac
                    m1 = ic(i)
                    m2 = m1 - me%ncon
                    m3 = 0
                    if (m2 > 0) m3 = abs(ms1(m2))
                    if (m2 <= 0) write (6, 5200) m1
                    if (m2 > 0) write (6, 5300) m3
                    do k = 1, me%ndv, 6
                        m1 = min(me%ndv, k + 5)
                        write (6, 5400) k, (a(j, i), j=k, m1)
                    end do
                    write (6, 8900)
                end do
            end if
        end if
        ! ------------------------------------------------------------------
        ! ******************  DETERMINE SEARCH DIRECTION *******************
        ! ------------------------------------------------------------------
        me%alp = 1.0e+20_wp
        if (me%nac > 0) go to 340
        ! ------------------------------------------------------------------
        !                    UNCONSTRAINED FUNCTION
        ! ------------------------------------------------------------------
        ! FIND DIRECTION OF STEEPEST DESCENT OR CONJUGATE DIRECTION.
!
!  S. N. 575 ADDED ON 2/25/81
!
330     me%nvc = 0
        me%nfeas = 0
        me%kcount = me%kcount + 1
        ! IF KCOUNT>ICNDIR  RESTART CONJUGATE DIRECTION ALGORITHM.
        if (me%kcount > me%icndir .or. me%iobj == 2) me%kcount = 1
        if (me%kcount == 1) me%jdir = 0
        ! IF JDIR = 0 FIND DIRECTION OF STEEPEST DESCENT.
        call me%cnmn02(me%jdir, me%slope, me%dftdf1, df, s)
        go to 380
        ! ------------------------------------------------------------------
        !                      CONSTRAINED FUNCTION
        ! ------------------------------------------------------------------
        ! FIND USABLE-FEASIBLE DIRECTION.
340     me%kcount = 0
        me%jdir = 0
        me%phi = 10.0_wp*me%phi
        if (me%phi > 1000.0_wp) me%phi = 1000.0_wp

        ! THE FOLLOWING LINE OF CODE WAS COMMENTED OUT ON 3/5/81
!
!       IF (NFEAS==1) PHI=5.
!       CALCULATE DIRECTION, S.
        call me%cnmn05(g, df, a, s, b, c, me%slope, me%phi, isc, ic, ms1, me%nvc, n1, n2, n3, n4, n5)

        ! THE FOLLOWING LINE WAS ADDED ON 2/25/81
        if (me%nac == 0) go to 330

        ! THE FOLLOWING FIVE LINES WERE COMMENTED OUT ON 3/5/81
        ! REASON : THEY WERE NOT IN G. VANDERPLAATS LISTING
!
!       IF THIS DESIGN IS FEASIBLE AND LAST ITERATION WAS INFEASIBLE,
!       SET ABOBJ1=.05 (5 PERCENT).
!       IF (NVC==0 .AND. NFEAS>1) ABOBJ1=.05
!       IF (NVC==0) NFEAS=0
        if (me%iprint >= 3) then
            write (6, 9000)
            do i = 1, me%nac, 6
                m1 = min(me%nac, i + 5)
                write (6, 5400) i, (a(ndv1, j), j=i, m1)
            end do
            write (6, 7400) s(ndv1)
        end if
        ! ------------------------------------------------------------------
        ! ****************** ONE-DIMENSIONAL SEARCH ************************
        ! ------------------------------------------------------------------
        if (s(ndv1) < 1.0e-6_wp .and. me%nvc == 0) go to 450
        ! ------------------------------------------------------------------
        !             FIND ALPHA TO OBTAIN A FEASIBLE DESIGN
        ! ------------------------------------------------------------------
        if (me%nvc /= 0) then
            me%alp = -1.0_wp
            do i = 1, me%nac
                nci = ic(i)
                c1 = g(nci)
                ctc = me%ctam
                if (isc(nci) > 0) ctc = me%ctbm
                if (c1 > ctc) then
                    alp1 = dot_product(s(1:me%ndv), a(1:me%ndv, i))
                    alp1 = alp1*a(ndv2, i)
                    if (abs(alp1) >= 1.0e-20_wp) then
                        alp1 = -c1/alp1
                        if (alp1 > me%alp) me%alp = alp1
                    end if
                end if
            end do
        end if
        ! ------------------------------------------------------------------
        !                   LIMIT CHANCE TO ABOBJ1*OBJ
        ! ------------------------------------------------------------------
380     alp1 = 1.0e+20_wp
        si = abs(me%obj)
        if (si < 0.01_wp) si = 0.01_wp
        if (abs(me%slope) > 1.0e-20_wp) alp1 = me%abobj1*si/me%slope
        alp1 = abs(alp1)
        if (me%nvc > 0) alp1 = 10.0_wp*alp1
        if (alp1 < me%alp) me%alp = alp1
        ! ------------------------------------------------------------------
        !               LIMIT CHANGE IN VARIABLE TO ALPHAX
        ! ------------------------------------------------------------------
        alp11 = 1.0e+20_wp
        do i = 1, me%ndv
            si = abs(s(i))
            me%xi = abs(x(i))
            if (si >= 1.0e-10_wp .and. me%xi >= 0.1_wp) then
                alp1 = me%alphax*me%xi/si
                if (alp1 < alp11) alp11 = alp1
            end if
        end do
        if (me%nvc > 0) alp11 = 10.0_wp*alp11
        if (alp11 < me%alp) me%alp = alp11
        if (me%alp > 1.0e+20_wp) me%alp = 1.0e+20_wp
        if (me%alp <= 1.0e-20_wp) me%alp = 1.0e-20_wp
        if (me%iprint >= 3) then
            write (6, 9100)
            do i = 1, me%ndv, 6
                m1 = min(me%ndv, i + 5)
                write (6, 5400) i, (s(j), j=i, m1)
            end do
            write (6, 6400) me%slope, me%alp
        end if
        if (me%ncon > 0 .or. me%nside > 0) go to 420
        ! ------------------------------------------------------------------
        !       DO ONE-DIMENSIONAL SEARCH FOR UNCONSTRAINED FUNCTION
        ! ------------------------------------------------------------------
        me%jgoto = 0
410     call me%cnmn03(x, s, me%slope, me%alp, me%fff, me%a1, me%a2, me%a3, me%a4, &
                       me%f1, me%f2, me%f3, me%f4, me%app, me%ncal, me%kount, me%jgoto)
        me%igoto = 4
        if (me%jgoto > 0) go to 580
        me%jdir = 1
        ! PROCEED TO CONVERGENCE CHECK.
        go to 450
        ! ------------------------------------------------------------------
        !   SOLVE ONE-DIMENSIONAL SEARCH PROBLEM FOR CONSTRAINED FUNCTION
        ! ------------------------------------------------------------------
420     me%jgoto = 0

430     call me%cnmn06(x, vlb, vub, g, scal, df, s, g1, g2, me%ctam, me%ctbm, me%slope, me%alp, me%a2, me%a3, me%a4, &
                       me%f1, me%f2, me%f3, me%cv1, me%cv2, me%cv3, me%cv4, me%alpca, me%alpfes, me%alpln, me%alpmin, me%alpnc, &
                       me%alpsav, me%alpsid, me%alptot, isc, me%ncal, me%nvc, me%icount, me%igood1, &
                       me%igood2, me%igood3, me%igood4, me%ibest, me%iii, me%nlnc, me%jgoto)
        me%igoto = 5
        if (me%jgoto > 0) go to 580
        if (me%nac == 0) me%jdir = 1
        ! ------------------------------------------------------------------
        ! *******************     UPDATE ALPHAX   **************************
        ! ------------------------------------------------------------------
450     if (me%alp > 1.0e+19_wp) me%alp = 0.0_wp
        ! UPDATE ALPHAX TO BE AVERAGE OF MAXIMUM CHANGE IN X(I) AND ALHPAX.
        alp11 = 0.0_wp
        do i = 1, me%ndv
            si = abs(s(i))
            me%xi = abs(x(i))
            if (me%xi >= 1.0e-10_wp) then
                alp1 = me%alp*si/me%xi
                if (alp1 > alp11) alp11 = alp1
            end if
        end do
        alp11 = 0.5_wp*(alp11 + me%alphax)
        alp12 = 5.0_wp*me%alphax
        if (alp11 > alp12) alp11 = alp12
        me%alphax = alp11
        me%ncobj = me%ncobj + 1
        ! ABSOLUTE CHANGE IN OBJECTIVE.
        objd = me%obj1 - me%obj
        objb = abs(objd)
        if (objb < 1.0e-10_wp) objb = 0.0_wp
        if (me%nac == 0 .or. objb > 0.0_wp) me%ncobj = 0
        if (me%ncobj > 1) me%ncobj = 0
        ! ------------------------------------------------------------------
        !                              PRINT
        ! ------------------------------------------------------------------
        ! PRINT MOVE PARAMETER, NEW X-VECTOR AND CONSTRAINTS.
        if (me%iprint >= 3) then
            write (6, 9200) me%alp
        end if
        if (me%iprint >= 2) then
            if (objb <= 0.) then
                if (me%iprint == 2) write (6, 9300) me%iter, me%obj
                if (me%iprint > 2) write (6, 9400) me%obj
            else
                if (me%iprint /= 2) then
                    write (6, 9500) me%obj
                else
                    write (6, 9600) me%iter, me%obj
                end if
            end if
            write (6, 9800)
            do i = 1, me%ndv
                ff1 = 1.
                if (me%nscal /= 0) ff1 = scal(i)
                g1(i) = ff1*x(i)
            end do
            do i = 1, me%ndv, 6
                m1 = min(me%ndv, i + 5)
                write (6, 5400) i, (g1(j), j=i, m1)
            end do
            if (me%ncon /= 0) then
                write (6, 10000)
                do i = 1, me%ncon, 6
                    m1 = min(me%ncon, i + 5)
                    write (6, 5400) i, (g(j), j=i, m1)
                end do
            end if
        end if

        !  THE FOLLOWING CODE WAS ADDED ON 3/5/81
        !
        !  IT HAD NOT BEEN REPORTED AS A FIX TO MAOB
        !  BUT WAS SENT TO JEFF STROUD A YEAR AGO
        !  SEE OTHER COMMENTS IN CONMIN SUBROUTINE FOR DELETIONS OF CODE
        !  ON 3/5/81 PERTAINING TO THIS FIX

        ! CHECK FEASIBILITY
        if (me%ncon > 0) then
            nfeasct = 10
            ! added by slp 11/17/94
            do i = 1, me%ncon
                c1 = me%ctam
                if (isc(i) > 0) c1 = me%ctbm
                if (g(i) > c1) then
                    me%nfeas = me%nfeas + 1
                    go to 510
                end if
            end do
            if (me%nfeas > 0) me%abobj1 = 0.05_wp
!cc
            me%nfeas = 0
            me%phi = 5.0_wp
510         if (me%nfeas >= nfeasct) go to 520
        end if

        ! END OF INSERTED FIX

        ! ------------------------------------------------------------------
        !                      CHECK CONVERGENCE
        ! ------------------------------------------------------------------
        ! STOP IF ITER EQUALS ITMAX.
        if (me%iter < me%itmax) then
            ! ------------------------------------------------------------------
            !                 ABSOLUTE CHANGE IN OBJECTIVE
            ! ------------------------------------------------------------------
            objb = abs(objd)
            me%kobj = me%kobj + 1
            if (objb >= me%dabfun .or. me%nfeas > 0) me%kobj = 0
            ! ------------------------------------------------------------------
            !                 RELATIVE CHANGE IN OBJECTIVE
            ! ------------------------------------------------------------------
            if (abs(me%obj1) > 1.0e-10_wp) objd = objd/abs(me%obj1)
            me%abobj1 = .5*(abs(me%abobj) + abs(objd))
            me%abobj = abs(objd)
            me%iobj = me%iobj + 1
            if (me%nvc > 0 .or. objd >= me%delfun) me%iobj = 0
            if (me%iobj < me%itrm .and. me%kobj < me%itrm) then
                me%obj1 = me%obj
                ! ------------------------------------------------------------------
                !       REDUCE CT IF OBJECTIVE FUNCTION IS CHANGING SLOWLY
                ! ------------------------------------------------------------------
                if (me%iobj < 1 .or. me%nac == 0) go to 130
                me%ct = me%dct*me%ct
                me%ctl = me%ctl*me%dctl
                if (abs(me%ct) < me%ctmin) me%ct = -me%ctmin
                if (abs(me%ctl) < me%ctlmin) me%ctl = -me%ctlmin
                go to 130
            end if
        end if

520     if (me%nac >= n3) write (6, 10200)
        ! ------------------------------------------------------------------
        ! ****************  FINAL FUNCTION INFORMATION  ********************
        ! ------------------------------------------------------------------
        if (me%nscal /= 0) then
            ! UN-SCALE THE DESIGN VARIABLES.
            do i = 1, me%ndv
                me%xi = scal(i)
                if (me%nside /= 0) then
                    vlb(i) = me%xi*vlb(i)
                    vub(i) = me%xi*vub(i)
                end if
                x(i) = me%xi*x(i)
            end do
        end if
        ! ------------------------------------------------------------------
        !                       PRINT FINAL RESULTS
        ! ------------------------------------------------------------------
        if (me%iprint /= 0 .and. me%nac < n3) then
            write (6, 10300)
            write (6, 9500) me%obj
            write (6, 9800)
            do i = 1, me%ndv, 6
                m1 = min(me%ndv, i + 5)
                write (6, 5400) i, (x(j), j=i, m1)
            end do
            if (me%ncon /= 0) then
                write (6, 10000)
                do i = 1, me%ncon, 6
                    m1 = min(me%ncon, i + 5)
                    write (6, 5400) i, (g(j), j=i, m1)
                end do
                ! DETERMINE WHICH CONSTRAINTS ARE ACTIVE AND PRINT.
                me%nac = 0
                me%nvc = 0
                do i = 1, me%ncon
                    me%cta = me%ctam
                    if (isc(i) > 0) me%cta = me%ctbm
                    gi = g(i)
                    if (gi <= me%cta) then
                        if (gi < me%ct .and. isc(i) == 0) cycle
                        if (gi < me%ctl .and. isc(i) > 0) cycle
                        me%nac = me%nac + 1
                        ic(me%nac) = i
                    else
                        me%nvc = me%nvc + 1
                        ms1(me%nvc) = i
                    end if
                end do
                write (6, 5900) me%nac
                if (me%nac /= 0) then
                    write (6, 6000)
                    write (6, 10100) ic(1:me%nac)
                end if
                write (6, 6100) me%nvc
                if (me%nvc /= 0) then
                    write (6, 6000)
                    write (6, 10100) ms1(1:me%nvc)
                end if
            end if
            if (me%nside /= 0) then
                ! DETERMINE WHICH SIDE CONSTRAINTS ARE ACTIVE AND PRINT.
                me%nac = 0
                do i = 1, me%ndv
                    me%xi = x(i)
                    xid = vlb(i)
                    x12 = abs(xid)
                    if (x12 < 1.0_wp) x12 = 1.0_wp
                    gi = (xid - me%xi)/x12
                    if (gi >= -1.0e-6_wp) then
                        me%nac = me%nac + 1
                        ms1(me%nac) = -i
                    end if
                    xid = vub(i)
                    x12 = abs(xid)
                    if (x12 < 1.0_wp) x12 = 1.0_wp
                    gi = (me%xi - xid)/x12
                    if (gi >= -1.0e-6_wp) then
                        me%nac = me%nac + 1
                        ms1(me%nac) = i
                    end if
                end do
                write (6, 6200) me%nac
                if (me%nac /= 0) then
                    write (6, 6300)
                    write (6, 10100) (ms1(j), j=1, me%nac)
                end if
            end if
            write (6, 6800)
            if (me%iter >= me%itmax) write (6, 6900)
            if (me%nfeas >= nfeasct) write (6, 7000)
            if (me%iobj >= me%itrm) write (6, 7100) me%itrm
            if (me%kobj >= me%itrm) write (6, 7200) me%itrm
            write (6, 7300) me%iter
            write (6, 10400) me%ncal(1)
            if (me%ncon > 0) write (6, 10500) me%ncal(1)
            if (me%nfdg /= 0) write (6, 10600) me%ncal(2)
            if (me%ncon > 0 .and. me%nfdg == 1) write (6, 10700) me%ncal(2)
        end if
        ! ------------------------------------------------------------------
        !               RE-SET BASIC PARAMETERS TO INPUT VALUES
        ! ------------------------------------------------------------------
        me%itrm = me%idm1
        me%itmax = me%idm2
        me%icndir = me%idm3
        me%delfun = me%dm1
        me%dabfun = me%dm2
        me%ct = me%dm3
        me%ctmin = me%dm4
        me%ctl = me%dm5
        me%ctlmin = me%dm6
        me%theta = me%dm7
        me%phi = me%dm8
        me%fdch = me%dm9
        me%fdchm = me%dm10
        me%abobj1 = me%dm11
        me%alphax = me%dm12
        me%igoto = 0

580     if (me%nscal == 0 .or. me%igoto == 0) return
        ! UN-SCALE VARIABLES.
        do i = 1, me%ndv
            c(i) = x(i)
            x(i) = x(i)*scal(i)
        end do
        return

        !  ------------------------------------------------------------------
        !                             FORMATS
        !  ------------------------------------------------------------------

5000    format(//t6, &
                'A COMPLETELY UNCONSTRAINED FUNCTION WITH A LINEAR OBJECTIVE IS SPECIFIED'// &
                t11, 'LINOBJ =', i5/t11, 'NCON   =', i5/ &
                t11, 'NSIDE  =', i5//t6, 'CONTROL RETURNED TO CALLING PROGRAM')
5100    format(//t6, &
                'CONMIN HAS ACHIEVED A SOLUTION OF OBJ LESS THAN -1.0E+40'/ &
                t6, 'SOLUTION APPEARS TO BE UNBOUNDED'/t6, 'OPTIMIZATION IS TERMINATED')
5200    format(t6, 'CONSTRAINT NUMBER', i5)
5300    format(t6, 'SIDE CONSTRAINT ON VARIABLE', i5)
5400    format(t4, i5, ')  ', 6e13.5)
5500    format(/t6, 'LINEAR CONSTRAINT IDENTIFIERS (ISC)'/ &
                t6, 'NON-ZERO INDICATES LINEAR CONSTRAINT')
5600    format(t4, i5, ')  ', 15i5)
5700    format(/t6, 'ALL CONSTRAINTS ARE LINEAR')
5800    format(/t6, 'ALL CONSTRAINTS ARE NON-LINEAR')
5900    format(/t6, 'THERE ARE', i5, ' ACTIVE CONSTRAINTS')
6000    format(t6, 'CONSTRAINT NUMBERS ARE')
6100    format(/t6, 'THERE ARE', i5, ' VIOLATED CONSTRAINTS')
6200    format(/t6, 'THERE ARE', i5, ' ACTIVE SIDE CONSTRAINTS')
6300    format(t6, 'DECISION VARIABLES AT LOWER OR UPPER BOUNDS', &
               ' (MINUS INDICATES LOWER BOUND)')
6400    format(/t6, 'ONE-DIMENSIONAL SEARCH'/t6, 'INITIAL SLOPE =', e12.4, &
                '  PROPOSED ALPHA =', e12.4)
6500    format(//t6, '* * CONMIN DETECTS VLB(I)>VUB(I)'/ &
                t6, 'FIX IS SET X(I)=VLB(I)=VUB(I) = .5*(VLB(I)+VUB(I) FOR I =', i5)
6600    format(//t6, '* * CONMIN DETECTS INITIAL X(I)<VLB(I)'/t6, &
                'X(I) =', e12.4, '  VLB(I) =', e12.4/t6, &
                'X(I) IS SET EQUAL TO VLB(I) FOR I =', i5)
6700    format(//t6, '* * CONMIN DETECTS INITIAL X(I)>VUB(I)'/t6, &
                'X(I) =', e12.4, '  VUB(I) =', e12.4/t6, &
                'X(I) IS SET EQUAL TO VUB(I) FOR I =', i5)
6800    format(/t6, 'TERMINATION CRITERION')
6900    format(t11, 'ITER EQUALS ITMAX')
7000    format(t11, &
               'NFEASCT CONSECUTIVE ITERATIONS FAILED TO PRODUCE A FEASIBLE DESIGN')
7100    format(t11, 'ABS(1-OBJ(I-1)/OBJ(I)) LESS THAN DELFUN FOR', i3, &
               ' ITERATIONS')
7200    format(t11, 'ABS(OBJ(I)-OBJ(I-1))   LESS THAN DABFUN FOR', i3, &
               ' ITERATIONS')
7300    format(/t6, 'NUMBER OF ITERATIONS =', i5)
7400    format(/t6, 'CONSTRAINT PARAMETER, BETA =', e14.5)

7500    format(//t13, '* * * * * * * * * * * * * * * * * * * * * * * * * * * '/t13,&
                '*', t65, '*'/t13, '*', t34, &
                'C O N M I N', t65, '*'/t13, '*', t65, '*'/t13, '*', t29, &
                ' FORTRAN PROGRAM FOR ', t65, '*'/t13, '*', t65, '*'/t13, '*', t23, &
                'CONSTRAINED FUNCTION MINIMIZATION', t65, '*'/t13, '*', t65, '*'/ &
                t13, '* * * * * * * * * * * * * * * * * * * * * * * * * * * ')
7600    format(//t6, 'CONSTRAINED FUNCTION MINIMIZATION'// &
                t6, 'CONTROL PARAMETERS')
7700    format(/t6, 'IPRINT  NDV    ITMAX    NCON    NSIDE  ICNDIR   NSCAL NFDG'/ &
                8i8//t6, 'LINOBJ  ITRM     N1      N2      N3      N4      N5'/8i8)
7800    format(/t10, 'FDCH', t26, 'FDCHM', t42, 'ALPHAX', t58, 'ABOBJ1'/ &
                ' ', *('  ', e14.5))
7900    format(/t10, 'CT', t26, 'CTMIN', t42, 'CTL', t58, 'CTLMIN'/ &
                ' ', *('  ', e14.5)// &
                t10, 'THETA', t26, 'PHI', t42, 'DELFUN', t58, 'DABFUN'/ &
                ' ', *('  ', e14.5))
8000    format(/t6, 'LOWER BOUNDS ON DECISION VARIABLES (VLB)')
8100    format(/t6, 'UPPER BOUNDS ON DECISION VARIABLES (VUB)')
8200    format(//t6, 'UNCONSTRAINED FUNCTION MINIMIZATION'//t6, &
                'CONTROL PARAMETERS')
8300    format(/t6, 'SCALING VECTOR (SCAL)')
8400    format(//t6, 'BEGIN ITERATION NUMBER', i5)
8500    format(/t6, 'CT =', e14.5, '     CTL =', e14.5, '     PHI =', e14.5)
8600    format(/t6, 'NEW SCALING VECTOR (SCAL)')
8700    format(/t6, 'GRADIENT OF OBJ')
8800    format(/t6, 'GRADIENTS OF ACTIVE AND VIOLATED CONSTRAINTS')
8900    format(' ')
9000    format(/t6, 'PUSH-OFF FACTORS, (THETA(I), I=1,NAC)')
9100    format(/t6, 'SEARCH DIRECTION (S-VECTOR)')
9200    format(/t6, 'CALCULATED ALPHA =', e14.5)
9300    format(//t6, 'ITER =', i5, '     OBJ =', e14.5, '     NO CHANGE IN OBJ')
9400    format(/t6, 'OBJ =', e15.6, '     NO CHANGE ON OBJ')
9500    format(/t6, 'OBJ =', e15.6)
9600    format(//t6, 'ITER =', i5, '     OBJ =', e14.5)
9700    format(//t6, 'INITIAL FUNCTION INFORMATION'//t6, 'OBJ =', e15.6)
9800    format(/t6, 'DECISION VARIABLES (X-VECTOR)')
9900    format(t4, 7e13.4)
10000   format(/t6, 'CONSTRAINT VALUES (G-VECTOR)')
10100   format(t6, 15i5)
10200   format(/t6, 'THE NUMBER OF ACTIVE AND VIOLATED CONSTRAINTS EXCEEDS N3-1.'/ &
                t6, 'DIMENSIONED SIZE OF MATRICES A AND B AND VECTOR IC IS INSUFFICIENT'/ &
                t6, 'OPTIMIZATION TERMINATED AND CONTROL RETURNED TO MAIN PROGRAM.')
10300   format(///'    FINAL OPTIMIZATION INFORMATION')
10400   format(/t6, 'OBJECTIVE FUNCTION WAS EVALUATED        ', i5, '  TIMES')
10500   format(/t6, 'CONSTRAINT FUNCTIONS WERE EVALUATED', i10, '  TIMES')
10600   format(/t6, 'GRADIENT OF OBJECTIVE WAS CALCULATED', i9, '  TIMES')
10700   format(/t6, 'GRADIENTS OF CONSTRAINTS WERE CALCULATED', i5, '  TIMES')
    end subroutine conmin

    subroutine cnmn01(me, jgoto, x, df, g, isc, ic, a, g1, vub, scal, ncal, dx, &
                      dx1, fi, xi, iii, n1, n2, n3)

        !!  Routine to calculate gradient information by finite difference.
        !!
        !!  BY G. N. VANDERPLAATS, JUNE, 1972.

        class(conmin_class), intent(inout) :: me
        integer, intent(inout)    :: jgoto
        real(wp), intent(in)      :: vub(:), scal(:)
        integer, intent(in)       :: n1, n2, n3
        real(wp), intent(inout)   :: x(:), df(:), g(n2), a(n1, n3), g1(n2)
        integer, intent(inout)    :: ic(n3), ncal(2)
        integer, intent(in)       :: isc(n2)
        real(wp), intent(out)     :: dx, dx1, fi, xi
        integer, intent(out)      :: iii

        real(wp)  :: fdch1, x1
        integer   :: i, i1, inf, j

        if (jgoto /= 1) then
            if (jgoto == 2) go to 40
            me%infog = 0
            inf = me%info
            me%nac = 0
            if (me%linobj == 0 .or. me%iter <= 1) then
                ! ------------------------------------------------------------------
                !                GRADIENT OF LINEAR OBJECTIVE
                ! ------------------------------------------------------------------
                if (me%nfdg == 2) jgoto = 1
                if (me%nfdg == 2) return
            end if
        end if
        jgoto = 0
        if (me%nfdg == 2 .and. me%ncon == 0) return
        if (me%ncon /= 0) then
            ! ------------------------------------------------------------------
            !   * * * DETERMINE WHICH CONSTRAINTS ARE ACTIVE OR VIOLATED * * *
            ! ------------------------------------------------------------------
            do i = 1, me%ncon
                if (g(i) >= me%ct) then
                    if (isc(i) <= 0 .or. g(i) >= me%ctl) then
                        me%nac = me%nac + 1
                        if (me%nac >= n3) return
                        ic(me%nac) = i
                    end if
                end if
            end do
            if (me%nfdg == 2 .and. me%nac == 0) return
            if (me%linobj > 0 .and. me%iter > 1 .and. me%nac == 0) return
            ! ------------------------------------------------------------------
            !              STORE VALUES OF CONSTRAINTS IN G1
            ! ------------------------------------------------------------------
            g1(1:me%ncon) = g(1:me%ncon)
        end if
        jgoto = 0
        if (me%nac == 0 .and. me%nfdg == 2) return
        ! ------------------------------------------------------------------
        !                        CALCULATE GRADIENTS
        ! ------------------------------------------------------------------
        me%infog = 1
        me%info = 1
        fi = me%obj
        iii = 0

30      iii = iii + 1
        xi = x(iii)
        dx = me%fdch*xi
        dx = abs(dx)
        fdch1 = me%fdchm
        if (me%nscal /= 0) fdch1 = me%fdchm/scal(iii)
        if (dx < fdch1) dx = fdch1
        x1 = xi + dx
        if (me%nside /= 0) then
            if (x1 > vub(iii)) dx = -dx
        end if
        dx1 = 1.0_wp/dx
        x(iii) = xi + dx
        ncal(1) = ncal(1) + 1
        !  ------------------------------------------------------------------
        !                      FUNCTION EVALUATION
        !  ------------------------------------------------------------------
        jgoto = 2
        return

40      x(iii) = xi
        if (me%nfdg == 0) df(iii) = dx1*(me%obj - fi)
        if (me%nac /= 0) then
            ! ------------------------------------------------------------------
            !         DETERMINE GRADIENT COMPONENTS OF ACTIVE CONSTRAINTS
            ! ------------------------------------------------------------------
            do j = 1, me%nac
                i1 = ic(j)
                a(iii, j) = dx1*(g(i1) - g1(i1))
            end do
        end if
        if (iii < me%ndv) go to 30
        me%infog = 0
        me%info = inf
        jgoto = 0
        me%obj = fi
        if (me%ncon == 0) return
        ! ------------------------------------------------------------------
        !         STORE CURRENT CONSTRAINT VALUES BACK IN G-VECTOR
        ! ------------------------------------------------------------------
        g(1:me%ncon) = g1(1:me%ncon)

    end subroutine cnmn01

    subroutine cnmn02(me, ncalc, slope, dftdf1, df, s)

        !!  Routine to determine conjugate direction vector or direction
        !!  of steepest descent for unconstrained function minimization.
        !!
        !!  BY G. N. VANDERPLAATS, APRIL, 1972.
        !!
        !!  Conjugate direction is found by fletcher-reeves algorithm.

        class(conmin_class), intent(inout) :: me
        integer, intent(inout)     :: ncalc !! NCALC = CALCULATION CONTROL.
                                            !!
                                            !!  * NCALC = 0,     S = STEEPEST DESCENT.
                                            !!  * NCALC = 1,     S = CONJUGATE DIRECTION.
        real(wp), intent(out)     :: slope
        real(wp), intent(inout)   :: dftdf1
        real(wp), intent(in)      :: df(:)
        real(wp), intent(inout)   :: s(:)

        integer  :: i
        real(wp) :: beta, dfi, dftdf, s1, s2, si
        logical  :: fletcher_reeves !! if the fletcher-reeves conjugate direction was computed

        dftdf = sum(df(1:me%ndv)**2) ! CALCULATE NORM OF GRADIENT VECTOR

        ! FIND DIRECTION S
        fletcher_reeves = .false.
        if (ncalc == 1) then
            if (dftdf1 >= 1.0e-20_wp) then
                ! FIND FLETCHER-REEVES CONJUGATE DIRECTION
                beta = dftdf/dftdf1
                slope = 0.0_wp
                do i = 1, me%ndv
                    dfi = df(i)
                    si = beta*s(i) - dfi
                    slope = slope + si*dfi
                    s(i) = si
                end do
                fletcher_reeves = .true.
            end if
        end if
        if (.not. fletcher_reeves) then
            ncalc = 0
            ! CALCULATE DIRECTION OF STEEPEST DESCENT
            s(1:me%ndv) = -df(1:me%ndv)
            slope = -dftdf
        end if
        ! NORMALIZE S TO MAX ABS VALUE OF UNITY
        s1 = 0.0_wp
        do i = 1, me%ndv
            s2 = abs(s(i))
            if (s2 > s1) s1 = s2
        end do
        if (s1 < 1.0e-20_wp) s1 = 1.0e-20_wp
        s1 = 1.0_wp/s1
        dftdf1 = dftdf*s1
        s(1:me%ndv) = s1*s(1:me%ndv)
        slope = s1*slope

    end subroutine cnmn02

    subroutine cnmn03(me, x, s, slope, alp, fff, a1, a2, a3, a4, f1, f2, f3, f4, app, &
                      ncal, kount, jgoto)

        !!  Routine to solve one-dimensional search in unconstrained
        !!  minimization using 2-point quadratic interpolation, 3-point
        !!  cubic interpolation and 4-point cubic interpolation.
        !!
        !!  BY G. N. VANDERPLAATS, APRIL, 1972.
        !!
        !!  OBJ = INITIAL FUNCTION VALUE.

        class(conmin_class), intent(inout) :: me
        real(wp), intent(inout)    :: x(:), s(:)
        real(wp), intent(inout)    :: slope !! SLOPE = INITIAL FUNCTION SLOPE = S-TRANSPOSE TIMES DF.
                                            !! SLOPE MUST BE NEGATIVE.
        real(wp), intent(inout)    :: alp !! PROPOSED MOVE PARAMETER
        real(wp), intent(inout)    :: fff, a1, a2, a3, a4, f1, f2, f3, f4
        real(wp), intent(inout)    :: app
        integer, intent(inout)     :: ncal(2)
        integer, intent(out)       :: kount
        integer, intent(inout)     :: jgoto

        real(wp)  :: aa, ab, ab2, ab3, ap, ff
        integer :: i, ii
        real(wp), parameter :: zro = 0.0_wp

        if (jgoto /= 0) then
            select case (jgoto)
            case (1); go to 30
            case (2); go to 50
            case (3); go to 80
            case (4); go to 110
            case (5); go to 150
            case (6); go to 190
            case (7); go to 230
            end select
        end if
        ! ------------------------------------------------------------------
        !                 INITIAL INFORMATION  (ALPHA=0)
        ! ------------------------------------------------------------------
        if (slope >= 0.0_wp) then
            alp = 0.0_wp
            return
        end if
        if (me%iprint > 4) write (6, 5000)
        fff = me%obj
        a1 = 0.0_wp
        f1 = me%obj
        a2 = alp
        a3 = 0.0_wp
        f3 = 0.0_wp
        ap = a2
        kount = 0
        ! ------------------------------------------------------------------
        !        MOVE A DISTANCE AP*S AND UPDATE FUNCTION VALUE
        ! ------------------------------------------------------------------
10      kount = kount + 1
        do i = 1, me%ndv
            x(i) = x(i) + ap*s(i)
        end do
        if (me%iprint > 4) write (6, 5100) ap
        if (me%iprint > 4) write (6, 5200) x(1:me%ndv)
        ncal(1) = ncal(1) + 1
        jgoto = 1
        return

30      f2 = me%obj
        if (me%iprint > 4) write (6, 5300) f2
        if (f2 < f1) go to 90
        ! ------------------------------------------------------------------
        !                 CHECK FOR ILL-CONDITIONING
        ! ------------------------------------------------------------------
        if (kount <= 5) then
            ff = 2.0_wp*abs(f1)
            if (f2 < ff) go to 60
            ff = 5.0_wp*abs(f1)
            if (f2 >= ff) then
                a2 = 0.5_wp*a2
                ap = -a2
                alp = a2
                go to 10
            end if
        end if
        f3 = f2
        a3 = a2
        a2 = 0.5_wp*a2
        ! ------------------------------------------------------------------
        !             UPDATE DESIGN VECTOR AND FUNCTION VALUE
        ! ------------------------------------------------------------------
        ap = a2 - alp
        alp = a2
        do i = 1, me%ndv
            x(i) = x(i) + ap*s(i)
        end do
        if (me%iprint > 4) write (6, 5100) a2
        if (me%iprint > 4) write (6, 5200) x(1:me%ndv)
        ncal(1) = ncal(1) + 1
        jgoto = 2
        return

50      f2 = me%obj
        if (me%iprint > 4) write (6, 5300) f2
        ! PROCEED TO CUBIC INTERPOLATION.
        go to 130
        ! ------------------------------------------------------------------
        ! **********        2-POINT QUADRATIC INTERPOLATION       **********
        ! ------------------------------------------------------------------
60      ii = 1
        call cnmn04(ii, app, zro, a1, f1, slope, a2, f2, zro, zro, zro, zro)
        if (app < zro .or. app > a2) go to 90
        f3 = f2
        a3 = a2
        a2 = app
        ! ------------------------------------------------------------------
        !              UPDATE DESIGN VECTOR AND FUNCTION VALUE
        ! ------------------------------------------------------------------
        ap = a2 - alp
        alp = a2
        do i = 1, me%ndv
            x(i) = x(i) + ap*s(i)
        end do
        if (me%iprint > 4) write (6, 5100) a2
        if (me%iprint > 4) write (6, 5200) x(1:me%ndv)
        ncal(1) = ncal(1) + 1
        jgoto = 3
        return

80      f2 = me%obj
        if (me%iprint > 4) write (6, 5300) f2
        go to 120

90      a3 = 2.*a2
        ! ------------------------------------------------------------------
        !              UPDATE DESIGN VECTOR AND FUNCTION VALUE
        ! ------------------------------------------------------------------
        ap = a3 - alp
        alp = a3
        do i = 1, me%ndv
            x(i) = x(i) + ap*s(i)
        end do
        if (me%iprint > 4) write (6, 5100) a3
        if (me%iprint > 4) write (6, 5200) x(1:me%ndv)
        ncal(1) = ncal(1) + 1
        jgoto = 4
        return

110     f3 = me%obj
        if (me%iprint > 4) write (6, 5300) f3

120     if (f3 < f2) go to 170

        ! ------------------------------------------------------------------
        ! **********       3-POINT CUBIC INTERPOLATION      **********
        ! ------------------------------------------------------------------
130     ii = 3
        call cnmn04(ii, app, zro, a1, f1, slope, a2, f2, a3, f3, zro, zro)
        if (app < zro .or. app > a3) go to 170
        ! ------------------------------------------------------------------
        ! UPDATE DESIGN VECTOR AND FUNCTION VALUE.
        ! ------------------------------------------------------------------
        ap = app - alp
        alp = app
        x(1:me%ndv) = x(1:me%ndv) + ap*s(1:me%ndv)
        if (me%iprint > 4) write (6, 5100) alp
        if (me%iprint > 4) write (6, 5200) x(1:me%ndv)
        ncal(1) = ncal(1) + 1
        jgoto = 5
        return

150     if (me%iprint > 4) write (6, 5300) me%obj
        ! ------------------------------------------------------------------
        !                     CHECK CONVERGENCE
        ! ------------------------------------------------------------------
        aa = 1.-app/a2
        ab2 = abs(f2)
        ab3 = abs(me%obj)
        ab = ab2
        if (ab3 > ab) ab = ab3
        if (ab < 1.0e-15_wp) ab = 1.0e-15_wp
        ab = (ab2 - ab3)/ab
        if (abs(ab) < 1.0e-15_wp .and. abs(aa) < 0.001_wp) go to 260
        a4 = a3
        f4 = f3
        a3 = app
        f3 = me%obj
        if (a3 > a2) go to 200
        a3 = a2
        f3 = f2
        a2 = app
        f2 = me%obj
        go to 200

        ! ------------------------------------------------------------------
        ! **********        4-POINT CUBIC INTERPOLATION       **********
        ! ------------------------------------------------------------------
170     a4 = 2.*a3
        ! UPDATE DESIGN VECTOR AND FUNCTION VALUE.
        ap = a4 - alp
        alp = a4
        x(1:me%ndv) = x(1:me%ndv) + ap*s(1:me%ndv)
        if (me%iprint > 4) write (6, 5100) alp
        if (me%iprint > 4) write (6, 5200) (x(i), i=1, me%ndv)
        ncal(1) = ncal(1) + 1
        jgoto = 6
        return

190     f4 = me%obj
        if (me%iprint > 4) write (6, 5300) f4
        if (f4 <= f3) then
            a1 = a2
            f1 = f2
            a2 = a3
            f2 = f3
            a3 = a4
            f3 = f4
            go to 170
        end if

200     ii = 4
        call cnmn04(ii, app, a1, a1, f1, slope, a2, f2, a3, f3, a4, f4)
        if (app <= a1) then
            ap = a1 - alp
            alp = a1
            me%obj = f1
            do i = 1, me%ndv
                x(i) = x(i) + ap*s(i)
            end do
            go to 240
        end if
        !  ------------------------------------------------------------------
        !              UPDATE DESIGN VECTOR AND FUNCTION VALUE
        !  ------------------------------------------------------------------
        ap = app - alp
        alp = app
        x(1:me%ndv) = x(1:me%ndv) + ap*s(1:me%ndv)
        if (me%iprint > 4) write (6, 5100) alp
        if (me%iprint > 4) write (6, 5200) x(1:me%ndv)
        ncal(1) = ncal(1) + 1
        jgoto = 7
        return

230     if (me%iprint > 4) write (6, 5300) me%obj

        ! ------------------------------------------------------------------
        !                CHECK FOR ILL-CONDITIONING
        ! ------------------------------------------------------------------
240     if (me%obj <= f2 .and. me%obj <= f3) then
            if (me%obj <= f1) go to 260
            ap = a1 - alp
            alp = a1
            me%obj = f1
        else
            if (f2 >= f3) then
                me%obj = f3
                ap = a3 - alp
                alp = a3
            else
                me%obj = f2
                ap = a2 - alp
                alp = a2
            end if
        end if
        ! ------------------------------------------------------------------
        !                   UPDATE DESIGN VECTOR
        ! ------------------------------------------------------------------
        x(1:me%ndv) = x(1:me%ndv) + ap*s(1:me%ndv)

        ! ------------------------------------------------------------------
        !                 CHECK FOR MULTIPLE MINIMA
        ! ------------------------------------------------------------------
260     if (me%obj > fff) then
            ! INITIAL FUNCTION IS MINIMUM.
            x(1:me%ndv) = x(1:me%ndv) - alp*s(1:me%ndv)
            alp = 0.0_wp
            me%obj = fff
        end if
        jgoto = 0
        return

        ! ------------------------------------------------------------------
        !                             FORMATS
        ! ------------------------------------------------------------------

5000    format(///t6, &
                '* * * UNCONSTRAINED ONE-DIMENSIONAL SEARCH INFORMATION * * *')
5100    format(/t6, 'ALPHA =', e14.5/t6, 'X-VECTOR')
5200    format(t6, 6e13.5)
5300    format(/t6, 'OBJ =', e14.5)

    end subroutine cnmn03

    subroutine cnmn04(ii, xbar, eps, x1, y1, slope, x2, y2, x3, y3, x4, y4)

        !!  Routine to find first `xbar>=eps` corresponding to a minimum
        !!  of a one-dimensional real function by polynomial interpolation.
        !!
        !!  BY G. N. VANDERPLAATS, APRIL, 1972.
        !!
        !!  IF REQUIRED MINIMUM ON Y DOES NOT EXITS, OR THE FUNCTION IS
        !!  ILL-CONDITIONED, XBAR = EPS-1.0 WILL BE RETURNED AS AN ERROR INDICATOR.
        !!  IF DESIRED INTERPOLATION IS ILL-CONDITIONED, A LOWER ORDER
        !!  INTERPOLATION, CONSISTANT WITH INPUT DATA, WILL BE ATTEMPTED,
        !!  AND II WILL BE CHANGED ACCORDINGLY.

        integer, intent(inout)   :: ii !! CALCULATION CONTROL:
                                       !!
                                       !!  1.  2-POINT QUADRATIC INTERPOLATION, GIVEN X1, Y1, SLOPE, X2 AND Y2.
                                       !!  2.  3-POINT QUADRATIC INTERPOLATION, GIVEN X1, Y1, X2, Y2, X3 AND Y3.
                                       !!  3.  3-POINT CUBIC INTERPOLATION, GIVEN X1, Y1, SLOPE, X2, Y2,
                                       !!      X3 AND Y3.
                                       !!  4.  4-POINT CUBIC INTERPOLATION, GIVEN X1, Y1, X2, Y2, X3,
                                       !!      Y3, X4 AND Y4.
        real(wp), intent(out)   :: xbar
        real(wp), intent(in)    :: eps !! may be negative
        real(wp), intent(in)    :: x1, y1, slope, x2, y2, x3, y3, x4, y4

        real(wp)  :: aa, bac, bb, cc, dnom, dx, q1, q2, q3, q4, q5, q6, qq, &
                     x11, x111, x21, x22, x222, x31, x32, x33, x41, x42, x44, xbar1
        integer   :: nslop

        xbar1 = eps - 1.0_wp
        xbar = xbar1
        x21 = x2 - x1
        if (abs(x21) < 1.0e-20_wp) return
        nslop = mod(ii, 2)
        select case (ii)
        case (1); go to 10
        case (2); go to 20
        case (3); go to 30
        case (4); go to 40
        end select

        ! ------------------------------------------------------------------
        !             II=1: 2-POINT QUADRATIC INTERPOLATION
        ! ------------------------------------------------------------------
10      ii = 1
        dx = x1 - x2
        if (abs(dx) < 1.0e-20_wp) return
        aa = (slope + (y2 - y1)/dx)/dx
        if (aa < 1.0e-20_wp) return
        bb = slope - 2.0_wp*aa*x1
        xbar = -0.5_wp*bb/aa
        if (xbar < eps) xbar = xbar1
        return

        ! ------------------------------------------------------------------
        !             II=2: 3-POINT QUADRATIC INTERPOLATION
        ! ------------------------------------------------------------------
20      ii = 2
        x21 = x2 - x1
        x31 = x3 - x1
        x32 = x3 - x2
        qq = x21*x31*x32
        if (abs(qq) < 1.0e-20_wp) return
        aa = (y1*x32 - y2*x31 + y3*x21)/qq
        if (aa >= 1.0e-20_wp) then
            bb = (y2 - y1)/x21 - aa*(x1 + x2)
            xbar = -0.5_wp*bb/aa
            if (xbar < eps) xbar = xbar1
            return
        end if
        if (nslop == 0) return
        go to 10

        ! ------------------------------------------------------------------
        !               II=3: 3-POINT CUBIC INTERPOLATION
        ! ------------------------------------------------------------------
30      ii = 3
        x21 = x2 - x1
        x31 = x3 - x1
        x32 = x3 - x2
        qq = x21*x31*x32
        if (abs(qq) < 1.0e-20_wp) return
        x11 = x1*x1
        dnom = x2*x2*x31 - x11*x32 - x3*x3*x21
        if (abs(dnom) < 1.0e-20_wp) go to 20
        aa = ((x31*x31*(y2 - y1) - x21*x21*(y3 - y1))/(x31*x21) - slope*x32)/dnom
        if (abs(aa) < 1.0e-20_wp) go to 20
        bb = ((y2 - y1)/x21 - slope - aa*(x2*x2 + x1*x2 - 2.0_wp*x11))/x21
        cc = slope - 3.0_wp*aa*x11 - 2.0_wp*bb*x1
        bac = bb*bb - 3.0_wp*aa*cc
        if (bac < 0.0_wp) go to 20
        bac = sqrt(bac)
        xbar = (bac - bb)/(3.0_wp*aa)
        if (xbar < eps) xbar = eps
        return

        ! ------------------------------------------------------------------
        !                II=4: 4-POINT CUBIC INTERPOLATION
        ! ------------------------------------------------------------------
40      x21 = x2 - x1
        x31 = x3 - x1
        x41 = x4 - x1
        x32 = x3 - x2
        x42 = x4 - x2
        x11 = x1*x1
        x22 = x2*x2
        x33 = x3*x3
        x44 = x4*x4
        x111 = x1*x11
        x222 = x2*x22
        q2 = x31*x21*x32
        if (abs(q2) < 1.0e-30_wp) return
        q1 = x111*x32 - x222*x31 + x3*x33*x21
        q4 = x111*x42 - x222*x41 + x4*x44*x21
        q5 = x41*x21*x42
        dnom = q2*q4 - q1*q5
        if (abs(dnom) >= 1.0e-30_wp) then
            q3 = y3*x21 - y2*x31 + y1*x32
            q6 = y4*x21 - y2*x41 + y1*x42
            aa = (q2*q6 - q3*q5)/dnom
            bb = (q3 - q1*aa)/q2
            cc = (y2 - y1 - aa*(x222 - x111))/x21 - bb*(x1 + x2)
            bac = bb*bb - 3.0_wp*aa*cc
            if (abs(aa) >= 1.0e-20_wp .and. bac >= 0.0_wp) then
                bac = sqrt(bac)
                xbar = (bac - bb)/(3.0_wp*aa)
                if (xbar < eps) xbar = xbar1
                return
            end if
        end if
        if (nslop == 1) go to 30
        go to 20

    end subroutine cnmn04

    subroutine cnmn05(me, g, df, a, s, b, c, slope, phi, isc, ic, ms1, nvc, n1, n2, n3, n4, n5)

        !!  Routine to solve direction finding problem in modified method of
        !!  feasible directions.
        !!
        !!  BY G. N. VANDERPLAATS, MAY, 1972.
        !!
        !!  NORM OF S VECTOR USED HERE IS S-TRANSPOSE TIMES S<=1.
        !!  IF NVC = 0 FIND DIRECTION BY ZOUTENDIJK'S METHOD.  OTHERWISE
        !!  FIND MODIFIED DIRECTION.

        class(conmin_class), intent(inout) :: me
        integer, intent(in)       :: n1, n2, n3, n4, n5
        real(wp), intent(inout)   :: df(:), g(n2), a(n1, n3), s(:), c(n4), b(n3, n3)
        real(wp), intent(inout)   :: slope, phi
        integer, intent(inout)    :: isc(n2), ic(n3), ms1(n5), nvc

        real(wp)  :: a1, c1, ct1, ct2, cta, ctam, ctb, ctbm, ctc, ctd, gg, &
                     s1, sg, thmax, tht
        integer   :: i, j, j1, k, nac1, nci, ncj, ndb, ndv1, ndv2, ner

        ! ------------------------------------------------------------------
        ! ***  NORMALIZE GRADIENTS, CALCULATE THETA'S AND DETERMINE NVC  ***
        ! ------------------------------------------------------------------
        ndv1 = me%ndv + 1
        ndv2 = me%ndv + 2
        nac1 = me%nac + 1
        nvc = 0
        thmax = 0.0_wp
        cta = abs(me%ct)
        ct1 = 1.0_wp/cta
        ctam = abs(me%ctmin)
        ctb = abs(me%ctl)
        ct2 = 1.0_wp/ctb
        ctbm = abs(me%ctlmin)
        a1 = 1.0_wp
        do i = 1, me%nac
            ! CALCULATE THETA
            nci = ic(i)
            ncj = 1
            if (nci <= me%ncon) ncj = isc(nci)
            c1 = g(nci)
            ctd = ct1
            ctc = ctam
            if (ncj > 0) then
                ctc = ctbm
                ctd = ct2
            end if
            if (c1 > ctc) nvc = nvc + 1
            tht = 0.
            gg = 1.+ctd*c1
            if (ncj == 0 .or. c1 > ctc) tht = me%theta*gg*gg
            if (tht > 50.0_wp) tht = 50.0_wp
            if (tht > thmax) thmax = tht
            a(ndv1, i) = tht
            ! ------------------------------------------------------------------
            !                NORMALIZE GRADIENTS OF CONSTRAINTS
            ! ------------------------------------------------------------------
            a(ndv2, i) = 1.0_wp
            if (nci <= me%ncon) then
                a1 = 0.
                do j = 1, me%ndv
                    a1 = a1 + a(j, i)**2
                end do
                if (a1 < 1.0e-20_wp) a1 = 1.0e-20_wp
                a1 = sqrt(a1)
                a(ndv2, i) = a1
                a1 = 1.0_wp/a1
                do j = 1, me%ndv
                    a(j, i) = a1*a(j, i)
                end do
            end if
        end do
        ! ------------------------------------------------------------------
        ! CHECK FOR ZERO GRADIENT.  PROGRAM CHANGE-FEB, 1981, GV.
        ! ------------------------------------------------------------------
        i = 0
40      i = i + 1
50      if (a(ndv2, i) <= 1.0e-6_wp) then
            ! ZERO GRADIENT IS FOUND.  WRITE ERROR MESSAGE.
            if (me%iprint >= 2) write (6, 5000) ic(i)
            ! REDUCE NAC BY ONE.
            me%nac = me%nac - 1
            ! SHIFT COLUMNS OF A AND ROWS OF IC IF I<=NAC.
            if (i > me%nac) go to 80
            ! SHIFT.
            do j = i, me%nac
                j1 = j + 1
                ic(j) = ic(j1)
                do k = 1, ndv2
                    a(k, j) = a(k, j1)
                end do
            end do
            if (i <= me%nac) go to 50
        end if
        if (i < me%nac) go to 40

80      if (me%nac <= 0) return
        nac1 = me%nac + 1
        ! DETERMINE IF CONSTRAINTS ARE VIOLATED.
        nvc = 0
        do i = 1, me%nac
            nci = ic(i)
            ncj = 1
            if (nci <= me%ncon) ncj = isc(nci)
            ctc = ctam
            if (ncj > 0) ctc = ctbm
            if (g(nci) > ctc) nvc = nvc + 1
        end do
        ! ------------------------------------------------------------------
        ! NORMALIZE GRADIENT OF OBJECTIVE FUNCTION AND STORE IN NAC+1
        ! COLUMN OF A
        ! ------------------------------------------------------------------
        a1 = 0.
        do i = 1, me%ndv
            a1 = a1 + df(i)**2
        end do
        if (a1 < 1.0e-20_wp) a1 = 1.0e-20_wp
        a1 = sqrt(a1)
        a1 = 1.0_wp/a1
        do i = 1, me%ndv
            a(i, nac1) = a1*df(i)
        end do
        ! BUILD C VECTOR.
        if (nvc <= 0) then
            ! ------------------------------------------------------------------
            !             BUILD C FOR CLASSICAL METHOD
            ! ------------------------------------------------------------------
            ndb = nac1
            a(ndv1, ndb) = 1.0_wp
            do i = 1, ndb
                c(i) = -a(ndv1, i)
            end do
        else
            ! ------------------------------------------------------------------
            !               BUILD C FOR MODIFIED METHOD
            ! ------------------------------------------------------------------
            ndb = me%nac
            a(ndv1, nac1) = -phi
            ! ------------------------------------------------------------------
            !       SCALE THETA'S SO THAT MAXIMUM THETA IS UNITY
            ! ------------------------------------------------------------------
            if (thmax > 0.00001_wp) thmax = 1./thmax
            do i = 1, ndb
                a(ndv1, i) = a(ndv1, i)*thmax
            end do
            do i = 1, ndb
                c(i) = 0.0_wp
                do j = 1, ndv1
                    c(i) = c(i) + a(j, i)*a(j, nac1)
                end do
            end do
        end if
        ! ------------------------------------------------------------------
        !                  BUILD B MATRIX
        ! ------------------------------------------------------------------
        do i = 1, ndb
            do j = 1, ndb
                b(i, j) = 0.0_wp
                do k = 1, ndv1
                    b(i, j) = b(i, j) - a(k, i)*a(k, j)
                end do
            end do
        end do
        ! ------------------------------------------------------------------
        !                SOLVE SPECIAL L. P. PROBLEM
        ! ------------------------------------------------------------------
        call cnmn08(ndb, ner, c, ms1, b, n3, n4)
        if (me%iprint > 1 .and. ner > 0) write (6, 5200)
        ! CALCULATE RESULTING DIRECTION VECTOR, S.
        slope = 0.0_wp
        ! ------------------------------------------------------------------
        !              USABLE-FEASIBLE DIRECTION
        ! ------------------------------------------------------------------
        do i = 1, me%ndv
            s1 = 0.0_wp
            if (nvc > 0) s1 = -a(i, nac1)
            do j = 1, ndb
                s1 = s1 - a(i, j)*c(j)
            end do
            slope = slope + s1*df(i)
            s(i) = s1
        end do
        s(ndv1) = 1.0_wp
        if (nvc > 0) s(ndv1) = -a(ndv1, nac1)
        do j = 1, ndb
            s(ndv1) = s(ndv1) - a(ndv1, j)*c(j)
        end do
        ! ------------------------------------------------------------------
        ! CHECK TO INSURE THE S-VECTOR IS FEASIBLE.
        ! PROGRAM MOD-FEB, 1981, GV.
        ! ------------------------------------------------------------------
        do j = 1, me%nac
            ! S DOT DEL(G).
            sg = dot_product(s(1:me%ndv), a(1:me%ndv, j))
!           IF(SG>0.) GO TO 176
!
!  THIS CHANGE MADE ON 4/8/81 FOR G. VANDERPLAATS

            if (sg > 1.0e-04_wp) go to 240
            ! FEASIBLE FOR THIS CONSTRAINT.  CONTINUE.
        end do
        go to 250

        ! S-VECTOR IS NOT FEASIBLE DUE TO SOME NUMERICAL PROBLEM.
240     if (me%iprint >= 2) write (6, 5100)
        s(ndv1) = 0.0_wp
        nvc = 0
        return

        ! ------------------------------------------------------------------
        !              NORMALIZE S TO MAX ABS OF UNITY
        ! ------------------------------------------------------------------
250     s1 = 0.0_wp
        do i = 1, me%ndv
            a1 = abs(s(i))
            if (a1 > s1) s1 = a1
        end do
!       IF (S1<1.0E-10) RETURN
!
!  E-10 CHANGED TO E-04 ON 1/12/81

        if (s1 < 1.0e-04_wp) return
        s1 = 1.0_wp/s1
        do i = 1, me%ndv
            s(i) = s1*s(i)
        end do
        slope = s1*slope
        s(ndv1) = s1*s(ndv1)
        return

! ------------------------------------------------------------------
!                       FORMATS
! ------------------------------------------------------------------

5000    format(t6, '** CONSTRAINT', i5, ' HAS ZERO GRADIENT'/t6, &
               'DELETED FROM ACTIVE SET')
5100    format(t6, '** CALCULATED S-VECTOR IS NOT FEASIBLE'/t6, &
               'BETA IS SET TO ZERO')
5200    format(//t6, '* * DIRECTION FINDING PROCESS DID NOT CONVERGE'/t6, &
                '* * S-VECTOR MAY NOT BE VALID')

    end subroutine cnmn05

    subroutine cnmn06(me, x, vlb, vub, g, scal, df, s, g1, g2, ctam, ctbm, slope, alp, a2, a3, a4, &
                      f1, f2, f3, cv1, cv2, cv3, cv4, alpca, alpfes, alpln, alpmin, alpnc, &
                      alpsav, alpsid, alptot, isc, ncal, nvc, icount, igood1, &
                      igood2, igood3, igood4, ibest, iii, nlnc, jgoto)

        !!  Routine to solve one-dimensional search problem for constrained
        !!  function minimization.
        !!
        !!  BY G. N. VANDERPLAATS, AUG., 1974.
        !!
        !!  * OBJ = INITIAL AND FINAL FUNCTION VALUE.
        !!  * ALP = MOVE PARAMETER.
        !!  * SLOPE = INITIAL SLOPE.
        !!  * ALPSID = MOVE TO SIDE CONSTRAINT.
        !!  * ALPFES = MOVE TO FEASIBLE REGION.
        !!  * ALPNC = MOVE TO NEW NON-LINEAR CONSTRAINT.
        !!  * ALPLN = MOVE TO LINEAR CONSTRAINT.
        !!  * ALPCA = MOVE TO RE-ENCOUNTER CURRENTLY ACTIVE CONSTRAINT.
        !!  * ALPMIN = MOVE TO MINIMIZE FUNCTION.
        !!  * ALPTOT = TOTAL MOVE PARAMETER.

        class(conmin_class), intent(inout) :: me
        real(wp), intent(inout)   :: x(:), vlb(:), vub(:), g(:), scal(:), df(:), &
                                     s(:), g1(:), g2(:), ctam, ctbm, slope, alp, &
                                     a2, a3, a4, f1, f2, f3, cv1, cv2, cv3, cv4, &
                                     alpca, alpfes, alpln, alpmin, alpnc, alpsav, &
                                     alpsid, alptot
        integer, intent(inout)     :: isc(:), ncal(2), nvc, icount, igood1, igood2, &
                                      igood3, igood4, ibest, iii, nlnc, jgoto

        real(wp)  :: alpa, alpb, c1, c2, c3, cc, f4, gi, si, xi, xi1, xi2
        integer   :: i, ii, jbest, ksid, nvc1
        real(wp), parameter :: zro = 0.0_wp

        if (jgoto /= 0) then
            select case (jgoto)
            case (1); go to 70
            case (2); go to 140
            case (3); go to 230
            end select
        end if
        if (me%iprint >= 5) write (6, 5100)
        alpsav = alp
        icount = 0
        alptot = 0.0_wp
        ! TOLERANCES.
        ctam = abs(me%ctmin)
        ctbm = abs(me%ctlmin)
        ! PROPOSED MOVE.
        ! ------------------------------------------------------------------
        ! *****  BEGIN SEARCH OR IMPOSE SIDE CONSTRAINT MODIFICATION  *****
        ! ------------------------------------------------------------------
10      a2 = alpsav
        icount = icount + 1
        alpsid = 1.0e+20_wp
        ! INITIAL ALPHA AND OBJ.
        alp = 0.0_wp
        f1 = me%obj
        ksid = 0
        if (me%nside /= 0) then
            ! ------------------------------------------------------------------
            ! FIND MOVE TO SIDE CONSTRAINT AND INSURE AGAINST VIOLATION OF
            ! SIDE CONSTRAINTS
            ! ------------------------------------------------------------------
            do i = 1, me%ndv
                si = s(i)
                if (abs(si) <= 1.0e-20_wp) then
                    ! ITH COMPONENT OF S IS SMALL.  SET TO ZERO.
                    s(i) = 0.0_wp
                    slope = slope - si*df(i)
                else
                    xi = x(i)
                    si = 1.0_wp/si
                    if (si <= 0.0_wp) then
                        ! LOWER BOUND.
                        xi2 = vlb(i)
                        xi1 = abs(xi2)
                        if (xi1 < 1.0_wp) xi1 = 1.0_wp
                        ! CONSTRAINT VALUE.
                        gi = (xi2 - xi)/xi1
                        if (gi > -1.0e-6_wp) go to 20
                        ! PROPOSED MOVE TO LOWER BOUND.
                        alpa = (xi2 - xi)*si
                        if (alpa < alpsid) alpsid = alpa
                        cycle
                    end if
                    ! UPPER BOUND.
                    xi2 = vub(i)
                    xi1 = abs(xi2)
                    if (xi1 < 1.0_wp) xi1 = 1.0_wp
                    ! CONSTRAINT VALUE.
                    gi = (xi - xi2)/xi1
                    if (gi <= -1.0e-6_wp) then
                        ! PROPOSED MOVE TO UPPER BOUND.
                        alpa = (xi2 - xi)*si
                        if (alpa < alpsid) alpsid = alpa
                        cycle
                    end if

                    ! MOVE WILL VIOLATE SIDE CONSTRAINT.  SET S(I)=0.
20                  slope = slope - s(i)*df(i)
                    s(i) = 0.0_wp
                    ksid = ksid + 1
                end if
            end do
            ! ALPSID IS UPPER BOUND ON ALPHA.
            if (a2 > alpsid) a2 = alpsid
        end if
        ! ------------------------------------------------------------------
        !           CHECK ILL-CONDITIONING
        ! ------------------------------------------------------------------
        if (ksid == me%ndv .or. icount > 10) go to 340
        if (nvc == 0 .and. slope > 0.0_wp) go to 340
        alpfes = -1.0_wp
        alpmin = -1.0_wp
        alpln = 1.1_wp*alpsid
        alpnc = alpsid
        alpca = alpsid
        if (me%ncon /= 0) then
            ! STORE CONSTRAINT VALUES IN G1.
            do i = 1, me%ncon
                g1(i) = g(i)
            end do
        end if
        !  ------------------------------------------------------------------
        !               MOVE A DISTANCE A2*S
        !  ------------------------------------------------------------------
        alptot = alptot + a2
        do i = 1, me%ndv
            x(i) = x(i) + a2*s(i)
        end do
        if (me%iprint >= 5) then
            write (6, 5200) a2
            if (me%nscal /= 0) then
                do i = 1, me%ndv
                    g(i) = scal(i)*x(i)
                end do
                write (6, 5300) (g(i), i=1, me%ndv)
            else
                write (6, 5300) (x(i), i=1, me%ndv)
            end if
        end if
        ! ------------------------------------------------------------------
        !               UPDATE FUNCTION AND CONSTRAINT VALUES
        ! ------------------------------------------------------------------
        ncal(1) = ncal(1) + 1
        jgoto = 1
        return

70      f2 = me%obj
        if (me%iprint >= 5) write (6, 5400) f2
        if (me%iprint >= 5 .and. me%ncon /= 0) then
            write (6, 5500)
            write (6, 5300) (g(i), i=1, me%ncon)
        end if
        ! ------------------------------------------------------------------
        !           IDENTIFY ACCAPTABILITY OF DESIGNS F1 AND F2
        ! ------------------------------------------------------------------
        ! IGOOD = 0 IS ACCAPTABLE.
        ! CV = MAXIMUM CONSTRAINT VIOLATION.
        igood1 = 0
        igood2 = 0
        cv1 = 0.0_wp
        cv2 = 0.0_wp
        nvc1 = 0
        if (me%ncon /= 0) then
            do i = 1, me%ncon
                cc = ctam
                if (isc(i) > 0) cc = ctbm
                c1 = g1(i) - cc
                c2 = g(i) - cc
                if (c2 > 0.0_wp) nvc1 = nvc1 + 1
                if (c1 > cv1) cv1 = c1
                if (c2 > cv2) cv2 = c2
            end do
            if (cv1 > 0.0_wp) igood1 = 1
            if (cv2 > 0.0_wp) igood2 = 1
        end if
        alp = a2
        me%obj = f2
        ! ------------------------------------------------------------------
        ! IF F2 VIOLATES FEWER CONSTRAINTS THAN F1 BUT STILL HAS CONSTRAINT
        ! VIOLATIONS RETURN
        ! ------------------------------------------------------------------
        if (nvc1 < nvc .and. nvc1 > 0) go to 340
        ! ------------------------------------------------------------------
        !         IDENTIFY BEST OF DESIGNS F1 ANF F2
        ! ------------------------------------------------------------------
        ! IBEST CORRESPONDS TO MINIMUM VALUE DESIGN.
        ! IF CONSTRAINTS ARE VIOLATED, IBEST CORRESPONDS TO MINIMUM
        ! CONSTRAINT VIOLATION.
        if (igood1 /= 0 .or. igood2 /= 0) then
            ! VIOLATED CONSTRAINTS.  PICK MINIMUM VIOLATION.
            ibest = 1
            if (cv1 >= cv2) ibest = 2
        else
            ! NO CONSTRAINT VIOLATION.  PICK MINIMUM F.
            ibest = 1
            if (f2 <= f1) ibest = 2
        end if
        ii = 1
        ! ------------------------------------------------------------------
        ! IF CV2 IS GREATER THAN CV1, SET MOVE LIMITS TO A2.
        ! PROGRAM MOD-FEB, 1981, GV.
        ! ------------------------------------------------------------------
        if (cv2 > cv1) then
            alpln = a2
            alpnc = a2
            alpca = a2
        end if
        if (me%ncon /= 0) then
            ! ------------------------------------------------------------------
            ! *****                 2 - POINT INTERPOLATION                *****
            ! ------------------------------------------------------------------
            iii = 0
90          iii = iii + 1
            c1 = g1(iii)
            c2 = g(iii)
            if (isc(iii) /= 0) then
                ! ------------------------------------------------------------------
                !                    LINEAR CONSTRAINT
                ! ------------------------------------------------------------------
                if (c1 >= 1.0e-5_wp .and. c1 <= ctbm) go to 100
                call cnmn07(ii, alp, zro, zro, c1, a2, c2, zro, zro)
                if (alp <= 0.0_wp) go to 100
                if (c1 > ctbm .and. alp > alpfes) alpfes = alp
                if (c1 < me%ctl .and. alp < alpln) alpln = alp
            else
                ! ------------------------------------------------------------------
                !                 NON-LINEAR CONSTRAINT
                ! ------------------------------------------------------------------
                if (c1 < 1.0e-5_wp .or. c1 > ctam) then
                    call cnmn07(ii, alp, zro, zro, c1, a2, c2, zro, zro)
                    if (alp > 0.0_wp) then
                        if (c1 > ctam .and. alp > alpfes) alpfes = alp
                        if (c1 < me%ct .and. alp < alpnc) alpnc = alp
                    end if
                end if
            end if

100         if (iii < me%ncon) go to 90
        end if
        if (me%linobj <= 0 .and. slope < 0.0_wp) then
            ! CALCULATE ALPHA TO MINIMIZE FUNCTION.
            call cnmn04(ii, alpmin, zro, zro, f1, slope, a2, f2, zro, zro, zro, zro)
        end if
        ! ------------------------------------------------------------------
        !                     PROPOSED MOVE
        ! ------------------------------------------------------------------
        ! MOVE AT LEAST FAR ENOUGH TO OVERCOME CONSTRAINT VIOLATIONS.
        a3 = alpfes
        ! MOVE TO MINIMIZE FUNCTION.
        if (alpmin > a3) a3 = alpmin
        ! IF A3<=0, SET A3 = ALPSID.
        if (a3 <= 0.0_wp) a3 = alpsid
        ! LIMIT MOVE TO NEW CONSTRAINT ENCOUNTER.
        if (a3 > alpnc) a3 = alpnc
        if (a3 > alpln) a3 = alpln
        ! MAKE A3 NON-ZERO.
        if (a3 <= 1.0e-20_wp) a3 = 1.0e-20_wp
        ! IF A3=A2=ALPSID AND F2 IS BEST, GO INVOKE SIDE CONSTRAINT
        ! MODIFICATION.
        alpb = 1.0_wp-a2/a3
        alpa = 1.0_wp-alpsid/a3
        jbest = 0
        if (abs(alpb) < 1.0e-10_wp .and. abs(alpa) < 1.0e-10_wp) jbest = 1
        if (jbest == 1 .and. ibest == 2) go to 10
        ! SIDE CONSTRAINT CHECK NOT SATISFIED.
        if (me%ncon /= 0) then
            ! STORE CONSTRAINT VALUES IN G2.
            do i = 1, me%ncon
                g2(i) = g(i)
            end do
        end if
        ! IF A3=A2, SET A3=.9*A2.
        if (abs(alpb) < 1.0e-10_wp) a3 = 0.9_wp*a2
        ! MOVE AT LEAST .01*A2.
        if (a3 < (0.01_wp*a2)) a3 = 0.01_wp*a2
        ! LIMIT MOVE TO 5.*A2.
        if (a3 > (5.0_wp*a2)) a3 = 5.0_wp*a2
        ! LIMIT MOVE TO ALPSID.
        if (a3 > alpsid) a3 = alpsid
        ! MOVE A DISTANCE A3*S.
        alp = a3 - a2
        alptot = alptot + alp
        do i = 1, me%ndv
            x(i) = x(i) + alp*s(i)
        end do
        if (me%iprint >= 5) then
            write (6, 5600)
            write (6, 5200) a3
            if (me%nscal /= 0) then
                g(1:me%ndv) = scal(1:me%ndv)*x(1:me%ndv)
                write (6, 5300) g(1:me%ndv)
            else
                write (6, 5300) x(1:me%ndv)
            end if
        end if
        ! ------------------------------------------------------------------
        !          UPDATE FUNCTION AND CONSTRAINT VALUES
        ! ------------------------------------------------------------------
        ncal(1) = ncal(1) + 1
        jgoto = 2
        return

140     f3 = me%obj
        if (me%iprint >= 5) write (6, 5400) f3
        if (me%iprint >= 5 .and. me%ncon /= 0) then
            write (6, 5500)
            write (6, 5300) (g(i), i=1, me%ncon)
        end if
        ! ------------------------------------------------------------------
        !   CALCULATE MAXIMUM CONSTRAINT VIOLATION AND PICK BEST DESIGN
        ! ------------------------------------------------------------------
        cv3 = 0.0_wp
        igood3 = 0
        nvc1 = 0
        if (me%ncon /= 0) then
            do i = 1, me%ncon
                cc = ctam
                if (isc(i) > 0) cc = ctbm
                c1 = g(i) - cc
                if (c1 > cv3) cv3 = c1
                if (c1 > 0.0_wp) nvc1 = nvc1 + 1
            end do
            if (cv3 > 0.0_wp) igood3 = 1
        end if
        ! DETERMINE BEST DESIGN.
        if (ibest /= 2) then
            ! CHOOSE BETWEEN F1 AND F3.
            if (igood1 /= 0 .or. igood3 /= 0) then
                if (cv1 >= cv3) ibest = 3
                go to 160
            end if
            if (f3 <= f1) ibest = 3
        else
            ! CHOOSE BETWEEN F2 AND F3.
            if (igood2 /= 0 .or. igood3 /= 0) then
                if (cv2 >= cv3) ibest = 3
            else
                if (f3 <= f2) ibest = 3
            end if
        end if

160     alp = a3
        me%obj = f3
        ! IF F3 VIOLATES FEWER CONSTRAINTS THAN F1 RETURN.
        if (nvc1 < nvc) go to 340
        ! IF OBJECTIVE AND ALL CONSTRAINTS ARE LINEAR, RETURN.
        if (me%linobj /= 0 .and. nlnc == me%ncon) go to 340
        ! IF A3 = ALPLN AND F3 IS BOTH GOOD AND BEST RETURN.
        alpb = 1.0_wp-alpln/a3
        if (abs(alpb) < 1.0e-20_wp .and. ibest == 3 .and. igood3 == 0) go to 340
        ! IF A3 = ALPSID AND F3 IS BEST, GO INVOKE SIDE CONSTRAINT MODIFICATION.
        alpa = 1.0_wp-alpsid/a3
        if (abs(alpa) < 1.0e-20_wp .and. ibest == 3) go to 10
        ! ------------------------------------------------------------------
        ! **********            3 - POINT INTERPOLATION            *********
        ! ------------------------------------------------------------------
        alpnc = alpsid
        alpca = alpsid
        alpfes = -1.0_wp
        alpmin = -1.0_wp
        ! ------------------------------------------------------------------
        ! IF A3 IS GREATER THAN A2 AND CV3 IS GREATER THAN CV2, SET
        ! MOVE LIMITS TO A3.  PROGRAM MOD-FEB, 1981, GV.
        ! ------------------------------------------------------------------
        if (a3 > a2 .and. cv3 > cv2) then
            alpln = a3
            alpnc = a3
            alpca = a3
        end if
        if (me%ncon /= 0) then
            iii = 0
170         iii = iii + 1
            c1 = g1(iii)
            c2 = g2(iii)
            c3 = g(iii)
            if (isc(iii) /= 0) then
                ! ------------------------------------------------------------------
                ! LINEAR CONSTRAINT.  FIND ALPFES ONLY.  ALPLN SAME AS BEFORE.
                ! ------------------------------------------------------------------
                if (c1 <= ctbm) go to 190
                ii = 1
                call cnmn07(ii, alp, zro, zro, c1, a3, c3, zro, zro)
                if (alp > alpfes) alpfes = alp
            else
                ! ------------------------------------------------------------------
                !                 NON-LINEAR CONSTRAINT
                ! ------------------------------------------------------------------
                ii = 2
                call cnmn07(ii, alp, zro, zro, c1, a2, c2, a3, c3)
                if (alp > zro) then
                    if (c1 < me%ct .or. c1 > 0.0_wp) then
                        if (c1 > ctam .or. c1 < 0.0_wp) go to 180
                    end if
                    ! ALP IS MINIMUM MOVE.  UPDATE FOR NEXT CONSTRAINT ENCOUNTER.
                    alpa = alp
                    call cnmn07(ii, alp, alpa, zro, c1, a2, c2, a3, c3)
                    if (alp < alpca .and. alp >= alpa) alpca = alp
                    go to 190

180                 if (alp > alpfes .and. c1 > ctam) alpfes = alp
                    if (alp < alpnc .and. c1 < 0.0_wp) alpnc = alp
                end if
            end if

190         if (iii < me%ncon) go to 170
        end if
        if (me%linobj <= 0 .and. slope <= 0.) then
            ! ------------------------------------------------------------------
            !          CALCULATE ALPHA TO MINIMIZE FUNCTION
            ! ------------------------------------------------------------------
            ii = 3
            if (a2 > a3 .and. (igood2 == 0 .and. ibest == 2)) ii = 2
            call cnmn04(ii, alpmin, zro, zro, f1, slope, a2, f2, a3, f3, zro, zro)
        end if
        ! ------------------------------------------------------------------
        !                   PROPOSED MOVE
        ! ------------------------------------------------------------------
        ! MOVE AT LEAST ENOUGH TO OVERCOME CONSTRAINT VIOLATIONS.
        a4 = alpfes
        ! MOVE TO MINIMIZE FUNCTION.
        if (alpmin > a4) a4 = alpmin
        ! IF A4<=0, SET A4 = ALPSID.
        if (a4 <= 0.0_wp) a4 = alpsid
        ! LIMIT MOVE TO NEW CONSTRAINT ENCOUNTER.
        if (a4 > alpln) a4 = alpln
        if (a4 > alpnc) a4 = alpnc
        ! LIMIT MOVE TO RE-ENCOUNTER CURRENTLY ACTIVE CONSTRAINT.
        if (a4 > alpca) a4 = alpca
        ! LIMIT A4 TO 5.*A3.
        if (a4 > (5.0_wp*a3)) a4 = 5.0_wp*a3
        ! UPDATE DESIGN.
        if (ibest == 3 .and. me%ncon /= 0) then
            ! STORE CONSTRAINT VALUES IN G2.  F3 IS BEST.  F2 IS NOT.
            do i = 1, me%ncon
                g2(i) = g(i)
            end do
        end if
        ! IF A4=A3 AND IGOOD1=0 AND IGOOD3=1, SET A4=.9*A3.
        alp = a4 - a3
        if (igood1 == 0 .and. igood3 == 1 .and. abs(alp) < 1.0e-20_wp) a4 = 0.9_wp*a3
        ! ------------------------------------------------------------------
        !               MOVE A DISTANCE A4*S
        ! ------------------------------------------------------------------
        alp = a4 - a3
        alptot = alptot + alp
        do i = 1, me%ndv
            x(i) = x(i) + alp*s(i)
        end do
        if (me%iprint >= 5) then
            write (6, 5000)
            write (6, 5200) a4
            if (me%nscal /= 0) then
                g(1:me%ndv) = scal(1:me%ndv)*x(1:me%ndv)
                write (6, 5300) g(1:me%ndv)
            else
                write (6, 5300) x(1:me%ndv)
            end if
        end if
        ! ------------------------------------------------------------------
        !          UPDATE FUNCTION AND CONSTRAINT VALUES
        ! ------------------------------------------------------------------
        ncal(1) = ncal(1) + 1
        jgoto = 3
        return

230     f4 = me%obj
        if (me%iprint >= 5) write (6, 5400) f4
        if (me%iprint >= 5 .and. me%ncon /= 0) then
            write (6, 5500)
            write (6, 5300) g(1:me%ncon)
        end if
        ! DETERMINE ACCAPTABILITY OF F4.
        igood4 = 0
        cv4 = 0.0_wp
        if (me%ncon /= 0) then
            do i = 1, me%ncon
                cc = ctam
                if (isc(i) > 0) cc = ctbm
                c1 = g(i) - cc
                if (c1 > cv4) cv4 = c1
            end do
            if (cv4 > 0.0_wp) igood4 = 1
        end if
        alp = a4
        me%obj = f4
        ! ------------------------------------------------------------------
        !                 DETERMINE BEST DESIGN
        ! ------------------------------------------------------------------
        select case (ibest)
        case (1)
            ! CHOOSE BETWEEN F1 AND F4.
            if (igood1 /= 0 .or. igood4 /= 0) then
                if (cv1 > cv4) go to 340
            else
                if (f4 <= f1) go to 340
            end if
            ! F1 IS BEST.
            alptot = alptot - a4
            me%obj = f1
            x(1:me%ndv) = x(1:me%ndv) - a4*s(1:me%ndv)
            if (me%ncon == 0) go to 340
            g(1:me%ncon) = g1(1:me%ncon)

        case (2)
            ! CHOOSE BETWEEN F2 AND F4.
            if (igood2 /= 0 .or. igood4 /= 0) then
                if (cv2 > cv4) go to 340
            else
                if (f4 <= f2) go to 340
            end if
            ! F2 IS BEST.
            me%obj = f2
            a2 = a4 - a2
            alptot = alptot - a2
            x(1:me%ndv) = x(1:me%ndv) - a2*s(1:me%ndv)
            if (me%ncon == 0) go to 340
            g(1:me%ncon) = g2(1:me%ncon)

        case (3)
            ! CHOOSE BETWEEN F3 AND F4.
            if (igood3 /= 0 .or. igood4 /= 0) then
                if (cv3 > cv4) go to 340
            else
                if (f4 <= f3) go to 340
            end if
            ! F3 IS BEST.
            me%obj = f3
            a3 = a4 - a3
            alptot = alptot - a3
            x(1:me%ndv) = x(1:me%ndv) - a3*s(1:me%ndv)
            if (me%ncon /= 0) then
                g(1:me%ncon) = g2(1:me%ncon)
            end if

        end select

340     alp = alptot
        if (me%iprint >= 5) write (6, 5700)
        jgoto = 0
        return

        ! ------------------------------------------------------------------
        !                              FORMATS
        ! ------------------------------------------------------------------

5000    format(/t6, 'THREE-POINT INTERPOLATION')
5100    format(///'* * * CONSTRAINED ONE-DIMENSIONAL SEARCH INFORMATION * * *')
5200    format(//t6, 'PROPOSED DESIGN'/t6, 'ALPHA =', e12.5/t6, 'X-VECTOR')
5300    format(' ', 8e12.4)
5400    format(/t6, 'OBJ =', e13.5)
5500    format(/t6, 'CONSTRAINT VALUES')
5600    format(/t6, 'TWO-POINT INTERPOLATION')
5700    format(/t6, '* * * END OF ONE-DIMENSIONAL SEARCH')

    end subroutine cnmn06

    subroutine cnmn07(ii, xbar, eps, x1, y1, x2, y2, x3, y3)

        !!  Routine to find first `xbar>=eps `corresponding to a real zero
        !!  of a one-dimensional function by polynomial interpolation.
        !!
        !!  BY G. N. VANDERPLAATS, APRIL, 1972.
        !!
        !!  If required zero on `y` does not exits, or the function is
        !!  ill-conditioned, `xbar = eps-1.0` will be returned as an error indicator.
        !!  if desired interpolation is ill-conditioned, a lower order
        !!  interpolation, consistant with input data, will be attempted and
        !!  ii will be changed accordingly.

        integer, intent(inout)     :: ii  !! CALCULATION CONTROL:
                                          !!
                                          !! 1.  2-POINT LINEAR INTERPOLATION, GIVEN X1, Y1, X2 AND Y2.
                                          !! 2.  3-POINT QUADRATIC INTERPOLATION, GIVEN X1, Y1, X2, Y2, X3 AND Y3.
        real(wp), intent(out)     :: xbar
        real(wp), intent(in)      :: eps !! may be negative
        real(wp), intent(in)      :: x1, y1, x2, y2, x3, y3

        real(wp)  :: aa, bac, bb, cc, dy, qq, x21, x31, x32, xb2, xbar1, yy
        integer   :: jj

        xbar1 = eps - 1.0_wp
        xbar = xbar1
        jj = 0
        x21 = x2 - x1
        if (abs(x21) < 1.0e-20_wp) return
        if (ii == 2) then
            ! ------------------------------------------------------------------
            !                 ii=2: 3-point quadratic interpolation
            ! ------------------------------------------------------------------
            jj = 1
            x31 = x3 - x1
            x32 = x3 - x2
            qq = x21*x31*x32
            if (abs(qq) < 1.0e-20_wp) return
            aa = (y1*x32 - y2*x31 + y3*x21)/qq
            if (abs(aa) >= 1.0e-20_wp) then
                bb = (y2 - y1)/x21 - aa*(x1 + x2)
                cc = y1 - x1*(aa*x1 + bb)
                bac = bb*bb - 4.0_wp*aa*cc
                if (bac >= 0.0_wp) then
                    bac = sqrt(bac)
                    aa = 0.5_wp/aa
                    xbar = aa*(bac - bb)
                    xb2 = -aa*(bac + bb)
                    if (xbar < eps) xbar = xb2
                    if (xb2 < xbar .and. xb2 > eps) xbar = xb2
                    if (xbar < eps) xbar = xbar1
                    return
                end if
            end if
        end if

        ! ------------------------------------------------------------------
        !                  ii=1: 2-point linear interpolation
        ! ------------------------------------------------------------------
        ii = 1
        yy = y1*y2
        if (jj /= 0 .and. yy >= 0.0_wp) then
            ! interpolate between x2 and x3.
            dy = y3 - y2
            if (abs(dy) >= 1.0e-20_wp) then
                xbar = x2 + y2*(x2 - x3)/dy
                if (xbar < eps) xbar = xbar1
                return
            end if
        end if
        dy = y2 - y1
        ! interpolate between x1 and x2.
        if (abs(dy) < 1.0e-20_wp) return
        xbar = x1 + y1*(x1 - x2)/dy
        if (xbar < eps) xbar = xbar1

    end subroutine cnmn07

    subroutine cnmn08(ndb, ner, c, ms1, b, n3, n4)

        !!  Routine to solve special linear problem for imposing s-transpose
        !!  times s<=1 bounds in the modified method of feasible directions.
        !!
        !!  FORM OF L. P. IS BX=C WHERE 1ST NDB COMPONENTS OF X CONTAIN VECTOR
        !!  U AND LAST NDB COMPONENTS CONTAIN VECTOR V.
        !!  CONSTRAINTS ARE U>=0, V>=0, AND U-TRANSPOSE TIMES V = 0.
        !!
        !!  BY G. N. VANDERPLAATS, APRIL, 1972.
        !!
        !!### Reference
        !!  * "[Structural optimization by methods of feasible directions](https://www.sciencedirect.com/science/article/abs/pii/0045794973900552)",
        !!    G. N. Vanderplaats and F. Moses, Journal of computers
        !!    and structures, vol 3, pp 739-755, 1973.

        !  ------------------------------------------------------------------
        !  CHOOSE INITIAL BASIC VARIABLES AS V, AND INITIALIZE VECTOR MS1
        !  ------------------------------------------------------------------

        integer, intent(in)        :: ndb, n3, n4
        integer, intent(out)       :: ner !! NER = ERROR FLAG.  IF NER/=0 ON RETURN, PROCESS HAS NOT
                                          !! CONVERGED IN 5*NDB ITERATIONS.
        integer, intent(out)       :: ms1(:) !! VECTOR MS1 IDENTIFIES THE SET OF BASIC VARIABLES.
        real(wp), intent(inout)   :: c(n4), b(n3, n3)

        integer   :: i, ichk, iter1, j, jj, kk, m2, nmax
        real(wp)  :: bb, bb1, bi, c1, cb, cbmax, cbmin, eps

        ner = 1
        m2 = 2*ndb
        ! CALCULATE CBMIN AND EPS AND INITIALIZE MS1.
        eps = -1.0e+10_wp
        cbmin = 0.0_wp
        do i = 1, ndb
            bi = b(i, i)
            cbmax = 0.
            if (bi < -1.0e-6_wp) cbmax = c(i)/bi
            if (bi > eps) eps = bi
            if (cbmax > cbmin) cbmin = cbmax
            ms1(i) = 0
        end do
        eps = 0.0001_wp*eps
!       IF (EPS<-1.0E-10) EPS=-1.0E-10
!
!  E-10 CHANGED TO E-03 ON 1/12/81
!
        if (eps < -1.0e-03_wp) eps = -1.0e-03_wp
        if (eps > -0.0001_wp) eps = -0.0001_wp
        cbmin = cbmin*1.0e-6_wp
!       IF (CBMIN<1.0e-10_wp) CBMIN=1.0e-10_wp
!
!  E-10 CHANGED TO E-05 ON 1/12/81
!
        if (cbmin < 1.0e-05_wp) cbmin = 1.0e-05_wp
        iter1 = 0
        nmax = 5*ndb
        main: do
            ! ------------------------------------------------------------------
            ! **********             BEGIN NEW ITERATION              **********
            ! ------------------------------------------------------------------
            iter1 = iter1 + 1
            if (iter1 > nmax) return
            ! FIND MAX. C(I)/B(I,I) FOR I=1,NDB.
            cbmax = 0.9_wp*cbmin
            ichk = 0
            do i = 1, ndb
                c1 = c(i)
                bi = b(i, i)
                !     IF (BI>EPS .OR. C1>0.) GO TO 30
                if (bi <= eps .and. c1 <= -1.0e-05_wp) then
                    !  0. CHANGED TO -1.0E-05 ON 1/12/81
                    cb = c1/bi
                    if (cb > cbmax) then
                        ichk = i
                        cbmax = cb
                    end if
                end if
            end do
            if (cbmax >= cbmin) then
                if (ichk /= 0) then
                    ! UPDATE VECTOR MS1.
                    jj = ichk
                    if (ms1(jj) == 0) jj = ichk + ndb
                    kk = jj + ndb
                    if (kk > m2) kk = jj - ndb
                    ms1(kk) = ichk
                    ms1(jj) = 0
                    ! ------------------------------------------------------------------
                    !                 PIVOT OF B(ICHK,ICHK)
                    ! ------------------------------------------------------------------
                    bb = 1.0_wp/b(ichk, ichk)
                    do j = 1, ndb
                        b(ichk, j) = bb*b(ichk, j)
                    end do
                    c(ichk) = cbmax
                    b(ichk, ichk) = bb
                    ! ELIMINATE COEFICIENTS ON VARIABLE ENTERING BASIS AND STORE
                    ! COEFICIENTS ON VARIABLE LEAVING BASIS IN THEIR PLACE.
                    do i = 1, ndb
                        if (i /= ichk) then
                            bb1 = b(i, ichk)
                            b(i, ichk) = 0.0_wp
                            do j = 1, ndb
                                b(i, j) = b(i, j) - bb1*b(ichk, j)
                            end do
                            c(i) = c(i) - bb1*cbmax
                        end if
                    end do
                    cycle main
                end if
            end if
            exit main
        end do main
        ner = 0
        ! ------------------------------------------------------------------
        ! STORE ONLY COMPONENTS OF U-VECTOR IN 'C'.  USE B(I,1) FOR
        ! TEMPORARY STORAGE
        ! ------------------------------------------------------------------
        b(1:ndb, 1) = c(1:ndb)
        do i = 1, ndb
            c(i) = 0.0_wp
            j = ms1(i)
            if (j > 0) c(i) = b(j, 1)
            if (c(i) < 0.0_wp) c(i) = 0.0_wp
        end do

    end subroutine cnmn08

end module conmin_module

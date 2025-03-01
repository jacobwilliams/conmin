

subroutine cnmn00(x,Vlb,Vub,g,Scal,Df,a,s,g1,g2,b,c,Isc,Ic,Ms1,n1,n2,n3,n4,n5)
   implicit none

   double precision a , a1 , a2 , a3 , a4 , Abobj , Abobj1 , Alp , alp1 , alp11 , alp12 , Alpca , Alpfes , Alphax , Alpln ,        &
                  & Alpmin , Alpnc , Alpsav , Alpsid , Alptot
   double precision App , b , c , c1 , Ct , ct1 , Cta , Ctam , Ctbm , ctc , Ctl , Ctlmin , Ctmin , Cv1 , Cv2 , Cv3 , Cv4 , Dabfun ,&
                  & Dct , Dctl
   double precision Delfun , Df , Dftdf1 , Dm1 , Dm10 , Dm11 , Dm12 , Dm2 , Dm3 , Dm4 , Dm5 , Dm6 , Dm7 , Dm8 , Dm9 , Dx , Dx1 ,   &
                  & f1 , f2 , f3
   double precision f4 , Fdch , Fdchm , ff1 , Fff , Fi , g , g1 , g2 , gi , Obj , Obj1 , objb , objd , Phi , Rspace , s , Scal ,   &
                  & scj , si
   double precision sib , Slope , Theta , Vlb , Vub , x , x1 , x12 , Xi , xid , xx
   integer i , Ibest , Ic , Icndir , Icount , Idm1 , Idm2 , Idm3 , Igood1 , Igood2 , Igood3 , Igood4 , Igoto , ii , Iii , Info ,   &
         & Infog , Iobj , Iout , Iprint
   integer Isc , Ispace , Iter , Itmax , Itrm , j , Jdir , Jgoto , k , Kcount , Kobj , Kount , Linobj , m1 , m2 , m3 , mcn1 , Ms1 ,&
         & Mscal , n1
   integer n2 , n3 , n4 , n5 , Nac , Ncal , nci , Ncobj , Ncon , Ndv , ndv1 , ndv2 , Nfdg , Nfeas , Nfeasct , nic , Nlnc , nnac ,  &
         & Nscal , Nside
   integer Nvc

   common /cnmn1 / Delfun , Dabfun , Fdch , Fdchm , Ct , Ctmin , Ctl , Ctlmin , Alphax , Abobj1 , Theta , Obj , Ndv , Ncon ,       &
                 & Nside , Iprint , Nfdg , Nscal , Linobj , Itmax , Itrm , Icndir , Igoto , Nac , Info , Infog , Iter , Nfeasct
   common /output/ Iout
!
!    NFEASCT ADDED TO COMMON BLOCK BY KCYOUNG ON 4/14/92 TO ALLOW MORE
!    THAN 10 ITERATION ATTEMPTS.  NFEASCT BECOMES AN INPUT VALUE
!
   dimension x(n1) , Vlb(n1) , Vub(n1) , g(n2) , Scal(n1) , Df(n1) , a(n1,n3) , s(n1) , g1(n2) , g2(n2) , b(n3,n3) , c(n4) ,       &
           & Isc(n2) , Ic(n3) , Ms1(n5)
   common /consav/ Dm1 , Dm2 , Dm3 , Dm4 , Dm5 , Dm6 , Dm7 , Dm8 , Dm9 , Dm10 , Dm11 , Dm12 , Dct , Dctl , Phi , Abobj , Cta ,     &
                 & Ctam , Ctbm , Obj1 , Slope , Dx , Dx1 , Fi , Xi , Dftdf1 , Alp , Fff , a1 , a2 , a3 , a4 , f1 , f2 , f3 , f4 ,  &
                 & Cv1 , Cv2 , Cv3 , Cv4 , App , Alpca , Alpfes , Alpln , Alpmin , Alpnc , Alpsav , Alpsid , Alptot , Rspace ,     &
                 & Idm1 , Idm2 , Idm3 , Jdir , Iobj , Kobj , Kcount , Ncal(2) , Nfeas , Mscal , Ncobj , Nvc , Kount , Icount ,     &
                 & Igood1 , Igood2 , Igood3 , Igood4 , Ibest , Iii , Nlnc , Jgoto , Ispace(2)
   integer :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: do
      select case (spag_nextblock_1)
      case (1)
!     ROUTINE TO SOLVE CONSTRAINED OR UNCONSTRAINED FUNCTION
!     MINIMIZATION.
!     BY G. N. VANDERPLAATS                          APRIL, 1972.
!     * * * * * * * * * * *   JUNE, 1979 VERSION   * * * * * * * * * * *
!     NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF.
!     REFERENCE;  CONMIN - A FORTRAN PROGRAM FOR CONSTRAINED FUNCTION
!         MINIMIZATION:  USER'S MANUAL,  BY G. N. VANDERPLAATS,
!         NASA TM X-62,282, AUGUST, 1973.
!     STORAGE REQUIREMENTS:
!         PROGRAM - 7000 DECIMAL WORDS (CDC COMPUTER)
!         ARRAYS  - APPROX. 2*(NDV**2)+26*NDV+4*NCON,
!               WHERE N3 = NDV+2.
!     RE-SCALE VARIABLES IF REQUIRED.
         if ( Nscal/=0 .and. Igoto/=0 ) then
            do i = 1 , Ndv
               x(i) = c(i)
            enddo
         endif
!     CONSTANTS.
         ndv1 = Ndv + 1
         ndv2 = Ndv + 2
         if ( Igoto/=0 ) then
!     ------------------------------------------------------------------
!                     CHECK FOR UNBOUNDED SOLUTION
!     ------------------------------------------------------------------
!     STOP IF OBJ IS LESS THAN -1.0E+20
            if ( Obj<=-1.0e+20 ) then
               write (Iout,99002)
99002          format (///5x,'CONMIN HAS ACHIEVED A SOLUTION OF OBJ LESS THAN -1.0E+40'/5x,'SOLUTION APPEARS TO BE UNBOUNDED'/5x,  &
                      &'OPTIMIZATION IS TERMINATED')
               spag_nextblock_1 = 10
               cycle SPAG_DispatchLoop_1
            elseif ( Igoto==1 ) then
               Obj1 = Obj
               if ( Dabfun<=0. ) Dabfun = .001*abs(Obj)
               if ( Dabfun<1.0e-10 ) Dabfun = 1.0e-10
               if ( Iprint>0 ) then
!     ------------------------------------------------------------------
!                    PRINT INITIAL DESIGN INFORMATION
!     ------------------------------------------------------------------
                  if ( Iprint>1 ) then
                     if ( Nside==0 .and. Ncon==0 ) write (Iout,99033)
99033                format (////5x,'UNCONSTRAINED FUNCTION MINIMIZATION'//5x,'CONTROL PARAMETERS')
                     if ( Nside/=0 .or. Ncon>0 ) write (Iout,99027)
99027                format (////5x,'CONSTRAINED FUNCTION MINIMIZATION'//5x,'CONTROL PARAMETERS')
                     write (Iout,99028) Iprint , Ndv , Itmax , Ncon , Nside , Icndir , Nscal , Nfdg , Linobj , Itrm , n1 , n2 ,    &
                                      & n3 , n4 , n5
99028                format (/5x,'IPRINT  NDV    ITMAX    NCON    NSIDE  ICNDIR   NSCAL   NFDG'/8i8//5x,'LINOBJ  ITRM',5x,'N1',6x, &
                            &'N2',6x,'N3',6x,'N4',6x,'N5'/8i8)
                     write (Iout,99030) Ct , Ctmin , Ctl , Ctlmin , Theta , Phi , Delfun , Dabfun
99030                format (/9x,'CT',14x,'CTMIN',11x,'CTL',13x,'CTLMIN'/1x,4(2x,e14.5)//9x,'THETA',11x,'PHI',13x,'DELFUN',10x,    &
                            &'DABFUN'/1x,4(2x,e14.5))
                     write (Iout,99029) Fdch , Fdchm , Alphax , Abobj1
99029                format (/9x,'FDCH',12x,'FDCHM',11x,'ALPHAX',10x,'ABOBJ1'/1x,4(2x,e14.5))
                     if ( Nside/=0 ) then
                        write (Iout,99031)
99031                   format (/5x,'LOWER BOUNDS ON DECISION VARIABLES (VLB)')
                        do i = 1 , Ndv , 6
                           m1 = min0(Ndv,i+5)
                           write (Iout,99005) i , (Vlb(j),j=i,m1)
                        enddo
                        write (Iout,99032)
99032                   format (/5x,'UPPER BOUNDS ON DECISION VARIABLES (VUB)')
                        do i = 1 , Ndv , 6
                           m1 = min0(Ndv,i+5)
                           write (Iout,99005) i , (Vub(j),j=i,m1)
                        enddo
                     endif
                     if ( Nscal<0 ) then
                        write (Iout,99034)
99034                   format (/5x,'SCALING VECTOR (SCAL)')
                        write (Iout,99050) (Scal(i),i=1,Ndv)
                     endif
                     if ( Ncon/=0 ) then
                        if ( Nlnc==0 .or. Nlnc==Ncon ) then
                           if ( Nlnc==Ncon ) write (Iout,99008)
99008                      format (/5x,'ALL CONSTRAINTS ARE LINEAR')
                           if ( Nlnc==0 ) write (Iout,99009)
99009                      format (/5x,'ALL CONSTRAINTS ARE NON-LINEAR')
                        else
                           write (Iout,99006)
99006                      format (/5x,'LINEAR CONSTRAINT IDENTIFIERS (ISC)'/5x,'NON-ZERO INDICATES LINEAR CONSTRAINT')
                           do i = 1 , Ncon , 15
                              m1 = min0(Ncon,i+14)
                              write (Iout,99007) i , (Isc(j),j=i,m1)
99007                         format (3x,i5,')',2x,15i5)
                           enddo
                        endif
                     endif
                  endif
                  write (Iout,99048) Obj
99048             format (//5x,'INITIAL FUNCTION INFORMATION'//5x,'OBJ =',e15.6)
                  write (Iout,99049)
                  do i = 1 , Ndv
                     x1 = 1.
                     if ( Nscal/=0 ) x1 = Scal(i)
                     g1(i) = x(i)*x1
                  enddo
                  do i = 1 , Ndv , 6
                     m1 = min0(Ndv,i+5)
                     write (Iout,99005) i , (g1(j),j=i,m1)
                  enddo
                  if ( Ncon/=0 ) then
                     write (Iout,99051)
                     do i = 1 , Ncon , 6
                        m1 = min0(Ncon,i+5)
                        write (Iout,99005) i , (g(j),j=i,m1)
                     enddo
                  endif
               endif
               if ( Iprint>1 ) write (Iout,99040)
               spag_nextblock_1 = 2
               cycle SPAG_DispatchLoop_1
            elseif ( Igoto==2 ) then
               spag_nextblock_1 = 4
               cycle SPAG_DispatchLoop_1
            elseif ( Igoto==3 ) then
               spag_nextblock_1 = 3
               cycle SPAG_DispatchLoop_1
            elseif ( Igoto==4 ) then
               spag_nextblock_1 = 7
               cycle SPAG_DispatchLoop_1
            elseif ( Igoto==5 ) then
               spag_nextblock_1 = 8
               cycle SPAG_DispatchLoop_1
            endif
         endif
!     ------------------------------------------------------------------
!                      SAVE INPUT CONTROL PARAMETERS
!     ------------------------------------------------------------------
         if ( Iprint>0 ) write (Iout,99026)
99026    format ('1',////12x,27('* ')/12x,'*',51x,'*'/12x,'*',20x,'C O N M I N',20x,'*'/12x,'*',51x,'*'/12x,'*',15x,               &
                &' FORTRAN PROGRAM FOR ',15x,'*'/12x,'*',51x,'*'/12x,'*',9x,'CONSTRAINED FUNCTION MINIMIZATION',9x,'*'/12x,'*',51x,&
                &'*'/12x,27('* '))
         if ( Linobj==0 .or. (Ncon>0 .or. Nside>0) ) then
            Idm1 = Itrm
            Idm2 = Itmax
            Idm3 = Icndir
            Dm1 = Delfun
            Dm2 = Dabfun
            Dm3 = Ct
            Dm4 = Ctmin
            Dm5 = Ctl
            Dm6 = Ctlmin
            Dm7 = Theta
            Dm8 = Phi
            Dm9 = Fdch
            Dm10 = Fdchm
            Dm11 = Abobj1
            Dm12 = Alphax
!     ------------------------------------------------------------------
!                                DEFAULTS
!     ------------------------------------------------------------------
            if ( Itrm<=0 ) Itrm = 3
            if ( Itmax<=0 ) Itmax = 20
            ndv1 = Ndv + 1
            if ( Icndir==0 ) Icndir = ndv1
            if ( Delfun<=0. ) Delfun = .0001
            Ct = -abs(Ct)
            if ( Ct>=0. ) Ct = -.1
            Ctmin = abs(Ctmin)
            if ( Ctmin<=0. ) Ctmin = .004
            Ctl = -abs(Ctl)
            if ( Ctl>=0. ) Ctl = -0.01
            Ctlmin = abs(Ctlmin)
            if ( Ctlmin<=0. ) Ctlmin = .001
            if ( Theta<=0. ) Theta = 1.
            if ( Abobj1<=0. ) Abobj1 = .1
            if ( Alphax<=0. ) Alphax = .1
            if ( Fdch<=0. ) Fdch = .01
            if ( Fdchm<=0. ) Fdchm = .01
!     ------------------------------------------------------------------
!                     INITIALIZE INTERNAL PARAMETERS
!     ------------------------------------------------------------------
            Infog = 0
            Iter = 0
            Jdir = 0
            Iobj = 0
            Kobj = 0
            ndv2 = Ndv + 2
            Kcount = 0
            Ncal(1) = 0
            Ncal(2) = 0
            Nac = 0
            Nfeas = 0
            Mscal = Nscal
            ct1 = Itrm
            ct1 = 1./ct1
            Dct = (Ctmin/abs(Ct))**ct1
            Dctl = (Ctlmin/abs(Ctl))**ct1
            Phi = 5.
            Abobj = Abobj1
            Ncobj = 0
            Ctam = abs(Ctmin)
            Ctbm = abs(Ctlmin)
!     CALCULATE NUMBER OF LINEAR CONSTRAINTS, NLNC.
            Nlnc = 0
            if ( Ncon/=0 ) then
               do i = 1 , Ncon
                  if ( Isc(i)>0 ) Nlnc = Nlnc + 1
               enddo
            endif
!     ------------------------------------------------------------------
!          CHECK TO BE SURE THAT SIDE CONSTRAINTS ARE SATISFIED
!     ------------------------------------------------------------------
            if ( Nside/=0 ) then
               do i = 1 , Ndv
                  if ( Vlb(i)>Vub(i) ) then
                     xx = .5*(Vlb(i)+Vub(i))
                     x(i) = xx
                     Vlb(i) = xx
                     Vub(i) = xx
                     write (Iout,99016) i
99016                format (///5x,'* * CONMIN DETECTS VLB(I).GT.VUB(I)'/5x,                                                       &
                            &'FIX IS SET X(I)=VLB(I)=VUB(I) = .5*(VLB(I)+VUB(I) FOR I =',i5)
                  endif
                  xx = x(i) - Vlb(i)
                  if ( xx>=0. ) then
                     xx = Vub(i) - x(i)
                     if ( xx<0. ) then
                        write (Iout,99018) x(i) , Vub(i) , i
99018                   format (///5x,'* * CONMIN DETECTS INITIAL X(I).GT.VUB(I)'/5x,'X(I) =',e12.4,2x,'VUB(I) =',e12.4/5x,        &
                               &'X(I) IS SET EQUAL TO VUB(I) FOR I =',i5)
                        x(i) = Vub(i)
                     endif
                  else
!     LOWER BOUND VIOLATED.
                     write (Iout,99017) x(i) , Vlb(i) , i
99017                format (///5x,'* * CONMIN DETECTS INITIAL X(I).LT.VLB(I)'/5x,'X(I) =',e12.4,2x,'VLB(I) =',e12.4/5x,           &
                            &'X(I) IS SET EQUAL TO VLB(I) FOR I =',i5)
                     x(i) = Vlb(i)
                  endif
               enddo
            endif
!     ------------------------------------------------------------------
!                        INITIALIZE SCALING VECTOR, SCAL
!     ------------------------------------------------------------------
            if ( Nscal/=0 ) then
               if ( Nscal<0 ) then
                  do i = 1 , Ndv
                     si = abs(Scal(i))
                     if ( si<1.0e-20 ) si = 1.0e-5
                     Scal(i) = si
                     si = 1./si
                     x(i) = x(i)*si
                     if ( Nside/=0 ) then
                        Vlb(i) = Vlb(i)*si
                        Vub(i) = Vub(i)*si
                     endif
                  enddo
               else
                  do i = 1 , Ndv
                     Scal(i) = 1.
                  enddo
               endif
            endif
!     ------------------------------------------------------------------
!     ***** CALCULATE INITIAL FUNCTION AND CONSTRAINT VALUES  *****
!     ------------------------------------------------------------------
            Info = 1
            Ncal(1) = 1
            Igoto = 1
            spag_nextblock_1 = 11
            cycle SPAG_DispatchLoop_1
         else
!     TOTALLY UNCONSTRAINED FUNCTION WITH LINEAR OBJECTIVE.
!     SOLUTION IS UNBOUNDED.
            write (Iout,99001) Linobj , Ncon , Nside
!     ------------------------------------------------------------------
!                                FORMATS
!     ------------------------------------------------------------------
!
!
99001       format (///5x,'A COMPLETELY UNCONSTRAINED FUNCTION WITH A LINEAR OBJECTIVE IS SPECIFIED'//10x,'LINOBJ =',i5/10x,       &
                   &'NCON   =',i5/10x,'NSIDE  =',i5//5x,'CONTROL RETURNED TO CALLING PROGRAM')
            return
         endif
      case (2)
!     ------------------------------------------------------------------
!     ********************  BEGIN MINIMIZATION  ************************
!     ------------------------------------------------------------------
         Iter = Iter + 1
         if ( Abobj1<.0001 ) Abobj1 = .0001
         if ( Abobj1>.2 ) Abobj1 = .2
         if ( Alphax>1. ) Alphax = 1.
         if ( Alphax<.001 ) Alphax = .001
!
!  THE FOLLOWING TWO LINES OF CODE WERE COMMENTED OUT ON 3/5/81
!
!     NFEAS=NFEAS+1
!     IF (NFEAS.GT.10) GO TO 810
         if ( Iprint>2 ) write (Iout,99035) Iter
99035    format (////5x,'BEGIN ITERATION NUMBER',i5)
         if ( Iprint>3 .and. Ncon>0 ) write (Iout,99036) Ct , Ctl , Phi
99036    format (/5x,'CT =',e14.5,5x,'CTL =',e14.5,5x,'PHI =',e14.5)
         Cta = abs(Ct)
         if ( Ncobj==0 ) then
            if ( Mscal>=Nscal .and. Nscal/=0 ) then
               if ( Nscal>=0 .or. Kcount>=Icndir ) then
                  Mscal = 0
                  Kcount = 0
!     ------------------------------------------------------------------
!                          SCALE VARIABLES
!     ------------------------------------------------------------------
                  do i = 1 , Ndv
                     si = Scal(i)
                     Xi = si*x(i)
                     sib = si
                     if ( Nscal>0 ) si = abs(Xi)
                     if ( si>=1.0e-10 ) then
                        Scal(i) = si
                        si = 1./si
                        x(i) = Xi*si
                        if ( Nside/=0 ) then
                           Vlb(i) = sib*si*Vlb(i)
                           Vub(i) = sib*si*Vub(i)
                        endif
                     endif
                  enddo
                  if ( .not.(Iprint<4 .or. (Nscal<0 .and. Iter>1)) ) then
                     write (Iout,99037)
99037                format (/5x,'NEW SCALING VECTOR (SCAL)')
                     write (Iout,99050) (Scal(i),i=1,Ndv)
                  endif
               endif
            endif
            Mscal = Mscal + 1
            Nac = 0
!     ------------------------------------------------------------------
!          OBTAIN GRADIENTS OF OBJECTIVE AND ACTIVE CONSTRAINTS
!     ------------------------------------------------------------------
            Info = 2
            Ncal(2) = Ncal(2) + 1
            if ( Nfdg/=1 ) then
               Jgoto = 0
            else
               Igoto = 2
               spag_nextblock_1 = 11
               cycle SPAG_DispatchLoop_1
            endif
         else
!     ------------------------------------------------------------------
!     NO MOVE ON LAST ITERATION.  DELETE CONSTRAINTS THAT ARE NO
!     LONGER ACTIVE.
!     ------------------------------------------------------------------
            nnac = Nac
            do i = 1 , nnac
               if ( Ic(i)>Ncon ) Nac = Nac - 1
            enddo
            if ( Nac>0 ) then
               nnac = Nac
               SPAG_Loop_1_2: do i = 1 , nnac
                  SPAG_Loop_2_1: do
                     nic = Ic(i)
                     ct1 = Ct
                     if ( Isc(nic)>0 ) ct1 = Ctl
                     if ( g(nic)>ct1 ) exit SPAG_Loop_2_1
                     Nac = Nac - 1
                     if ( i>Nac ) exit SPAG_Loop_1_2
                     do k = i , Nac
                        ii = k + 1
                        do j = 1 , ndv2
                           a(j,k) = a(j,ii)
                        enddo
                        Ic(k) = Ic(ii)
                     enddo
                  enddo SPAG_Loop_2_1
               enddo SPAG_Loop_1_2
            endif
            spag_nextblock_1 = 5
            cycle SPAG_DispatchLoop_1
         endif
         spag_nextblock_1 = 3
      case (3)
         call cnmn01(Jgoto,x,Df,g,Isc,Ic,a,g1,Vlb,Vub,Scal,c,Ncal,Dx,Dx1,Fi,Xi,Iii,n1,n2,n3,n4)
         Igoto = 3
         if ( Jgoto>0 ) then
            spag_nextblock_1 = 11
            cycle SPAG_DispatchLoop_1
         endif
         spag_nextblock_1 = 4
      case (4)
         Info = 1
         if ( Nac>=n3 ) then
            spag_nextblock_1 = 10
            cycle SPAG_DispatchLoop_1
         endif
         if ( Nscal/=0 .and. Nfdg/=0 ) then
!     ------------------------------------------------------------------
!                              SCALE GRADIENTS
!     ------------------------------------------------------------------
!     SCALE GRADIENT OF OBJECTIVE FUNCTION.
            do i = 1 , Ndv
               Df(i) = Df(i)*Scal(i)
            enddo
            if ( Nfdg/=2 .and. Nac/=0 ) then
!     SCALE GRADIENTS OF ACTIVE CONSTRAINTS.
               do j = 1 , Ndv
                  scj = Scal(j)
                  do i = 1 , Nac
                     a(j,i) = a(j,i)*scj
                  enddo
               enddo
            endif
         endif
         spag_nextblock_1 = 5
      case (5)
         if ( Iprint>=3 .and. Ncon/=0 ) then
!     ------------------------------------------------------------------
!                                   PRINT
!     ------------------------------------------------------------------
!     PRINT ACTIVE AND VIOLATED CONSTRAINT NUMBERS.
            m1 = 0
            m2 = n3
            if ( Nac/=0 ) then
               do i = 1 , Nac
                  j = Ic(i)
                  if ( j<=Ncon ) then
                     gi = g(j)
                     c1 = Ctam
                     if ( Isc(j)>0 ) c1 = Ctbm
                     gi = gi - c1
                     if ( gi>0. ) then
                        m2 = m2 + 1
!     VIOLATED CONSTRAINT.
                        Ms1(m2) = j
                     else
!     ACTIVE CONSTRAINT.
                        m1 = m1 + 1
                        Ms1(m1) = j
                     endif
                  endif
               enddo
            endif
            m3 = m2 - n3
            write (Iout,99010) m1
            if ( m1/=0 ) then
               write (Iout,99011)
               write (Iout,99052) (Ms1(i),i=1,m1)
            endif
            write (Iout,99012) m3
            if ( m3/=0 ) then
               write (Iout,99011)
               m3 = n3 + 1
               write (Iout,99052) (Ms1(i),i=m3,m2)
            endif
         endif
!     ------------------------------------------------------------------
!            CALCULATE GRADIENTS OF ACTIVE SIDE CONSTRAINTS
!     ------------------------------------------------------------------
         if ( Nside/=0 ) then
            mcn1 = Ncon
            m1 = 0
            do i = 1 , Ndv
!     LOWER BOUND.
               Xi = x(i)
               xid = Vlb(i)
               x12 = abs(xid)
               if ( x12<1. ) x12 = 1.
               gi = (xid-Xi)/x12
               if ( gi>=-1.0e-6 ) then
                  m1 = m1 + 1
                  Ms1(m1) = -i
                  Nac = Nac + 1
                  if ( Nac>=n3 ) then
                     spag_nextblock_1 = 10
                     cycle SPAG_DispatchLoop_1
                  endif
                  mcn1 = mcn1 + 1
                  do j = 1 , Ndv
                     a(j,Nac) = 0.
                  enddo
                  a(i,Nac) = -1.
                  Ic(Nac) = mcn1
                  g(mcn1) = gi
                  Isc(mcn1) = 1
               endif
!     UPPER BOUND.
               xid = Vub(i)
               x12 = abs(xid)
               if ( x12<1. ) x12 = 1.
               gi = (Xi-xid)/x12
               if ( gi>=-1.0e-6 ) then
                  m1 = m1 + 1
                  Ms1(m1) = i
                  Nac = Nac + 1
                  if ( Nac>=n3 ) then
                     spag_nextblock_1 = 10
                     cycle SPAG_DispatchLoop_1
                  endif
                  mcn1 = mcn1 + 1
                  do j = 1 , Ndv
                     a(j,Nac) = 0.
                  enddo
                  a(i,Nac) = 1.
                  Ic(Nac) = mcn1
                  g(mcn1) = gi
                  Isc(mcn1) = 1
               endif
            enddo
!     ------------------------------------------------------------------
!                                  PRINT
!     ------------------------------------------------------------------
!     PRINT ACTIVE SIDE CONSTRAINT NUMBERS.
            if ( Iprint>=3 ) then
               write (Iout,99013) m1
               if ( m1/=0 ) then
                  write (Iout,99014)
                  write (6,99052) (Ms1(j),j=1,m1)
               endif
            endif
         endif
!     PRINT GRADIENTS OF ACTIVE AND VIOLATED CONSTRAINTS.
         if ( Iprint>=4 ) then
            write (Iout,99038)
99038       format (/5x,'GRADIENT OF OBJ')
            do i = 1 , Ndv , 6
               m1 = min0(Ndv,i+5)
               write (Iout,99005) i , (Df(j),j=i,m1)
            enddo
            if ( Nac/=0 ) then
               write (Iout,99039)
99039          format (/5x,'GRADIENTS OF ACTIVE AND VIOLATED CONSTRAINTS')
               do i = 1 , Nac
                  m1 = Ic(i)
                  m2 = m1 - Ncon
                  m3 = 0
                  if ( m2>0 ) m3 = iabs(Ms1(m2))
                  if ( m2<=0 ) write (Iout,99003) m1
99003             format (5x,'CONSTRAINT NUMBER',i5)
                  if ( m2>0 ) write (Iout,99004) m3
99004             format (5x,'SIDE CONSTRAINT ON VARIABLE',i5)
                  do k = 1 , Ndv , 6
                     m1 = min0(Ndv,k+5)
                     write (Iout,99005) k , (a(j,i),j=k,m1)
                  enddo
                  write (Iout,99040)
               enddo
            endif
         endif
!     ------------------------------------------------------------------
!     ******************  DETERMINE SEARCH DIRECTION *******************
!     ------------------------------------------------------------------
         Alp = 1.0e+20
         if ( Nac>0 ) then
!     ------------------------------------------------------------------
!                          CONSTRAINED FUNCTION
!     ------------------------------------------------------------------
!     FIND USABLE-FEASIBLE DIRECTION.
            Kcount = 0
            Jdir = 0
            Phi = 10.*Phi
            if ( Phi>1000. ) Phi = 1000.
!
!  THE FOLLOWING LINE OF CODE WAS COMMENTED OUT ON 3/5/81
!
!     IF (NFEAS.EQ.1) PHI=5.
!     CALCULATE DIRECTION, S.
            call cnmn05(g,Df,a,s,b,c,Slope,Phi,Isc,Ic,Ms1,Nvc,n1,n2,n3,n4,n5)
!
!  THE FOLLOWING LINE WAS ADDED ON 2/25/81
!
            if ( Nac/=0 ) then
!
!  THE FOLLOWING FIVE LINES WERE COMMENTED OUT ON 3/5/81
!  REASON : THEY WERE NOT IN G. VANDERPLAATS LISTING
!
!     IF THIS DESIGN IS FEASIBLE AND LAST ITERATION WAS INFEASIBLE,
!     SET ABOBJ1=.05 (5 PERCENT).
!     IF (NVC.EQ.0.AND.NFEAS.GT.1) ABOBJ1=.05
!     IF (NVC.EQ.0) NFEAS=0
               if ( Iprint>=3 ) then
                  write (Iout,99041)
99041             format (/5x,'PUSH-OFF FACTORS, (THETA(I), I=1,NAC)')
                  do i = 1 , Nac , 6
                     m1 = min0(Nac,i+5)
                     write (Iout,99005) i , (a(ndv1,j),j=i,m1)
                  enddo
                  write (Iout,99025) s(ndv1)
99025             format (/5x,'CONSTRAINT PARAMETER, BETA =',e14.5)
               endif
!     ------------------------------------------------------------------
!     ****************** ONE-DIMENSIONAL SEARCH ************************
!     ------------------------------------------------------------------
               if ( s(ndv1)<1.0e-6 .and. Nvc==0 ) then
                  spag_nextblock_1 = 9
                  cycle SPAG_DispatchLoop_1
               endif
!     ------------------------------------------------------------------
!                 FIND ALPHA TO OBTAIN A FEASIBLE DESIGN
!     ------------------------------------------------------------------
               if ( Nvc/=0 ) then
                  Alp = -1.
                  do i = 1 , Nac
                     nci = Ic(i)
                     c1 = g(nci)
                     ctc = Ctam
                     if ( Isc(nci)>0 ) ctc = Ctbm
                     if ( c1>ctc ) then
                        alp1 = 0.
                        do j = 1 , Ndv
                           alp1 = alp1 + s(j)*a(j,i)
                        enddo
                        alp1 = alp1*a(ndv2,i)
                        if ( abs(alp1)>=1.0e-20 ) then
                           alp1 = -c1/alp1
                           if ( alp1>Alp ) Alp = alp1
                        endif
                     endif
                  enddo
               endif
               spag_nextblock_1 = 6
               cycle SPAG_DispatchLoop_1
            endif
         endif
!     ------------------------------------------------------------------
!                        UNCONSTRAINED FUNCTION
!     ------------------------------------------------------------------
!     FIND DIRECTION OF STEEPEST DESCENT OR CONJUGATE DIRECTION.
!
!  S. N. 575 ADDED ON 2/25/81
!
         Nvc = 0
         Nfeas = 0
         Kcount = Kcount + 1
!     IF KCOUNT.GT.ICNDIR  RESTART CONJUGATE DIRECTION ALGORITHM.
         if ( Kcount>Icndir .or. Iobj==2 ) Kcount = 1
         if ( Kcount==1 ) Jdir = 0
!     IF JDIR = 0 FIND DIRECTION OF STEEPEST DESCENT.
         call cnmn02(Jdir,Slope,Dftdf1,Df,s,n1)
         spag_nextblock_1 = 6
      case (6)
!     ------------------------------------------------------------------
!                       LIMIT CHANCE TO ABOBJ1*OBJ
!     ------------------------------------------------------------------
         alp1 = 1.0e+20
         si = abs(Obj)
         if ( si<.01 ) si = .01
         if ( abs(Slope)>1.0e-20 ) alp1 = Abobj1*si/Slope
         alp1 = abs(alp1)
         if ( Nvc>0 ) alp1 = 10.*alp1
         if ( alp1<Alp ) Alp = alp1
!     ------------------------------------------------------------------
!                   LIMIT CHANGE IN VARIABLE TO ALPHAX
!     ------------------------------------------------------------------
         alp11 = 1.0e+20
         do i = 1 , Ndv
            si = abs(s(i))
            Xi = abs(x(i))
            if ( si>=1.0e-10 .and. Xi>=0.1 ) then
               alp1 = Alphax*Xi/si
               if ( alp1<alp11 ) alp11 = alp1
            endif
         enddo
         if ( Nvc>0 ) alp11 = 10.*alp11
         if ( alp11<Alp ) Alp = alp11
         if ( Alp>1.0e+20 ) Alp = 1.0e+20
         if ( Alp<=1.0e-20 ) Alp = 1.0e-20
         if ( Iprint>=3 ) then
            write (Iout,99042)
99042       format (/5x,'SEARCH DIRECTION (S-VECTOR)')
            do i = 1 , Ndv , 6
               m1 = min0(Ndv,i+5)
               write (Iout,99005) i , (s(j),j=i,m1)
            enddo
            write (Iout,99015) Slope , Alp
99015       format (/5x,'ONE-DIMENSIONAL SEARCH'/5x,'INITIAL SLOPE =',e12.4,2x,'PROPOSED ALPHA =',e12.4)
         endif
         if ( Ncon>0 .or. Nside>0 ) then
!     ------------------------------------------------------------------
!       SOLVE ONE-DIMENSIONAL SEARCH PROBLEM FOR CONSTRAINED FUNCTION
!     ------------------------------------------------------------------
            Jgoto = 0
            spag_nextblock_1 = 8
            cycle SPAG_DispatchLoop_1
         else
!     ------------------------------------------------------------------
!           DO ONE-DIMENSIONAL SEARCH FOR UNCONSTRAINED FUNCTION
!     ------------------------------------------------------------------
            Jgoto = 0
         endif
         spag_nextblock_1 = 7
      case (7)
         call cnmn03(x,s,Slope,Alp,Fff,a1,a2,a3,a4,f1,f2,f3,f4,App,n1,Ncal,Kount,Jgoto)
         Igoto = 4
         if ( Jgoto>0 ) then
            spag_nextblock_1 = 11
            cycle SPAG_DispatchLoop_1
         endif
!     PROCEED TO CONVERGENCE CHECK.
         Jdir = 1
         spag_nextblock_1 = 9
         cycle SPAG_DispatchLoop_1
      case (8)
         call cnmn06(x,Vlb,Vub,g,Scal,Df,s,g1,g2,Ctam,Ctbm,Slope,Alp,a2,a3,a4,f1,f2,f3,Cv1,Cv2,Cv3,Cv4,Alpca,Alpfes,Alpln,Alpmin,  &
                   & Alpnc,Alpsav,Alpsid,Alptot,Isc,n1,n2,Ncal,Nvc,Icount,Igood1,Igood2,Igood3,Igood4,Ibest,Iii,Nlnc,Jgoto)
         Igoto = 5
         if ( Jgoto>0 ) then
            spag_nextblock_1 = 11
            cycle SPAG_DispatchLoop_1
         endif
         if ( Nac==0 ) Jdir = 1
         spag_nextblock_1 = 9
      case (9)
!     ------------------------------------------------------------------
!     *******************     UPDATE ALPHAX   **************************
!     ------------------------------------------------------------------
         if ( Alp>1.0e+19 ) Alp = 0.
!     UPDATE ALPHAX TO BE AVERAGE OF MAXIMUM CHANGE IN X(I)
!     AND ALHPAX.
         alp11 = 0.
         do i = 1 , Ndv
            si = abs(s(i))
            Xi = abs(x(i))
            if ( Xi>=1.0e-10 ) then
               alp1 = Alp*si/Xi
               if ( alp1>alp11 ) alp11 = alp1
            endif
         enddo
         alp11 = .5*(alp11+Alphax)
         alp12 = 5.*Alphax
         if ( alp11>alp12 ) alp11 = alp12
         Alphax = alp11
         Ncobj = Ncobj + 1
!     ABSOLUTE CHANGE IN OBJECTIVE.
         objd = Obj1 - Obj
         objb = abs(objd)
         if ( objb<1.0e-10 ) objb = 0.
         if ( Nac==0 .or. objb>0. ) Ncobj = 0
         if ( Ncobj>1 ) Ncobj = 0
!     ------------------------------------------------------------------
!                                  PRINT
!     ------------------------------------------------------------------
!     PRINT MOVE PARAMETER, NEW X-VECTOR AND CONSTRAINTS.
         if ( Iprint>=3 ) then
            write (Iout,99043) Alp
99043       format (/5x,'CALCULATED ALPHA =',e14.5)
         endif
         if ( Iprint>=2 ) then
            if ( objb<=0. ) then
               if ( Iprint==2 ) write (Iout,99044) Iter , Obj
99044          format (////5x,'ITER =',i5,5x,'OBJ =',e14.5,5x,'NO CHANGE IN OBJ')
               if ( Iprint>2 ) write (Iout,99045) Obj
99045          format (/5x,'OBJ =',e15.6,5x,'NO CHANGE ON OBJ')
            elseif ( Iprint==2 ) then
               write (Iout,99047) Iter , Obj
99047          format (////5x,'ITER =',i5,5x,'OBJ =',e14.5)
            else
               write (Iout,99046) Obj
            endif
            write (Iout,99049)
            do i = 1 , Ndv
               ff1 = 1.
               if ( Nscal/=0 ) ff1 = Scal(i)
               g1(i) = ff1*x(i)
            enddo
            do i = 1 , Ndv , 6
               m1 = min0(Ndv,i+5)
               write (Iout,99005) i , (g1(j),j=i,m1)
            enddo
            if ( Ncon/=0 ) then
               write (Iout,99051)
               do i = 1 , Ncon , 6
                  m1 = min0(Ncon,i+5)
                  write (Iout,99005) i , (g(j),j=i,m1)
               enddo
            endif
         endif
!
!  THE FOLLOWING CODE WAS ADDED ON 3/5/81
!
!  IT HAD NOT BEEN REPORTED AS A FIX TO MAOB
!  BUT WAS SENT TO JEFF STROUD A YEAR AGO
!  SEE OTHER COMMENTS IN CONMIN SUBROUTINE FOR DELETIONS OF CODE
!  ON 3/5/81 PERTAINING TO THIS FIX
!
!
!                   CHECK FEASIBILITY
!
         if ( Ncon>0 ) then
            Nfeasct = 10
            do i = 1 , Ncon
               c1 = Ctam
               if ( Isc(i)>0 ) c1 = Ctbm
               if ( g(i)>c1 ) then
                  Nfeas = Nfeas + 1
                  goto 10
               endif
            enddo
            if ( Nfeas>0 ) Abobj1 = .05
            Nfeas = 0
            Phi = 5.
 10         if ( Nfeas>=Nfeasct ) then
               spag_nextblock_1 = 10
               cycle SPAG_DispatchLoop_1
            endif
         endif
!
!  END OF INSERTED FIX
!
!     ------------------------------------------------------------------
!                          CHECK CONVERGENCE
!     ------------------------------------------------------------------
!     STOP IF ITER EQUALS ITMAX.
         if ( Iter<Itmax ) then
!     ------------------------------------------------------------------
!                     ABSOLUTE CHANGE IN OBJECTIVE
!     ------------------------------------------------------------------
            objb = abs(objd)
            Kobj = Kobj + 1
            if ( objb>=Dabfun .or. Nfeas>0 ) Kobj = 0
!     ------------------------------------------------------------------
!                     RELATIVE CHANGE IN OBJECTIVE
!     ------------------------------------------------------------------
            if ( abs(Obj1)>1.0e-10 ) objd = objd/abs(Obj1)
            Abobj1 = .5*(abs(Abobj)+abs(objd))
            Abobj = abs(objd)
            Iobj = Iobj + 1
            if ( Nvc>0 .or. objd>=Delfun ) Iobj = 0
            if ( Iobj<Itrm .and. Kobj<Itrm ) then
               Obj1 = Obj
!     ------------------------------------------------------------------
!           REDUCE CT IF OBJECTIVE FUNCTION IS CHANGING SLOWLY
!     ------------------------------------------------------------------
               if ( Iobj>=1 .and. Nac/=0 ) then
                  Ct = Dct*Ct
                  Ctl = Ctl*Dctl
                  if ( abs(Ct)<Ctmin ) Ct = -Ctmin
                  if ( abs(Ctl)<Ctlmin ) Ctl = -Ctlmin
               endif
               spag_nextblock_1 = 2
               cycle SPAG_DispatchLoop_1
            endif
         endif
         spag_nextblock_1 = 10
      case (10)
         if ( Nac>=n3 ) write (Iout,99053)
99053    format (/5x,'THE NUMBER OF ACTIVE AND VIOLATED CONSTRAINTS EXCEEDS N3-1.'/5x,                                             &
                &'DIMENSIONED SIZE OF MATRICES A AND B AND VECTOR IC IS INSUFFICIENT'/5x,                                          &
                &'OPTIMIZATION TERMINATED AND CONTROL RETURNED TO MAIN PROGRAM.')
!     ------------------------------------------------------------------
!     ****************  FINAL FUNCTION INFORMATION  ********************
!     ------------------------------------------------------------------
         if ( Nscal/=0 ) then
!     UN-SCALE THE DESIGN VARIABLES.
            do i = 1 , Ndv
               Xi = Scal(i)
               if ( Nside/=0 ) then
                  Vlb(i) = Xi*Vlb(i)
                  Vub(i) = Xi*Vub(i)
               endif
               x(i) = Xi*x(i)
            enddo
         endif
!     ------------------------------------------------------------------
!                           PRINT FINAL RESULTS
!     ------------------------------------------------------------------
         if ( Iprint/=0 .and. Nac<n3 ) then
            write (Iout,99054)
99054       format ('1',////4x,'FINAL OPTIMIZATION INFORMATION')
            write (Iout,99046) Obj
            write (Iout,99049)
            do i = 1 , Ndv , 6
               m1 = min0(Ndv,i+5)
               write (Iout,99005) i , (x(j),j=i,m1)
            enddo
            if ( Ncon/=0 ) then
               write (Iout,99051)
               do i = 1 , Ncon , 6
                  m1 = min0(Ncon,i+5)
                  write (Iout,99005) i , (g(j),j=i,m1)
               enddo
!     DETERMINE WHICH CONSTRAINTS ARE ACTIVE AND PRINT.
               Nac = 0
               Nvc = 0
               do i = 1 , Ncon
                  Cta = Ctam
                  if ( Isc(i)>0 ) Cta = Ctbm
                  gi = g(i)
                  if ( gi>Cta ) then
                     Nvc = Nvc + 1
                     Ms1(Nvc) = i
                  elseif ( gi>=Ct .or. Isc(i)/=0 ) then
                     if ( gi>=Ctl .or. Isc(i)<=0 ) then
                        Nac = Nac + 1
                        Ic(Nac) = i
                     endif
                  endif
               enddo
               write (Iout,99010) Nac
               if ( Nac/=0 ) then
                  write (Iout,99011)
                  write (Iout,99052) (Ic(j),j=1,Nac)
               endif
               write (Iout,99012) Nvc
               if ( Nvc/=0 ) then
                  write (Iout,99011)
                  write (Iout,99052) (Ms1(j),j=1,Nvc)
               endif
            endif
            if ( Nside/=0 ) then
!     DETERMINE WHICH SIDE CONSTRAINTS ARE ACTIVE AND PRINT.
               Nac = 0
               do i = 1 , Ndv
                  Xi = x(i)
                  xid = Vlb(i)
                  x12 = abs(xid)
                  if ( x12<1. ) x12 = 1.
                  gi = (xid-Xi)/x12
                  if ( gi>=-1.0e-6 ) then
                     Nac = Nac + 1
                     Ms1(Nac) = -i
                  endif
                  xid = Vub(i)
                  x12 = abs(xid)
                  if ( x12<1. ) x12 = 1.
                  gi = (Xi-xid)/x12
                  if ( gi>=-1.0e-6 ) then
                     Nac = Nac + 1
                     Ms1(Nac) = i
                  endif
               enddo
               write (Iout,99013) Nac
               if ( Nac/=0 ) then
                  write (Iout,99014)
                  write (Iout,99052) (Ms1(j),j=1,Nac)
               endif
            endif
            write (Iout,99019)
99019       format (/5x,'TERMINATION CRITERION')
            if ( Iter>=Itmax ) write (Iout,99020)
99020       format (10x,'ITER EQUALS ITMAX')
            if ( Nfeas>=Nfeasct ) write (Iout,99021)
99021       format (10x,'NFEASCT CONSECUTIVE ITERATIONS FAILED TO PRODUCE A   FEASIBLE DESIGN')
            if ( Iobj>=Itrm ) write (Iout,99022) Itrm
99022       format (10x,'ABS(1-OBJ(I-1)/OBJ(I)) LESS THAN DELFUN FOR',i3,' ITERATIONS')
            if ( Kobj>=Itrm ) write (Iout,99023) Itrm
99023       format (10x,'ABS(OBJ(I)-OBJ(I-1))   LESS THAN DABFUN FOR',i3,' ITERATIONS')
            write (Iout,99024) Iter
99024       format (/5x,'NUMBER OF ITERATIONS =',i5)
            write (Iout,99055) Ncal(1)
99055       format (/5x,'OBJECTIVE FUNCTION WAS EVALUATED',8x,i5,2x,'TIMES')
            if ( Ncon>0 ) write (Iout,99056) Ncal(1)
99056       format (/5x,'CONSTRAINT FUNCTIONS WERE EVALUATED',i10,2x,'TIMES')
            if ( Nfdg/=0 ) write (Iout,99057) Ncal(2)
99057       format (/5x,'GRADIENT OF OBJECTIVE WAS CALCULATED',i9,2x,'TIMES')
            if ( Ncon>0 .and. Nfdg==1 ) write (Iout,99058) Ncal(2)
99058       format (/5x,'GRADIENTS OF CONSTRAINTS WERE CALCULATED',i5,2x,'TIMES')
         endif
!     ------------------------------------------------------------------
!                   RE-SET BASIC PARAMETERS TO INPUT VALUES
!     ------------------------------------------------------------------
         Itrm = Idm1
         Itmax = Idm2
         Icndir = Idm3
         Delfun = Dm1
         Dabfun = Dm2
         Ct = Dm3
         Ctmin = Dm4
         Ctl = Dm5
         Ctlmin = Dm6
         Theta = Dm7
         Phi = Dm8
         Fdch = Dm9
         Fdchm = Dm10
         Abobj1 = Dm11
         Alphax = Dm12
         Igoto = 0
         spag_nextblock_1 = 11
      case (11)
         if ( Nscal==0 .or. Igoto==0 ) return
!     UN-SCALE VARIABLES.
         do i = 1 , Ndv
            c(i) = x(i)
            x(i) = x(i)*Scal(i)
         enddo
         return
      end select
   enddo SPAG_DispatchLoop_1
99005 format (3x,i5,')',2x,6e13.5)
99010 format (/5x,'THERE ARE',i5,' ACTIVE CONSTRAINTS')
99011 format (5x,'CONSTRAINT NUMBERS ARE')
99012 format (/5x,'THERE ARE',i5,' VIOLATED CONSTRAINTS')
99013 format (/5x,'THERE ARE',i5,' ACTIVE SIDE CONSTRAINTS')
99014 format (5x,'DECISION VARIABLES AT LOWER OR UPPER BOUNDS',' (MINUS INDICATES LOWER BOUND)')
99040 format (' ')
99046 format (/5x,'OBJ =',e15.6)
99049 format (/5x,'DECISION VARIABLES (X-VECTOR)')
99050 format (3x,7e13.4)
99051 format (/5x,'CONSTRAINT VALUES (G-VECTOR)')
99052 format (5x,15i5)
end subroutine cnmn00


subroutine cnmn01(Jgoto,x,Df,g,Isc,Ic,a,g1,Vlb,Vub,Scal,c,Ncal,Dx,Dx1,Fi,Xi,Iii,n1,n2,n3,n4)
   implicit none

   double precision a , Abobj1 , Alphax , c , Ct , Ctl , Ctlmin , Ctmin , Dabfun , Delfun , Df , Dx , Dx1 , Fdch , fdch1 , Fdchm , &
                  & Fi , g , g1 , Obj
   double precision Scal , Theta , Vlb , Vub , x , x1 , Xi
   integer i , i1 , Ic , Icndir , Igoto , Iii , inf , Info , Infog , Iprint , Isc , Iter , Itmax , Itrm , j , Jgoto , Linobj , n1 ,&
         & n2 , n3
   integer n4 , Nac , Ncal , Ncon , Ndv , Nfdg , Nfeasct , Nscal , Nside

   common /cnmn1 / Delfun , Dabfun , Fdch , Fdchm , Ct , Ctmin , Ctl , Ctlmin , Alphax , Abobj1 , Theta , Obj , Ndv , Ncon ,       &
                 & Nside , Iprint , Nfdg , Nscal , Linobj , Itmax , Itrm , Icndir , Igoto , Nac , Info , Infog , Iter , Nfeasct
   dimension x(n1) , Df(n1) , g(n2) , Isc(n2) , Ic(n3) , a(n1,n3) , g1(n2) , Vlb(n1) , Vub(n1) , Scal(n1) , Ncal(2) , c(n4)
!     ROUTINE TO CALCULATE GRADIENT INFORMATION BY FINITE DIFFERENCE.
!     BY G. N. VANDERPLAATS                         JUNE, 1972.
!     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF.
   if ( Jgoto/=1 ) then
      if ( Jgoto==2 ) then
         x(Iii) = Xi
         if ( Nfdg==0 ) Df(Iii) = Dx1*(Obj-Fi)
         if ( Nac/=0 ) then
!     ------------------------------------------------------------------
!             DETERMINE GRADIENT COMPONENTS OF ACTIVE CONSTRAINTS
!     ------------------------------------------------------------------
            do j = 1 , Nac
               i1 = Ic(j)
               a(Iii,j) = Dx1*(g(i1)-g1(i1))
            enddo
         endif
         if ( Iii<Ndv ) then
            call spag_block_1
            return
         endif
         Infog = 0
         Info = inf
         Jgoto = 0
         Obj = Fi
         if ( Ncon==0 ) return
!     ------------------------------------------------------------------
!             STORE CURRENT CONSTRAINT VALUES BACK IN G-VECTOR
!     ------------------------------------------------------------------
         do i = 1 , Ncon
            g(i) = g1(i)
         enddo
         return
      else
         Infog = 0
         inf = Info
         Nac = 0
         if ( Linobj==0 .or. Iter<=1 ) then
!     ------------------------------------------------------------------
!                    GRADIENT OF LINEAR OBJECTIVE
!     ------------------------------------------------------------------
            if ( Nfdg==2 ) Jgoto = 1
            if ( Nfdg==2 ) return
         endif
      endif
   endif
   Jgoto = 0
   if ( Nfdg==2 .and. Ncon==0 ) return
   if ( Ncon/=0 ) then
!     ------------------------------------------------------------------
!       * * * DETERMINE WHICH CONSTRAINTS ARE ACTIVE OR VIOLATED * * *
!     ------------------------------------------------------------------
      do i = 1 , Ncon
         if ( g(i)>=Ct ) then
            if ( Isc(i)<=0 .or. g(i)>=Ctl ) then
               Nac = Nac + 1
               if ( Nac>=n3 ) return
               Ic(Nac) = i
            endif
         endif
      enddo
      if ( Nfdg==2 .and. Nac==0 ) return
      if ( (Linobj>0 .and. Iter>1) .and. Nac==0 ) return
!     ------------------------------------------------------------------
!                  STORE VALUES OF CONSTRAINTS IN G1
!     ------------------------------------------------------------------
      do i = 1 , Ncon
         g1(i) = g(i)
      enddo
   endif
   Jgoto = 0
   if ( Nac==0 .and. Nfdg==2 ) return
!     ------------------------------------------------------------------
!                            CALCULATE GRADIENTS
!     ------------------------------------------------------------------
   Infog = 1
   Info = 1
   Fi = Obj
   Iii = 0
   call spag_block_1
contains
   subroutine spag_block_1
      Iii = Iii + 1
      Xi = x(Iii)
      Dx = Fdch*Xi
      Dx = abs(Dx)
      fdch1 = Fdchm
      if ( Nscal/=0 ) fdch1 = Fdchm/Scal(Iii)
      if ( Dx<fdch1 ) Dx = fdch1
      x1 = Xi + Dx
      if ( Nside/=0 ) then
         if ( x1>Vub(Iii) ) Dx = -Dx
      endif
      Dx1 = 1./Dx
      x(Iii) = Xi + Dx
      Ncal(1) = Ncal(1) + 1
!     ------------------------------------------------------------------
!                         FUNCTION EVALUATION
!     ------------------------------------------------------------------
      Jgoto = 2
      return
   end subroutine spag_block_1
end subroutine cnmn01


subroutine cnmn02(Ncalc,Slope,Dftdf1,Df,s,n1)
   implicit none

   double precision Abobj1 , Alphax , beta , Ct , Ctl , Ctlmin , Ctmin , Dabfun , Delfun , Df , dfi , dftdf , Dftdf1 , Fdch ,      &
                  & Fdchm , Obj , s , s1 , s2 , si
   double precision Slope , Theta
   integer i , Icndir , Igoto , Info , Infog , Iprint , Iter , Itmax , Itrm , Linobj , n1 , Nac , Ncalc , Ncon , Ndv , Nfdg ,      &
         & Nfeasct , Nscal , Nside

   common /cnmn1 / Delfun , Dabfun , Fdch , Fdchm , Ct , Ctmin , Ctl , Ctlmin , Alphax , Abobj1 , Theta , Obj , Ndv , Ncon ,       &
                 & Nside , Iprint , Nfdg , Nscal , Linobj , Itmax , Itrm , Icndir , Igoto , Nac , Info , Infog , Iter , Nfeasct
   dimension Df(n1) , s(n1)
!     ROUTINE TO DETERMINE CONJUGATE DIRECTION VECTOR OR DIRECTION
!     OF STEEPEST DESCENT FOR UNCONSTRAINED FUNCTION MINIMIZATION.
!     BY G. N. VANDERPLAATS                       APRIL, 1972.
!     NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF.
!     NCALC = CALCULATION CONTROL.
!         NCALC = 0,     S = STEEPEST DESCENT.
!         NCALC = 1,     S = CONJUGATE DIRECTION.
!     CONJUGATE DIRECTION IS FOUND BY FLETCHER-REEVES ALGORITHM.
!     ------------------------------------------------------------------
!                   CALCULATE NORM OF GRADIENT VECTOR
!     ------------------------------------------------------------------
   dftdf = 0.
   do i = 1 , Ndv
      dfi = Df(i)
      dftdf = dftdf + dfi*dfi
   enddo
!     ------------------------------------------------------------------
!     **********                FIND DIRECTION S              **********
!     ------------------------------------------------------------------
   if ( Ncalc==1 ) then
      if ( Dftdf1>=1.0e-20 ) then
!     ------------------------------------------------------------------
!                 FIND FLETCHER-REEVES CONJUGATE DIRECTION
!     ------------------------------------------------------------------
         beta = dftdf/Dftdf1
         Slope = 0.
         do i = 1 , Ndv
            dfi = Df(i)
            si = beta*s(i) - dfi
            Slope = Slope + si*dfi
            s(i) = si
         enddo
         call spag_block_1
         return
      endif
   endif
   Ncalc = 0
!     ------------------------------------------------------------------
!                  CALCULATE DIRECTION OF STEEPEST DESCENT
!     ------------------------------------------------------------------
   do i = 1 , Ndv
      s(i) = -Df(i)
   enddo
   Slope = -dftdf
   call spag_block_1
contains
   subroutine spag_block_1
!     ------------------------------------------------------------------
!                  NORMALIZE S TO MAX ABS VALUE OF UNITY
!     ------------------------------------------------------------------
      s1 = 0.
      do i = 1 , Ndv
         s2 = abs(s(i))
         if ( s2>s1 ) s1 = s2
      enddo
      if ( s1<1.0e-20 ) s1 = 1.0e-20
      s1 = 1./s1
      Dftdf1 = dftdf*s1
      do i = 1 , Ndv
         s(i) = s1*s(i)
      enddo
      Slope = s1*Slope
   end subroutine spag_block_1
end subroutine cnmn02


subroutine cnmn03(x,s,Slope,Alp,Fff,a1,a2,a3,a4,f1,f2,f3,f4,App,n1,Ncal,Kount,Jgoto)
   implicit none

   double precision a1 , a2 , a3 , a4 , aa , ab , ab2 , ab3 , Abobj1 , Alp , Alphax , ap , ap1 , App , Ct , Ctl , Ctlmin , Ctmin , &
                  & Dabfun , Delfun
   double precision f1 , f2 , f3 , f4 , Fdch , Fdchm , ff , Fff , Obj , s , Slope , Theta , x , zro
   integer i , Icndir , Igoto , ii , Info , Infog , Iout , Iprint , Iter , Itmax , Itrm , Jgoto , jj , Kount , Linobj , n1 , Nac , &
         & Ncal , Ncon , Ndv
   integer Nfdg , Nfeasct , Nscal , Nside

   common /cnmn1 / Delfun , Dabfun , Fdch , Fdchm , Ct , Ctmin , Ctl , Ctlmin , Alphax , Abobj1 , Theta , Obj , Ndv , Ncon ,       &
                 & Nside , Iprint , Nfdg , Nscal , Linobj , Itmax , Itrm , Icndir , Igoto , Nac , Info , Infog , Iter , Nfeasct
   common /output/ Iout
   dimension x(n1) , s(n1) , Ncal(2)
   integer :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: do
      select case (spag_nextblock_1)
      case (1)
!     ROUTINE TO SOLVE ONE-DIMENSIONAL SEARCH IN UNCONSTRAINED
!     MINIMIZATION USING 2-POINT QUADRATIC INTERPOLATION, 3-POINT
!     CUBIC INTERPOLATION AND 4-POINT CUBIC INTERPOLATION.
!     BY G. N. VANDERPLAATS                         APRIL, 1972.
!     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF.
!     ALP = PROPOSED MOVE PARAMETER.
!     SLOPE = INITIAL FUNCTION SLOPE = S-TRANSPOSE TIMES DF.
!     SLOPE MUST BE NEGATIVE.
!     OBJ = INITIAL FUNCTION VALUE.
         zro = 0.
         if ( Jgoto/=0 ) then
            if ( Jgoto==1 ) then
               f2 = Obj
               if ( Iprint>4 ) write (Iout,99004) f2
               if ( f2<f1 ) then
                  spag_nextblock_1 = 4
                  cycle SPAG_DispatchLoop_1
               endif
!     ------------------------------------------------------------------
!                     CHECK FOR ILL-CONDITIONING
!     ------------------------------------------------------------------
               if ( Kount>5 ) then
                  spag_nextblock_1 = 3
                  cycle SPAG_DispatchLoop_1
               endif
               ff = 2.*abs(f1)
               if ( f2<ff ) then
!     ------------------------------------------------------------------
!     **********        2-POINT QUADRATIC INTERPOLATION       **********
!     ------------------------------------------------------------------
                  jj = 1
                  ii = 1
                  call cnmn04(ii,App,zro,a1,f1,Slope,a2,f2,zro,zro,zro,zro)
                  if ( App<zro .or. App>a2 ) then
                     spag_nextblock_1 = 4
                     cycle SPAG_DispatchLoop_1
                  endif
                  f3 = f2
                  a3 = a2
                  a2 = App
                  jj = 0
!     ------------------------------------------------------------------
!                  UPDATE DESIGN VECTOR AND FUNCTION VALUE
!     ------------------------------------------------------------------
                  ap = a2 - Alp
                  Alp = a2
                  do i = 1 , Ndv
                     x(i) = x(i) + ap*s(i)
                  enddo
                  if ( Iprint>4 ) write (Iout,99002) a2
                  if ( Iprint>4 ) write (Iout,99003) (x(i),i=1,Ndv)
                  Ncal(1) = Ncal(1) + 1
                  Jgoto = 3
                  return
               else
                  ff = 5.*abs(f1)
                  if ( f2<ff ) then
                     spag_nextblock_1 = 3
                     cycle SPAG_DispatchLoop_1
                  endif
                  a2 = .5*a2
                  ap = -a2
                  Alp = a2
                  spag_nextblock_1 = 2
                  cycle SPAG_DispatchLoop_1
               endif
            elseif ( Jgoto==2 ) then
               f2 = Obj
!     PROCEED TO CUBIC INTERPOLATION.
               if ( Iprint>4 ) write (Iout,99004) f2
               spag_nextblock_1 = 6
               cycle SPAG_DispatchLoop_1
            elseif ( Jgoto==3 ) then
               f2 = Obj
               if ( Iprint>4 ) write (Iout,99004) f2
               spag_nextblock_1 = 5
               cycle SPAG_DispatchLoop_1
            elseif ( Jgoto==4 ) then
               f3 = Obj
               if ( Iprint>4 ) write (Iout,99004) f3
               spag_nextblock_1 = 5
               cycle SPAG_DispatchLoop_1
            elseif ( Jgoto==5 ) then
               if ( Iprint>4 ) write (Iout,99004) Obj
!     ------------------------------------------------------------------
!                         CHECK CONVERGENCE
!     ------------------------------------------------------------------
               aa = 1. - App/a2
               ab2 = abs(f2)
               ab3 = abs(Obj)
               ab = ab2
               if ( ab3>ab ) ab = ab3
               if ( ab<1.0e-15 ) ab = 1.0e-15
               ab = (ab2-ab3)/ab
               if ( abs(ab)<1.0e-15 .and. abs(aa)<.001 ) then
                  spag_nextblock_1 = 10
                  cycle SPAG_DispatchLoop_1
               endif
               a4 = a3
               f4 = f3
               a3 = App
               f3 = Obj
               if ( a3<=a2 ) then
                  a3 = a2
                  f3 = f2
                  a2 = App
                  f2 = Obj
               endif
               spag_nextblock_1 = 8
               cycle SPAG_DispatchLoop_1
            elseif ( Jgoto==6 ) then
               f4 = Obj
               if ( Iprint>4 ) write (Iout,99004) f4
               if ( f4>f3 ) then
                  spag_nextblock_1 = 8
                  cycle SPAG_DispatchLoop_1
               endif
               a1 = a2
               f1 = f2
               a2 = a3
               f2 = f3
               a3 = a4
               f3 = f4
               spag_nextblock_1 = 7
               cycle SPAG_DispatchLoop_1
            elseif ( Jgoto==7 ) then
               if ( Iprint>4 ) write (Iout,99004) Obj
               spag_nextblock_1 = 9
               cycle SPAG_DispatchLoop_1
            endif
         endif
!     ------------------------------------------------------------------
!                     INITIAL INFORMATION  (ALPHA=0)
!     ------------------------------------------------------------------
         if ( Slope<0. ) then
            if ( Iprint>4 ) write (Iout,99001)
!     ------------------------------------------------------------------
!                                 FORMATS
!     ------------------------------------------------------------------
!
!
99001       format (/////5x,'* * * UNCONSTRAINED ONE-DIMENSIONAL SEARCH INFORMATION * * *')
            Fff = Obj
            ap1 = 0.
            a1 = 0.
            f1 = Obj
            a2 = Alp
            a3 = 0.
            f3 = 0.
            ap = a2
            Kount = 0
         else
            Alp = 0.
            return
         endif
         spag_nextblock_1 = 2
      case (2)
!     ------------------------------------------------------------------
!            MOVE A DISTANCE AP*S AND UPDATE FUNCTION VALUE
!     ------------------------------------------------------------------
         Kount = Kount + 1
         do i = 1 , Ndv
            x(i) = x(i) + ap*s(i)
         enddo
         if ( Iprint>4 ) write (Iout,99002) ap
         if ( Iprint>4 ) write (Iout,99003) (x(i),i=1,Ndv)
         Ncal(1) = Ncal(1) + 1
         Jgoto = 1
         return
      case (3)
         f3 = f2
         a3 = a2
         a2 = .5*a2
!     ------------------------------------------------------------------
!                 UPDATE DESIGN VECTOR AND FUNCTION VALUE
!     ------------------------------------------------------------------
         ap = a2 - Alp
         Alp = a2
         do i = 1 , Ndv
            x(i) = x(i) + ap*s(i)
         enddo
         if ( Iprint>4 ) write (Iout,99002) a2
         if ( Iprint>4 ) write (Iout,99003) (x(i),i=1,Ndv)
         Ncal(1) = Ncal(1) + 1
         Jgoto = 2
         return
      case (4)
         a3 = 2.*a2
!     ------------------------------------------------------------------
!                  UPDATE DESIGN VECTOR AND FUNCTION VALUE
!     ------------------------------------------------------------------
         ap = a3 - Alp
         Alp = a3
         do i = 1 , Ndv
            x(i) = x(i) + ap*s(i)
         enddo
         if ( Iprint>4 ) write (Iout,99002) a3
         if ( Iprint>4 ) write (Iout,99003) (x(i),i=1,Ndv)
         Ncal(1) = Ncal(1) + 1
         Jgoto = 4
         return
      case (5)
         if ( f3<f2 ) then
            spag_nextblock_1 = 7
            cycle SPAG_DispatchLoop_1
         endif
         spag_nextblock_1 = 6
      case (6)
!     ------------------------------------------------------------------
!     **********       3-POINT CUBIC INTERPOLATION      **********
!     ------------------------------------------------------------------
         ii = 3
         call cnmn04(ii,App,zro,a1,f1,Slope,a2,f2,a3,f3,zro,zro)
         if ( App>=zro .and. App<=a3 ) then
!     ------------------------------------------------------------------
!     UPDATE DESIGN VECTOR AND FUNCTION VALUE.
!     ------------------------------------------------------------------
            ap1 = App
            ap = App - Alp
            Alp = App
            do i = 1 , Ndv
               x(i) = x(i) + ap*s(i)
            enddo
            if ( Iprint>4 ) write (Iout,99002) Alp
            if ( Iprint>4 ) write (Iout,99003) (x(i),i=1,Ndv)
            Ncal(1) = Ncal(1) + 1
            Jgoto = 5
            return
         endif
         spag_nextblock_1 = 7
      case (7)
!     ------------------------------------------------------------------
!     **********        4-POINT CUBIC INTERPOLATION       **********
!     ------------------------------------------------------------------
         a4 = 2.*a3
!     UPDATE DESIGN VECTOR AND FUNCTION VALUE.
         ap = a4 - Alp
         Alp = a4
         do i = 1 , Ndv
            x(i) = x(i) + ap*s(i)
         enddo
         if ( Iprint>4 ) write (Iout,99002) Alp
         if ( Iprint>4 ) write (Iout,99003) (x(i),i=1,Ndv)
         Ncal(1) = Ncal(1) + 1
         Jgoto = 6
         return
      case (8)
         ii = 4
         call cnmn04(ii,App,a1,a1,f1,Slope,a2,f2,a3,f3,a4,f4)
         if ( App>a1 ) then
!     ------------------------------------------------------------------
!                 UPDATE DESIGN VECTOR AND FUNCTION VALUE
!     ------------------------------------------------------------------
            ap = App - Alp
            Alp = App
            do i = 1 , Ndv
               x(i) = x(i) + ap*s(i)
            enddo
            if ( Iprint>4 ) write (Iout,99002) Alp
            if ( Iprint>4 ) write (Iout,99003) (x(i),i=1,Ndv)
            Ncal(1) = Ncal(1) + 1
            Jgoto = 7
            return
         else
            ap = a1 - Alp
            Alp = a1
            Obj = f1
            do i = 1 , Ndv
               x(i) = x(i) + ap*s(i)
            enddo
         endif
         spag_nextblock_1 = 9
      case (9)
!     ------------------------------------------------------------------
!                    CHECK FOR ILL-CONDITIONING
!     ------------------------------------------------------------------
         if ( Obj<=f2 .and. Obj<=f3 ) then
            if ( Obj<=f1 ) then
               spag_nextblock_1 = 10
               cycle SPAG_DispatchLoop_1
            endif
            ap = a1 - Alp
            Alp = a1
            Obj = f1
         elseif ( f2<f3 ) then
            Obj = f2
            ap = a2 - Alp
            Alp = a2
         else
            Obj = f3
            ap = a3 - Alp
            Alp = a3
         endif
!     ------------------------------------------------------------------
!                       UPDATE DESIGN VECTOR
!     ------------------------------------------------------------------
         do i = 1 , Ndv
            x(i) = x(i) + ap*s(i)
         enddo
         spag_nextblock_1 = 10
      case (10)
!     ------------------------------------------------------------------
!                     CHECK FOR MULTIPLE MINIMA
!     ------------------------------------------------------------------
         if ( Obj>Fff ) then
!     INITIAL FUNCTION IS MINIMUM.
            do i = 1 , Ndv
               x(i) = x(i) - Alp*s(i)
            enddo
            Alp = 0.
            Obj = Fff
         endif
         Jgoto = 0
         return
      end select
   enddo SPAG_DispatchLoop_1
99002 format (/5x,'ALPHA =',e14.5/5x,'X-VECTOR')
99003 format (5x,6e13.5)
99004 format (/5x,'OBJ =',e14.5)
end subroutine cnmn03


subroutine cnmn04(Ii,Xbar,Eps,x1,y1,Slope,x2,y2,x3,y3,x4,y4)
   implicit none

   double precision aa , bac , bb , cc , dnom , dx , Eps , q1 , q2 , q3 , q4 , q5 , q6 , qq , Slope , x1 , x11 , x111 , x2 , x21
   double precision x22 , x222 , x3 , x31 , x32 , x33 , x4 , x41 , x42 , x44 , Xbar , xbar1 , y1 , y2 , y3 , y4
   integer Ii , nslop

!     ROUTINE TO FIND FIRST XBAR.GE.EPS CORRESPONDING TO A MINIMUM
!     OF A ONE-DIMENSIONAL REAL FUNCTION BY POLYNOMIEL INTERPOLATION.
!     BY G. N. VANDERPLAATS                          APRIL, 1972.
!     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF.
!
!     II = CALCULATION CONTROL.
!          1:  2-POINT QUADRATIC INTERPOLATION, GIVEN X1, Y1, SLOPE,
!              X2 AND Y2.
!          2:  3-POINT QUADRATIC INTERPOLATION, GIVEN X1, Y1, X2, Y2,
!              X3 AND Y3.
!          3:  3-POINT CUBIC INTERPOLATION, GIVEN X1, Y1, SLOPE, X2, Y2,
!              X3 AND Y3.
!          4:  4-POINT CUBIC INTERPOLATION, GIVEN X1, Y1, X2, Y2, X3,
!              Y3, X4 AND Y4.
!     EPS MAY BE NEGATIVE.
!     IF REQUIRED MINIMUM ON Y DOES NOT EXITS, OR THE FUNCTION IS
!     ILL-CONDITIONED, XBAR = EPS-1.0 WILL BE RETURNED AS AN ERROR
!     INDICATOR.
!     IF DESIRED INTERPOLATION IS ILL-CONDITIONED, A LOWER ORDER
!     INTERPOLATION, CONSISTANT WITH INPUT DATA, WILL BE ATTEMPTED,
!     AND II WILL BE CHANGED ACCORDINGLY.
   xbar1 = Eps - 1.
   Xbar = xbar1
   x21 = x2 - x1
   if ( abs(x21)<1.0e-20 ) return
   nslop = mod(Ii,2)
   if ( Ii==2 ) then
      call spag_block_2
      return
   endif
   if ( Ii==3 ) then
      call spag_block_3
      return
   endif
   if ( Ii==4 ) then
!     ------------------------------------------------------------------
!                    II=4: 4-POINT CUBIC INTERPOLATION
!     ------------------------------------------------------------------
      x21 = x2 - x1
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
      if ( abs(q2)<1.0e-30 ) return
      q1 = x111*x32 - x222*x31 + x3*x33*x21
      q4 = x111*x42 - x222*x41 + x4*x44*x21
      q5 = x41*x21*x42
      dnom = q2*q4 - q1*q5
      if ( abs(dnom)<1.0e-30 ) then
         call spag_block_4
         return
      endif
      q3 = y3*x21 - y2*x31 + y1*x32
      q6 = y4*x21 - y2*x41 + y1*x42
      aa = (q2*q6-q3*q5)/dnom
      bb = (q3-q1*aa)/q2
      cc = (y2-y1-aa*(x222-x111))/x21 - bb*(x1+x2)
      bac = bb*bb - 3.*aa*cc
      if ( abs(aa)<1.0e-20 .or. bac<0. ) then
         call spag_block_4
         return
      endif
      bac = sqrt(bac)
      Xbar = (bac-bb)/(3.*aa)
      if ( Xbar<Eps ) Xbar = xbar1
      return
   endif
   call spag_block_1
contains
   subroutine spag_block_1
!     ------------------------------------------------------------------
!                 II=1: 2-POINT QUADRATIC INTERPOLATION
!     ------------------------------------------------------------------
      Ii = 1
      dx = x1 - x2
      if ( abs(dx)<1.0e-20 ) return
      aa = (Slope+(y2-y1)/dx)/dx
      if ( aa<1.0e-20 ) return
      bb = Slope - 2.*aa*x1
      Xbar = -.5*bb/aa
      if ( Xbar<Eps ) Xbar = xbar1
      return
   end subroutine spag_block_1
   subroutine spag_block_2
!     ------------------------------------------------------------------
!                 II=2: 3-POINT QUADRATIC INTERPOLATION
!     ------------------------------------------------------------------
      Ii = 2
      x21 = x2 - x1
      x31 = x3 - x1
      x32 = x3 - x2
      qq = x21*x31*x32
      if ( abs(qq)<1.0e-20 ) return
      aa = (y1*x32-y2*x31+y3*x21)/qq
      if ( aa<1.0e-20 ) then
         if ( nslop==0 ) return
         call spag_block_1
         return
      else
         bb = (y2-y1)/x21 - aa*(x1+x2)
         Xbar = -.5*bb/aa
         if ( Xbar<Eps ) Xbar = xbar1
         return
      endif
   end subroutine spag_block_2
   subroutine spag_block_3
!     ------------------------------------------------------------------
!                   II=3: 3-POINT CUBIC INTERPOLATION
!     ------------------------------------------------------------------
      Ii = 3
      x21 = x2 - x1
      x31 = x3 - x1
      x32 = x3 - x2
      qq = x21*x31*x32
      if ( abs(qq)<1.0e-20 ) return
      x11 = x1*x1
      dnom = x2*x2*x31 - x11*x32 - x3*x3*x21
      if ( abs(dnom)<1.0e-20 ) then
         call spag_block_2
         return
      endif
      aa = ((x31*x31*(y2-y1)-x21*x21*(y3-y1))/(x31*x21)-Slope*x32)/dnom
      if ( abs(aa)<1.0e-20 ) then
         call spag_block_2
         return
      endif
      bb = ((y2-y1)/x21-Slope-aa*(x2*x2+x1*x2-2.*x11))/x21
      cc = Slope - 3.*aa*x11 - 2.*bb*x1
      bac = bb*bb - 3.*aa*cc
      if ( bac<0. ) then
         call spag_block_2
         return
      endif
      bac = sqrt(bac)
      Xbar = (bac-bb)/(3.*aa)
      if ( Xbar<Eps ) Xbar = Eps
      return
   end subroutine spag_block_3
   subroutine spag_block_4
      if ( nslop/=1 ) then
         call spag_block_2
         return
      endif
      call spag_block_3
      return
   end subroutine spag_block_4
end subroutine cnmn04


subroutine cnmn05(g,Df,a,s,b,c,Slope,Phi,Isc,Ic,Ms1,Nvc,n1,n2,n3,n4,n5)
   implicit none

   double precision a , a1 , Abobj1 , Alphax , b , c , c1 , Ct , ct1 , ct2 , cta , ctam , ctb , ctbm , ctc , ctd , Ctl , Ctlmin ,  &
                  & Ctmin , Dabfun
   double precision Delfun , Df , Fdch , Fdchm , g , gg , Obj , Phi , s , s1 , sg , Slope , Theta , thmax , tht
   integer i , Ic , Icndir , Igoto , Info , Infog , Iout , Iprint , Isc , Iter , Itmax , Itrm , j , j1 , k , Linobj , Ms1 , n1 ,   &
         & n2 , n3
   integer n4 , n5 , Nac , nac1 , nci , ncj , Ncon , ndb , Ndv , ndv1 , ndv2 , ner , Nfdg , Nfeasct , Nscal , Nside , Nvc

   common /cnmn1 / Delfun , Dabfun , Fdch , Fdchm , Ct , Ctmin , Ctl , Ctlmin , Alphax , Abobj1 , Theta , Obj , Ndv , Ncon ,       &
                 & Nside , Iprint , Nfdg , Nscal , Linobj , Itmax , Itrm , Icndir , Igoto , Nac , Info , Infog , Iter , Nfeasct
   common /output/ Iout
   dimension Df(n1) , g(n2) , Isc(n2) , Ic(n3) , a(n1,n3) , s(n1) , c(n4) , Ms1(n5) , b(n3,n3)
   integer :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: do
      select case (spag_nextblock_1)
      case (1)
!     ROUTINE TO SOLVE DIRECTION FINDING PROBLEM IN MODIFIED METHOD OF
!     FEASIBLE DIRECTIONS.
!     BY G. N. VANDERPLAATS                            MAY, 1972.
!     NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF.
!     NORM OF S VECTOR USED HERE IS S-TRANSPOSE TIMES S.LE.1.
!     IF NVC = 0 FIND DIRECTION BY ZOUTENDIJK'S METHOD.  OTHERWISE
!     FIND MODIFIED DIRECTION.
!     ------------------------------------------------------------------
!     ***  NORMALIZE GRADIENTS, CALCULATE THETA'S AND DETERMINE NVC  ***
!     ------------------------------------------------------------------
         ndv1 = Ndv + 1
         ndv2 = Ndv + 2
         nac1 = Nac + 1
         Nvc = 0
         thmax = 0.
         cta = abs(Ct)
         ct1 = 1./cta
         ctam = abs(Ctmin)
         ctb = abs(Ctl)
         ct2 = 1./ctb
         ctbm = abs(Ctlmin)
         a1 = 1.
         do i = 1 , Nac
!     CALCULATE THETA
            nci = Ic(i)
            ncj = 1
            if ( nci<=Ncon ) ncj = Isc(nci)
            c1 = g(nci)
            ctd = ct1
            ctc = ctam
            if ( ncj>0 ) then
               ctc = ctbm
               ctd = ct2
            endif
            if ( c1>ctc ) Nvc = Nvc + 1
            tht = 0.
            gg = 1. + ctd*c1
            if ( ncj==0 .or. c1>ctc ) tht = Theta*gg*gg
            if ( tht>50. ) tht = 50.
            if ( tht>thmax ) thmax = tht
            a(ndv1,i) = tht
!     ------------------------------------------------------------------
!                    NORMALIZE GRADIENTS OF CONSTRAINTS
!     ------------------------------------------------------------------
            a(ndv2,i) = 1.
            if ( nci<=Ncon ) then
               a1 = 0.
               do j = 1 , Ndv
                  a1 = a1 + a(j,i)**2
               enddo
               if ( a1<1.0e-20 ) a1 = 1.0e-20
               a1 = sqrt(a1)
               a(ndv2,i) = a1
               a1 = 1./a1
               do j = 1 , Ndv
                  a(j,i) = a1*a(j,i)
               enddo
            endif
         enddo
!     ------------------------------------------------------------------
!     CHECK FOR ZERO GRADIENT.  PROGRAM CHANGE-FEB, 1981, GV.
!     ------------------------------------------------------------------
         i = 0
         spag_nextblock_1 = 2
      case (2)
         i = i + 1
         SPAG_Loop_1_1: do while ( a(ndv2,i)<=1.0e-6 )
!     ZERO GRADIENT IS FOUND.  WRITE ERROR MESSAGE.
            if ( Iprint>=2 ) write (Iout,99001) Ic(i)
99001       format (5x,'** CONSTRAINT',i5,' HAS ZERO GRADIENT'/5x,'DELETED FROM ACTIVE SET')
!     REDUCE NAC BY ONE.
            Nac = Nac - 1
!     SHIFT COLUMNS OF A AND ROWS OF IC IF I.LE.NAC.
            if ( i>Nac ) then
               spag_nextblock_1 = 3
               cycle SPAG_DispatchLoop_1
            endif
!     SHIFT.
            do j = i , Nac
               j1 = j + 1
               Ic(j) = Ic(j1)
               do k = 1 , ndv2
                  a(k,j) = a(k,j1)
               enddo
            enddo
            if ( i>Nac ) exit SPAG_Loop_1_1
         enddo SPAG_Loop_1_1
         if ( i<Nac ) then
            spag_nextblock_1 = 2
            cycle SPAG_DispatchLoop_1
         endif
         spag_nextblock_1 = 3
      case (3)
         if ( Nac<=0 ) return
         nac1 = Nac + 1
!     DETERMINE IF CONSTRAINTS ARE VIOLATED.
         Nvc = 0
         do i = 1 , Nac
            nci = Ic(i)
            ncj = 1
            if ( nci<=Ncon ) ncj = Isc(nci)
            ctc = ctam
            if ( ncj>0 ) ctc = ctbm
            if ( g(nci)>ctc ) Nvc = Nvc + 1
         enddo
!     ------------------------------------------------------------------
!     NORMALIZE GRADIENT OF OBJECTIVE FUNCTION AND STORE IN NAC+1
!     COLUMN OF A
!     ------------------------------------------------------------------
         a1 = 0.
         do i = 1 , Ndv
            a1 = a1 + Df(i)**2
         enddo
         if ( a1<1.0e-20 ) a1 = 1.0e-20
         a1 = sqrt(a1)
         a1 = 1./a1
         do i = 1 , Ndv
            a(i,nac1) = a1*Df(i)
         enddo
!     BUILD C VECTOR.
         if ( Nvc>0 ) then
!     ------------------------------------------------------------------
!                   BUILD C FOR MODIFIED METHOD
!     ------------------------------------------------------------------
            ndb = Nac
            a(ndv1,nac1) = -Phi
!     ------------------------------------------------------------------
!           SCALE THETA'S SO THAT MAXIMUM THETA IS UNITY
!     ------------------------------------------------------------------
            if ( thmax>0.00001 ) thmax = 1./thmax
            do i = 1 , ndb
               a(ndv1,i) = a(ndv1,i)*thmax
            enddo
            do i = 1 , ndb
               c(i) = 0.
               do j = 1 , ndv1
                  c(i) = c(i) + a(j,i)*a(j,nac1)
               enddo
            enddo
         else
!     ------------------------------------------------------------------
!                 BUILD C FOR CLASSICAL METHOD
!     ------------------------------------------------------------------
            ndb = nac1
            a(ndv1,ndb) = 1.
            do i = 1 , ndb
               c(i) = -a(ndv1,i)
            enddo
         endif
!     ------------------------------------------------------------------
!                      BUILD B MATRIX
!     ------------------------------------------------------------------
         do i = 1 , ndb
            do j = 1 , ndb
               b(i,j) = 0.
               do k = 1 , ndv1
                  b(i,j) = b(i,j) - a(k,i)*a(k,j)
               enddo
            enddo
         enddo
!     ------------------------------------------------------------------
!                    SOLVE SPECIAL L. P. PROBLEM
!     ------------------------------------------------------------------
         call cnmn08(ndb,ner,c,Ms1,b,n3,n4,n5)
         if ( Iprint>1 .and. ner>0 ) write (Iout,99002)
!     ------------------------------------------------------------------
!                           FORMATS
!     ------------------------------------------------------------------
!
!
99002    format (//5x,'* * DIRECTION FINDING PROCESS DID NOT CONVERGE'/5x,'* * S-VECTOR MAY NOT BE VALID')
!     CALCULATE RESULTING DIRECTION VECTOR, S.
         Slope = 0.
!     ------------------------------------------------------------------
!                  USABLE-FEASIBLE DIRECTION
!     ------------------------------------------------------------------
         do i = 1 , Ndv
            s1 = 0.
            if ( Nvc>0 ) s1 = -a(i,nac1)
            do j = 1 , ndb
               s1 = s1 - a(i,j)*c(j)
            enddo
            Slope = Slope + s1*Df(i)
            s(i) = s1
         enddo
         s(ndv1) = 1.
         if ( Nvc>0 ) s(ndv1) = -a(ndv1,nac1)
         do j = 1 , ndb
            s(ndv1) = s(ndv1) - a(ndv1,j)*c(j)
         enddo
!     ------------------------------------------------------------------
!     CHECK TO INSURE THE S-VECTOR IS FEASIBLE.
!     PROGRAM MOD-FEB, 1981, GV.
!     ------------------------------------------------------------------
         do j = 1 , Nac
!     S DOT DEL(G).
            sg = 0.
            do i = 1 , Ndv
               sg = sg + s(i)*a(i,j)
            enddo
!     IF(SG.GT.0.) GO TO 176
!
!  THIS CHANGE MADE ON 4/8/81 FOR G. VANDERPLAATS
!
            if ( sg>1.0e-04 ) then
               spag_nextblock_1 = 4
               cycle SPAG_DispatchLoop_1
            endif
!     FEASIBLE FOR THIS CONSTRAINT.  CONTINUE.
         enddo
!     ------------------------------------------------------------------
!                  NORMALIZE S TO MAX ABS OF UNITY
!     ------------------------------------------------------------------
         s1 = 0.
         do i = 1 , Ndv
            a1 = abs(s(i))
            if ( a1>s1 ) s1 = a1
         enddo
!     IF (S1.LT.1.0E-10) RETURN
!
!  E-10 CHANGED TO E-04 ON 1/12/81
!
         if ( s1<1.0e-04 ) return
         s1 = 1./s1
         do i = 1 , Ndv
            s(i) = s1*s(i)
         enddo
         Slope = s1*Slope
         s(ndv1) = s1*s(ndv1)
         return
      case (4)
!     S-VECTOR IS NOT FEASIBLE DUE TO SOME NUMERICAL PROBLEM.
         if ( Iprint>=2 ) write (Iout,99003)
99003    format (5x,'** CALCULATED S-VECTOR IS NOT FEASIBLE'/5x,'BETA IS SET TO ZERO')
         s(ndv1) = 0.
         Nvc = 0
         return
      end select
   enddo SPAG_DispatchLoop_1
end subroutine cnmn05


subroutine cnmn06(x,Vlb,Vub,g,Scal,Df,s,g1,g2,Ctam,Ctbm,Slope,Alp,a2,a3,a4,f1,f2,f3,Cv1,Cv2,Cv3,Cv4,Alpca,Alpfes,Alpln,Alpmin,     &
                & Alpnc,Alpsav,Alpsid,Alptot,Isc,n1,n2,Ncal,Nvc,Icount,Igood1,Igood2,Igood3,Igood4,Ibest,Iii,Nlnc,Jgoto)
   implicit none

   double precision a2 , a3 , a4 , Abobj1 , Alp , alpa , alpb , Alpca , Alpfes , Alphax , Alpln , Alpmin , Alpnc , Alpsav ,        &
                  & Alpsid , Alptot , c1 , c2 , c3 , cc
   double precision Ct , Ctam , Ctbm , Ctl , Ctlmin , Ctmin , Cv1 , Cv2 , Cv3 , Cv4 , Dabfun , Delfun , Df , f1 , f2 , f3 , f4 ,   &
                  & Fdch , Fdchm , g
   double precision g1 , g2 , gi , Obj , s , Scal , si , Slope , Theta , Vlb , Vub , x , xi , xi1 , xi2 , zro
   integer i , Ibest , Icndir , Icount , Igood1 , Igood2 , Igood3 , Igood4 , Igoto , ii , Iii , Info , Infog , Iout , Iprint ,     &
         & Isc , Iter , Itmax , Itrm , jbest
   integer Jgoto , ksid , Linobj , n1 , n2 , Nac , Ncal , Ncon , Ndv , Nfdg , Nfeasct , Nlnc , Nscal , Nside , Nvc , nvc1

   common /cnmn1 / Delfun , Dabfun , Fdch , Fdchm , Ct , Ctmin , Ctl , Ctlmin , Alphax , Abobj1 , Theta , Obj , Ndv , Ncon ,       &
                 & Nside , Iprint , Nfdg , Nscal , Linobj , Itmax , Itrm , Icndir , Igoto , Nac , Info , Infog , Iter , Nfeasct
   common /output/ Iout
   dimension x(n1) , Vlb(n1) , Vub(n1) , g(n2) , Scal(n1) , Df(n1) , s(n1) , g1(n2) , g2(n2) , Isc(n2) , Ncal(2)
   integer :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: do
      select case (spag_nextblock_1)
      case (1)
!     ROUTINE TO SOLVE ONE-DIMENSIONAL SEARCH PROBLEM FOR CONSTRAINED
!     FUNCTION MINIMIZATION.
!     BY G. N. VANDERPLAATS                           AUG., 1974.
!     NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF.
!     OBJ = INITIAL AND FINAL FUNCTION VALUE.
!     ALP = MOVE PARAMETER.
!     SLOPE = INITIAL SLOPE.
!
!     ALPSID = MOVE TO SIDE CONSTRAINT.
!     ALPFES = MOVE TO FEASIBLE REGION.
!     ALPNC = MOVE TO NEW NON-LINEAR CONSTRAINT.
!     ALPLN = MOVE TO LINEAR CONSTRAINT.
!     ALPCA = MOVE TO RE-ENCOUNTER CURRENTLY ACTIVE CONSTRAINT.
!     ALPMIN = MOVE TO MINIMIZE FUNCTION.
!     ALPTOT = TOTAL MOVE PARAMETER.
         zro = 0.
         if ( Jgoto/=0 ) then
            if ( Jgoto==1 ) then
               f2 = Obj
               if ( Iprint>=5 ) write (Iout,99005) f2
               if ( Iprint>=5 .and. Ncon/=0 ) then
                  write (Iout,99006)
                  write (Iout,99004) (g(i),i=1,Ncon)
               endif
!     ------------------------------------------------------------------
!               IDENTIFY ACCAPTABILITY OF DESIGNS F1 AND F2
!     ------------------------------------------------------------------
!     IGOOD = 0 IS ACCAPTABLE.
!     CV = MAXIMUM CONSTRAINT VIOLATION.
               Igood1 = 0
               Igood2 = 0
               Cv1 = 0.
               Cv2 = 0.
               nvc1 = 0
               if ( Ncon/=0 ) then
                  do i = 1 , Ncon
                     cc = Ctam
                     if ( Isc(i)>0 ) cc = Ctbm
                     c1 = g1(i) - cc
                     c2 = g(i) - cc
                     if ( c2>0. ) nvc1 = nvc1 + 1
                     if ( c1>Cv1 ) Cv1 = c1
                     if ( c2>Cv2 ) Cv2 = c2
                  enddo
                  if ( Cv1>0. ) Igood1 = 1
                  if ( Cv2>0. ) Igood2 = 1
               endif
               Alp = a2
               Obj = f2
!     ------------------------------------------------------------------
!     IF F2 VIOLATES FEWER CONSTRAINTS THAN F1 BUT STILL HAS CONSTRAINT
!     VIOLATIONS RETURN
!     ------------------------------------------------------------------
               if ( nvc1<Nvc .and. nvc1>0 ) then
                  spag_nextblock_1 = 9
                  cycle SPAG_DispatchLoop_1
               endif
!     ------------------------------------------------------------------
!             IDENTIFY BEST OF DESIGNS F1 ANF F2
!     ------------------------------------------------------------------
!     IBEST CORRESPONDS TO MINIMUM VALUE DESIGN.
!     IF CONSTRAINTS ARE VIOLATED, IBEST CORRESPONDS TO MINIMUM
!     CONSTRAINT VIOLATION.
               if ( Igood1==0 .and. Igood2==0 ) then
!     NO CONSTRAINT VIOLATION.  PICK MINIMUM F.
                  Ibest = 1
                  if ( f2<=f1 ) Ibest = 2
               else
!     VIOLATED CONSTRAINTS.  PICK MINIMUM VIOLATION.
                  Ibest = 1
                  if ( Cv1>=Cv2 ) Ibest = 2
               endif
               ii = 1
!     ------------------------------------------------------------------
!     IF CV2 IS GREATER THAN CV1, SET MOVE LIMITS TO A2.
!     PROGRAM MOD-FEB, 1981, GV.
!     ------------------------------------------------------------------
               if ( Cv2>Cv1 ) then
                  Alpln = a2
                  Alpnc = a2
                  Alpca = a2
               endif
               if ( Ncon==0 ) then
                  spag_nextblock_1 = 3
                  cycle SPAG_DispatchLoop_1
               endif
!     ------------------------------------------------------------------
!     *****                 2 - POINT INTERPOLATION                *****
!     ------------------------------------------------------------------
               Iii = 0
               do
                  Iii = Iii + 1
                  c1 = g1(Iii)
                  c2 = g(Iii)
                  if ( Isc(Iii)==0 ) then
!     ------------------------------------------------------------------
!                     NON-LINEAR CONSTRAINT
!     ------------------------------------------------------------------
                     if ( c1<1.0e-5 .or. c1>Ctam ) then
                        call cnmn07(ii,Alp,zro,zro,c1,a2,c2,zro,zro)
                        if ( Alp>0. ) then
                           if ( c1>Ctam .and. Alp>Alpfes ) Alpfes = Alp
                           if ( c1<Ct .and. Alp<Alpnc ) Alpnc = Alp
                        endif
                     endif
!     ------------------------------------------------------------------
!                        LINEAR CONSTRAINT
!     ------------------------------------------------------------------
                  elseif ( c1<1.0e-5 .or. c1>Ctbm ) then
                     call cnmn07(ii,Alp,zro,zro,c1,a2,c2,zro,zro)
                     if ( Alp>0. ) then
                        if ( c1>Ctbm .and. Alp>Alpfes ) Alpfes = Alp
                        if ( c1<Ctl .and. Alp<Alpln ) Alpln = Alp
                     endif
                  endif
                  if ( Iii>=Ncon ) then
                     spag_nextblock_1 = 3
                     cycle SPAG_DispatchLoop_1
                  endif
               enddo
            elseif ( Jgoto==2 ) then
               f3 = Obj
               if ( Iprint>=5 ) write (Iout,99005) f3
               if ( Iprint>=5 .and. Ncon/=0 ) then
                  write (Iout,99006)
                  write (Iout,99004) (g(i),i=1,Ncon)
               endif
!     ------------------------------------------------------------------
!       CALCULATE MAXIMUM CONSTRAINT VIOLATION AND PICK BEST DESIGN
!     ------------------------------------------------------------------
               Cv3 = 0.
               Igood3 = 0
               nvc1 = 0
               if ( Ncon/=0 ) then
                  do i = 1 , Ncon
                     cc = Ctam
                     if ( Isc(i)>0 ) cc = Ctbm
                     c1 = g(i) - cc
                     if ( c1>Cv3 ) Cv3 = c1
                     if ( c1>0. ) nvc1 = nvc1 + 1
                  enddo
                  if ( Cv3>0. ) Igood3 = 1
               endif
!     DETERMINE BEST DESIGN.
               if ( Ibest==2 ) then
!     CHOOSE BETWEEN F2 AND F3.
                  if ( Igood2==0 .and. Igood3==0 ) then
                     if ( f3<=f2 ) Ibest = 3
                  else
                     if ( Cv2>=Cv3 ) Ibest = 3
                  endif
!     CHOOSE BETWEEN F1 AND F3.
               elseif ( Igood1==0 .and. Igood3==0 ) then
                  if ( f3<=f1 ) Ibest = 3
               else
                  if ( Cv1>=Cv3 ) Ibest = 3
               endif
               Alp = a3
               Obj = f3
!     IF F3 VIOLATES FEWER CONSTRAINTS THAN F1 RETURN.
               if ( nvc1<Nvc ) then
                  spag_nextblock_1 = 9
                  cycle SPAG_DispatchLoop_1
               endif
!     IF OBJECTIVE AND ALL CONSTRAINTS ARE LINEAR, RETURN.
               if ( Linobj/=0 .and. Nlnc==Ncon ) then
                  spag_nextblock_1 = 9
                  cycle SPAG_DispatchLoop_1
               endif
!     IF A3 = ALPLN AND F3 IS BOTH GOOD AND BEST RETURN.
               alpb = 1. - Alpln/a3
               if ( (abs(alpb)<1.0e-20 .and. Ibest==3) .and. (Igood3==0) ) then
                  spag_nextblock_1 = 9
                  cycle SPAG_DispatchLoop_1
               endif
!     IF A3 = ALPSID AND F3 IS BEST, GO INVOKE SIDE CONSTRAINT
!     MODIFICATION.
               alpa = 1. - Alpsid/a3
               if ( abs(alpa)<1.0e-20 .and. Ibest==3 ) then
                  spag_nextblock_1 = 2
                  cycle SPAG_DispatchLoop_1
               endif
!     ------------------------------------------------------------------
!     **********            3 - POINT INTERPOLATION            *********
!     ------------------------------------------------------------------
               Alpnc = Alpsid
               Alpca = Alpsid
               Alpfes = -1.
               Alpmin = -1.
!     ------------------------------------------------------------------
!     IF A3 IS GREATER THAN A2 AND CV3 IS GREATER THAN CV2, SET
!     MOVE LIMITS TO A3.  PROGRAM MOD-FEB, 1981, GV.
!     ------------------------------------------------------------------
               if ( a3>a2 .and. Cv3>Cv2 ) then
                  Alpln = a3
                  Alpnc = a3
                  Alpca = a3
               endif
               if ( Ncon==0 ) then
                  spag_nextblock_1 = 6
                  cycle SPAG_DispatchLoop_1
               endif
               Iii = 0
               spag_nextblock_1 = 4
               cycle SPAG_DispatchLoop_1
            elseif ( Jgoto==3 ) then
               f4 = Obj
               if ( Iprint>=5 ) write (Iout,99005) f4
               if ( Iprint>=5 .and. Ncon/=0 ) then
                  write (Iout,99006)
                  write (Iout,99004) (g(i),i=1,Ncon)
               endif
!     DETERMINE ACCAPTABILITY OF F4.
               Igood4 = 0
               Cv4 = 0.
               if ( Ncon/=0 ) then
                  do i = 1 , Ncon
                     cc = Ctam
                     if ( Isc(i)>0 ) cc = Ctbm
                     c1 = g(i) - cc
                     if ( c1>Cv4 ) Cv4 = c1
                  enddo
                  if ( Cv4>0. ) Igood4 = 1
               endif
               Alp = a4
               Obj = f4
!     ------------------------------------------------------------------
!                     DETERMINE BEST DESIGN
!     ------------------------------------------------------------------
               if ( Ibest==2 ) then
!     CHOOSE BETWEEN F2 AND F4.
                  if ( Igood2==0 .and. Igood4==0 ) then
                     if ( f4>f2 ) then
                        spag_nextblock_1 = 7
                        cycle SPAG_DispatchLoop_1
                     endif
                     spag_nextblock_1 = 9
                     cycle SPAG_DispatchLoop_1
                  else
                     if ( Cv2<=Cv4 ) then
                        spag_nextblock_1 = 7
                        cycle SPAG_DispatchLoop_1
                     endif
                     spag_nextblock_1 = 9
                     cycle SPAG_DispatchLoop_1
                  endif
               elseif ( Ibest==3 ) then
!     CHOOSE BETWEEN F3 AND F4.
                  if ( Igood3==0 .and. Igood4==0 ) then
                     if ( f4>f3 ) then
                        spag_nextblock_1 = 8
                        cycle SPAG_DispatchLoop_1
                     endif
                     spag_nextblock_1 = 9
                     cycle SPAG_DispatchLoop_1
                  else
                     if ( Cv3<=Cv4 ) then
                        spag_nextblock_1 = 8
                        cycle SPAG_DispatchLoop_1
                     endif
                     spag_nextblock_1 = 9
                     cycle SPAG_DispatchLoop_1
                  endif
               else
!     CHOOSE BETWEEN F1 AND F4.
                  if ( Igood1==0 .and. Igood4==0 ) then
                     if ( f4<=f1 ) then
                        spag_nextblock_1 = 9
                        cycle SPAG_DispatchLoop_1
                     endif
                  elseif ( Cv1>Cv4 ) then
                     spag_nextblock_1 = 9
                     cycle SPAG_DispatchLoop_1
                  endif
!     F1 IS BEST.
                  Alptot = Alptot - a4
                  Obj = f1
                  do i = 1 , Ndv
                     x(i) = x(i) - a4*s(i)
                  enddo
                  if ( Ncon/=0 ) then
                     do i = 1 , Ncon
                        g(i) = g1(i)
                     enddo
                  endif
                  spag_nextblock_1 = 9
                  cycle SPAG_DispatchLoop_1
               endif
            endif
         endif
         if ( Iprint>=5 ) write (Iout,99002)
99002    format (/////'* * * CONSTRAINED ONE-DIMENSIONAL SEARCH INFORMATION * * *')
         Alpsav = Alp
         Icount = 0
         Alptot = 0.
!     TOLERANCES.
         Ctam = abs(Ctmin)
         Ctbm = abs(Ctlmin)
         spag_nextblock_1 = 2
      case (2)
!     PROPOSED MOVE.
!     ------------------------------------------------------------------
!     *****  BEGIN SEARCH OR IMPOSE SIDE CONSTRAINT MODIFICATION  *****
!     ------------------------------------------------------------------
         a2 = Alpsav
         Icount = Icount + 1
         Alpsid = 1.0e+20
!     INITIAL ALPHA AND OBJ.
         Alp = 0.
         f1 = Obj
         ksid = 0
         if ( Nside/=0 ) then
!     ------------------------------------------------------------------
!     FIND MOVE TO SIDE CONSTRAINT AND INSURE AGAINST VIOLATION OF
!     SIDE CONSTRAINTS
!     ------------------------------------------------------------------
            do i = 1 , Ndv
               si = s(i)
               if ( abs(si)>1.0e-20 ) then
                  xi = x(i)
                  si = 1./si
                  if ( si>0. ) then
!     UPPER BOUND.
                     xi2 = Vub(i)
                     xi1 = abs(xi2)
                     if ( xi1<1. ) xi1 = 1.
!     CONSTRAINT VALUE.
                     gi = (xi-xi2)/xi1
                     if ( gi<=-1.0e-6 ) then
!     PROPOSED MOVE TO UPPER BOUND.
                        alpa = (xi2-xi)*si
                        if ( alpa<Alpsid ) Alpsid = alpa
                        cycle
                     endif
                  else
!     LOWER BOUND.
                     xi2 = Vlb(i)
                     xi1 = abs(xi2)
                     if ( xi1<1. ) xi1 = 1.
!     CONSTRAINT VALUE.
                     gi = (xi2-xi)/xi1
                     if ( gi<=-1.0e-6 ) then
!     PROPOSED MOVE TO LOWER BOUND.
                        alpa = (xi2-xi)*si
                        if ( alpa<Alpsid ) Alpsid = alpa
                        cycle
                     endif
                  endif
!     MOVE WILL VIOLATE SIDE CONSTRAINT.  SET S(I)=0.
                  Slope = Slope - s(i)*Df(i)
                  s(i) = 0.
                  ksid = ksid + 1
               else
!     ITH COMPONENT OF S IS SMALL.  SET TO ZERO.
                  s(i) = 0.
                  Slope = Slope - si*Df(i)
               endif
            enddo
!     ALPSID IS UPPER BOUND ON ALPHA.
            if ( a2>Alpsid ) a2 = Alpsid
         endif
!     ------------------------------------------------------------------
!               CHECK ILL-CONDITIONING
!     ------------------------------------------------------------------
         if ( ksid==Ndv .or. Icount>10 ) then
            spag_nextblock_1 = 9
            cycle SPAG_DispatchLoop_1
         endif
         if ( Nvc==0 .and. Slope>0. ) then
            spag_nextblock_1 = 9
            cycle SPAG_DispatchLoop_1
         endif
         Alpfes = -1.
         Alpmin = -1.
         Alpln = 1.1*Alpsid
         Alpnc = Alpsid
         Alpca = Alpsid
         if ( Ncon/=0 ) then
!     STORE CONSTRAINT VALUES IN G1.
            do i = 1 , Ncon
               g1(i) = g(i)
            enddo
         endif
!     ------------------------------------------------------------------
!                  MOVE A DISTANCE A2*S
!     ------------------------------------------------------------------
         Alptot = Alptot + a2
         do i = 1 , Ndv
            x(i) = x(i) + a2*s(i)
         enddo
         if ( Iprint>=5 ) then
            write (Iout,99003) a2
            if ( Nscal==0 ) then
               write (Iout,99004) (x(i),i=1,Ndv)
            else
               do i = 1 , Ndv
                  g(i) = Scal(i)*x(i)
               enddo
               write (Iout,99004) (g(i),i=1,Ndv)
            endif
         endif
!     ------------------------------------------------------------------
!                   UPDATE FUNCTION AND CONSTRAINT VALUES
!     ------------------------------------------------------------------
         Ncal(1) = Ncal(1) + 1
         Jgoto = 1
         return
      case (3)
!     CALCULATE ALPHA TO MINIMIZE FUNCTION.
         if ( Linobj<=0 .and. Slope<0. ) call cnmn04(ii,Alpmin,zro,zro,f1,Slope,a2,f2,zro,zro,zro,zro)
!     ------------------------------------------------------------------
!                         PROPOSED MOVE
!     ------------------------------------------------------------------
!     MOVE AT LEAST FAR ENOUGH TO OVERCOME CONSTRAINT VIOLATIONS.
         a3 = Alpfes
!     MOVE TO MINIMIZE FUNCTION.
         if ( Alpmin>a3 ) a3 = Alpmin
!     IF A3.LE.0, SET A3 = ALPSID.
         if ( a3<=0. ) a3 = Alpsid
!     LIMIT MOVE TO NEW CONSTRAINT ENCOUNTER.
         if ( a3>Alpnc ) a3 = Alpnc
         if ( a3>Alpln ) a3 = Alpln
!     MAKE A3 NON-ZERO.
         if ( a3<=1.0e-20 ) a3 = 1.0e-20
!     IF A3=A2=ALPSID AND F2 IS BEST, GO INVOKE SIDE CONSTRAINT
!     MODIFICATION.
         alpb = 1. - a2/a3
         alpa = 1. - Alpsid/a3
         jbest = 0
         if ( abs(alpb)<1.0e-10 .and. abs(alpa)<1.0e-10 ) jbest = 1
         if ( jbest==1 .and. Ibest==2 ) then
            spag_nextblock_1 = 2
            cycle SPAG_DispatchLoop_1
         endif
!     SIDE CONSTRAINT CHECK NOT SATISFIED.
         if ( Ncon/=0 ) then
!     STORE CONSTRAINT VALUES IN G2.
            do i = 1 , Ncon
               g2(i) = g(i)
            enddo
         endif
!     IF A3=A2, SET A3=.9*A2.
         if ( abs(alpb)<1.0e-10 ) a3 = .9*a2
!     MOVE AT LEAST .01*A2.
         if ( a3<(.01*a2) ) a3 = .01*a2
!     LIMIT MOVE TO 5.*A2.
         if ( a3>(5.*a2) ) a3 = 5.*a2
!     LIMIT MOVE TO ALPSID.
         if ( a3>Alpsid ) a3 = Alpsid
!     MOVE A DISTANCE A3*S.
         Alp = a3 - a2
         Alptot = Alptot + Alp
         do i = 1 , Ndv
            x(i) = x(i) + Alp*s(i)
         enddo
         if ( Iprint>=5 ) then
            write (Iout,99007)
99007       format (/5x,'TWO-POINT INTERPOLATION')
            write (Iout,99003) a3
            if ( Nscal==0 ) then
               write (Iout,99004) (x(i),i=1,Ndv)
            else
               do i = 1 , Ndv
                  g(i) = Scal(i)*x(i)
               enddo
               write (Iout,99004) (g(i),i=1,Ndv)
            endif
         endif
!     ------------------------------------------------------------------
!              UPDATE FUNCTION AND CONSTRAINT VALUES
!     ------------------------------------------------------------------
         Ncal(1) = Ncal(1) + 1
         Jgoto = 2
         return
      case (4)
         Iii = Iii + 1
         c1 = g1(Iii)
         c2 = g2(Iii)
         c3 = g(Iii)
         if ( Isc(Iii)==0 ) then
!     ------------------------------------------------------------------
!                     NON-LINEAR CONSTRAINT
!     ------------------------------------------------------------------
            ii = 2
            call cnmn07(ii,Alp,zro,zro,c1,a2,c2,a3,c3)
            if ( Alp>zro ) then
               if ( c1<Ct .or. c1>0. ) then
                  if ( c1>Ctam .or. c1<0. ) then
                     if ( Alp>Alpfes .and. c1>Ctam ) Alpfes = Alp
                     if ( Alp<Alpnc .and. c1<0. ) Alpnc = Alp
                     spag_nextblock_1 = 5
                     cycle SPAG_DispatchLoop_1
                  endif
               endif
!     ALP IS MINIMUM MOVE.  UPDATE FOR NEXT CONSTRAINT ENCOUNTER.
               alpa = Alp
               call cnmn07(ii,Alp,alpa,zro,c1,a2,c2,a3,c3)
               if ( Alp<Alpca .and. Alp>=alpa ) Alpca = Alp
            endif
!     ------------------------------------------------------------------
!     LINEAR CONSTRAINT.  FIND ALPFES ONLY.  ALPLN SAME AS BEFORE.
!     ------------------------------------------------------------------
         elseif ( c1>Ctbm ) then
            ii = 1
            call cnmn07(ii,Alp,zro,zro,c1,a3,c3,zro,zro)
            if ( Alp>Alpfes ) Alpfes = Alp
         endif
         spag_nextblock_1 = 5
      case (5)
         if ( Iii<Ncon ) then
            spag_nextblock_1 = 4
            cycle SPAG_DispatchLoop_1
         endif
         spag_nextblock_1 = 6
      case (6)
         if ( Linobj<=0 .and. Slope<=0. ) then
!     ------------------------------------------------------------------
!              CALCULATE ALPHA TO MINIMIZE FUNCTION
!     ------------------------------------------------------------------
            ii = 3
            if ( a2>a3 .and. (Igood2==0 .and. Ibest==2) ) ii = 2
            call cnmn04(ii,Alpmin,zro,zro,f1,Slope,a2,f2,a3,f3,zro,zro)
         endif
!     ------------------------------------------------------------------
!                       PROPOSED MOVE
!     ------------------------------------------------------------------
!     MOVE AT LEAST ENOUGH TO OVERCOME CONSTRAINT VIOLATIONS.
         a4 = Alpfes
!     MOVE TO MINIMIZE FUNCTION.
         if ( Alpmin>a4 ) a4 = Alpmin
!     IF A4.LE.0, SET A4 = ALPSID.
         if ( a4<=0. ) a4 = Alpsid
!     LIMIT MOVE TO NEW CONSTRAINT ENCOUNTER.
         if ( a4>Alpln ) a4 = Alpln
         if ( a4>Alpnc ) a4 = Alpnc
!     LIMIT MOVE TO RE-ENCOUNTER CURRENTLY ACTIVE CONSTRAINT.
         if ( a4>Alpca ) a4 = Alpca
!     LIMIT A4 TO 5.*A3.
         if ( a4>(5.*a3) ) a4 = 5.*a3
!     UPDATE DESIGN.
         if ( Ibest==3 .and. Ncon/=0 ) then
!     STORE CONSTRAINT VALUES IN G2.  F3 IS BEST.  F2 IS NOT.
            do i = 1 , Ncon
               g2(i) = g(i)
            enddo
         endif
!     IF A4=A3 AND IGOOD1=0 AND IGOOD3=1, SET A4=.9*A3.
         Alp = a4 - a3
         if ( (Igood1==0 .and. Igood3==1) .and. abs(Alp)<1.0e-20 ) a4 = .9*a3
!     ------------------------------------------------------------------
!                   MOVE A DISTANCE A4*S
!     ------------------------------------------------------------------
         Alp = a4 - a3
         Alptot = Alptot + Alp
         do i = 1 , Ndv
            x(i) = x(i) + Alp*s(i)
         enddo
         if ( Iprint>=5 ) then
            write (Iout,99001)
!     ------------------------------------------------------------------
!                                  FORMATS
!     ------------------------------------------------------------------
!
!
99001       format (/5x,'THREE-POINT INTERPOLATION')
            write (Iout,99003) a4
            if ( Nscal==0 ) then
               write (Iout,99004) (x(i),i=1,Ndv)
            else
               do i = 1 , Ndv
                  g(i) = Scal(i)*x(i)
               enddo
               write (Iout,99004) (g(i),i=1,Ndv)
            endif
         endif
!     ------------------------------------------------------------------
!              UPDATE FUNCTION AND CONSTRAINT VALUES
!     ------------------------------------------------------------------
         Ncal(1) = Ncal(1) + 1
         Jgoto = 3
         return
      case (7)
!     F2 IS BEST.
         Obj = f2
         a2 = a4 - a2
         Alptot = Alptot - a2
         do i = 1 , Ndv
            x(i) = x(i) - a2*s(i)
         enddo
         if ( Ncon/=0 ) then
            do i = 1 , Ncon
               g(i) = g2(i)
            enddo
         endif
         spag_nextblock_1 = 9
         cycle SPAG_DispatchLoop_1
      case (8)
!     F3 IS BEST.
         Obj = f3
         a3 = a4 - a3
         Alptot = Alptot - a3
         do i = 1 , Ndv
            x(i) = x(i) - a3*s(i)
         enddo
         if ( Ncon/=0 ) then
            do i = 1 , Ncon
               g(i) = g2(i)
            enddo
         endif
         spag_nextblock_1 = 9
      case (9)
         Alp = Alptot
         if ( Iprint>=5 ) write (Iout,99008)
99008    format (/5x,'* * * END OF ONE-DIMENSIONAL SEARCH')
         Jgoto = 0
         return
      end select
   enddo SPAG_DispatchLoop_1
99003 format (//5x,'PROPOSED DESIGN'/5x,'ALPHA =',e12.5/5x,'X-VECTOR')
99004 format (1x,8e12.4)
99005 format (/5x,'OBJ =',e13.5)
99006 format (/5x,'CONSTRAINT VALUES')
end subroutine cnmn06


subroutine cnmn07(Ii,Xbar,Eps,x1,y1,x2,y2,x3,y3)
   implicit none

   double precision aa , bac , bb , cc , dy , Eps , qq , x1 , x2 , x21 , x3 , x31 , x32 , xb2 , Xbar , xbar1 , y1 , y2 , y3 , yy
   integer Ii , jj

!     ROUTINE TO FIND FIRST XBAR.GE.EPS CORRESPONDING TO A REAL ZERO
!     OF A ONE-DIMENSIONAL FUNCTION BY POLYNOMIEL INTERPOLATION.
!     BY G. N. VANDERPLAATS                          APRIL, 1972.
!     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF.
!     II = CALCULATION CONTROL.
!          1:  2-POINT LINEAR INTERPOLATION, GIVEN X1, Y1, X2 AND Y2.
!          2:  3-POINT QUADRATIC INTERPOLATION, GIVEN X1, Y1, X2, Y2,
!              X3 AND Y3.
!     EPS MAY BE NEGATIVE.
!     IF REQUIRED ZERO ON Y DOES NOT EXITS, OR THE FUNCTION IS
!     ILL-CONDITIONED, XBAR = EPS-1.0 WILL BE RETURNED AS AN ERROR
!     INDICATOR.
!     IF DESIRED INTERPOLATION IS ILL-CONDITIONED, A LOWER ORDER
!     INTERPOLATION, CONSISTANT WITH INPUT DATA, WILL BE ATTEMPTED AND
!     II WILL BE CHANGED ACCORDINGLY.
   xbar1 = Eps - 1.
   Xbar = xbar1
   jj = 0
   x21 = x2 - x1
   if ( abs(x21)<1.0e-20 ) return
   if ( Ii==2 ) then
!     ------------------------------------------------------------------
!                 II=2: 3-POINT QUADRATIC INTERPOLATION
!     ------------------------------------------------------------------
      jj = 1
      x31 = x3 - x1
      x32 = x3 - x2
      qq = x21*x31*x32
      if ( abs(qq)<1.0e-20 ) return
      aa = (y1*x32-y2*x31+y3*x21)/qq
      if ( abs(aa)>=1.0e-20 ) then
         bb = (y2-y1)/x21 - aa*(x1+x2)
         cc = y1 - x1*(aa*x1+bb)
         bac = bb*bb - 4.*aa*cc
         if ( bac>=0. ) then
            bac = sqrt(bac)
            aa = .5/aa
            Xbar = aa*(bac-bb)
            xb2 = -aa*(bac+bb)
            if ( Xbar<Eps ) Xbar = xb2
            if ( xb2<Xbar .and. xb2>Eps ) Xbar = xb2
            if ( Xbar<Eps ) Xbar = xbar1
            return
         endif
      endif
   endif
!
!     ------------------------------------------------------------------
!                  II=1: 2-POINT LINEAR INTERPOLATION
!     ------------------------------------------------------------------
   Ii = 1
   yy = y1*y2
   if ( jj/=0 .and. yy>=0. ) then
!     INTERPOLATE BETWEEN X2 AND X3.
      dy = y3 - y2
      if ( abs(dy)>=1.0e-20 ) then
         Xbar = x2 + y2*(x2-x3)/dy
         if ( Xbar<Eps ) Xbar = xbar1
         return
      endif
   endif
   dy = y2 - y1
!     INTERPOLATE BETWEEN X1 AND X2.
   if ( abs(dy)<1.0e-20 ) return
   Xbar = x1 + y1*(x1-x2)/dy
   if ( Xbar<Eps ) Xbar = xbar1
   return
end subroutine cnmn07


subroutine cnmn08(Ndb,Ner,c,Ms1,b,n3,n4,n5)
   implicit none

   double precision b , bb , bb1 , bi , c , c1 , cb , cbmax , cbmin , eps
   integer i , ichk , iter1 , j , jj , kk , m2 , Ms1 , n3 , n4 , n5 , Ndb , Ner , nmax

   dimension c(n4) , b(n3,n3) , Ms1(n5)
!     ROUTINE TO SOLVE SPECIAL LINEAR PROBLEM FOR IMPOSING S-TRANSPOSE
!     TIMES S.LE.1 BOUNDS IN THE MODIFIED METHOD OF FEASIBLE DIRECTIONS.
!     BY G. N. VANDERPLAATS                             APRIL, 1972.
!     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF.
!     REF.  'STRUCTURAL OPTIMIZATION BY METHODS OF FEASIBLE DIRECTIONS',
!     G. N. VANDERPLAATS AND F. MOSES, JOURNAL OF COMPUTERS
!     AND STRUCTURES, VOL 3, PP 739-755, 1973.
!     FORM OF L. P. IS BX=C WHERE 1ST NDB COMPONENTS OF X CONTAIN VECTOR
!     U AND LAST NDB COMPONENTS CONTAIN VECTOR V.  CONSTRAINTS ARE
!     U.GE.0, V.GE.0, AND U-TRANSPOSE TIMES V = 0.
!     NER = ERROR FLAG.  IF NER.NE.0 ON RETURN, PROCESS HAS NOT
!     CONVERGED IN 5*NDB ITERATIONS.
!     VECTOR MS1 IDENTIFIES THE SET OF BASIC VARIABLES.
!     ------------------------------------------------------------------
!     CHOOSE INITIAL BASIC VARIABLES AS V, AND INITIALIZE VECTOR MS1
!     ------------------------------------------------------------------
   Ner = 1
   m2 = 2*Ndb
!     CALCULATE CBMIN AND EPS AND INITIALIZE MS1.
   eps = -1.0e+10
   cbmin = 0.
   do i = 1 , Ndb
      bi = b(i,i)
      cbmax = 0.
      if ( bi<-1.0e-6 ) cbmax = c(i)/bi
      if ( bi>eps ) eps = bi
      if ( cbmax>cbmin ) cbmin = cbmax
      Ms1(i) = 0
   enddo
   eps = .0001*eps
!     IF (EPS.LT.-1.0E-10) EPS=-1.0E-10
!
!  E-10 CHANGED TO E-03 ON 1/12/81
!
   if ( eps<-1.0e-03 ) eps = -1.0e-03
   if ( eps>-.0001 ) eps = -.0001
   cbmin = cbmin*1.0e-6
!     IF (CBMIN.LT.1.0E-10) CBMIN=1.0E-10
!
!  E-10 CHANGED TO E-05 ON 1/12/81
!
   if ( cbmin<1.0e-05 ) cbmin = 1.0e-05
   iter1 = 0
   nmax = 5*Ndb
   SPAG_Loop_1_1: do
!     ------------------------------------------------------------------
!     **********             BEGIN NEW ITERATION              **********
!     ------------------------------------------------------------------
      iter1 = iter1 + 1
      if ( iter1>nmax ) return
!     FIND MAX. C(I)/B(I,I) FOR I=1,NDB.
      cbmax = .9*cbmin
      ichk = 0
      do i = 1 , Ndb
         c1 = c(i)
         bi = b(i,i)
!     IF (BI.GT.EPS.OR.C1.GT.0.) GO TO 30
         if ( bi<=eps .and. c1<=-1.0e-05 ) then
!
!  0. CHANGED TO -1.0E-05 ON 1/12/81
!
            cb = c1/bi
            if ( cb>cbmax ) then
               ichk = i
               cbmax = cb
            endif
         endif
      enddo
      if ( cbmax<cbmin ) exit SPAG_Loop_1_1
      if ( ichk==0 ) exit SPAG_Loop_1_1
!     UPDATE VECTOR MS1.
      jj = ichk
      if ( Ms1(jj)==0 ) jj = ichk + Ndb
      kk = jj + Ndb
      if ( kk>m2 ) kk = jj - Ndb
      Ms1(kk) = ichk
      Ms1(jj) = 0
!     ------------------------------------------------------------------
!                     PIVOT OF B(ICHK,ICHK)
!     ------------------------------------------------------------------
      bb = 1./b(ichk,ichk)
      do j = 1 , Ndb
         b(ichk,j) = bb*b(ichk,j)
      enddo
      c(ichk) = cbmax
      b(ichk,ichk) = bb
!     ELIMINATE COEFICIENTS ON VARIABLE ENTERING BASIS AND STORE
!     COEFICIENTS ON VARIABLE LEAVING BASIS IN THEIR PLACE.
      do i = 1 , Ndb
         if ( i/=ichk ) then
            bb1 = b(i,ichk)
            b(i,ichk) = 0.
            do j = 1 , Ndb
               b(i,j) = b(i,j) - bb1*b(ichk,j)
            enddo
            c(i) = c(i) - bb1*cbmax
         endif
      enddo
   enddo SPAG_Loop_1_1
   Ner = 0
!     ------------------------------------------------------------------
!     STORE ONLY COMPONENTS OF U-VECTOR IN 'C'.  USE B(I,1) FOR
!     TEMPORARY STORAGE
!     ------------------------------------------------------------------
   do i = 1 , Ndb
      b(i,1) = c(i)
   enddo
   do i = 1 , Ndb
      c(i) = 0.
      j = Ms1(i)
      if ( j>0 ) c(i) = b(j,1)
      if ( c(i)<0. ) c(i) = 0.
   enddo
end subroutine cnmn08


subroutine cnmn09(cnmnfun,cnmngrd,x,g,Ic,Df,a,n1,n2,n3,n4,n5)

   implicit none

   double precision a , Abobj1 , Alphax , Aobj , Ct , Ctl , Ctlmin , Ctmin , Dabfun , Delfun , Df , Fdch , Fdchm , g , Obj ,       &
                  & Theta , x
   integer Ic , Icndir , Igoto , Info , Infog , Iprint , Iter , Itmax , Itrm , Linobj , n1 , n2 , n3 , n4 , n5 , Nac , Ncon , Ndv ,&
         & Nfdg , Nfun
   integer Ngrd , Nscal , Nside


   common /cnmn1 / Delfun , Dabfun , Fdch , Fdchm , Ct , Ctmin , Ctl , Ctlmin , Alphax , Abobj1 , Theta , Obj , Ndv , Ncon ,       &
                 & Nside , Iprint , Nfdg , Nscal , Linobj , Itmax , Itrm , Icndir , Igoto , Nac , Info , Infog , Iter

   common /varable/ Aobj
   common /fevals/ Nfun , Ngrd

   dimension x(n1) , g(n2) , Ic(n2) , Df(n1) , a(n1,n2)

   external cnmnfun , cnmngrd

!
!
   if ( Info>=2 ) then

!
!
!	GRADIENT INFORMATION
!
      call cnmngrd(n1,n2,x,Aobj,g,Ct,Df,a,Ic,Nac)
      Ngrd = Ngrd + 1
   else
!
!  OBJECTIVE FUNCTION & CONSTRAINTS
!
      call cnmnfun(n1,n2,x,Aobj,g)
      Nfun = Nfun + 1
   endif


end subroutine cnmn09


subroutine conmin(Ndv_,Ncon_,x_,Vlb_,Vub_,Obj_,g_,n1,n2,n3,n4,n5,Iprint_,Iout_,Ifile,Itmax_,Delfun_,Dabfun_,Itrm_,Nfeasct_,Nfdg_,  &
                & Nfun_,Ngrd_,Cnmnfun,Cnmngrd)

   implicit none

   double precision a , Abobj1 , Alphax , Aobj , b , c , Ct , Ctl , Ctlmin , Ctmin , Dabfun , Dabfun_ ,        &
                  & Delfun , Delfun_ , df , Fdch , Fdchm , g
   double precision g1 , g2 , g_ , Obj , Obj_ , phi , s , scal , Theta , vlb , Vlb_ , vub , Vub_ , x , x_
   integer i , ic , Icndir , Igoto , Info , Infog , Iout , Iout_ , Iprint , Iprint_ , isc , Iter , Itmax , Itmax_ , Itrm , Itrm_ , &
         & j , Linobj , loopcnt , ms1
   integer n1 , n2 , n3 , n4 , n5 , Nac , Ncon , Ncon_ , Ndv , Ndv_ , Nfdg , Nfdg_ , Nfeasct , Nfeasct_ , Nfun , Nfun_ , Ngrd ,    &
         & Ngrd_ , nlim , Nscal
   integer Nside


   external Cnmnfun , Cnmngrd

   character*(*) Ifile

   dimension x_(Ndv_) , Vlb_(Ndv_) , Vub_(Ndv_) , g_(Ncon_)

   dimension x(n1) , vlb(n1) , vub(n1) , scal(n1) , s(n1) , df(n1)
   dimension g(n2) , g1(n2) , g2(n2) , isc(n2)
   dimension ic(n3) , b(n3,n3)
   dimension c(n4)
   dimension ms1(n5)
   dimension a(n1,n3)

   common /cnmn1 / Delfun , Dabfun , Fdch , Fdchm , Ct , Ctmin , Ctl , Ctlmin , Alphax , Abobj1 , Theta , Obj , Ndv , Ncon ,       &
                 & Nside , Iprint , Nfdg , Nscal , Linobj , Itmax , Itrm , Icndir , Igoto , Nac , Info , Infog , Iter , Nfeasct

   common /varable/ Aobj
   common /output/ Iout
   common /fevals/ Nfun , Ngrd
!
!  INITIALIZE
!
   Infog = 0
   Info = 0
   Ndv = Ndv_
   Ncon = Ncon_
   do i = 1 , Ndv
      x(i) = x_(i)
      vlb(i) = Vlb_(i)
      vub(i) = Vub_(i)
   enddo
   Iprint = Iprint_
   Itmax = Itmax_
   do j = 1 , Ncon
      isc(j) = 0
   enddo
   Nfdg = Nfdg_
   Nside = 1
   Icndir = 0
   Nscal = 0
   Nfeasct = Nfeasct_
   Fdch = 0.0
   Fdchm = 0.0
   Ct = 0.0
   Ctmin = 0.0
   Ctl = 0.0
   Ctlmin = 0.0
   Theta = 0.0
   phi = 0.0
   Delfun = Delfun_
   Dabfun = Dabfun_
   Linobj = 0.0
   Itrm = Itrm_
   Alphax = 0.0
   Abobj1 = 0.0
!
   Nfun = 0
   Ngrd = 0
!
!     OPEN WRITE FILE
!
   Iout = Iout_
   if ( Iprint/=0 ) open (unit=Iout,file=Ifile(1:len_trim(Ifile)),status='UNKNOWN')
!
!     MAXIMUM NUMBER OF ITERATIONS
!
   nlim = Itmax*(Ndv+5)
!
!     NON-ITERATIVE PART OF ANALYSIS
!
   Igoto = 0
!
!     ITERATIVE PART OF ANALYSIS
!
   SPAG_Loop_1_1: do i = 1 , nlim
!
      loopcnt = i
!
!       CALL THE OPTIMIZATION ROUTINE CONMIN
!
      call cnmn00(x,vlb,vub,g,scal,df,a,s,g1,g2,b,c,isc,ic,ms1,n1,n2,n3,n4,n5)
!
!     CHECK TERMINATION CRITERIA
!
      if ( Igoto==0 ) loopcnt = -999
!
!       ANALYSIS MODULE
!
      call cnmn09(Cnmnfun,Cnmngrd,x,g,ic,df,a,n1,n2,n3,n4,n5)
      Obj = Aobj
!
      if ( Igoto==0 ) exit SPAG_Loop_1_1
   enddo SPAG_Loop_1_1
!
!
!  PRINT FINAL RESULTS
!
   if ( Iprint/=0 ) then
      write (6,99001) Nfun - 1

!  ------------------------------------------------------------------
!  FORMATS
!  ------------------------------------------------------------------
99001 format (//8x,'NUMBER OF FUNC-CALLS:  NFUN =',i5)
      write (6,99002) Ngrd
99002 format (8x,'NUMBER OF GRAD-CALLS:  NGRD =',i5)
   endif
!
!
!     OUTPUT HANDLING
!
   do i = 1 , Ndv
      x_(i) = x(i)
   enddo
   Obj_ = Obj
   do j = 1 , Ncon
      g_(j) = g(j)
   enddo
   Nfun_ = Nfun - 1
   Ngrd_ = Ngrd

   return

end subroutine conmin

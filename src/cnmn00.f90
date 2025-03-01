!*==cnmn00.f90 processed by SPAG 8.01MH 10:24  1 Mar 2025
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE cnmn00(X,Vlb,Vub,G,Scal,Df,A,S,G1,G2,B,C,Isc,Ic,Ms1,N1,N2,N3,N4,N5)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   DOUBLE PRECISION A , A1 , A2 , A3 , A4 , Abobj , Abobj1 , Alp , alp1 , alp11 , alp12 , Alpca , Alpfes , Alphax , Alpln ,        &
                  & Alpmin , Alpnc , Alpsav , Alpsid , Alptot
   DOUBLE PRECISION App , B , C , c1 , Ct , ct1 , Cta , Ctam , Ctbm , ctc , Ctl , Ctlmin , Ctmin , Cv1 , Cv2 , Cv3 , Cv4 , Dabfun ,&
                  & Dct , Dctl
   DOUBLE PRECISION Delfun , Df , Dftdf1 , Dm1 , Dm10 , Dm11 , Dm12 , Dm2 , Dm3 , Dm4 , Dm5 , Dm6 , Dm7 , Dm8 , Dm9 , Dx , Dx1 ,   &
                  & F1 , F2 , F3
   DOUBLE PRECISION F4 , Fdch , Fdchm , ff1 , Fff , Fi , G , G1 , G2 , gi , Obj , Obj1 , objb , objd , Phi , Rspace , S , Scal ,   &
                  & scj , si
   DOUBLE PRECISION sib , Slope , Theta , Vlb , Vub , X , x1 , x12 , Xi , xid , xx
   INTEGER i , Ibest , Ic , Icndir , Icount , Idm1 , Idm2 , Idm3 , Igood1 , Igood2 , Igood3 , Igood4 , Igoto , ii , Iii , Info ,   &
         & Infog , Iobj , Iout , Iprint
   INTEGER Isc , Ispace , Iter , Itmax , Itrm , j , Jdir , Jgoto , k , Kcount , Kobj , Kount , Linobj , m1 , m2 , m3 , mcn1 , Ms1 ,&
         & Mscal , N1
   INTEGER N2 , N3 , N4 , N5 , Nac , Ncal , nci , Ncobj , Ncon , Ndv , ndv1 , ndv2 , Nfdg , Nfeas , Nfeasct , nic , Nlnc , nnac ,  &
         & Nscal , Nside
   INTEGER Nvc
!*** End of declarations inserted by SPAG
   COMMON /cnmn1 / Delfun , Dabfun , Fdch , Fdchm , Ct , Ctmin , Ctl , Ctlmin , Alphax , Abobj1 , Theta , Obj , Ndv , Ncon ,       &
                 & Nside , Iprint , Nfdg , Nscal , Linobj , Itmax , Itrm , Icndir , Igoto , Nac , Info , Infog , Iter , Nfeasct
   COMMON /output/ Iout
!
!    NFEASCT ADDED TO COMMON BLOCK BY KCYOUNG ON 4/14/92 TO ALLOW MORE
!    THAN 10 ITERATION ATTEMPTS.  NFEASCT BECOMES AN INPUT VALUE
!
   DIMENSION X(N1) , Vlb(N1) , Vub(N1) , G(N2) , Scal(N1) , Df(N1) , A(N1,N3) , S(N1) , G1(N2) , G2(N2) , B(N3,N3) , C(N4) ,       &
           & Isc(N2) , Ic(N3) , Ms1(N5)
   COMMON /consav/ Dm1 , Dm2 , Dm3 , Dm4 , Dm5 , Dm6 , Dm7 , Dm8 , Dm9 , Dm10 , Dm11 , Dm12 , Dct , Dctl , Phi , Abobj , Cta ,     &
                 & Ctam , Ctbm , Obj1 , Slope , Dx , Dx1 , Fi , Xi , Dftdf1 , Alp , Fff , A1 , A2 , A3 , A4 , F1 , F2 , F3 , F4 ,  &
                 & Cv1 , Cv2 , Cv3 , Cv4 , App , Alpca , Alpfes , Alpln , Alpmin , Alpnc , Alpsav , Alpsid , Alptot , Rspace ,     &
                 & Idm1 , Idm2 , Idm3 , Jdir , Iobj , Kobj , Kcount , Ncal(2) , Nfeas , Mscal , Ncobj , Nvc , Kount , Icount ,     &
                 & Igood1 , Igood2 , Igood3 , Igood4 , Ibest , Iii , Nlnc , Jgoto , Ispace(2)
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
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
         IF ( Nscal/=0 .AND. Igoto/=0 ) THEN
            DO i = 1 , Ndv
               X(i) = C(i)
            ENDDO
         ENDIF
!     CONSTANTS.
         ndv1 = Ndv + 1
         ndv2 = Ndv + 2
         IF ( Igoto/=0 ) THEN
!     ------------------------------------------------------------------
!                     CHECK FOR UNBOUNDED SOLUTION
!     ------------------------------------------------------------------
!     STOP IF OBJ IS LESS THAN -1.0E+20
            IF ( Obj<=-1.0E+20 ) THEN
               WRITE (Iout,99002)
99002          FORMAT (///5X,'CONMIN HAS ACHIEVED A SOLUTION OF OBJ LESS THAN -1.0E+40'/5X,'SOLUTION APPEARS TO BE UNBOUNDED'/5X,  &
                      &'OPTIMIZATION IS TERMINATED')
               spag_nextblock_1 = 10
               CYCLE SPAG_DispatchLoop_1
            ELSEIF ( Igoto==1 ) THEN
               Obj1 = Obj
               IF ( Dabfun<=0. ) Dabfun = .001*abs(Obj)
               IF ( Dabfun<1.0E-10 ) Dabfun = 1.0E-10
               IF ( Iprint>0 ) THEN
!     ------------------------------------------------------------------
!                    PRINT INITIAL DESIGN INFORMATION
!     ------------------------------------------------------------------
                  IF ( Iprint>1 ) THEN
                     IF ( Nside==0 .AND. Ncon==0 ) WRITE (Iout,99033)
99033                FORMAT (////5X,'UNCONSTRAINED FUNCTION MINIMIZATION'//5X,'CONTROL PARAMETERS')
                     IF ( Nside/=0 .OR. Ncon>0 ) WRITE (Iout,99027)
99027                FORMAT (////5X,'CONSTRAINED FUNCTION MINIMIZATION'//5X,'CONTROL PARAMETERS')
                     WRITE (Iout,99028) Iprint , Ndv , Itmax , Ncon , Nside , Icndir , Nscal , Nfdg , Linobj , Itrm , N1 , N2 ,    &
                                      & N3 , N4 , N5
99028                FORMAT (/5X,'IPRINT  NDV    ITMAX    NCON    NSIDE  ICNDIR   NSCAL   NFDG'/8I8//5X,'LINOBJ  ITRM',5X,'N1',6X, &
                            &'N2',6X,'N3',6X,'N4',6X,'N5'/8I8)
                     WRITE (Iout,99030) Ct , Ctmin , Ctl , Ctlmin , Theta , Phi , Delfun , Dabfun
99030                FORMAT (/9X,'CT',14X,'CTMIN',11X,'CTL',13X,'CTLMIN'/1X,4(2X,E14.5)//9X,'THETA',11X,'PHI',13X,'DELFUN',10X,    &
                            &'DABFUN'/1X,4(2X,E14.5))
                     WRITE (Iout,99029) Fdch , Fdchm , Alphax , Abobj1
99029                FORMAT (/9X,'FDCH',12X,'FDCHM',11X,'ALPHAX',10X,'ABOBJ1'/1X,4(2X,E14.5))
                     IF ( Nside/=0 ) THEN
                        WRITE (Iout,99031)
99031                   FORMAT (/5X,'LOWER BOUNDS ON DECISION VARIABLES (VLB)')
                        DO i = 1 , Ndv , 6
                           m1 = min0(Ndv,i+5)
                           WRITE (Iout,99005) i , (Vlb(j),j=i,m1)
                        ENDDO
                        WRITE (Iout,99032)
99032                   FORMAT (/5X,'UPPER BOUNDS ON DECISION VARIABLES (VUB)')
                        DO i = 1 , Ndv , 6
                           m1 = min0(Ndv,i+5)
                           WRITE (Iout,99005) i , (Vub(j),j=i,m1)
                        ENDDO
                     ENDIF
                     IF ( Nscal<0 ) THEN
                        WRITE (Iout,99034)
99034                   FORMAT (/5X,'SCALING VECTOR (SCAL)')
                        WRITE (Iout,99050) (Scal(i),i=1,Ndv)
                     ENDIF
                     IF ( Ncon/=0 ) THEN
                        IF ( Nlnc==0 .OR. Nlnc==Ncon ) THEN
                           IF ( Nlnc==Ncon ) WRITE (Iout,99008)
99008                      FORMAT (/5X,'ALL CONSTRAINTS ARE LINEAR')
                           IF ( Nlnc==0 ) WRITE (Iout,99009)
99009                      FORMAT (/5X,'ALL CONSTRAINTS ARE NON-LINEAR')
                        ELSE
                           WRITE (Iout,99006)
99006                      FORMAT (/5X,'LINEAR CONSTRAINT IDENTIFIERS (ISC)'/5X,'NON-ZERO INDICATES LINEAR CONSTRAINT')
                           DO i = 1 , Ncon , 15
                              m1 = min0(Ncon,i+14)
                              WRITE (Iout,99007) i , (Isc(j),j=i,m1)
99007                         FORMAT (3X,I5,')',2X,15I5)
                           ENDDO
                        ENDIF
                     ENDIF
                  ENDIF
                  WRITE (Iout,99048) Obj
99048             FORMAT (//5X,'INITIAL FUNCTION INFORMATION'//5X,'OBJ =',E15.6)
                  WRITE (Iout,99049)
                  DO i = 1 , Ndv
                     x1 = 1.
                     IF ( Nscal/=0 ) x1 = Scal(i)
                     G1(i) = X(i)*x1
                  ENDDO
                  DO i = 1 , Ndv , 6
                     m1 = min0(Ndv,i+5)
                     WRITE (Iout,99005) i , (G1(j),j=i,m1)
                  ENDDO
                  IF ( Ncon/=0 ) THEN
                     WRITE (Iout,99051)
                     DO i = 1 , Ncon , 6
                        m1 = min0(Ncon,i+5)
                        WRITE (Iout,99005) i , (G(j),j=i,m1)
                     ENDDO
                  ENDIF
               ENDIF
               IF ( Iprint>1 ) WRITE (Iout,99040)
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
            ELSEIF ( Igoto==2 ) THEN
               spag_nextblock_1 = 4
               CYCLE SPAG_DispatchLoop_1
            ELSEIF ( Igoto==3 ) THEN
               spag_nextblock_1 = 3
               CYCLE SPAG_DispatchLoop_1
            ELSEIF ( Igoto==4 ) THEN
               spag_nextblock_1 = 7
               CYCLE SPAG_DispatchLoop_1
            ELSEIF ( Igoto==5 ) THEN
               spag_nextblock_1 = 8
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDIF
!     ------------------------------------------------------------------
!                      SAVE INPUT CONTROL PARAMETERS
!     ------------------------------------------------------------------
         IF ( Iprint>0 ) WRITE (Iout,99026)
99026    FORMAT ('1',////12X,27('* ')/12X,'*',51X,'*'/12X,'*',20X,'C O N M I N',20X,'*'/12X,'*',51X,'*'/12X,'*',15X,               &
                &' FORTRAN PROGRAM FOR ',15X,'*'/12X,'*',51X,'*'/12X,'*',9X,'CONSTRAINED FUNCTION MINIMIZATION',9X,'*'/12X,'*',51X,&
                &'*'/12X,27('* '))
         IF ( Linobj==0 .OR. (Ncon>0 .OR. Nside>0) ) THEN
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
            IF ( Itrm<=0 ) Itrm = 3
            IF ( Itmax<=0 ) Itmax = 20
            ndv1 = Ndv + 1
            IF ( Icndir==0 ) Icndir = ndv1
            IF ( Delfun<=0. ) Delfun = .0001
            Ct = -abs(Ct)
            IF ( Ct>=0. ) Ct = -.1
            Ctmin = abs(Ctmin)
            IF ( Ctmin<=0. ) Ctmin = .004
            Ctl = -abs(Ctl)
            IF ( Ctl>=0. ) Ctl = -0.01
            Ctlmin = abs(Ctlmin)
            IF ( Ctlmin<=0. ) Ctlmin = .001
            IF ( Theta<=0. ) Theta = 1.
            IF ( Abobj1<=0. ) Abobj1 = .1
            IF ( Alphax<=0. ) Alphax = .1
            IF ( Fdch<=0. ) Fdch = .01
            IF ( Fdchm<=0. ) Fdchm = .01
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
            IF ( Ncon/=0 ) THEN
               DO i = 1 , Ncon
                  IF ( Isc(i)>0 ) Nlnc = Nlnc + 1
               ENDDO
            ENDIF
!     ------------------------------------------------------------------
!          CHECK TO BE SURE THAT SIDE CONSTRAINTS ARE SATISFIED
!     ------------------------------------------------------------------
            IF ( Nside/=0 ) THEN
               DO i = 1 , Ndv
                  IF ( Vlb(i)>Vub(i) ) THEN
                     xx = .5*(Vlb(i)+Vub(i))
                     X(i) = xx
                     Vlb(i) = xx
                     Vub(i) = xx
                     WRITE (Iout,99016) i
99016                FORMAT (///5X,'* * CONMIN DETECTS VLB(I).GT.VUB(I)'/5X,                                                       &
                            &'FIX IS SET X(I)=VLB(I)=VUB(I) = .5*(VLB(I)+VUB(I) FOR I =',I5)
                  ENDIF
                  xx = X(i) - Vlb(i)
                  IF ( xx>=0. ) THEN
                     xx = Vub(i) - X(i)
                     IF ( xx<0. ) THEN
                        WRITE (Iout,99018) X(i) , Vub(i) , i
99018                   FORMAT (///5X,'* * CONMIN DETECTS INITIAL X(I).GT.VUB(I)'/5X,'X(I) =',E12.4,2X,'VUB(I) =',E12.4/5X,        &
                               &'X(I) IS SET EQUAL TO VUB(I) FOR I =',I5)
                        X(i) = Vub(i)
                     ENDIF
                  ELSE
!     LOWER BOUND VIOLATED.
                     WRITE (Iout,99017) X(i) , Vlb(i) , i
99017                FORMAT (///5X,'* * CONMIN DETECTS INITIAL X(I).LT.VLB(I)'/5X,'X(I) =',E12.4,2X,'VLB(I) =',E12.4/5X,           &
                            &'X(I) IS SET EQUAL TO VLB(I) FOR I =',I5)
                     X(i) = Vlb(i)
                  ENDIF
               ENDDO
            ENDIF
!     ------------------------------------------------------------------
!                        INITIALIZE SCALING VECTOR, SCAL
!     ------------------------------------------------------------------
            IF ( Nscal/=0 ) THEN
               IF ( Nscal<0 ) THEN
                  DO i = 1 , Ndv
                     si = abs(Scal(i))
                     IF ( si<1.0E-20 ) si = 1.0E-5
                     Scal(i) = si
                     si = 1./si
                     X(i) = X(i)*si
                     IF ( Nside/=0 ) THEN
                        Vlb(i) = Vlb(i)*si
                        Vub(i) = Vub(i)*si
                     ENDIF
                  ENDDO
               ELSE
                  DO i = 1 , Ndv
                     Scal(i) = 1.
                  ENDDO
               ENDIF
            ENDIF
!     ------------------------------------------------------------------
!     ***** CALCULATE INITIAL FUNCTION AND CONSTRAINT VALUES  *****
!     ------------------------------------------------------------------
            Info = 1
            Ncal(1) = 1
            Igoto = 1
            spag_nextblock_1 = 11
            CYCLE SPAG_DispatchLoop_1
         ELSE
!     TOTALLY UNCONSTRAINED FUNCTION WITH LINEAR OBJECTIVE.
!     SOLUTION IS UNBOUNDED.
            WRITE (Iout,99001) Linobj , Ncon , Nside
!     ------------------------------------------------------------------
!                                FORMATS
!     ------------------------------------------------------------------
!
!
99001       FORMAT (///5X,'A COMPLETELY UNCONSTRAINED FUNCTION WITH A LINEAR OBJECTIVE IS SPECIFIED'//10X,'LINOBJ =',I5/10X,       &
                   &'NCON   =',I5/10X,'NSIDE  =',I5//5X,'CONTROL RETURNED TO CALLING PROGRAM')
            RETURN
         ENDIF
      CASE (2)
!     ------------------------------------------------------------------
!     ********************  BEGIN MINIMIZATION  ************************
!     ------------------------------------------------------------------
         Iter = Iter + 1
         IF ( Abobj1<.0001 ) Abobj1 = .0001
         IF ( Abobj1>.2 ) Abobj1 = .2
         IF ( Alphax>1. ) Alphax = 1.
         IF ( Alphax<.001 ) Alphax = .001
!
!  THE FOLLOWING TWO LINES OF CODE WERE COMMENTED OUT ON 3/5/81
!
!     NFEAS=NFEAS+1
!     IF (NFEAS.GT.10) GO TO 810
         IF ( Iprint>2 ) WRITE (Iout,99035) Iter
99035    FORMAT (////5X,'BEGIN ITERATION NUMBER',I5)
         IF ( Iprint>3 .AND. Ncon>0 ) WRITE (Iout,99036) Ct , Ctl , Phi
99036    FORMAT (/5X,'CT =',E14.5,5X,'CTL =',E14.5,5X,'PHI =',E14.5)
         Cta = abs(Ct)
         IF ( Ncobj==0 ) THEN
            IF ( Mscal>=Nscal .AND. Nscal/=0 ) THEN
               IF ( Nscal>=0 .OR. Kcount>=Icndir ) THEN
                  Mscal = 0
                  Kcount = 0
!     ------------------------------------------------------------------
!                          SCALE VARIABLES
!     ------------------------------------------------------------------
                  DO i = 1 , Ndv
                     si = Scal(i)
                     Xi = si*X(i)
                     sib = si
                     IF ( Nscal>0 ) si = abs(Xi)
                     IF ( si>=1.0E-10 ) THEN
                        Scal(i) = si
                        si = 1./si
                        X(i) = Xi*si
                        IF ( Nside/=0 ) THEN
                           Vlb(i) = sib*si*Vlb(i)
                           Vub(i) = sib*si*Vub(i)
                        ENDIF
                     ENDIF
                  ENDDO
                  IF ( .NOT.(Iprint<4 .OR. (Nscal<0 .AND. Iter>1)) ) THEN
                     WRITE (Iout,99037)
99037                FORMAT (/5X,'NEW SCALING VECTOR (SCAL)')
                     WRITE (Iout,99050) (Scal(i),i=1,Ndv)
                  ENDIF
               ENDIF
            ENDIF
            Mscal = Mscal + 1
            Nac = 0
!     ------------------------------------------------------------------
!          OBTAIN GRADIENTS OF OBJECTIVE AND ACTIVE CONSTRAINTS
!     ------------------------------------------------------------------
            Info = 2
            Ncal(2) = Ncal(2) + 1
            IF ( Nfdg/=1 ) THEN
               Jgoto = 0
            ELSE
               Igoto = 2
               spag_nextblock_1 = 11
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ELSE
!     ------------------------------------------------------------------
!     NO MOVE ON LAST ITERATION.  DELETE CONSTRAINTS THAT ARE NO
!     LONGER ACTIVE.
!     ------------------------------------------------------------------
            nnac = Nac
            DO i = 1 , nnac
               IF ( Ic(i)>Ncon ) Nac = Nac - 1
            ENDDO
            IF ( Nac>0 ) THEN
               nnac = Nac
               SPAG_Loop_1_2: DO i = 1 , nnac
                  SPAG_Loop_2_1: DO
                     nic = Ic(i)
                     ct1 = Ct
                     IF ( Isc(nic)>0 ) ct1 = Ctl
                     IF ( G(nic)>ct1 ) EXIT SPAG_Loop_2_1
                     Nac = Nac - 1
                     IF ( i>Nac ) EXIT SPAG_Loop_1_2
                     DO k = i , Nac
                        ii = k + 1
                        DO j = 1 , ndv2
                           A(j,k) = A(j,ii)
                        ENDDO
                        Ic(k) = Ic(ii)
                     ENDDO
                  ENDDO SPAG_Loop_2_1
               ENDDO SPAG_Loop_1_2
            ENDIF
            spag_nextblock_1 = 5
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         spag_nextblock_1 = 3
      CASE (3)
         CALL cnmn01(Jgoto,X,Df,G,Isc,Ic,A,G1,Vlb,Vub,Scal,C,Ncal,Dx,Dx1,Fi,Xi,Iii,N1,N2,N3,N4)
         Igoto = 3
         IF ( Jgoto>0 ) THEN
            spag_nextblock_1 = 11
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         spag_nextblock_1 = 4
      CASE (4)
         Info = 1
         IF ( Nac>=N3 ) THEN
            spag_nextblock_1 = 10
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         IF ( Nscal/=0 .AND. Nfdg/=0 ) THEN
!     ------------------------------------------------------------------
!                              SCALE GRADIENTS
!     ------------------------------------------------------------------
!     SCALE GRADIENT OF OBJECTIVE FUNCTION.
            DO i = 1 , Ndv
               Df(i) = Df(i)*Scal(i)
            ENDDO
            IF ( Nfdg/=2 .AND. Nac/=0 ) THEN
!     SCALE GRADIENTS OF ACTIVE CONSTRAINTS.
               DO j = 1 , Ndv
                  scj = Scal(j)
                  DO i = 1 , Nac
                     A(j,i) = A(j,i)*scj
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
         spag_nextblock_1 = 5
      CASE (5)
         IF ( Iprint>=3 .AND. Ncon/=0 ) THEN
!     ------------------------------------------------------------------
!                                   PRINT
!     ------------------------------------------------------------------
!     PRINT ACTIVE AND VIOLATED CONSTRAINT NUMBERS.
            m1 = 0
            m2 = N3
            IF ( Nac/=0 ) THEN
               DO i = 1 , Nac
                  j = Ic(i)
                  IF ( j<=Ncon ) THEN
                     gi = G(j)
                     c1 = Ctam
                     IF ( Isc(j)>0 ) c1 = Ctbm
                     gi = gi - c1
                     IF ( gi>0. ) THEN
                        m2 = m2 + 1
!     VIOLATED CONSTRAINT.
                        Ms1(m2) = j
                     ELSE
!     ACTIVE CONSTRAINT.
                        m1 = m1 + 1
                        Ms1(m1) = j
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
            m3 = m2 - N3
            WRITE (Iout,99010) m1
            IF ( m1/=0 ) THEN
               WRITE (Iout,99011)
               WRITE (Iout,99052) (Ms1(i),i=1,m1)
            ENDIF
            WRITE (Iout,99012) m3
            IF ( m3/=0 ) THEN
               WRITE (Iout,99011)
               m3 = N3 + 1
               WRITE (Iout,99052) (Ms1(i),i=m3,m2)
            ENDIF
         ENDIF
!     ------------------------------------------------------------------
!            CALCULATE GRADIENTS OF ACTIVE SIDE CONSTRAINTS
!     ------------------------------------------------------------------
         IF ( Nside/=0 ) THEN
            mcn1 = Ncon
            m1 = 0
            DO i = 1 , Ndv
!     LOWER BOUND.
               Xi = X(i)
               xid = Vlb(i)
               x12 = abs(xid)
               IF ( x12<1. ) x12 = 1.
               gi = (xid-Xi)/x12
               IF ( gi>=-1.0E-6 ) THEN
                  m1 = m1 + 1
                  Ms1(m1) = -i
                  Nac = Nac + 1
                  IF ( Nac>=N3 ) THEN
                     spag_nextblock_1 = 10
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
                  mcn1 = mcn1 + 1
                  DO j = 1 , Ndv
                     A(j,Nac) = 0.
                  ENDDO
                  A(i,Nac) = -1.
                  Ic(Nac) = mcn1
                  G(mcn1) = gi
                  Isc(mcn1) = 1
               ENDIF
!     UPPER BOUND.
               xid = Vub(i)
               x12 = abs(xid)
               IF ( x12<1. ) x12 = 1.
               gi = (Xi-xid)/x12
               IF ( gi>=-1.0E-6 ) THEN
                  m1 = m1 + 1
                  Ms1(m1) = i
                  Nac = Nac + 1
                  IF ( Nac>=N3 ) THEN
                     spag_nextblock_1 = 10
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
                  mcn1 = mcn1 + 1
                  DO j = 1 , Ndv
                     A(j,Nac) = 0.
                  ENDDO
                  A(i,Nac) = 1.
                  Ic(Nac) = mcn1
                  G(mcn1) = gi
                  Isc(mcn1) = 1
               ENDIF
            ENDDO
!     ------------------------------------------------------------------
!                                  PRINT
!     ------------------------------------------------------------------
!     PRINT ACTIVE SIDE CONSTRAINT NUMBERS.
            IF ( Iprint>=3 ) THEN
               WRITE (Iout,99013) m1
               IF ( m1/=0 ) THEN
                  WRITE (Iout,99014)
                  WRITE (6,99052) (Ms1(j),j=1,m1)
               ENDIF
            ENDIF
         ENDIF
!     PRINT GRADIENTS OF ACTIVE AND VIOLATED CONSTRAINTS.
         IF ( Iprint>=4 ) THEN
            WRITE (Iout,99038)
99038       FORMAT (/5X,'GRADIENT OF OBJ')
            DO i = 1 , Ndv , 6
               m1 = min0(Ndv,i+5)
               WRITE (Iout,99005) i , (Df(j),j=i,m1)
            ENDDO
            IF ( Nac/=0 ) THEN
               WRITE (Iout,99039)
99039          FORMAT (/5X,'GRADIENTS OF ACTIVE AND VIOLATED CONSTRAINTS')
               DO i = 1 , Nac
                  m1 = Ic(i)
                  m2 = m1 - Ncon
                  m3 = 0
                  IF ( m2>0 ) m3 = iabs(Ms1(m2))
                  IF ( m2<=0 ) WRITE (Iout,99003) m1
99003             FORMAT (5X,'CONSTRAINT NUMBER',I5)
                  IF ( m2>0 ) WRITE (Iout,99004) m3
99004             FORMAT (5X,'SIDE CONSTRAINT ON VARIABLE',I5)
                  DO k = 1 , Ndv , 6
                     m1 = min0(Ndv,k+5)
                     WRITE (Iout,99005) k , (A(j,i),j=k,m1)
                  ENDDO
                  WRITE (Iout,99040)
               ENDDO
            ENDIF
         ENDIF
!     ------------------------------------------------------------------
!     ******************  DETERMINE SEARCH DIRECTION *******************
!     ------------------------------------------------------------------
         Alp = 1.0E+20
         IF ( Nac>0 ) THEN
!     ------------------------------------------------------------------
!                          CONSTRAINED FUNCTION
!     ------------------------------------------------------------------
!     FIND USABLE-FEASIBLE DIRECTION.
            Kcount = 0
            Jdir = 0
            Phi = 10.*Phi
            IF ( Phi>1000. ) Phi = 1000.
!
!  THE FOLLOWING LINE OF CODE WAS COMMENTED OUT ON 3/5/81
!
!     IF (NFEAS.EQ.1) PHI=5.
!     CALCULATE DIRECTION, S.
            CALL cnmn05(G,Df,A,S,B,C,Slope,Phi,Isc,Ic,Ms1,Nvc,N1,N2,N3,N4,N5)
!
!  THE FOLLOWING LINE WAS ADDED ON 2/25/81
!
            IF ( Nac/=0 ) THEN
!
!  THE FOLLOWING FIVE LINES WERE COMMENTED OUT ON 3/5/81
!  REASON : THEY WERE NOT IN G. VANDERPLAATS LISTING
!
!     IF THIS DESIGN IS FEASIBLE AND LAST ITERATION WAS INFEASIBLE,
!     SET ABOBJ1=.05 (5 PERCENT).
!     IF (NVC.EQ.0.AND.NFEAS.GT.1) ABOBJ1=.05
!     IF (NVC.EQ.0) NFEAS=0
               IF ( Iprint>=3 ) THEN
                  WRITE (Iout,99041)
99041             FORMAT (/5X,'PUSH-OFF FACTORS, (THETA(I), I=1,NAC)')
                  DO i = 1 , Nac , 6
                     m1 = min0(Nac,i+5)
                     WRITE (Iout,99005) i , (A(ndv1,j),j=i,m1)
                  ENDDO
                  WRITE (Iout,99025) S(ndv1)
99025             FORMAT (/5X,'CONSTRAINT PARAMETER, BETA =',E14.5)
               ENDIF
!     ------------------------------------------------------------------
!     ****************** ONE-DIMENSIONAL SEARCH ************************
!     ------------------------------------------------------------------
               IF ( S(ndv1)<1.0E-6 .AND. Nvc==0 ) THEN
                  spag_nextblock_1 = 9
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
!     ------------------------------------------------------------------
!                 FIND ALPHA TO OBTAIN A FEASIBLE DESIGN
!     ------------------------------------------------------------------
               IF ( Nvc/=0 ) THEN
                  Alp = -1.
                  DO i = 1 , Nac
                     nci = Ic(i)
                     c1 = G(nci)
                     ctc = Ctam
                     IF ( Isc(nci)>0 ) ctc = Ctbm
                     IF ( c1>ctc ) THEN
                        alp1 = 0.
                        DO j = 1 , Ndv
                           alp1 = alp1 + S(j)*A(j,i)
                        ENDDO
                        alp1 = alp1*A(ndv2,i)
                        IF ( abs(alp1)>=1.0E-20 ) THEN
                           alp1 = -c1/alp1
                           IF ( alp1>Alp ) Alp = alp1
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
               spag_nextblock_1 = 6
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDIF
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
         IF ( Kcount>Icndir .OR. Iobj==2 ) Kcount = 1
         IF ( Kcount==1 ) Jdir = 0
!     IF JDIR = 0 FIND DIRECTION OF STEEPEST DESCENT.
         CALL cnmn02(Jdir,Slope,Dftdf1,Df,S,N1)
         spag_nextblock_1 = 6
      CASE (6)
!     ------------------------------------------------------------------
!                       LIMIT CHANCE TO ABOBJ1*OBJ
!     ------------------------------------------------------------------
         alp1 = 1.0E+20
         si = abs(Obj)
         IF ( si<.01 ) si = .01
         IF ( abs(Slope)>1.0E-20 ) alp1 = Abobj1*si/Slope
         alp1 = abs(alp1)
         IF ( Nvc>0 ) alp1 = 10.*alp1
         IF ( alp1<Alp ) Alp = alp1
!     ------------------------------------------------------------------
!                   LIMIT CHANGE IN VARIABLE TO ALPHAX
!     ------------------------------------------------------------------
         alp11 = 1.0E+20
         DO i = 1 , Ndv
            si = abs(S(i))
            Xi = abs(X(i))
            IF ( si>=1.0E-10 .AND. Xi>=0.1 ) THEN
               alp1 = Alphax*Xi/si
               IF ( alp1<alp11 ) alp11 = alp1
            ENDIF
         ENDDO
         IF ( Nvc>0 ) alp11 = 10.*alp11
         IF ( alp11<Alp ) Alp = alp11
         IF ( Alp>1.0E+20 ) Alp = 1.0E+20
         IF ( Alp<=1.0E-20 ) Alp = 1.0E-20
         IF ( Iprint>=3 ) THEN
            WRITE (Iout,99042)
99042       FORMAT (/5X,'SEARCH DIRECTION (S-VECTOR)')
            DO i = 1 , Ndv , 6
               m1 = min0(Ndv,i+5)
               WRITE (Iout,99005) i , (S(j),j=i,m1)
            ENDDO
            WRITE (Iout,99015) Slope , Alp
99015       FORMAT (/5X,'ONE-DIMENSIONAL SEARCH'/5X,'INITIAL SLOPE =',E12.4,2X,'PROPOSED ALPHA =',E12.4)
         ENDIF
         IF ( Ncon>0 .OR. Nside>0 ) THEN
!     ------------------------------------------------------------------
!       SOLVE ONE-DIMENSIONAL SEARCH PROBLEM FOR CONSTRAINED FUNCTION
!     ------------------------------------------------------------------
            Jgoto = 0
            spag_nextblock_1 = 8
            CYCLE SPAG_DispatchLoop_1
         ELSE
!     ------------------------------------------------------------------
!           DO ONE-DIMENSIONAL SEARCH FOR UNCONSTRAINED FUNCTION
!     ------------------------------------------------------------------
            Jgoto = 0
         ENDIF
         spag_nextblock_1 = 7
      CASE (7)
         CALL cnmn03(X,S,Slope,Alp,Fff,A1,A2,A3,A4,F1,F2,F3,F4,App,N1,Ncal,Kount,Jgoto)
         Igoto = 4
         IF ( Jgoto>0 ) THEN
            spag_nextblock_1 = 11
            CYCLE SPAG_DispatchLoop_1
         ENDIF
!     PROCEED TO CONVERGENCE CHECK.
         Jdir = 1
         spag_nextblock_1 = 9
         CYCLE SPAG_DispatchLoop_1
      CASE (8)
         CALL cnmn06(X,Vlb,Vub,G,Scal,Df,S,G1,G2,Ctam,Ctbm,Slope,Alp,A2,A3,A4,F1,F2,F3,Cv1,Cv2,Cv3,Cv4,Alpca,Alpfes,Alpln,Alpmin,  &
                   & Alpnc,Alpsav,Alpsid,Alptot,Isc,N1,N2,Ncal,Nvc,Icount,Igood1,Igood2,Igood3,Igood4,Ibest,Iii,Nlnc,Jgoto)
         Igoto = 5
         IF ( Jgoto>0 ) THEN
            spag_nextblock_1 = 11
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         IF ( Nac==0 ) Jdir = 1
         spag_nextblock_1 = 9
      CASE (9)
!     ------------------------------------------------------------------
!     *******************     UPDATE ALPHAX   **************************
!     ------------------------------------------------------------------
         IF ( Alp>1.0E+19 ) Alp = 0.
!     UPDATE ALPHAX TO BE AVERAGE OF MAXIMUM CHANGE IN X(I)
!     AND ALHPAX.
         alp11 = 0.
         DO i = 1 , Ndv
            si = abs(S(i))
            Xi = abs(X(i))
            IF ( Xi>=1.0E-10 ) THEN
               alp1 = Alp*si/Xi
               IF ( alp1>alp11 ) alp11 = alp1
            ENDIF
         ENDDO
         alp11 = .5*(alp11+Alphax)
         alp12 = 5.*Alphax
         IF ( alp11>alp12 ) alp11 = alp12
         Alphax = alp11
         Ncobj = Ncobj + 1
!     ABSOLUTE CHANGE IN OBJECTIVE.
         objd = Obj1 - Obj
         objb = abs(objd)
         IF ( objb<1.0E-10 ) objb = 0.
         IF ( Nac==0 .OR. objb>0. ) Ncobj = 0
         IF ( Ncobj>1 ) Ncobj = 0
!     ------------------------------------------------------------------
!                                  PRINT
!     ------------------------------------------------------------------
!     PRINT MOVE PARAMETER, NEW X-VECTOR AND CONSTRAINTS.
         IF ( Iprint>=3 ) THEN
            WRITE (Iout,99043) Alp
99043       FORMAT (/5X,'CALCULATED ALPHA =',E14.5)
         ENDIF
         IF ( Iprint>=2 ) THEN
            IF ( objb<=0. ) THEN
               IF ( Iprint==2 ) WRITE (Iout,99044) Iter , Obj
99044          FORMAT (////5X,'ITER =',I5,5X,'OBJ =',E14.5,5X,'NO CHANGE IN OBJ')
               IF ( Iprint>2 ) WRITE (Iout,99045) Obj
99045          FORMAT (/5X,'OBJ =',E15.6,5X,'NO CHANGE ON OBJ')
            ELSEIF ( Iprint==2 ) THEN
               WRITE (Iout,99047) Iter , Obj
99047          FORMAT (////5X,'ITER =',I5,5X,'OBJ =',E14.5)
            ELSE
               WRITE (Iout,99046) Obj
            ENDIF
            WRITE (Iout,99049)
            DO i = 1 , Ndv
               ff1 = 1.
               IF ( Nscal/=0 ) ff1 = Scal(i)
               G1(i) = ff1*X(i)
            ENDDO
            DO i = 1 , Ndv , 6
               m1 = min0(Ndv,i+5)
               WRITE (Iout,99005) i , (G1(j),j=i,m1)
            ENDDO
            IF ( Ncon/=0 ) THEN
               WRITE (Iout,99051)
               DO i = 1 , Ncon , 6
                  m1 = min0(Ncon,i+5)
                  WRITE (Iout,99005) i , (G(j),j=i,m1)
               ENDDO
            ENDIF
         ENDIF
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
         IF ( Ncon>0 ) THEN
            Nfeasct = 10
            DO i = 1 , Ncon
               c1 = Ctam
               IF ( Isc(i)>0 ) c1 = Ctbm
               IF ( G(i)>c1 ) THEN
                  Nfeas = Nfeas + 1
                  GOTO 10
               ENDIF
            ENDDO
            IF ( Nfeas>0 ) Abobj1 = .05
            Nfeas = 0
            Phi = 5.
 10         IF ( Nfeas>=Nfeasct ) THEN
               spag_nextblock_1 = 10
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDIF
!
!  END OF INSERTED FIX
!
!     ------------------------------------------------------------------
!                          CHECK CONVERGENCE
!     ------------------------------------------------------------------
!     STOP IF ITER EQUALS ITMAX.
         IF ( Iter<Itmax ) THEN
!     ------------------------------------------------------------------
!                     ABSOLUTE CHANGE IN OBJECTIVE
!     ------------------------------------------------------------------
            objb = abs(objd)
            Kobj = Kobj + 1
            IF ( objb>=Dabfun .OR. Nfeas>0 ) Kobj = 0
!     ------------------------------------------------------------------
!                     RELATIVE CHANGE IN OBJECTIVE
!     ------------------------------------------------------------------
            IF ( abs(Obj1)>1.0E-10 ) objd = objd/abs(Obj1)
            Abobj1 = .5*(abs(Abobj)+abs(objd))
            Abobj = abs(objd)
            Iobj = Iobj + 1
            IF ( Nvc>0 .OR. objd>=Delfun ) Iobj = 0
            IF ( Iobj<Itrm .AND. Kobj<Itrm ) THEN
               Obj1 = Obj
!     ------------------------------------------------------------------
!           REDUCE CT IF OBJECTIVE FUNCTION IS CHANGING SLOWLY
!     ------------------------------------------------------------------
               IF ( Iobj>=1 .AND. Nac/=0 ) THEN
                  Ct = Dct*Ct
                  Ctl = Ctl*Dctl
                  IF ( abs(Ct)<Ctmin ) Ct = -Ctmin
                  IF ( abs(Ctl)<Ctlmin ) Ctl = -Ctlmin
               ENDIF
               spag_nextblock_1 = 2
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDIF
         spag_nextblock_1 = 10
      CASE (10)
         IF ( Nac>=N3 ) WRITE (Iout,99053)
99053    FORMAT (/5X,'THE NUMBER OF ACTIVE AND VIOLATED CONSTRAINTS EXCEEDS N3-1.'/5X,                                             &
                &'DIMENSIONED SIZE OF MATRICES A AND B AND VECTOR IC IS INSUFFICIENT'/5X,                                          &
                &'OPTIMIZATION TERMINATED AND CONTROL RETURNED TO MAIN PROGRAM.')
!     ------------------------------------------------------------------
!     ****************  FINAL FUNCTION INFORMATION  ********************
!     ------------------------------------------------------------------
         IF ( Nscal/=0 ) THEN
!     UN-SCALE THE DESIGN VARIABLES.
            DO i = 1 , Ndv
               Xi = Scal(i)
               IF ( Nside/=0 ) THEN
                  Vlb(i) = Xi*Vlb(i)
                  Vub(i) = Xi*Vub(i)
               ENDIF
               X(i) = Xi*X(i)
            ENDDO
         ENDIF
!     ------------------------------------------------------------------
!                           PRINT FINAL RESULTS
!     ------------------------------------------------------------------
         IF ( Iprint/=0 .AND. Nac<N3 ) THEN
            WRITE (Iout,99054)
99054       FORMAT ('1',////4X,'FINAL OPTIMIZATION INFORMATION')
            WRITE (Iout,99046) Obj
            WRITE (Iout,99049)
            DO i = 1 , Ndv , 6
               m1 = min0(Ndv,i+5)
               WRITE (Iout,99005) i , (X(j),j=i,m1)
            ENDDO
            IF ( Ncon/=0 ) THEN
               WRITE (Iout,99051)
               DO i = 1 , Ncon , 6
                  m1 = min0(Ncon,i+5)
                  WRITE (Iout,99005) i , (G(j),j=i,m1)
               ENDDO
!     DETERMINE WHICH CONSTRAINTS ARE ACTIVE AND PRINT.
               Nac = 0
               Nvc = 0
               DO i = 1 , Ncon
                  Cta = Ctam
                  IF ( Isc(i)>0 ) Cta = Ctbm
                  gi = G(i)
                  IF ( gi>Cta ) THEN
                     Nvc = Nvc + 1
                     Ms1(Nvc) = i
                  ELSEIF ( gi>=Ct .OR. Isc(i)/=0 ) THEN
                     IF ( gi>=Ctl .OR. Isc(i)<=0 ) THEN
                        Nac = Nac + 1
                        Ic(Nac) = i
                     ENDIF
                  ENDIF
               ENDDO
               WRITE (Iout,99010) Nac
               IF ( Nac/=0 ) THEN
                  WRITE (Iout,99011)
                  WRITE (Iout,99052) (Ic(j),j=1,Nac)
               ENDIF
               WRITE (Iout,99012) Nvc
               IF ( Nvc/=0 ) THEN
                  WRITE (Iout,99011)
                  WRITE (Iout,99052) (Ms1(j),j=1,Nvc)
               ENDIF
            ENDIF
            IF ( Nside/=0 ) THEN
!     DETERMINE WHICH SIDE CONSTRAINTS ARE ACTIVE AND PRINT.
               Nac = 0
               DO i = 1 , Ndv
                  Xi = X(i)
                  xid = Vlb(i)
                  x12 = abs(xid)
                  IF ( x12<1. ) x12 = 1.
                  gi = (xid-Xi)/x12
                  IF ( gi>=-1.0E-6 ) THEN
                     Nac = Nac + 1
                     Ms1(Nac) = -i
                  ENDIF
                  xid = Vub(i)
                  x12 = abs(xid)
                  IF ( x12<1. ) x12 = 1.
                  gi = (Xi-xid)/x12
                  IF ( gi>=-1.0E-6 ) THEN
                     Nac = Nac + 1
                     Ms1(Nac) = i
                  ENDIF
               ENDDO
               WRITE (Iout,99013) Nac
               IF ( Nac/=0 ) THEN
                  WRITE (Iout,99014)
                  WRITE (Iout,99052) (Ms1(j),j=1,Nac)
               ENDIF
            ENDIF
            WRITE (Iout,99019)
99019       FORMAT (/5X,'TERMINATION CRITERION')
            IF ( Iter>=Itmax ) WRITE (Iout,99020)
99020       FORMAT (10X,'ITER EQUALS ITMAX')
            IF ( Nfeas>=Nfeasct ) WRITE (Iout,99021)
99021       FORMAT (10X,'NFEASCT CONSECUTIVE ITERATIONS FAILED TO PRODUCE A   FEASIBLE DESIGN')
            IF ( Iobj>=Itrm ) WRITE (Iout,99022) Itrm
99022       FORMAT (10X,'ABS(1-OBJ(I-1)/OBJ(I)) LESS THAN DELFUN FOR',I3,' ITERATIONS')
            IF ( Kobj>=Itrm ) WRITE (Iout,99023) Itrm
99023       FORMAT (10X,'ABS(OBJ(I)-OBJ(I-1))   LESS THAN DABFUN FOR',I3,' ITERATIONS')
            WRITE (Iout,99024) Iter
99024       FORMAT (/5X,'NUMBER OF ITERATIONS =',I5)
            WRITE (Iout,99055) Ncal(1)
99055       FORMAT (/5X,'OBJECTIVE FUNCTION WAS EVALUATED',8X,I5,2X,'TIMES')
            IF ( Ncon>0 ) WRITE (Iout,99056) Ncal(1)
99056       FORMAT (/5X,'CONSTRAINT FUNCTIONS WERE EVALUATED',I10,2X,'TIMES')
            IF ( Nfdg/=0 ) WRITE (Iout,99057) Ncal(2)
99057       FORMAT (/5X,'GRADIENT OF OBJECTIVE WAS CALCULATED',I9,2X,'TIMES')
            IF ( Ncon>0 .AND. Nfdg==1 ) WRITE (Iout,99058) Ncal(2)
99058       FORMAT (/5X,'GRADIENTS OF CONSTRAINTS WERE CALCULATED',I5,2X,'TIMES')
         ENDIF
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
      CASE (11)
         IF ( Nscal==0 .OR. Igoto==0 ) RETURN
!     UN-SCALE VARIABLES.
         DO i = 1 , Ndv
            C(i) = X(i)
            X(i) = X(i)*Scal(i)
         ENDDO
         RETURN
      END SELECT
   ENDDO SPAG_DispatchLoop_1
99005 FORMAT (3X,I5,')',2X,6E13.5)
99010 FORMAT (/5X,'THERE ARE',I5,' ACTIVE CONSTRAINTS')
99011 FORMAT (5X,'CONSTRAINT NUMBERS ARE')
99012 FORMAT (/5X,'THERE ARE',I5,' VIOLATED CONSTRAINTS')
99013 FORMAT (/5X,'THERE ARE',I5,' ACTIVE SIDE CONSTRAINTS')
99014 FORMAT (5X,'DECISION VARIABLES AT LOWER OR UPPER BOUNDS',' (MINUS INDICATES LOWER BOUND)')
99040 FORMAT (' ')
99046 FORMAT (/5X,'OBJ =',E15.6)
99049 FORMAT (/5X,'DECISION VARIABLES (X-VECTOR)')
99050 FORMAT (3X,7E13.4)
99051 FORMAT (/5X,'CONSTRAINT VALUES (G-VECTOR)')
99052 FORMAT (5X,15I5)
END SUBROUTINE cnmn00

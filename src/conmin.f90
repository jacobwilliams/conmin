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
!*==cnmn01.f90 processed by SPAG 8.01MH 10:24  1 Mar 2025
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE cnmn01(Jgoto,X,Df,G,Isc,Ic,A,G1,Vlb,Vub,Scal,C,Ncal,Dx,Dx1,Fi,Xi,Iii,N1,N2,N3,N4)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   DOUBLE PRECISION A , Abobj1 , Alphax , C , Ct , Ctl , Ctlmin , Ctmin , Dabfun , Delfun , Df , Dx , Dx1 , Fdch , fdch1 , Fdchm , &
                  & Fi , G , G1 , Obj
   DOUBLE PRECISION Scal , Theta , Vlb , Vub , X , x1 , Xi
   INTEGER i , i1 , Ic , Icndir , Igoto , Iii , inf , Info , Infog , Iprint , Isc , Iter , Itmax , Itrm , j , Jgoto , Linobj , N1 ,&
         & N2 , N3
   INTEGER N4 , Nac , Ncal , Ncon , Ndv , Nfdg , Nfeasct , Nscal , Nside
!*** End of declarations inserted by SPAG
   COMMON /cnmn1 / Delfun , Dabfun , Fdch , Fdchm , Ct , Ctmin , Ctl , Ctlmin , Alphax , Abobj1 , Theta , Obj , Ndv , Ncon ,       &
                 & Nside , Iprint , Nfdg , Nscal , Linobj , Itmax , Itrm , Icndir , Igoto , Nac , Info , Infog , Iter , Nfeasct
   DIMENSION X(N1) , Df(N1) , G(N2) , Isc(N2) , Ic(N3) , A(N1,N3) , G1(N2) , Vlb(N1) , Vub(N1) , Scal(N1) , Ncal(2) , C(N4)
!     ROUTINE TO CALCULATE GRADIENT INFORMATION BY FINITE DIFFERENCE.
!     BY G. N. VANDERPLAATS                         JUNE, 1972.
!     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF.
   IF ( Jgoto/=1 ) THEN
      IF ( Jgoto==2 ) THEN
         X(Iii) = Xi
         IF ( Nfdg==0 ) Df(Iii) = Dx1*(Obj-Fi)
         IF ( Nac/=0 ) THEN
!     ------------------------------------------------------------------
!             DETERMINE GRADIENT COMPONENTS OF ACTIVE CONSTRAINTS
!     ------------------------------------------------------------------
            DO j = 1 , Nac
               i1 = Ic(j)
               A(Iii,j) = Dx1*(G(i1)-G1(i1))
            ENDDO
         ENDIF
         IF ( Iii<Ndv ) THEN
            CALL spag_block_1
            RETURN
         ENDIF
         Infog = 0
         Info = inf
         Jgoto = 0
         Obj = Fi
         IF ( Ncon==0 ) RETURN
!     ------------------------------------------------------------------
!             STORE CURRENT CONSTRAINT VALUES BACK IN G-VECTOR
!     ------------------------------------------------------------------
         DO i = 1 , Ncon
            G(i) = G1(i)
         ENDDO
         RETURN
      ELSE
         Infog = 0
         inf = Info
         Nac = 0
         IF ( Linobj==0 .OR. Iter<=1 ) THEN
!     ------------------------------------------------------------------
!                    GRADIENT OF LINEAR OBJECTIVE
!     ------------------------------------------------------------------
            IF ( Nfdg==2 ) Jgoto = 1
            IF ( Nfdg==2 ) RETURN
         ENDIF
      ENDIF
   ENDIF
   Jgoto = 0
   IF ( Nfdg==2 .AND. Ncon==0 ) RETURN
   IF ( Ncon/=0 ) THEN
!     ------------------------------------------------------------------
!       * * * DETERMINE WHICH CONSTRAINTS ARE ACTIVE OR VIOLATED * * *
!     ------------------------------------------------------------------
      DO i = 1 , Ncon
         IF ( G(i)>=Ct ) THEN
            IF ( Isc(i)<=0 .OR. G(i)>=Ctl ) THEN
               Nac = Nac + 1
               IF ( Nac>=N3 ) RETURN
               Ic(Nac) = i
            ENDIF
         ENDIF
      ENDDO
      IF ( Nfdg==2 .AND. Nac==0 ) RETURN
      IF ( (Linobj>0 .AND. Iter>1) .AND. Nac==0 ) RETURN
!     ------------------------------------------------------------------
!                  STORE VALUES OF CONSTRAINTS IN G1
!     ------------------------------------------------------------------
      DO i = 1 , Ncon
         G1(i) = G(i)
      ENDDO
   ENDIF
   Jgoto = 0
   IF ( Nac==0 .AND. Nfdg==2 ) RETURN
!     ------------------------------------------------------------------
!                            CALCULATE GRADIENTS
!     ------------------------------------------------------------------
   Infog = 1
   Info = 1
   Fi = Obj
   Iii = 0
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
      Iii = Iii + 1
      Xi = X(Iii)
      Dx = Fdch*Xi
      Dx = abs(Dx)
      fdch1 = Fdchm
      IF ( Nscal/=0 ) fdch1 = Fdchm/Scal(Iii)
      IF ( Dx<fdch1 ) Dx = fdch1
      x1 = Xi + Dx
      IF ( Nside/=0 ) THEN
         IF ( x1>Vub(Iii) ) Dx = -Dx
      ENDIF
      Dx1 = 1./Dx
      X(Iii) = Xi + Dx
      Ncal(1) = Ncal(1) + 1
!     ------------------------------------------------------------------
!                         FUNCTION EVALUATION
!     ------------------------------------------------------------------
      Jgoto = 2
      RETURN
   END SUBROUTINE spag_block_1
END SUBROUTINE cnmn01
!*==cnmn02.f90 processed by SPAG 8.01MH 10:24  1 Mar 2025
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE cnmn02(Ncalc,Slope,Dftdf1,Df,S,N1)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   DOUBLE PRECISION Abobj1 , Alphax , beta , Ct , Ctl , Ctlmin , Ctmin , Dabfun , Delfun , Df , dfi , dftdf , Dftdf1 , Fdch ,      &
                  & Fdchm , Obj , S , s1 , s2 , si
   DOUBLE PRECISION Slope , Theta
   INTEGER i , Icndir , Igoto , Info , Infog , Iprint , Iter , Itmax , Itrm , Linobj , N1 , Nac , Ncalc , Ncon , Ndv , Nfdg ,      &
         & Nfeasct , Nscal , Nside
!*** End of declarations inserted by SPAG
   COMMON /cnmn1 / Delfun , Dabfun , Fdch , Fdchm , Ct , Ctmin , Ctl , Ctlmin , Alphax , Abobj1 , Theta , Obj , Ndv , Ncon ,       &
                 & Nside , Iprint , Nfdg , Nscal , Linobj , Itmax , Itrm , Icndir , Igoto , Nac , Info , Infog , Iter , Nfeasct
   DIMENSION Df(N1) , S(N1)
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
   DO i = 1 , Ndv
      dfi = Df(i)
      dftdf = dftdf + dfi*dfi
   ENDDO
!     ------------------------------------------------------------------
!     **********                FIND DIRECTION S              **********
!     ------------------------------------------------------------------
   IF ( Ncalc==1 ) THEN
      IF ( Dftdf1>=1.0E-20 ) THEN
!     ------------------------------------------------------------------
!                 FIND FLETCHER-REEVES CONJUGATE DIRECTION
!     ------------------------------------------------------------------
         beta = dftdf/Dftdf1
         Slope = 0.
         DO i = 1 , Ndv
            dfi = Df(i)
            si = beta*S(i) - dfi
            Slope = Slope + si*dfi
            S(i) = si
         ENDDO
         CALL spag_block_1
         RETURN
      ENDIF
   ENDIF
   Ncalc = 0
!     ------------------------------------------------------------------
!                  CALCULATE DIRECTION OF STEEPEST DESCENT
!     ------------------------------------------------------------------
   DO i = 1 , Ndv
      S(i) = -Df(i)
   ENDDO
   Slope = -dftdf
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
!     ------------------------------------------------------------------
!                  NORMALIZE S TO MAX ABS VALUE OF UNITY
!     ------------------------------------------------------------------
      s1 = 0.
      DO i = 1 , Ndv
         s2 = abs(S(i))
         IF ( s2>s1 ) s1 = s2
      ENDDO
      IF ( s1<1.0E-20 ) s1 = 1.0E-20
      s1 = 1./s1
      Dftdf1 = dftdf*s1
      DO i = 1 , Ndv
         S(i) = s1*S(i)
      ENDDO
      Slope = s1*Slope
   END SUBROUTINE spag_block_1
END SUBROUTINE cnmn02
!*==cnmn03.f90 processed by SPAG 8.01MH 10:24  1 Mar 2025
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE cnmn03(X,S,Slope,Alp,Fff,A1,A2,A3,A4,F1,F2,F3,F4,App,N1,Ncal,Kount,Jgoto)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   DOUBLE PRECISION A1 , A2 , A3 , A4 , aa , ab , ab2 , ab3 , Abobj1 , Alp , Alphax , ap , ap1 , App , Ct , Ctl , Ctlmin , Ctmin , &
                  & Dabfun , Delfun
   DOUBLE PRECISION F1 , F2 , F3 , F4 , Fdch , Fdchm , ff , Fff , Obj , S , Slope , Theta , X , zro
   INTEGER i , Icndir , Igoto , ii , Info , Infog , Iout , Iprint , Iter , Itmax , Itrm , Jgoto , jj , Kount , Linobj , N1 , Nac , &
         & Ncal , Ncon , Ndv
   INTEGER Nfdg , Nfeasct , Nscal , Nside
!*** End of declarations inserted by SPAG
   COMMON /cnmn1 / Delfun , Dabfun , Fdch , Fdchm , Ct , Ctmin , Ctl , Ctlmin , Alphax , Abobj1 , Theta , Obj , Ndv , Ncon ,       &
                 & Nside , Iprint , Nfdg , Nscal , Linobj , Itmax , Itrm , Icndir , Igoto , Nac , Info , Infog , Iter , Nfeasct
   COMMON /output/ Iout
   DIMENSION X(N1) , S(N1) , Ncal(2)
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
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
         IF ( Jgoto/=0 ) THEN
            IF ( Jgoto==1 ) THEN
               F2 = Obj
               IF ( Iprint>4 ) WRITE (Iout,99004) F2
               IF ( F2<F1 ) THEN
                  spag_nextblock_1 = 4
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
!     ------------------------------------------------------------------
!                     CHECK FOR ILL-CONDITIONING
!     ------------------------------------------------------------------
               IF ( Kount>5 ) THEN
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               ff = 2.*abs(F1)
               IF ( F2<ff ) THEN
!     ------------------------------------------------------------------
!     **********        2-POINT QUADRATIC INTERPOLATION       **********
!     ------------------------------------------------------------------
                  jj = 1
                  ii = 1
                  CALL cnmn04(ii,App,zro,A1,F1,Slope,A2,F2,zro,zro,zro,zro)
                  IF ( App<zro .OR. App>A2 ) THEN
                     spag_nextblock_1 = 4
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
                  F3 = F2
                  A3 = A2
                  A2 = App
                  jj = 0
!     ------------------------------------------------------------------
!                  UPDATE DESIGN VECTOR AND FUNCTION VALUE
!     ------------------------------------------------------------------
                  ap = A2 - Alp
                  Alp = A2
                  DO i = 1 , Ndv
                     X(i) = X(i) + ap*S(i)
                  ENDDO
                  IF ( Iprint>4 ) WRITE (Iout,99002) A2
                  IF ( Iprint>4 ) WRITE (Iout,99003) (X(i),i=1,Ndv)
                  Ncal(1) = Ncal(1) + 1
                  Jgoto = 3
                  RETURN
               ELSE
                  ff = 5.*abs(F1)
                  IF ( F2<ff ) THEN
                     spag_nextblock_1 = 3
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
                  A2 = .5*A2
                  ap = -A2
                  Alp = A2
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
            ELSEIF ( Jgoto==2 ) THEN
               F2 = Obj
!     PROCEED TO CUBIC INTERPOLATION.
               IF ( Iprint>4 ) WRITE (Iout,99004) F2
               spag_nextblock_1 = 6
               CYCLE SPAG_DispatchLoop_1
            ELSEIF ( Jgoto==3 ) THEN
               F2 = Obj
               IF ( Iprint>4 ) WRITE (Iout,99004) F2
               spag_nextblock_1 = 5
               CYCLE SPAG_DispatchLoop_1
            ELSEIF ( Jgoto==4 ) THEN
               F3 = Obj
               IF ( Iprint>4 ) WRITE (Iout,99004) F3
               spag_nextblock_1 = 5
               CYCLE SPAG_DispatchLoop_1
            ELSEIF ( Jgoto==5 ) THEN
               IF ( Iprint>4 ) WRITE (Iout,99004) Obj
!     ------------------------------------------------------------------
!                         CHECK CONVERGENCE
!     ------------------------------------------------------------------
               aa = 1. - App/A2
               ab2 = abs(F2)
               ab3 = abs(Obj)
               ab = ab2
               IF ( ab3>ab ) ab = ab3
               IF ( ab<1.0E-15 ) ab = 1.0E-15
               ab = (ab2-ab3)/ab
               IF ( abs(ab)<1.0E-15 .AND. abs(aa)<.001 ) THEN
                  spag_nextblock_1 = 10
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               A4 = A3
               F4 = F3
               A3 = App
               F3 = Obj
               IF ( A3<=A2 ) THEN
                  A3 = A2
                  F3 = F2
                  A2 = App
                  F2 = Obj
               ENDIF
               spag_nextblock_1 = 8
               CYCLE SPAG_DispatchLoop_1
            ELSEIF ( Jgoto==6 ) THEN
               F4 = Obj
               IF ( Iprint>4 ) WRITE (Iout,99004) F4
               IF ( F4>F3 ) THEN
                  spag_nextblock_1 = 8
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               A1 = A2
               F1 = F2
               A2 = A3
               F2 = F3
               A3 = A4
               F3 = F4
               spag_nextblock_1 = 7
               CYCLE SPAG_DispatchLoop_1
            ELSEIF ( Jgoto==7 ) THEN
               IF ( Iprint>4 ) WRITE (Iout,99004) Obj
               spag_nextblock_1 = 9
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDIF
!     ------------------------------------------------------------------
!                     INITIAL INFORMATION  (ALPHA=0)
!     ------------------------------------------------------------------
         IF ( Slope<0. ) THEN
            IF ( Iprint>4 ) WRITE (Iout,99001)
!     ------------------------------------------------------------------
!                                 FORMATS
!     ------------------------------------------------------------------
!
!
99001       FORMAT (/////5X,'* * * UNCONSTRAINED ONE-DIMENSIONAL SEARCH INFORMATION * * *')
            Fff = Obj
            ap1 = 0.
            A1 = 0.
            F1 = Obj
            A2 = Alp
            A3 = 0.
            F3 = 0.
            ap = A2
            Kount = 0
         ELSE
            Alp = 0.
            RETURN
         ENDIF
         spag_nextblock_1 = 2
      CASE (2)
!     ------------------------------------------------------------------
!            MOVE A DISTANCE AP*S AND UPDATE FUNCTION VALUE
!     ------------------------------------------------------------------
         Kount = Kount + 1
         DO i = 1 , Ndv
            X(i) = X(i) + ap*S(i)
         ENDDO
         IF ( Iprint>4 ) WRITE (Iout,99002) ap
         IF ( Iprint>4 ) WRITE (Iout,99003) (X(i),i=1,Ndv)
         Ncal(1) = Ncal(1) + 1
         Jgoto = 1
         RETURN
      CASE (3)
         F3 = F2
         A3 = A2
         A2 = .5*A2
!     ------------------------------------------------------------------
!                 UPDATE DESIGN VECTOR AND FUNCTION VALUE
!     ------------------------------------------------------------------
         ap = A2 - Alp
         Alp = A2
         DO i = 1 , Ndv
            X(i) = X(i) + ap*S(i)
         ENDDO
         IF ( Iprint>4 ) WRITE (Iout,99002) A2
         IF ( Iprint>4 ) WRITE (Iout,99003) (X(i),i=1,Ndv)
         Ncal(1) = Ncal(1) + 1
         Jgoto = 2
         RETURN
      CASE (4)
         A3 = 2.*A2
!     ------------------------------------------------------------------
!                  UPDATE DESIGN VECTOR AND FUNCTION VALUE
!     ------------------------------------------------------------------
         ap = A3 - Alp
         Alp = A3
         DO i = 1 , Ndv
            X(i) = X(i) + ap*S(i)
         ENDDO
         IF ( Iprint>4 ) WRITE (Iout,99002) A3
         IF ( Iprint>4 ) WRITE (Iout,99003) (X(i),i=1,Ndv)
         Ncal(1) = Ncal(1) + 1
         Jgoto = 4
         RETURN
      CASE (5)
         IF ( F3<F2 ) THEN
            spag_nextblock_1 = 7
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         spag_nextblock_1 = 6
      CASE (6)
!     ------------------------------------------------------------------
!     **********       3-POINT CUBIC INTERPOLATION      **********
!     ------------------------------------------------------------------
         ii = 3
         CALL cnmn04(ii,App,zro,A1,F1,Slope,A2,F2,A3,F3,zro,zro)
         IF ( App>=zro .AND. App<=A3 ) THEN
!     ------------------------------------------------------------------
!     UPDATE DESIGN VECTOR AND FUNCTION VALUE.
!     ------------------------------------------------------------------
            ap1 = App
            ap = App - Alp
            Alp = App
            DO i = 1 , Ndv
               X(i) = X(i) + ap*S(i)
            ENDDO
            IF ( Iprint>4 ) WRITE (Iout,99002) Alp
            IF ( Iprint>4 ) WRITE (Iout,99003) (X(i),i=1,Ndv)
            Ncal(1) = Ncal(1) + 1
            Jgoto = 5
            RETURN
         ENDIF
         spag_nextblock_1 = 7
      CASE (7)
!     ------------------------------------------------------------------
!     **********        4-POINT CUBIC INTERPOLATION       **********
!     ------------------------------------------------------------------
         A4 = 2.*A3
!     UPDATE DESIGN VECTOR AND FUNCTION VALUE.
         ap = A4 - Alp
         Alp = A4
         DO i = 1 , Ndv
            X(i) = X(i) + ap*S(i)
         ENDDO
         IF ( Iprint>4 ) WRITE (Iout,99002) Alp
         IF ( Iprint>4 ) WRITE (Iout,99003) (X(i),i=1,Ndv)
         Ncal(1) = Ncal(1) + 1
         Jgoto = 6
         RETURN
      CASE (8)
         ii = 4
         CALL cnmn04(ii,App,A1,A1,F1,Slope,A2,F2,A3,F3,A4,F4)
         IF ( App>A1 ) THEN
!     ------------------------------------------------------------------
!                 UPDATE DESIGN VECTOR AND FUNCTION VALUE
!     ------------------------------------------------------------------
            ap = App - Alp
            Alp = App
            DO i = 1 , Ndv
               X(i) = X(i) + ap*S(i)
            ENDDO
            IF ( Iprint>4 ) WRITE (Iout,99002) Alp
            IF ( Iprint>4 ) WRITE (Iout,99003) (X(i),i=1,Ndv)
            Ncal(1) = Ncal(1) + 1
            Jgoto = 7
            RETURN
         ELSE
            ap = A1 - Alp
            Alp = A1
            Obj = F1
            DO i = 1 , Ndv
               X(i) = X(i) + ap*S(i)
            ENDDO
         ENDIF
         spag_nextblock_1 = 9
      CASE (9)
!     ------------------------------------------------------------------
!                    CHECK FOR ILL-CONDITIONING
!     ------------------------------------------------------------------
         IF ( Obj<=F2 .AND. Obj<=F3 ) THEN
            IF ( Obj<=F1 ) THEN
               spag_nextblock_1 = 10
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            ap = A1 - Alp
            Alp = A1
            Obj = F1
         ELSEIF ( F2<F3 ) THEN
            Obj = F2
            ap = A2 - Alp
            Alp = A2
         ELSE
            Obj = F3
            ap = A3 - Alp
            Alp = A3
         ENDIF
!     ------------------------------------------------------------------
!                       UPDATE DESIGN VECTOR
!     ------------------------------------------------------------------
         DO i = 1 , Ndv
            X(i) = X(i) + ap*S(i)
         ENDDO
         spag_nextblock_1 = 10
      CASE (10)
!     ------------------------------------------------------------------
!                     CHECK FOR MULTIPLE MINIMA
!     ------------------------------------------------------------------
         IF ( Obj>Fff ) THEN
!     INITIAL FUNCTION IS MINIMUM.
            DO i = 1 , Ndv
               X(i) = X(i) - Alp*S(i)
            ENDDO
            Alp = 0.
            Obj = Fff
         ENDIF
         Jgoto = 0
         RETURN
      END SELECT
   ENDDO SPAG_DispatchLoop_1
99002 FORMAT (/5X,'ALPHA =',E14.5/5X,'X-VECTOR')
99003 FORMAT (5X,6E13.5)
99004 FORMAT (/5X,'OBJ =',E14.5)
END SUBROUTINE cnmn03
!*==cnmn04.f90 processed by SPAG 8.01MH 10:24  1 Mar 2025
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE cnmn04(Ii,Xbar,Eps,X1,Y1,Slope,X2,Y2,X3,Y3,X4,Y4)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   DOUBLE PRECISION aa , bac , bb , cc , dnom , dx , Eps , q1 , q2 , q3 , q4 , q5 , q6 , qq , Slope , X1 , x11 , x111 , X2 , x21
   DOUBLE PRECISION x22 , x222 , X3 , x31 , x32 , x33 , X4 , x41 , x42 , x44 , Xbar , xbar1 , Y1 , Y2 , Y3 , Y4
   INTEGER Ii , nslop
!*** End of declarations inserted by SPAG
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
   x21 = X2 - X1
   IF ( abs(x21)<1.0E-20 ) RETURN
   nslop = mod(Ii,2)
   IF ( Ii==2 ) THEN
      CALL spag_block_2
      RETURN
   ENDIF
   IF ( Ii==3 ) THEN
      CALL spag_block_3
      RETURN
   ENDIF
   IF ( Ii==4 ) THEN
!     ------------------------------------------------------------------
!                    II=4: 4-POINT CUBIC INTERPOLATION
!     ------------------------------------------------------------------
      x21 = X2 - X1
      x31 = X3 - X1
      x41 = X4 - X1
      x32 = X3 - X2
      x42 = X4 - X2
      x11 = X1*X1
      x22 = X2*X2
      x33 = X3*X3
      x44 = X4*X4
      x111 = X1*x11
      x222 = X2*x22
      q2 = x31*x21*x32
      IF ( abs(q2)<1.0E-30 ) RETURN
      q1 = x111*x32 - x222*x31 + X3*x33*x21
      q4 = x111*x42 - x222*x41 + X4*x44*x21
      q5 = x41*x21*x42
      dnom = q2*q4 - q1*q5
      IF ( abs(dnom)<1.0E-30 ) THEN
         CALL spag_block_4
         RETURN
      ENDIF
      q3 = Y3*x21 - Y2*x31 + Y1*x32
      q6 = Y4*x21 - Y2*x41 + Y1*x42
      aa = (q2*q6-q3*q5)/dnom
      bb = (q3-q1*aa)/q2
      cc = (Y2-Y1-aa*(x222-x111))/x21 - bb*(X1+X2)
      bac = bb*bb - 3.*aa*cc
      IF ( abs(aa)<1.0E-20 .OR. bac<0. ) THEN
         CALL spag_block_4
         RETURN
      ENDIF
      bac = sqrt(bac)
      Xbar = (bac-bb)/(3.*aa)
      IF ( Xbar<Eps ) Xbar = xbar1
      RETURN
   ENDIF
   CALL spag_block_1
CONTAINS
   SUBROUTINE spag_block_1
!     ------------------------------------------------------------------
!                 II=1: 2-POINT QUADRATIC INTERPOLATION
!     ------------------------------------------------------------------
      Ii = 1
      dx = X1 - X2
      IF ( abs(dx)<1.0E-20 ) RETURN
      aa = (Slope+(Y2-Y1)/dx)/dx
      IF ( aa<1.0E-20 ) RETURN
      bb = Slope - 2.*aa*X1
      Xbar = -.5*bb/aa
      IF ( Xbar<Eps ) Xbar = xbar1
      RETURN
   END SUBROUTINE spag_block_1
   SUBROUTINE spag_block_2
!     ------------------------------------------------------------------
!                 II=2: 3-POINT QUADRATIC INTERPOLATION
!     ------------------------------------------------------------------
      Ii = 2
      x21 = X2 - X1
      x31 = X3 - X1
      x32 = X3 - X2
      qq = x21*x31*x32
      IF ( abs(qq)<1.0E-20 ) RETURN
      aa = (Y1*x32-Y2*x31+Y3*x21)/qq
      IF ( aa<1.0E-20 ) THEN
         IF ( nslop==0 ) RETURN
         CALL spag_block_1
         RETURN
      ELSE
         bb = (Y2-Y1)/x21 - aa*(X1+X2)
         Xbar = -.5*bb/aa
         IF ( Xbar<Eps ) Xbar = xbar1
         RETURN
      ENDIF
   END SUBROUTINE spag_block_2
   SUBROUTINE spag_block_3
!     ------------------------------------------------------------------
!                   II=3: 3-POINT CUBIC INTERPOLATION
!     ------------------------------------------------------------------
      Ii = 3
      x21 = X2 - X1
      x31 = X3 - X1
      x32 = X3 - X2
      qq = x21*x31*x32
      IF ( abs(qq)<1.0E-20 ) RETURN
      x11 = X1*X1
      dnom = X2*X2*x31 - x11*x32 - X3*X3*x21
      IF ( abs(dnom)<1.0E-20 ) THEN
         CALL spag_block_2
         RETURN
      ENDIF
      aa = ((x31*x31*(Y2-Y1)-x21*x21*(Y3-Y1))/(x31*x21)-Slope*x32)/dnom
      IF ( abs(aa)<1.0E-20 ) THEN
         CALL spag_block_2
         RETURN
      ENDIF
      bb = ((Y2-Y1)/x21-Slope-aa*(X2*X2+X1*X2-2.*x11))/x21
      cc = Slope - 3.*aa*x11 - 2.*bb*X1
      bac = bb*bb - 3.*aa*cc
      IF ( bac<0. ) THEN
         CALL spag_block_2
         RETURN
      ENDIF
      bac = sqrt(bac)
      Xbar = (bac-bb)/(3.*aa)
      IF ( Xbar<Eps ) Xbar = Eps
      RETURN
   END SUBROUTINE spag_block_3
   SUBROUTINE spag_block_4
      IF ( nslop/=1 ) THEN
         CALL spag_block_2
         RETURN
      ENDIF
      CALL spag_block_3
      RETURN
   END SUBROUTINE spag_block_4
END SUBROUTINE cnmn04
!*==cnmn05.f90 processed by SPAG 8.01MH 10:24  1 Mar 2025
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE cnmn05(G,Df,A,S,B,C,Slope,Phi,Isc,Ic,Ms1,Nvc,N1,N2,N3,N4,N5)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   DOUBLE PRECISION A , a1 , Abobj1 , Alphax , B , C , c1 , Ct , ct1 , ct2 , cta , ctam , ctb , ctbm , ctc , ctd , Ctl , Ctlmin ,  &
                  & Ctmin , Dabfun
   DOUBLE PRECISION Delfun , Df , Fdch , Fdchm , G , gg , Obj , Phi , S , s1 , sg , Slope , Theta , thmax , tht
   INTEGER i , Ic , Icndir , Igoto , Info , Infog , Iout , Iprint , Isc , Iter , Itmax , Itrm , j , j1 , k , Linobj , Ms1 , N1 ,   &
         & N2 , N3
   INTEGER N4 , N5 , Nac , nac1 , nci , ncj , Ncon , ndb , Ndv , ndv1 , ndv2 , ner , Nfdg , Nfeasct , Nscal , Nside , Nvc
!*** End of declarations inserted by SPAG
   COMMON /cnmn1 / Delfun , Dabfun , Fdch , Fdchm , Ct , Ctmin , Ctl , Ctlmin , Alphax , Abobj1 , Theta , Obj , Ndv , Ncon ,       &
                 & Nside , Iprint , Nfdg , Nscal , Linobj , Itmax , Itrm , Icndir , Igoto , Nac , Info , Infog , Iter , Nfeasct
   COMMON /output/ Iout
   DIMENSION Df(N1) , G(N2) , Isc(N2) , Ic(N3) , A(N1,N3) , S(N1) , C(N4) , Ms1(N5) , B(N3,N3)
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
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
         DO i = 1 , Nac
!     CALCULATE THETA
            nci = Ic(i)
            ncj = 1
            IF ( nci<=Ncon ) ncj = Isc(nci)
            c1 = G(nci)
            ctd = ct1
            ctc = ctam
            IF ( ncj>0 ) THEN
               ctc = ctbm
               ctd = ct2
            ENDIF
            IF ( c1>ctc ) Nvc = Nvc + 1
            tht = 0.
            gg = 1. + ctd*c1
            IF ( ncj==0 .OR. c1>ctc ) tht = Theta*gg*gg
            IF ( tht>50. ) tht = 50.
            IF ( tht>thmax ) thmax = tht
            A(ndv1,i) = tht
!     ------------------------------------------------------------------
!                    NORMALIZE GRADIENTS OF CONSTRAINTS
!     ------------------------------------------------------------------
            A(ndv2,i) = 1.
            IF ( nci<=Ncon ) THEN
               a1 = 0.
               DO j = 1 , Ndv
                  a1 = a1 + A(j,i)**2
               ENDDO
               IF ( a1<1.0E-20 ) a1 = 1.0E-20
               a1 = sqrt(a1)
               A(ndv2,i) = a1
               a1 = 1./a1
               DO j = 1 , Ndv
                  A(j,i) = a1*A(j,i)
               ENDDO
            ENDIF
         ENDDO
!     ------------------------------------------------------------------
!     CHECK FOR ZERO GRADIENT.  PROGRAM CHANGE-FEB, 1981, GV.
!     ------------------------------------------------------------------
         i = 0
         spag_nextblock_1 = 2
      CASE (2)
         i = i + 1
         SPAG_Loop_1_1: DO WHILE ( A(ndv2,i)<=1.0E-6 )
!     ZERO GRADIENT IS FOUND.  WRITE ERROR MESSAGE.
            IF ( Iprint>=2 ) WRITE (Iout,99001) Ic(i)
99001       FORMAT (5X,'** CONSTRAINT',I5,' HAS ZERO GRADIENT'/5X,'DELETED FROM ACTIVE SET')
!     REDUCE NAC BY ONE.
            Nac = Nac - 1
!     SHIFT COLUMNS OF A AND ROWS OF IC IF I.LE.NAC.
            IF ( i>Nac ) THEN
               spag_nextblock_1 = 3
               CYCLE SPAG_DispatchLoop_1
            ENDIF
!     SHIFT.
            DO j = i , Nac
               j1 = j + 1
               Ic(j) = Ic(j1)
               DO k = 1 , ndv2
                  A(k,j) = A(k,j1)
               ENDDO
            ENDDO
            IF ( i>Nac ) EXIT SPAG_Loop_1_1
         ENDDO SPAG_Loop_1_1
         IF ( i<Nac ) THEN
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         spag_nextblock_1 = 3
      CASE (3)
         IF ( Nac<=0 ) RETURN
         nac1 = Nac + 1
!     DETERMINE IF CONSTRAINTS ARE VIOLATED.
         Nvc = 0
         DO i = 1 , Nac
            nci = Ic(i)
            ncj = 1
            IF ( nci<=Ncon ) ncj = Isc(nci)
            ctc = ctam
            IF ( ncj>0 ) ctc = ctbm
            IF ( G(nci)>ctc ) Nvc = Nvc + 1
         ENDDO
!     ------------------------------------------------------------------
!     NORMALIZE GRADIENT OF OBJECTIVE FUNCTION AND STORE IN NAC+1
!     COLUMN OF A
!     ------------------------------------------------------------------
         a1 = 0.
         DO i = 1 , Ndv
            a1 = a1 + Df(i)**2
         ENDDO
         IF ( a1<1.0E-20 ) a1 = 1.0E-20
         a1 = sqrt(a1)
         a1 = 1./a1
         DO i = 1 , Ndv
            A(i,nac1) = a1*Df(i)
         ENDDO
!     BUILD C VECTOR.
         IF ( Nvc>0 ) THEN
!     ------------------------------------------------------------------
!                   BUILD C FOR MODIFIED METHOD
!     ------------------------------------------------------------------
            ndb = Nac
            A(ndv1,nac1) = -Phi
!     ------------------------------------------------------------------
!           SCALE THETA'S SO THAT MAXIMUM THETA IS UNITY
!     ------------------------------------------------------------------
            IF ( thmax>0.00001 ) thmax = 1./thmax
            DO i = 1 , ndb
               A(ndv1,i) = A(ndv1,i)*thmax
            ENDDO
            DO i = 1 , ndb
               C(i) = 0.
               DO j = 1 , ndv1
                  C(i) = C(i) + A(j,i)*A(j,nac1)
               ENDDO
            ENDDO
         ELSE
!     ------------------------------------------------------------------
!                 BUILD C FOR CLASSICAL METHOD
!     ------------------------------------------------------------------
            ndb = nac1
            A(ndv1,ndb) = 1.
            DO i = 1 , ndb
               C(i) = -A(ndv1,i)
            ENDDO
         ENDIF
!     ------------------------------------------------------------------
!                      BUILD B MATRIX
!     ------------------------------------------------------------------
         DO i = 1 , ndb
            DO j = 1 , ndb
               B(i,j) = 0.
               DO k = 1 , ndv1
                  B(i,j) = B(i,j) - A(k,i)*A(k,j)
               ENDDO
            ENDDO
         ENDDO
!     ------------------------------------------------------------------
!                    SOLVE SPECIAL L. P. PROBLEM
!     ------------------------------------------------------------------
         CALL cnmn08(ndb,ner,C,Ms1,B,N3,N4,N5)
         IF ( Iprint>1 .AND. ner>0 ) WRITE (Iout,99002)
!     ------------------------------------------------------------------
!                           FORMATS
!     ------------------------------------------------------------------
!
!
99002    FORMAT (//5X,'* * DIRECTION FINDING PROCESS DID NOT CONVERGE'/5X,'* * S-VECTOR MAY NOT BE VALID')
!     CALCULATE RESULTING DIRECTION VECTOR, S.
         Slope = 0.
!     ------------------------------------------------------------------
!                  USABLE-FEASIBLE DIRECTION
!     ------------------------------------------------------------------
         DO i = 1 , Ndv
            s1 = 0.
            IF ( Nvc>0 ) s1 = -A(i,nac1)
            DO j = 1 , ndb
               s1 = s1 - A(i,j)*C(j)
            ENDDO
            Slope = Slope + s1*Df(i)
            S(i) = s1
         ENDDO
         S(ndv1) = 1.
         IF ( Nvc>0 ) S(ndv1) = -A(ndv1,nac1)
         DO j = 1 , ndb
            S(ndv1) = S(ndv1) - A(ndv1,j)*C(j)
         ENDDO
!     ------------------------------------------------------------------
!     CHECK TO INSURE THE S-VECTOR IS FEASIBLE.
!     PROGRAM MOD-FEB, 1981, GV.
!     ------------------------------------------------------------------
         DO j = 1 , Nac
!     S DOT DEL(G).
            sg = 0.
            DO i = 1 , Ndv
               sg = sg + S(i)*A(i,j)
            ENDDO
!     IF(SG.GT.0.) GO TO 176
!
!  THIS CHANGE MADE ON 4/8/81 FOR G. VANDERPLAATS
!
            IF ( sg>1.0E-04 ) THEN
               spag_nextblock_1 = 4
               CYCLE SPAG_DispatchLoop_1
            ENDIF
!     FEASIBLE FOR THIS CONSTRAINT.  CONTINUE.
         ENDDO
!     ------------------------------------------------------------------
!                  NORMALIZE S TO MAX ABS OF UNITY
!     ------------------------------------------------------------------
         s1 = 0.
         DO i = 1 , Ndv
            a1 = abs(S(i))
            IF ( a1>s1 ) s1 = a1
         ENDDO
!     IF (S1.LT.1.0E-10) RETURN
!
!  E-10 CHANGED TO E-04 ON 1/12/81
!
         IF ( s1<1.0E-04 ) RETURN
         s1 = 1./s1
         DO i = 1 , Ndv
            S(i) = s1*S(i)
         ENDDO
         Slope = s1*Slope
         S(ndv1) = s1*S(ndv1)
         RETURN
      CASE (4)
!     S-VECTOR IS NOT FEASIBLE DUE TO SOME NUMERICAL PROBLEM.
         IF ( Iprint>=2 ) WRITE (Iout,99003)
99003    FORMAT (5X,'** CALCULATED S-VECTOR IS NOT FEASIBLE'/5X,'BETA IS SET TO ZERO')
         S(ndv1) = 0.
         Nvc = 0
         RETURN
      END SELECT
   ENDDO SPAG_DispatchLoop_1
END SUBROUTINE cnmn05
!*==cnmn06.f90 processed by SPAG 8.01MH 10:24  1 Mar 2025
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE cnmn06(X,Vlb,Vub,G,Scal,Df,S,G1,G2,Ctam,Ctbm,Slope,Alp,A2,A3,A4,F1,F2,F3,Cv1,Cv2,Cv3,Cv4,Alpca,Alpfes,Alpln,Alpmin,     &
                & Alpnc,Alpsav,Alpsid,Alptot,Isc,N1,N2,Ncal,Nvc,Icount,Igood1,Igood2,Igood3,Igood4,Ibest,Iii,Nlnc,Jgoto)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   DOUBLE PRECISION A2 , A3 , A4 , Abobj1 , Alp , alpa , alpb , Alpca , Alpfes , Alphax , Alpln , Alpmin , Alpnc , Alpsav ,        &
                  & Alpsid , Alptot , c1 , c2 , c3 , cc
   DOUBLE PRECISION Ct , Ctam , Ctbm , Ctl , Ctlmin , Ctmin , Cv1 , Cv2 , Cv3 , Cv4 , Dabfun , Delfun , Df , F1 , F2 , F3 , f4 ,   &
                  & Fdch , Fdchm , G
   DOUBLE PRECISION G1 , G2 , gi , Obj , S , Scal , si , Slope , Theta , Vlb , Vub , X , xi , xi1 , xi2 , zro
   INTEGER i , Ibest , Icndir , Icount , Igood1 , Igood2 , Igood3 , Igood4 , Igoto , ii , Iii , Info , Infog , Iout , Iprint ,     &
         & Isc , Iter , Itmax , Itrm , jbest
   INTEGER Jgoto , ksid , Linobj , N1 , N2 , Nac , Ncal , Ncon , Ndv , Nfdg , Nfeasct , Nlnc , Nscal , Nside , Nvc , nvc1
!*** End of declarations inserted by SPAG
   COMMON /cnmn1 / Delfun , Dabfun , Fdch , Fdchm , Ct , Ctmin , Ctl , Ctlmin , Alphax , Abobj1 , Theta , Obj , Ndv , Ncon ,       &
                 & Nside , Iprint , Nfdg , Nscal , Linobj , Itmax , Itrm , Icndir , Igoto , Nac , Info , Infog , Iter , Nfeasct
   COMMON /output/ Iout
   DIMENSION X(N1) , Vlb(N1) , Vub(N1) , G(N2) , Scal(N1) , Df(N1) , S(N1) , G1(N2) , G2(N2) , Isc(N2) , Ncal(2)
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
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
         IF ( Jgoto/=0 ) THEN
            IF ( Jgoto==1 ) THEN
               F2 = Obj
               IF ( Iprint>=5 ) WRITE (Iout,99005) F2
               IF ( Iprint>=5 .AND. Ncon/=0 ) THEN
                  WRITE (Iout,99006)
                  WRITE (Iout,99004) (G(i),i=1,Ncon)
               ENDIF
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
               IF ( Ncon/=0 ) THEN
                  DO i = 1 , Ncon
                     cc = Ctam
                     IF ( Isc(i)>0 ) cc = Ctbm
                     c1 = G1(i) - cc
                     c2 = G(i) - cc
                     IF ( c2>0. ) nvc1 = nvc1 + 1
                     IF ( c1>Cv1 ) Cv1 = c1
                     IF ( c2>Cv2 ) Cv2 = c2
                  ENDDO
                  IF ( Cv1>0. ) Igood1 = 1
                  IF ( Cv2>0. ) Igood2 = 1
               ENDIF
               Alp = A2
               Obj = F2
!     ------------------------------------------------------------------
!     IF F2 VIOLATES FEWER CONSTRAINTS THAN F1 BUT STILL HAS CONSTRAINT
!     VIOLATIONS RETURN
!     ------------------------------------------------------------------
               IF ( nvc1<Nvc .AND. nvc1>0 ) THEN
                  spag_nextblock_1 = 9
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
!     ------------------------------------------------------------------
!             IDENTIFY BEST OF DESIGNS F1 ANF F2
!     ------------------------------------------------------------------
!     IBEST CORRESPONDS TO MINIMUM VALUE DESIGN.
!     IF CONSTRAINTS ARE VIOLATED, IBEST CORRESPONDS TO MINIMUM
!     CONSTRAINT VIOLATION.
               IF ( Igood1==0 .AND. Igood2==0 ) THEN
!     NO CONSTRAINT VIOLATION.  PICK MINIMUM F.
                  Ibest = 1
                  IF ( F2<=F1 ) Ibest = 2
               ELSE
!     VIOLATED CONSTRAINTS.  PICK MINIMUM VIOLATION.
                  Ibest = 1
                  IF ( Cv1>=Cv2 ) Ibest = 2
               ENDIF
               ii = 1
!     ------------------------------------------------------------------
!     IF CV2 IS GREATER THAN CV1, SET MOVE LIMITS TO A2.
!     PROGRAM MOD-FEB, 1981, GV.
!     ------------------------------------------------------------------
               IF ( Cv2>Cv1 ) THEN
                  Alpln = A2
                  Alpnc = A2
                  Alpca = A2
               ENDIF
               IF ( Ncon==0 ) THEN
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
!     ------------------------------------------------------------------
!     *****                 2 - POINT INTERPOLATION                *****
!     ------------------------------------------------------------------
               Iii = 0
               DO
                  Iii = Iii + 1
                  c1 = G1(Iii)
                  c2 = G(Iii)
                  IF ( Isc(Iii)==0 ) THEN
!     ------------------------------------------------------------------
!                     NON-LINEAR CONSTRAINT
!     ------------------------------------------------------------------
                     IF ( c1<1.0E-5 .OR. c1>Ctam ) THEN
                        CALL cnmn07(ii,Alp,zro,zro,c1,A2,c2,zro,zro)
                        IF ( Alp>0. ) THEN
                           IF ( c1>Ctam .AND. Alp>Alpfes ) Alpfes = Alp
                           IF ( c1<Ct .AND. Alp<Alpnc ) Alpnc = Alp
                        ENDIF
                     ENDIF
!     ------------------------------------------------------------------
!                        LINEAR CONSTRAINT
!     ------------------------------------------------------------------
                  ELSEIF ( c1<1.0E-5 .OR. c1>Ctbm ) THEN
                     CALL cnmn07(ii,Alp,zro,zro,c1,A2,c2,zro,zro)
                     IF ( Alp>0. ) THEN
                        IF ( c1>Ctbm .AND. Alp>Alpfes ) Alpfes = Alp
                        IF ( c1<Ctl .AND. Alp<Alpln ) Alpln = Alp
                     ENDIF
                  ENDIF
                  IF ( Iii>=Ncon ) THEN
                     spag_nextblock_1 = 3
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
               ENDDO
            ELSEIF ( Jgoto==2 ) THEN
               F3 = Obj
               IF ( Iprint>=5 ) WRITE (Iout,99005) F3
               IF ( Iprint>=5 .AND. Ncon/=0 ) THEN
                  WRITE (Iout,99006)
                  WRITE (Iout,99004) (G(i),i=1,Ncon)
               ENDIF
!     ------------------------------------------------------------------
!       CALCULATE MAXIMUM CONSTRAINT VIOLATION AND PICK BEST DESIGN
!     ------------------------------------------------------------------
               Cv3 = 0.
               Igood3 = 0
               nvc1 = 0
               IF ( Ncon/=0 ) THEN
                  DO i = 1 , Ncon
                     cc = Ctam
                     IF ( Isc(i)>0 ) cc = Ctbm
                     c1 = G(i) - cc
                     IF ( c1>Cv3 ) Cv3 = c1
                     IF ( c1>0. ) nvc1 = nvc1 + 1
                  ENDDO
                  IF ( Cv3>0. ) Igood3 = 1
               ENDIF
!     DETERMINE BEST DESIGN.
               IF ( Ibest==2 ) THEN
!     CHOOSE BETWEEN F2 AND F3.
                  IF ( Igood2==0 .AND. Igood3==0 ) THEN
                     IF ( F3<=F2 ) Ibest = 3
                  ELSE
                     IF ( Cv2>=Cv3 ) Ibest = 3
                  ENDIF
!     CHOOSE BETWEEN F1 AND F3.
               ELSEIF ( Igood1==0 .AND. Igood3==0 ) THEN
                  IF ( F3<=F1 ) Ibest = 3
               ELSE
                  IF ( Cv1>=Cv3 ) Ibest = 3
               ENDIF
               Alp = A3
               Obj = F3
!     IF F3 VIOLATES FEWER CONSTRAINTS THAN F1 RETURN.
               IF ( nvc1<Nvc ) THEN
                  spag_nextblock_1 = 9
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
!     IF OBJECTIVE AND ALL CONSTRAINTS ARE LINEAR, RETURN.
               IF ( Linobj/=0 .AND. Nlnc==Ncon ) THEN
                  spag_nextblock_1 = 9
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
!     IF A3 = ALPLN AND F3 IS BOTH GOOD AND BEST RETURN.
               alpb = 1. - Alpln/A3
               IF ( (abs(alpb)<1.0E-20 .AND. Ibest==3) .AND. (Igood3==0) ) THEN
                  spag_nextblock_1 = 9
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
!     IF A3 = ALPSID AND F3 IS BEST, GO INVOKE SIDE CONSTRAINT
!     MODIFICATION.
               alpa = 1. - Alpsid/A3
               IF ( abs(alpa)<1.0E-20 .AND. Ibest==3 ) THEN
                  spag_nextblock_1 = 2
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
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
               IF ( A3>A2 .AND. Cv3>Cv2 ) THEN
                  Alpln = A3
                  Alpnc = A3
                  Alpca = A3
               ENDIF
               IF ( Ncon==0 ) THEN
                  spag_nextblock_1 = 6
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               Iii = 0
               spag_nextblock_1 = 4
               CYCLE SPAG_DispatchLoop_1
            ELSEIF ( Jgoto==3 ) THEN
               f4 = Obj
               IF ( Iprint>=5 ) WRITE (Iout,99005) f4
               IF ( Iprint>=5 .AND. Ncon/=0 ) THEN
                  WRITE (Iout,99006)
                  WRITE (Iout,99004) (G(i),i=1,Ncon)
               ENDIF
!     DETERMINE ACCAPTABILITY OF F4.
               Igood4 = 0
               Cv4 = 0.
               IF ( Ncon/=0 ) THEN
                  DO i = 1 , Ncon
                     cc = Ctam
                     IF ( Isc(i)>0 ) cc = Ctbm
                     c1 = G(i) - cc
                     IF ( c1>Cv4 ) Cv4 = c1
                  ENDDO
                  IF ( Cv4>0. ) Igood4 = 1
               ENDIF
               Alp = A4
               Obj = f4
!     ------------------------------------------------------------------
!                     DETERMINE BEST DESIGN
!     ------------------------------------------------------------------
               IF ( Ibest==2 ) THEN
!     CHOOSE BETWEEN F2 AND F4.
                  IF ( Igood2==0 .AND. Igood4==0 ) THEN
                     IF ( f4>F2 ) THEN
                        spag_nextblock_1 = 7
                        CYCLE SPAG_DispatchLoop_1
                     ENDIF
                     spag_nextblock_1 = 9
                     CYCLE SPAG_DispatchLoop_1
                  ELSE
                     IF ( Cv2<=Cv4 ) THEN
                        spag_nextblock_1 = 7
                        CYCLE SPAG_DispatchLoop_1
                     ENDIF
                     spag_nextblock_1 = 9
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
               ELSEIF ( Ibest==3 ) THEN
!     CHOOSE BETWEEN F3 AND F4.
                  IF ( Igood3==0 .AND. Igood4==0 ) THEN
                     IF ( f4>F3 ) THEN
                        spag_nextblock_1 = 8
                        CYCLE SPAG_DispatchLoop_1
                     ENDIF
                     spag_nextblock_1 = 9
                     CYCLE SPAG_DispatchLoop_1
                  ELSE
                     IF ( Cv3<=Cv4 ) THEN
                        spag_nextblock_1 = 8
                        CYCLE SPAG_DispatchLoop_1
                     ENDIF
                     spag_nextblock_1 = 9
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
               ELSE
!     CHOOSE BETWEEN F1 AND F4.
                  IF ( Igood1==0 .AND. Igood4==0 ) THEN
                     IF ( f4<=F1 ) THEN
                        spag_nextblock_1 = 9
                        CYCLE SPAG_DispatchLoop_1
                     ENDIF
                  ELSEIF ( Cv1>Cv4 ) THEN
                     spag_nextblock_1 = 9
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
!     F1 IS BEST.
                  Alptot = Alptot - A4
                  Obj = F1
                  DO i = 1 , Ndv
                     X(i) = X(i) - A4*S(i)
                  ENDDO
                  IF ( Ncon/=0 ) THEN
                     DO i = 1 , Ncon
                        G(i) = G1(i)
                     ENDDO
                  ENDIF
                  spag_nextblock_1 = 9
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
            ENDIF
         ENDIF
         IF ( Iprint>=5 ) WRITE (Iout,99002)
99002    FORMAT (/////'* * * CONSTRAINED ONE-DIMENSIONAL SEARCH INFORMATION * * *')
         Alpsav = Alp
         Icount = 0
         Alptot = 0.
!     TOLERANCES.
         Ctam = abs(Ctmin)
         Ctbm = abs(Ctlmin)
         spag_nextblock_1 = 2
      CASE (2)
!     PROPOSED MOVE.
!     ------------------------------------------------------------------
!     *****  BEGIN SEARCH OR IMPOSE SIDE CONSTRAINT MODIFICATION  *****
!     ------------------------------------------------------------------
         A2 = Alpsav
         Icount = Icount + 1
         Alpsid = 1.0E+20
!     INITIAL ALPHA AND OBJ.
         Alp = 0.
         F1 = Obj
         ksid = 0
         IF ( Nside/=0 ) THEN
!     ------------------------------------------------------------------
!     FIND MOVE TO SIDE CONSTRAINT AND INSURE AGAINST VIOLATION OF
!     SIDE CONSTRAINTS
!     ------------------------------------------------------------------
            DO i = 1 , Ndv
               si = S(i)
               IF ( abs(si)>1.0E-20 ) THEN
                  xi = X(i)
                  si = 1./si
                  IF ( si>0. ) THEN
!     UPPER BOUND.
                     xi2 = Vub(i)
                     xi1 = abs(xi2)
                     IF ( xi1<1. ) xi1 = 1.
!     CONSTRAINT VALUE.
                     gi = (xi-xi2)/xi1
                     IF ( gi<=-1.0E-6 ) THEN
!     PROPOSED MOVE TO UPPER BOUND.
                        alpa = (xi2-xi)*si
                        IF ( alpa<Alpsid ) Alpsid = alpa
                        CYCLE
                     ENDIF
                  ELSE
!     LOWER BOUND.
                     xi2 = Vlb(i)
                     xi1 = abs(xi2)
                     IF ( xi1<1. ) xi1 = 1.
!     CONSTRAINT VALUE.
                     gi = (xi2-xi)/xi1
                     IF ( gi<=-1.0E-6 ) THEN
!     PROPOSED MOVE TO LOWER BOUND.
                        alpa = (xi2-xi)*si
                        IF ( alpa<Alpsid ) Alpsid = alpa
                        CYCLE
                     ENDIF
                  ENDIF
!     MOVE WILL VIOLATE SIDE CONSTRAINT.  SET S(I)=0.
                  Slope = Slope - S(i)*Df(i)
                  S(i) = 0.
                  ksid = ksid + 1
               ELSE
!     ITH COMPONENT OF S IS SMALL.  SET TO ZERO.
                  S(i) = 0.
                  Slope = Slope - si*Df(i)
               ENDIF
            ENDDO
!     ALPSID IS UPPER BOUND ON ALPHA.
            IF ( A2>Alpsid ) A2 = Alpsid
         ENDIF
!     ------------------------------------------------------------------
!               CHECK ILL-CONDITIONING
!     ------------------------------------------------------------------
         IF ( ksid==Ndv .OR. Icount>10 ) THEN
            spag_nextblock_1 = 9
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         IF ( Nvc==0 .AND. Slope>0. ) THEN
            spag_nextblock_1 = 9
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         Alpfes = -1.
         Alpmin = -1.
         Alpln = 1.1*Alpsid
         Alpnc = Alpsid
         Alpca = Alpsid
         IF ( Ncon/=0 ) THEN
!     STORE CONSTRAINT VALUES IN G1.
            DO i = 1 , Ncon
               G1(i) = G(i)
            ENDDO
         ENDIF
!     ------------------------------------------------------------------
!                  MOVE A DISTANCE A2*S
!     ------------------------------------------------------------------
         Alptot = Alptot + A2
         DO i = 1 , Ndv
            X(i) = X(i) + A2*S(i)
         ENDDO
         IF ( Iprint>=5 ) THEN
            WRITE (Iout,99003) A2
            IF ( Nscal==0 ) THEN
               WRITE (Iout,99004) (X(i),i=1,Ndv)
            ELSE
               DO i = 1 , Ndv
                  G(i) = Scal(i)*X(i)
               ENDDO
               WRITE (Iout,99004) (G(i),i=1,Ndv)
            ENDIF
         ENDIF
!     ------------------------------------------------------------------
!                   UPDATE FUNCTION AND CONSTRAINT VALUES
!     ------------------------------------------------------------------
         Ncal(1) = Ncal(1) + 1
         Jgoto = 1
         RETURN
      CASE (3)
!     CALCULATE ALPHA TO MINIMIZE FUNCTION.
         IF ( Linobj<=0 .AND. Slope<0. ) CALL cnmn04(ii,Alpmin,zro,zro,F1,Slope,A2,F2,zro,zro,zro,zro)
!     ------------------------------------------------------------------
!                         PROPOSED MOVE
!     ------------------------------------------------------------------
!     MOVE AT LEAST FAR ENOUGH TO OVERCOME CONSTRAINT VIOLATIONS.
         A3 = Alpfes
!     MOVE TO MINIMIZE FUNCTION.
         IF ( Alpmin>A3 ) A3 = Alpmin
!     IF A3.LE.0, SET A3 = ALPSID.
         IF ( A3<=0. ) A3 = Alpsid
!     LIMIT MOVE TO NEW CONSTRAINT ENCOUNTER.
         IF ( A3>Alpnc ) A3 = Alpnc
         IF ( A3>Alpln ) A3 = Alpln
!     MAKE A3 NON-ZERO.
         IF ( A3<=1.0E-20 ) A3 = 1.0E-20
!     IF A3=A2=ALPSID AND F2 IS BEST, GO INVOKE SIDE CONSTRAINT
!     MODIFICATION.
         alpb = 1. - A2/A3
         alpa = 1. - Alpsid/A3
         jbest = 0
         IF ( abs(alpb)<1.0E-10 .AND. abs(alpa)<1.0E-10 ) jbest = 1
         IF ( jbest==1 .AND. Ibest==2 ) THEN
            spag_nextblock_1 = 2
            CYCLE SPAG_DispatchLoop_1
         ENDIF
!     SIDE CONSTRAINT CHECK NOT SATISFIED.
         IF ( Ncon/=0 ) THEN
!     STORE CONSTRAINT VALUES IN G2.
            DO i = 1 , Ncon
               G2(i) = G(i)
            ENDDO
         ENDIF
!     IF A3=A2, SET A3=.9*A2.
         IF ( abs(alpb)<1.0E-10 ) A3 = .9*A2
!     MOVE AT LEAST .01*A2.
         IF ( A3<(.01*A2) ) A3 = .01*A2
!     LIMIT MOVE TO 5.*A2.
         IF ( A3>(5.*A2) ) A3 = 5.*A2
!     LIMIT MOVE TO ALPSID.
         IF ( A3>Alpsid ) A3 = Alpsid
!     MOVE A DISTANCE A3*S.
         Alp = A3 - A2
         Alptot = Alptot + Alp
         DO i = 1 , Ndv
            X(i) = X(i) + Alp*S(i)
         ENDDO
         IF ( Iprint>=5 ) THEN
            WRITE (Iout,99007)
99007       FORMAT (/5X,'TWO-POINT INTERPOLATION')
            WRITE (Iout,99003) A3
            IF ( Nscal==0 ) THEN
               WRITE (Iout,99004) (X(i),i=1,Ndv)
            ELSE
               DO i = 1 , Ndv
                  G(i) = Scal(i)*X(i)
               ENDDO
               WRITE (Iout,99004) (G(i),i=1,Ndv)
            ENDIF
         ENDIF
!     ------------------------------------------------------------------
!              UPDATE FUNCTION AND CONSTRAINT VALUES
!     ------------------------------------------------------------------
         Ncal(1) = Ncal(1) + 1
         Jgoto = 2
         RETURN
      CASE (4)
         Iii = Iii + 1
         c1 = G1(Iii)
         c2 = G2(Iii)
         c3 = G(Iii)
         IF ( Isc(Iii)==0 ) THEN
!     ------------------------------------------------------------------
!                     NON-LINEAR CONSTRAINT
!     ------------------------------------------------------------------
            ii = 2
            CALL cnmn07(ii,Alp,zro,zro,c1,A2,c2,A3,c3)
            IF ( Alp>zro ) THEN
               IF ( c1<Ct .OR. c1>0. ) THEN
                  IF ( c1>Ctam .OR. c1<0. ) THEN
                     IF ( Alp>Alpfes .AND. c1>Ctam ) Alpfes = Alp
                     IF ( Alp<Alpnc .AND. c1<0. ) Alpnc = Alp
                     spag_nextblock_1 = 5
                     CYCLE SPAG_DispatchLoop_1
                  ENDIF
               ENDIF
!     ALP IS MINIMUM MOVE.  UPDATE FOR NEXT CONSTRAINT ENCOUNTER.
               alpa = Alp
               CALL cnmn07(ii,Alp,alpa,zro,c1,A2,c2,A3,c3)
               IF ( Alp<Alpca .AND. Alp>=alpa ) Alpca = Alp
            ENDIF
!     ------------------------------------------------------------------
!     LINEAR CONSTRAINT.  FIND ALPFES ONLY.  ALPLN SAME AS BEFORE.
!     ------------------------------------------------------------------
         ELSEIF ( c1>Ctbm ) THEN
            ii = 1
            CALL cnmn07(ii,Alp,zro,zro,c1,A3,c3,zro,zro)
            IF ( Alp>Alpfes ) Alpfes = Alp
         ENDIF
         spag_nextblock_1 = 5
      CASE (5)
         IF ( Iii<Ncon ) THEN
            spag_nextblock_1 = 4
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         spag_nextblock_1 = 6
      CASE (6)
         IF ( Linobj<=0 .AND. Slope<=0. ) THEN
!     ------------------------------------------------------------------
!              CALCULATE ALPHA TO MINIMIZE FUNCTION
!     ------------------------------------------------------------------
            ii = 3
            IF ( A2>A3 .AND. (Igood2==0 .AND. Ibest==2) ) ii = 2
            CALL cnmn04(ii,Alpmin,zro,zro,F1,Slope,A2,F2,A3,F3,zro,zro)
         ENDIF
!     ------------------------------------------------------------------
!                       PROPOSED MOVE
!     ------------------------------------------------------------------
!     MOVE AT LEAST ENOUGH TO OVERCOME CONSTRAINT VIOLATIONS.
         A4 = Alpfes
!     MOVE TO MINIMIZE FUNCTION.
         IF ( Alpmin>A4 ) A4 = Alpmin
!     IF A4.LE.0, SET A4 = ALPSID.
         IF ( A4<=0. ) A4 = Alpsid
!     LIMIT MOVE TO NEW CONSTRAINT ENCOUNTER.
         IF ( A4>Alpln ) A4 = Alpln
         IF ( A4>Alpnc ) A4 = Alpnc
!     LIMIT MOVE TO RE-ENCOUNTER CURRENTLY ACTIVE CONSTRAINT.
         IF ( A4>Alpca ) A4 = Alpca
!     LIMIT A4 TO 5.*A3.
         IF ( A4>(5.*A3) ) A4 = 5.*A3
!     UPDATE DESIGN.
         IF ( Ibest==3 .AND. Ncon/=0 ) THEN
!     STORE CONSTRAINT VALUES IN G2.  F3 IS BEST.  F2 IS NOT.
            DO i = 1 , Ncon
               G2(i) = G(i)
            ENDDO
         ENDIF
!     IF A4=A3 AND IGOOD1=0 AND IGOOD3=1, SET A4=.9*A3.
         Alp = A4 - A3
         IF ( (Igood1==0 .AND. Igood3==1) .AND. abs(Alp)<1.0E-20 ) A4 = .9*A3
!     ------------------------------------------------------------------
!                   MOVE A DISTANCE A4*S
!     ------------------------------------------------------------------
         Alp = A4 - A3
         Alptot = Alptot + Alp
         DO i = 1 , Ndv
            X(i) = X(i) + Alp*S(i)
         ENDDO
         IF ( Iprint>=5 ) THEN
            WRITE (Iout,99001)
!     ------------------------------------------------------------------
!                                  FORMATS
!     ------------------------------------------------------------------
!
!
99001       FORMAT (/5X,'THREE-POINT INTERPOLATION')
            WRITE (Iout,99003) A4
            IF ( Nscal==0 ) THEN
               WRITE (Iout,99004) (X(i),i=1,Ndv)
            ELSE
               DO i = 1 , Ndv
                  G(i) = Scal(i)*X(i)
               ENDDO
               WRITE (Iout,99004) (G(i),i=1,Ndv)
            ENDIF
         ENDIF
!     ------------------------------------------------------------------
!              UPDATE FUNCTION AND CONSTRAINT VALUES
!     ------------------------------------------------------------------
         Ncal(1) = Ncal(1) + 1
         Jgoto = 3
         RETURN
      CASE (7)
!     F2 IS BEST.
         Obj = F2
         A2 = A4 - A2
         Alptot = Alptot - A2
         DO i = 1 , Ndv
            X(i) = X(i) - A2*S(i)
         ENDDO
         IF ( Ncon/=0 ) THEN
            DO i = 1 , Ncon
               G(i) = G2(i)
            ENDDO
         ENDIF
         spag_nextblock_1 = 9
         CYCLE SPAG_DispatchLoop_1
      CASE (8)
!     F3 IS BEST.
         Obj = F3
         A3 = A4 - A3
         Alptot = Alptot - A3
         DO i = 1 , Ndv
            X(i) = X(i) - A3*S(i)
         ENDDO
         IF ( Ncon/=0 ) THEN
            DO i = 1 , Ncon
               G(i) = G2(i)
            ENDDO
         ENDIF
         spag_nextblock_1 = 9
      CASE (9)
         Alp = Alptot
         IF ( Iprint>=5 ) WRITE (Iout,99008)
99008    FORMAT (/5X,'* * * END OF ONE-DIMENSIONAL SEARCH')
         Jgoto = 0
         RETURN
      END SELECT
   ENDDO SPAG_DispatchLoop_1
99003 FORMAT (//5X,'PROPOSED DESIGN'/5X,'ALPHA =',E12.5/5X,'X-VECTOR')
99004 FORMAT (1X,8E12.4)
99005 FORMAT (/5X,'OBJ =',E13.5)
99006 FORMAT (/5X,'CONSTRAINT VALUES')
END SUBROUTINE cnmn06
!*==cnmn07.f90 processed by SPAG 8.01MH 10:24  1 Mar 2025
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE cnmn07(Ii,Xbar,Eps,X1,Y1,X2,Y2,X3,Y3)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   DOUBLE PRECISION aa , bac , bb , cc , dy , Eps , qq , X1 , X2 , x21 , X3 , x31 , x32 , xb2 , Xbar , xbar1 , Y1 , Y2 , Y3 , yy
   INTEGER Ii , jj
!*** End of declarations inserted by SPAG
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
   x21 = X2 - X1
   IF ( abs(x21)<1.0E-20 ) RETURN
   IF ( Ii==2 ) THEN
!     ------------------------------------------------------------------
!                 II=2: 3-POINT QUADRATIC INTERPOLATION
!     ------------------------------------------------------------------
      jj = 1
      x31 = X3 - X1
      x32 = X3 - X2
      qq = x21*x31*x32
      IF ( abs(qq)<1.0E-20 ) RETURN
      aa = (Y1*x32-Y2*x31+Y3*x21)/qq
      IF ( abs(aa)>=1.0E-20 ) THEN
         bb = (Y2-Y1)/x21 - aa*(X1+X2)
         cc = Y1 - X1*(aa*X1+bb)
         bac = bb*bb - 4.*aa*cc
         IF ( bac>=0. ) THEN
            bac = sqrt(bac)
            aa = .5/aa
            Xbar = aa*(bac-bb)
            xb2 = -aa*(bac+bb)
            IF ( Xbar<Eps ) Xbar = xb2
            IF ( xb2<Xbar .AND. xb2>Eps ) Xbar = xb2
            IF ( Xbar<Eps ) Xbar = xbar1
            RETURN
         ENDIF
      ENDIF
   ENDIF
!
!     ------------------------------------------------------------------
!                  II=1: 2-POINT LINEAR INTERPOLATION
!     ------------------------------------------------------------------
   Ii = 1
   yy = Y1*Y2
   IF ( jj/=0 .AND. yy>=0. ) THEN
!     INTERPOLATE BETWEEN X2 AND X3.
      dy = Y3 - Y2
      IF ( abs(dy)>=1.0E-20 ) THEN
         Xbar = X2 + Y2*(X2-X3)/dy
         IF ( Xbar<Eps ) Xbar = xbar1
         RETURN
      ENDIF
   ENDIF
   dy = Y2 - Y1
!     INTERPOLATE BETWEEN X1 AND X2.
   IF ( abs(dy)<1.0E-20 ) RETURN
   Xbar = X1 + Y1*(X1-X2)/dy
   IF ( Xbar<Eps ) Xbar = xbar1
   RETURN
END SUBROUTINE cnmn07
!*==cnmn08.f90 processed by SPAG 8.01MH 10:24  1 Mar 2025
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE cnmn08(Ndb,Ner,C,Ms1,B,N3,N4,N5)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   DOUBLE PRECISION B , bb , bb1 , bi , C , c1 , cb , cbmax , cbmin , eps
   INTEGER i , ichk , iter1 , j , jj , kk , m2 , Ms1 , N3 , N4 , N5 , Ndb , Ner , nmax
!*** End of declarations inserted by SPAG
   DIMENSION C(N4) , B(N3,N3) , Ms1(N5)
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
   eps = -1.0E+10
   cbmin = 0.
   DO i = 1 , Ndb
      bi = B(i,i)
      cbmax = 0.
      IF ( bi<-1.0E-6 ) cbmax = C(i)/bi
      IF ( bi>eps ) eps = bi
      IF ( cbmax>cbmin ) cbmin = cbmax
      Ms1(i) = 0
   ENDDO
   eps = .0001*eps
!     IF (EPS.LT.-1.0E-10) EPS=-1.0E-10
!
!  E-10 CHANGED TO E-03 ON 1/12/81
!
   IF ( eps<-1.0E-03 ) eps = -1.0E-03
   IF ( eps>-.0001 ) eps = -.0001
   cbmin = cbmin*1.0E-6
!     IF (CBMIN.LT.1.0E-10) CBMIN=1.0E-10
!
!  E-10 CHANGED TO E-05 ON 1/12/81
!
   IF ( cbmin<1.0E-05 ) cbmin = 1.0E-05
   iter1 = 0
   nmax = 5*Ndb
   SPAG_Loop_1_1: DO
!     ------------------------------------------------------------------
!     **********             BEGIN NEW ITERATION              **********
!     ------------------------------------------------------------------
      iter1 = iter1 + 1
      IF ( iter1>nmax ) RETURN
!     FIND MAX. C(I)/B(I,I) FOR I=1,NDB.
      cbmax = .9*cbmin
      ichk = 0
      DO i = 1 , Ndb
         c1 = C(i)
         bi = B(i,i)
!     IF (BI.GT.EPS.OR.C1.GT.0.) GO TO 30
         IF ( bi<=eps .AND. c1<=-1.0E-05 ) THEN
!
!  0. CHANGED TO -1.0E-05 ON 1/12/81
!
            cb = c1/bi
            IF ( cb>cbmax ) THEN
               ichk = i
               cbmax = cb
            ENDIF
         ENDIF
      ENDDO
      IF ( cbmax<cbmin ) EXIT SPAG_Loop_1_1
      IF ( ichk==0 ) EXIT SPAG_Loop_1_1
!     UPDATE VECTOR MS1.
      jj = ichk
      IF ( Ms1(jj)==0 ) jj = ichk + Ndb
      kk = jj + Ndb
      IF ( kk>m2 ) kk = jj - Ndb
      Ms1(kk) = ichk
      Ms1(jj) = 0
!     ------------------------------------------------------------------
!                     PIVOT OF B(ICHK,ICHK)
!     ------------------------------------------------------------------
      bb = 1./B(ichk,ichk)
      DO j = 1 , Ndb
         B(ichk,j) = bb*B(ichk,j)
      ENDDO
      C(ichk) = cbmax
      B(ichk,ichk) = bb
!     ELIMINATE COEFICIENTS ON VARIABLE ENTERING BASIS AND STORE
!     COEFICIENTS ON VARIABLE LEAVING BASIS IN THEIR PLACE.
      DO i = 1 , Ndb
         IF ( i/=ichk ) THEN
            bb1 = B(i,ichk)
            B(i,ichk) = 0.
            DO j = 1 , Ndb
               B(i,j) = B(i,j) - bb1*B(ichk,j)
            ENDDO
            C(i) = C(i) - bb1*cbmax
         ENDIF
      ENDDO
   ENDDO SPAG_Loop_1_1
   Ner = 0
!     ------------------------------------------------------------------
!     STORE ONLY COMPONENTS OF U-VECTOR IN 'C'.  USE B(I,1) FOR
!     TEMPORARY STORAGE
!     ------------------------------------------------------------------
   DO i = 1 , Ndb
      B(i,1) = C(i)
   ENDDO
   DO i = 1 , Ndb
      C(i) = 0.
      j = Ms1(i)
      IF ( j>0 ) C(i) = B(j,1)
      IF ( C(i)<0. ) C(i) = 0.
   ENDDO
END SUBROUTINE cnmn08
!*==cnmn09.f90 processed by SPAG 8.01MH 10:24  1 Mar 2025
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE cnmn09(cnmnfun,cnmngrd,X,G,Ic,Df,A,N1,N2,N3,N4,N5)

   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   DOUBLE PRECISION A , Abobj1 , Alphax , Aobj , Ct , Ctl , Ctlmin , Ctmin , Dabfun , Delfun , Df , Fdch , Fdchm , G , Obj ,       &
                  & Theta , X
   INTEGER Ic , Icndir , Igoto , Info , Infog , Iprint , Iter , Itmax , Itrm , Linobj , N1 , N2 , N3 , N4 , N5 , Nac , Ncon , Ndv ,&
         & Nfdg , Nfun
   INTEGER Ngrd , Nscal , Nside
!*** End of declarations inserted by SPAG

   COMMON /cnmn1 / Delfun , Dabfun , Fdch , Fdchm , Ct , Ctmin , Ctl , Ctlmin , Alphax , Abobj1 , Theta , Obj , Ndv , Ncon ,       &
                 & Nside , Iprint , Nfdg , Nscal , Linobj , Itmax , Itrm , Icndir , Igoto , Nac , Info , Infog , Iter

   COMMON /varable/ Aobj
   COMMON /fevals/ Nfun , Ngrd

   DIMENSION X(N1) , G(N2) , Ic(N2) , Df(N1) , A(N1,N2)

   EXTERNAL cnmnfun , cnmngrd

!
!
   IF ( Info>=2 ) THEN

!
!
!	GRADIENT INFORMATION
!
      CALL cnmngrd(N1,N2,X,Aobj,G,Ct,Df,A,Ic,Nac)
      Ngrd = Ngrd + 1
   ELSE
!
!  OBJECTIVE FUNCTION & CONSTRAINTS
!
      CALL cnmnfun(N1,N2,X,Aobj,G)
      Nfun = Nfun + 1
   ENDIF


END SUBROUTINE cnmn09
!*==conmin.f90 processed by SPAG 8.01MH 10:24  1 Mar 2025
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE conmin(Ndv_,Ncon_,X_,Vlb_,Vub_,Obj_,G_,N1,N2,N3,N4,N5,Iprint_,Iout_,Ifile,Itmax_,Delfun_,Dabfun_,Itrm_,Nfeasct_,Nfdg_,  &
                & Nfun_,Ngrd_,Cnmnfun,Cnmngrd)

   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   DOUBLE PRECISION a , Abobj1 , Alphax , Aobj , b , c , Ct , Ctl , Ctlmin , Ctmin , Dabfun , Dabfun_ ,        &
                  & Delfun , Delfun_ , df , Fdch , Fdchm , g
   DOUBLE PRECISION g1 , g2 , G_ , Obj , Obj_ , phi , s , scal , Theta , vlb , Vlb_ , vub , Vub_ , x , X_
   INTEGER i , ic , Icndir , Igoto , Info , Infog , Iout , Iout_ , Iprint , Iprint_ , isc , Iter , Itmax , Itmax_ , Itrm , Itrm_ , &
         & j , Linobj , loopcnt , ms1
   INTEGER N1 , N2 , N3 , N4 , N5 , Nac , Ncon , Ncon_ , Ndv , Ndv_ , Nfdg , Nfdg_ , Nfeasct , Nfeasct_ , Nfun , Nfun_ , Ngrd ,    &
         & Ngrd_ , nlim , Nscal
   INTEGER Nside
!*** End of declarations inserted by SPAG

   EXTERNAL Cnmnfun , Cnmngrd

   CHARACTER*(*) Ifile

   DIMENSION X_(Ndv_) , Vlb_(Ndv_) , Vub_(Ndv_) , G_(Ncon_)

   DIMENSION x(N1) , vlb(N1) , vub(N1) , scal(N1) , s(N1) , df(N1)
   DIMENSION g(N2) , g1(N2) , g2(N2) , isc(N2)
   DIMENSION ic(N3) , b(N3,N3)
   DIMENSION c(N4)
   DIMENSION ms1(N5)
   DIMENSION a(N1,N3)

   COMMON /cnmn1 / Delfun , Dabfun , Fdch , Fdchm , Ct , Ctmin , Ctl , Ctlmin , Alphax , Abobj1 , Theta , Obj , Ndv , Ncon ,       &
                 & Nside , Iprint , Nfdg , Nscal , Linobj , Itmax , Itrm , Icndir , Igoto , Nac , Info , Infog , Iter , Nfeasct

   COMMON /varable/ Aobj
   COMMON /output/ Iout
   COMMON /fevals/ Nfun , Ngrd
!
!  INITIALIZE
!
   Infog = 0
   Info = 0
   Ndv = Ndv_
   Ncon = Ncon_
   DO i = 1 , Ndv
      x(i) = X_(i)
      vlb(i) = Vlb_(i)
      vub(i) = Vub_(i)
   ENDDO
   Iprint = Iprint_
   Itmax = Itmax_
   DO j = 1 , Ncon
      isc(j) = 0
   ENDDO
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
   IF ( Iprint/=0 ) OPEN (UNIT=Iout,FILE=Ifile(1:len_trim(Ifile)),STATUS='UNKNOWN')
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
   SPAG_Loop_1_1: DO i = 1 , nlim
!
      loopcnt = i
!
!       CALL THE OPTIMIZATION ROUTINE CONMIN
!
      CALL cnmn00(x,vlb,vub,g,scal,df,a,s,g1,g2,b,c,isc,ic,ms1,N1,N2,N3,N4,N5)
!
!     CHECK TERMINATION CRITERIA
!
      IF ( Igoto==0 ) loopcnt = -999
!
!       ANALYSIS MODULE
!
      CALL cnmn09(Cnmnfun,Cnmngrd,x,g,ic,df,a,N1,N2,N3,N4,N5)
      Obj = Aobj
!
      IF ( Igoto==0 ) EXIT SPAG_Loop_1_1
   ENDDO SPAG_Loop_1_1
!
!
!  PRINT FINAL RESULTS
!
   IF ( Iprint/=0 ) THEN
      WRITE (6,99001) Nfun - 1

!  ------------------------------------------------------------------
!  FORMATS
!  ------------------------------------------------------------------
99001 FORMAT (//8X,'NUMBER OF FUNC-CALLS:  NFUN =',I5)
      WRITE (6,99002) Ngrd
99002 FORMAT (8X,'NUMBER OF GRAD-CALLS:  NGRD =',I5)
   ENDIF
!
!
!     OUTPUT HANDLING
!
   DO i = 1 , Ndv
      X_(i) = x(i)
   ENDDO
   Obj_ = Obj
   DO j = 1 , Ncon
      G_(j) = g(j)
   ENDDO
   Nfun_ = Nfun - 1
   Ngrd_ = Ngrd

   RETURN

END SUBROUTINE conmin

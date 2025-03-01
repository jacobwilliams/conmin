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

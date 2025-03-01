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

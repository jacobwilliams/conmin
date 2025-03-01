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

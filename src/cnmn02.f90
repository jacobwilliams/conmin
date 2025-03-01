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

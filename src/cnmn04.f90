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

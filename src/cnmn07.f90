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

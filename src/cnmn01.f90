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

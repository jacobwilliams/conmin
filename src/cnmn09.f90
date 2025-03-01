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

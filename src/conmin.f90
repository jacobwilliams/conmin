!*==conmin.f90 processed by SPAG 8.01MH 10:24  1 Mar 2025
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
SUBROUTINE conmin(Ndv_,Ncon_,X_,Vlb_,Vub_,Obj_,G_,N1,N2,N3,N4,N5,Iprint_,Iout_,Ifile,Itmax_,Delfun_,Dabfun_,Itrm_,Nfeasct_,Nfdg_,  &
                & Nfun_,Ngrd_,Cnmnfun,Cnmngrd)
 
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   DOUBLE PRECISION a , Abobj1 , Alphax , Aobj , b , c , Cnmnfun , Cnmngrd , Ct , Ctl , Ctlmin , Ctmin , Dabfun , Dabfun_ ,        &
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

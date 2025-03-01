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

*CMZ :  8.14/00 01/07/98  20.38.08  by  Stephen Burke
*CMZU:  8.13/00 24/04/98  16.03.21  by  C.Niebuhr
*CMZ :  8.04/00 27/06/96  20.27.36  by  Stephen Burke
*CMZ :  7.02/00 12/09/95  16.37.24  by  Gaby Raedel
*-- Author :  Volker Blobel
*
* The descriptions inserted below will appear asis with the
* initial asterisk removed in the WWW writeup. HTML commands
* (links etc.) may be included in the description if desired.
* The regions can be expanded to as many comment lines as
* required up to a maximum of 100 from *HTMLP to *HTMLE inclusive.
*
*HTMLP     : Describe the Purpose of the routine
*
*
*
*HTMLI     : Describe the Input variables to the routine
*
*
*
*HTMLO     : Describe the Output of the routine
*
*
*
*HTMLE     : Terminates the HTML documentation
*
      SUBROUTINE CIRPARW(X,Y,W,N,TR,CDF)
*
*     Circle parameter from N points (X,Y)
*
*     Input:
*            X(.), Y(.) = arrays of N data points in a plane
*            W(.)       = arrays of N weights (=1/sigma**2)
*
*     Result:
*     TR(1) =RNV = inverse radius
*     TR(2) =DCA = distance of closest approach to 0.0, 0.0
*     TR(3) =PHI = phi angle, all in H1 convention
*            CDF = chi square, divided by (N-3)
*
*-- Author :  Volker Blobel, Algorithm by Veikko Karimaeki
*     Veikko Karimaeki, Fast Code to fit circular Arcs, University of
*     Helsinki Report HU-SEFT-1991-10

      REAL X(*),Y(*),W(*),TR(3)
      DOUBLE PRECISION SW,SX,SY,SXX,SXY,SYY,SXR2,SYR2,SR2,SR4
      DOUBLE PRECISION CX,CY,CXX,CXY,CYY,CXR2,CYR2,CR2,CR4,R2
      DOUBLE PRECISION Q1,Q2,FKA,DLT,SQT,RNV,PHI,DCA,COSPHI,SINPHI
      DOUBLE PRECISION U
      DOUBLE PRECISION ZERO(10)
      EQUIVALENCE (ZERO( 1),SW  ),(ZERO( 2),SX  ),(ZERO( 3),SY ),
     +            (ZERO( 4),SXX ),(ZERO( 5),SXY ),(ZERO( 6),SYY),
     +            (ZERO( 7),SXR2),(ZERO( 8),SYR2),(ZERO( 9),SR2),
     +            (ZERO(10),SR4 )
*     ...
      CDF=0.0
*     ------------------------------------------------------------------
*     fit to shifted reference point (to reduce round-off errors)
*     ------------------------------------------------------------------
      IF(N.LT.3) GOTO 100
*     zero all sums
      DO 10 I=1,10
   10 ZERO(I)=0.0D0
*     shift to first point
      DO 20 I=1,N
      SW  =SW  +DBLE(W(I))
      SX  =SX  +DBLE(W(I))*(DBLE(X(I))-DBLE(X(1)))
   20 SY  =SY  +DBLE(W(I))*(DBLE(Y(I))-DBLE(Y(1)))
      CX=SX/SW
      CY=SY/SW
*     form sums ...
      DO 30 I=1,N
*     ... using shifted coordinates
      XI=DBLE(X(I))-DBLE(X(1))
      YI=DBLE(Y(I))-DBLE(Y(1))
      R2=XI*XI+YI*YI
      SXX =SXX +DBLE(W(I))*(XI-CX)*(XI-CX)
      SXY =SXY +DBLE(W(I))*(XI-CX)*(YI-CY)
      SYY =SYY +DBLE(W(I))*(YI-CY)*(YI-CY)
      SXR2=SXR2+DBLE(W(I))*XI*R2
      SYR2=SYR2+DBLE(W(I))*YI*R2
      SR2 =SR2 +DBLE(W(I))*R2
   30 SR4 =SR4 +DBLE(W(I))*R2*R2
*     calculate averages and moments
      CXX =SXX/SW
      CXY =SXY/SW
      CYY =SYY/SW
      SXX =SXX+SW*CX*CX
      SXY =SXY+SW*CX*CY
      SYY =SYY+SW*CY*CY
      CXR2=(SXR2-CX*SR2)/SW
      CYR2=(SYR2-CY*SR2)/SW
      CR2 =SR2/SW
      CR4 =(SR4-SR2*CR2)/SW
*     calculate circle parameter
      Q1=CR4*CXY-CXR2*CYR2
      Q2=CR4*(CXX-CYY)-CXR2*CXR2+CYR2*CYR2
      PHI=0.5D0*DATAN2(2.0D0*Q1,Q2)
          SINPHI=DSIN(PHI)
          COSPHI=DCOS(PHI)
*     compare PHI with initial track direction
C     IF(COSPHI*CX+SINPHI*CY.LT.0.0) THEN
      IF(COSPHI*(DBLE(X(1))+CX)+SINPHI*(DBLE(Y(1))+CY).LT.0.0) THEN
*         reverse direction
          IF(PHI.LT.0.0) THEN
             PHI=PHI+3.141592654
          ELSE
             PHI=PHI-3.141592654
          END IF
          COSPHI=-COSPHI
          SINPHI=-SINPHI
      END IF
      FKA=(SINPHI*CXR2-COSPHI*CYR2)/CR4
      DLT=-FKA*CR2+SINPHI*CX-COSPHI*CY
          SQT=DSQRT(1.0D0-4.0D0*DLT*FKA)
      RNV=2.0D0*FKA/SQT
      DCA=2.0D0*DLT/(1.0D0+SQT)
      U =1.0D0+RNV*DCA
*     ------------------------------------------------------------------
*     calculate estimate for standard deviation
*     ------------------------------------------------------------------
      IF(N.GT.3) THEN
*     estimate for chi square/(n-3) (note: fit is to unweighted data!)..
         CDF=SW*(1.0D0+RNV*DCA)**2
     +     *(SINPHI*SINPHI*CXX-2.0D0*SINPHI*COSPHI*CXY+COSPHI*COSPHI*CYY
     +       -FKA*FKA*CR4)/DFLOAT(N-3)
      END IF
*     ------------------------------------------------------------------
*     shift to reference point (0,0)
*     ------------------------------------------------------------------
      DTR=DBLE(X(1))*SINPHI-DBLE(Y(1))*COSPHI+DCA
      DLN=DBLE(X(1))*COSPHI+DBLE(Y(1))*SINPHI
      SPA=2.0*DTR+RNV*(DTR*DTR+DLN*DLN)
      SPB= RNV*DBLE(X(1))+U*SINPHI
      SPC=-RNV*DBLE(Y(1))+U*COSPHI
      SPT=SQRT(1.0+RNV*SPA)
*     shifted parameters
      TR(1)=RNV
      TR(2)=SPA/(1.0+SPT)
      TR(3)=ATAN2(SPB,SPC)
*     ------------------------------------------------------------------
*     change sign of RNV
*     ------------------------------------------------------------------
      TR(1)=-TR(1)
  100 RETURN
      END






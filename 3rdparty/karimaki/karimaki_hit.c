/***************************************************************************************************
Karimaki:

perform for all reconstructed track Karimaki (circle) fit in 2D
Add s-z fit subsequently

*            X(.), Y(.) = arrays of N data points in a plane
*            W(.)       = arrays of N weights (=1/sigma**2)
*
*     Result:
*     TR(1) =RNV = inverse radius
*     TR(2) =DCA = distance of closest approach to 0.0, 0.0
*     TR(3) =PHI = phi angle, all in H1 convention
*            CDF = chi square, divided by (N-3)
*
*     assume origin (vertex) at (0,0)!!!
*
***************************************************************************************************/

#include <stdlib.h> 
#include <time.h> 
#include <math.h> 
#include <stdio.h> 

#include "basics.h"
#include "geometry.h"
#include "sim_track.h"
#include "hits_include.h"
#include "pixrec.h"
#include "karimaki.h"

#include "c_wrappedcode.h"

//#define CALCWEIGHT
#define WEIGHTMAX 10.0

#define MINRESOLUTION 1E-6  // minimum resolution for spatial position

int karimaki_hit(KARITRACK &karires, int npoints, double xp[MAXLAYER], double yp[MAXLAYER],  double zp[MAXLAYER], double tres[MAXLAYER], double zres[MAXLAYER], double rres[MAXLAYER]) {

  int i;

  float ConstrX[MAXLAYER];
  float ConstrY[MAXLAYER];
  float ConstrZ[MAXLAYER];
  float ConstrS[MAXLAYER];
  float TWeight[MAXLAYER];
  float ZWeight[MAXLAYER];

  float TCh2dF=0;
  float ZCh2dF=0;

  int NumConstr=0;
  float TPar[3]; 
  float ZPar[2]; 
  float rad;
  float r_hit;
  float phi_hit;
  float arg;
  float alpha;
  float theta,r3d;

  int fitszw (float *,float *,float *,int ,float *,float *);

  /* number hits */
  NumConstr=0;

  if (npoints>=MAXLAYER) {
    printf("karimaki_hit: EROOR! npoints>=MAXLAYER \n");
    exit(0);
  }

  /* hit positions */
  for (i=0;i<npoints;i++) {

    ConstrX[i]=xp[i];
    ConstrY[i]=yp[i];
    ConstrZ[i]=zp[i];

    TWeight[i]=1./tres[i];  // standard weight
    printf("i=%d: x=%f   y=%f   z=%f TWeight[i]=%f\n",i, ConstrX[i],ConstrY[i],ConstrZ[i],TWeight[i]);
    NumConstr++;
  }

  // 2d track fit in transverse plane
  c_cirparw (ConstrX,ConstrY,TWeight,&NumConstr,TPar,&TCh2dF);
  printf("Karimaki transverse track parameter: %f %f %f  chi2norm=%f,\n",TPar[0],TPar[1],TPar[2],TCh2dF);

  // radius parameter
  rad=1.0/TPar[0];

  printf("karimaki_hit result: rad=%f, chi2=%f \n",rad,TCh2dF );

  //TODO arg mit sinus näherung, dann ist bogenlänge gleich radius

  // calculate arc length
  for (i=0;i<npoints;i++) {
    r_hit=sqrt(pow(ConstrX[i],2)+pow(ConstrY[i],2));
    phi_hit=atan2(ConstrY[NumConstr],ConstrX[NumConstr]);
    arg=0.5*r_hit/rad;

// hit reached ?
    if (fabs(arg)<1) {
      alpha=asin(arg);
      ConstrS[i]=2*rad*alpha;  // arc length in 2d --> TODO radius einsetzen falls arg klein
      ZWeight[i]=1.0/zres[i];
    } else {
      ConstrS[i]=0.0;
      ZWeight[i]=0.0;  // standard weight
    }
  }


  // first 2d track fit in longitudinal plane
  fitszw (ConstrS,ConstrZ,ZWeight,NumConstr,ZPar,&ZCh2dF);
  printf("Karimaki longitudinal track parameter: %f %f  chi2norm=%f,\n",ZPar[0],ZPar[1],ZCh2dF);

  // redo weight calculation of used hits
  theta=atan2(1.0,ZPar[1]);

  for (i=0;i<npoints;i++) {
    if (ZWeight[i]>0) {
      ZWeight[i]=1.0/sqrt(pow(zres[i]/sin(theta),2)+pow(rres[i]/cos(theta),2));
    }
  }

  // redo fit with correct weights
  fitszw (ConstrS,ConstrZ,ZWeight,NumConstr,ZPar,&ZCh2dF);
  printf("Karimaki longitudinal track parameter: %f %f  chi2norm=%f,\n",ZPar[0],ZPar[1],ZCh2dF);

  // store results
  theta=atan2(1.0,ZPar[1]);
  r3d=1.0/TPar[0]/sin(theta);

//    karires.rad=TPar[0];
    karires.rad=rad;
    karires.r3d=r3d;
    karires.dca=TPar[1];
    karires.phi=TPar[2];
    karires.z0=ZPar[0];
    karires.theta=theta;
    karires.tchi2n=TCh2dF;
    karires.zchi2n=ZCh2dF;

//  karires.rad=TPar[1];


  printf("karimaki_hit result: rad=%f,r3d=%f  theta=%f  chi2=%f %f\n",karires.rad,karires.r3d,karires.theta,TCh2dF,ZCh2dF );

  return 0;
}

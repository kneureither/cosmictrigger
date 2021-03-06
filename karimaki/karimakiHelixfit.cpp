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

#include "karimakiHelixfit.h"
#include "fitszw.h"
#include "wrappers.h"

const bool DEBUG_CONSOLE_LOG = false;

//#define CALCWEIGHT
#define WEIGHTMAX 10.0

#define MINRESOLUTION 1E-6  // minimum resolution for spatial position

int karimakiHelixfit(KARITRACK &karires, int npoints, double *xp, double *yp, double *zp, double *phip, double *thetas, double *tres, double *zres, double *rres) {

  if(DEBUG_CONSOLE_LOG) printf("\n\tKARIMAKI_HIT\n\t-----------\n\n");
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
  float alpha;
  float d_hit;
  float theta,r3d, z_off;

  /* number hits */
  NumConstr=0;

  if (npoints>=MAXLAYER) {
    printf("karimakiHelixfit: EROOR! npoints>=MAXLAYER \n");
    exit(0);
  }

  /* hit positions */
  for (i=0;i<npoints;i++) {

    ConstrX[i]=xp[i];
    ConstrY[i]=yp[i];
    ConstrZ[i]=zp[i];

    //Calculate Weights for the
    double phi_hit = atan2(yp[i], xp[i]);
    double dphi =  phi_hit - phip[i];
    TWeight[i]=1./pow(sin(dphi)*tres[i], 2);

    if(DEBUG_CONSOLE_LOG) printf("\tphi_hit=%f, phi_track=%f, dphi=%f, sin(dphi)=%f\t", phi_hit, phip[i], dphi, fabs(sin(dphi)));
    if(DEBUG_CONSOLE_LOG) printf("\ti=%d: x=%f   y=%f   z=%f TWeight[i]=%f\n",i, ConstrX[i],ConstrY[i],ConstrZ[i],TWeight[i]);
    NumConstr++;
  }

  // 2d track fit in transverse plane
  c_cirparw (ConstrX,ConstrY,TWeight,&NumConstr,TPar,&TCh2dF);
  if(DEBUG_CONSOLE_LOG) printf("\tKarimaki transverse track parameter: %f %f %f  chi2norm=%f,\n",TPar[0],TPar[1],TPar[2],TCh2dF);

  // radius parameter
  rad=1.0/TPar[0];
//  printf("karimakiHelixfit result: rad=%f, chi2=%f \n",rad,TCh2dF );


//// Do the helix fit in 3rd dimension

  for (i=0;i<npoints;i++) {
      //transform to first hit at (0,0) frame
      ConstrX[i]=xp[i] - xp[0];
      ConstrY[i]=yp[i] - yp[0];

      //get arc length in 1st ref frame
      d_hit = sqrt(pow(ConstrX[i], 2) + pow(ConstrY[i], 2));
      alpha = 2*asin(0.5*d_hit / rad);

      phi_hit=atan2(ConstrY[NumConstr],ConstrX[NumConstr]);

      ConstrS[i] = rad*alpha; //corrected ConstrS
      ZWeight[i]=1.0/zres[i];

// // hit reached ? -- OBSOLETE! TODO
//    if (fabs(arg)<1) {
//      ConstrS[i] = rad*alpha; //corrected ConstrS
// //        ConstrS[i]=2*rad*alpha;  // arc length in 2d
//        ZWeight[i]=1.0/zres[i];
//        if(DEBUG_CONSOLE_LOG) printf("\tZweight as in call (RMS)=%f theta from ms=%f \t alpha=%f arg=%f Constr[%d]=%f\n", ZWeight[i], thetas[i], alpha, arg, i, ConstrS[i]);
//    } else {
//      ConstrS[i]=0.0;
//      ZWeight[i]=0.0;  // standard weight
//      if(DEBUG_CONSOLE_LOG) printf("\tSomething went wrong!\n");
//    }
  }


  // first 2d track fit in longitudinal plane
  fitszw (ConstrS,ConstrZ,ZWeight,NumConstr,ZPar,&ZCh2dF);
  if(DEBUG_CONSOLE_LOG) printf("\tKarimaki longitudinal track parameter: %f %f  chi2norm=%f,\n",ZPar[0],ZPar[1],ZCh2dF);

  // redo weight calculation of used hits with fitted parameters
  theta=atan2(1.0,ZPar[1]);

  for (i=0;i<npoints;i++){
      if (ZWeight[i]>0) {
          theta = thetas[i];
          ZWeight[i]=1.0/abs(zres[i]*sin(theta));
      }
  }

  // redo fit with correct weights
  fitszw (ConstrS,ConstrZ,ZWeight,NumConstr,ZPar,&ZCh2dF);
  if(DEBUG_CONSOLE_LOG) printf("\tKarimaki longitudinal track parameter: %f %f  chi2norm=%f,\n",ZPar[0],ZPar[1],ZCh2dF);

  // store results
  theta=atan2(1.0,ZPar[1]);
  r3d=1.0/TPar[0]/sin(theta);
  z_off = sqrt(pow(xp[0], 2) + pow(yp[0], 2)) / tan(theta);


    //karires.rad=TPar[0];
    karires.rad=rad;
    karires.r3d=r3d;
    karires.dca=TPar[1];
    karires.phi=TPar[2];
    karires.z0=ZPar[0] + z_off;
    karires.theta=theta;
    karires.tchi2n=TCh2dF;
    karires.zchi2n=ZCh2dF;


  if(DEBUG_CONSOLE_LOG) printf("\tkarimakiHelixfit result: rad=%f,r3d=%f  theta=%f  chi2t=%f chi2zs=%f\n",karires.rad,karires.r3d,karires.theta,TCh2dF,ZCh2dF );

  return 0;
}

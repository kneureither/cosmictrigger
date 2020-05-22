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
* Note: this version takes onput from global variables!!!
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

#include "c_wrappedcode.h"

//#define CALCWEIGHT
#define WEIGHTMAX 10.0

#define MINRESOLUTION 1E-6  // minimum resolution for spatial position





int karimaki() {

  int i,k,itr;

  float ConstrX[MAXLAYER];
  float ConstrY[MAXLAYER];
  float ConstrZ[MAXLAYER];
  float ConstrS[MAXLAYER];
  float TWeight[MAXLAYER];
  float ZWeight[MAXLAYER];

  float TCh2dF=0;
  float ZCh2dF=0;

  int NumConstr;
  float TPar[3]; 
  float ZPar[2]; 
  float th1,rad;
  float r_hit;
  //  float phi_hit;
  float arg;
  float alpha;
  float theta,r3d;

  int fitszw (float *,float *,float *,int ,float *,float *);

  /* loop over reconstructed tracks */
  for (itr=0;itr<nrectrack;itr++) {

    /* simple estimate of weight by extrapolating track to hit radius and calculation of inclination angle*/
    th1=track_rec[itr].th1;
    rad=sin(th1)*track_rec[itr].r3d;

    /* number hits */
    NumConstr=0.0;

    /* hit positions */
    for (i=0;i<track_rec[itr].nfit;i++) {
      if (track_rec[itr].hits[i]>=0 && track_rec[itr].hits[i]<MAXHITS) {
        // layer
        k=track_rec[itr].layers[i];
        // barrel
        if (track_rec[itr].dettypes[i]==0) {
          ConstrX[NumConstr]=hit_mem_barr[track_rec[itr].hits[i]][k].x;
          ConstrY[NumConstr]=hit_mem_barr[track_rec[itr].hits[i]][k].y;
          ConstrZ[NumConstr]=hit_mem_barr[track_rec[itr].hits[i]][k].z;
          // rings
        } else if (track_rec[itr].dettypes[i]==1) {
          ConstrX[NumConstr]=hit_mem_ring[track_rec[itr].hits[i]][k].x;
          ConstrY[NumConstr]=hit_mem_ring[track_rec[itr].hits[i]][k].y;
          ConstrZ[NumConstr]=hit_mem_ring[track_rec[itr].hits[i]][k].z;
          // endcaps
        } else {
          ConstrX[NumConstr]=hit_mem_endc[track_rec[itr].hits[i]][k].x;
          ConstrY[NumConstr]=hit_mem_endc[track_rec[itr].hits[i]][k].y;
          ConstrZ[NumConstr]=hit_mem_endc[track_rec[itr].hits[i]][k].z;
        }

	/* simple estimate of weight by extrapolating track to hit radius and calculation of inclination angle*/
	r_hit=sqrt(pow(ConstrX[NumConstr],2)+pow(ConstrY[NumConstr],2));
	//	phi_hit=atan2(ConstrY[NumConstr],ConstrX[NumConstr]);
	arg=0.5*r_hit/rad;

// layer reached ?
	if (fabs(arg)<1) {
	  alpha=asin(arg);
	  ConstrS[NumConstr]=2*rad*alpha;  // arc length in 2d
	  TWeight[NumConstr]=1.0;  // standard weight
	  ZWeight[NumConstr]=1.0;  // standard weight

#ifdef CALCWEIGHT
	  TWeight[NumConstr]=1.0/fabs(cos(alpha));
	  if (TWeight[NumConstr]>WEIGHTMAX) {
	    TWeight[NumConstr]=WEIGHTMAX;
	  }
#endif
	} else {
	  ConstrS[NumConstr]=0.0;  // not defined
	  TWeight[NumConstr]=0.0;   //cannot reach this layer!
	  ZWeight[NumConstr]=0.0;   //cannot reach this layer!
	}
	//	printf("Karimaki Weight correction: hit=%d weight=%f \n",NumConstr,TWeight[NumConstr]);

	// increase hit counter
	NumConstr++;
      }
    }
    // 2d track fit in transverse plane
    c_cirparw (ConstrX,ConstrY,TWeight,&NumConstr,TPar,&TCh2dF);
//    printf("Karimaki transverse track parameter: %f %f %f  chi2norm=%f,\n",TPar[0],TPar[1],TPar[2],TCh2dF);

    // 2d track fit in longitudinal plane
    fitszw (ConstrS,ConstrZ,ZWeight,NumConstr,ZPar,&ZCh2dF);
//    printf("Karimaki longitudinal track parameter: %f %f  chi2norm=%f,\n",ZPar[0],ZPar[1],ZCh2dF);

    // store results
    theta=atan2(1.0,ZPar[1]);
    r3d=1.0/TPar[0]/sin(theta);

    track_rec[itr].kari_r3d=r3d;
    track_rec[itr].kari_dca=TPar[1];
    track_rec[itr].kari_phi=TPar[2];
    track_rec[itr].kari_z0=ZPar[0];
    track_rec[itr].kari_theta=theta;

#ifdef TRES_HIT
    track_rec[itr].kari_tchi2n=TCh2dF/pow(TRES_HIT,2);
#else
    track_rec[itr].kari_tchi2n=TCh2dF/pow(MINRESOLUTION,2);
#endif
#ifdef ZRES_HIT
    track_rec[itr].kari_zchi2n=ZCh2dF/pow(ZRES_HIT,2);
#else
    track_rec[itr].kari_zchi2n=ZCh2dF/pow(MINRESOLUTION,2);
#endif

  }
  return 0;
}



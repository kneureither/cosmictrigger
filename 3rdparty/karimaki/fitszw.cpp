#include <stdlib.h> 
#include <time.h> 
#include <math.h> 
#include <stdio.h>
#include "fitszw.h"

int fitszw (float *ConstrS,float *ConstrZ,float *Weight,int NumConstr,float *ZPar,float *ZCh2dF) {
/***************************************************************************************************
Simple linear fit with weights:  z=z0+b*s
ZPar[0]=z0
ZPar[1]=b
***************************************************************************************************/

//    if(DEBUG_CONSOLE_LOG) printf("\n\tFITSZW\n");
    int i;
    double Sum=0.0;
    double SumS=0.0;
    double SumS2=0.0;
    double SumSZ=0.0;
    double SumZ=0.0;
    double SumZ2=0.0;
    double w2;
    double denom;
    float chi2;

    for (i=0;i<NumConstr;i++) {
        w2=pow(Weight[i],2);
        Sum+=w2;
        SumS+=w2*ConstrS[i];
        SumZ+=w2*ConstrZ[i];
        SumS2+=w2*pow(ConstrS[i],2);
        SumZ2+=w2*pow(ConstrZ[i],2);
        SumSZ+=w2*ConstrS[i]*ConstrZ[i];
    }

    denom=SumS*SumS-Sum*SumS2;

    if (denom==0) {
        printf("Error in  fitszw(): denominator is zero!\n");
        exit(0);
    }

    ZPar[0]=(SumS*SumSZ-SumZ*SumS2)/denom;
    ZPar[1]=(SumS*SumZ-Sum*SumSZ)/denom;

    chi2=SumZ2+pow(ZPar[0],2)*Sum+pow(ZPar[1],2)*SumS2-2.0*ZPar[0]*SumZ-2.0*ZPar[1]*SumSZ+2.0*ZPar[0]*ZPar[1]*SumS;
    *ZCh2dF=chi2/(NumConstr-2.0);

    return 0;
}


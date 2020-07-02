#ifdef __cplusplus
extern "C" {
#endif

  void cirpar_ (float *,float *,int *,float *,float *);
  void cirparw_ (float *,float *,float *,int *,float *,float *);


#include "wrappers.h"
#include <string.h> /* strlen is need by wrappers! */

  void c_cirpar (float *ConstrX,float *ConstrY,int *NumConstr,float *TPar,float *Ch2dF)
  {
    /*    printf("Phi=%f Dca=%f Kappa=%f\n",TPar[2],TPar[1],TPar[0]); */
    cirpar_ (ConstrX,ConstrY,NumConstr,TPar,Ch2dF);
/*    printf("Phi=%f Dca=%f Kappa=%f\n",TPar[2],TPar[1],TPar[0]);   */
  }


  void c_cirparw (float *ConstrX,float *ConstrY,float *Weight,int *NumConstr,float *TPar,float *Ch2dF)
  {
    /*    printf("Phi=%f Dca=%f Kappa=%f\n",TPar[2],TPar[1],TPar[0]); */
    cirparw_ (ConstrX,ConstrY,Weight,NumConstr,TPar,Ch2dF);
/*    printf("Phi=%f Dca=%f Kappa=%f\n",TPar[2],TPar[1],TPar[0]);   */
  }


#ifdef __cplusplus
} //extern "C"
#endif

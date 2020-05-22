#ifndef C_WRAPPEDCODE_H
#define C_WRAPPEDCODE_H

#ifdef __cplusplus
extern "C" {
#endif

  void c_cirpar (float *ConstrX,float *ConstrY,int *NumConstr,float *TPar,float *CH2df);
  void c_cirparw (float *ConstrX,float *ConstrY,float *Weight,int *NumConstr,float *TPar,float *CH2df);

#ifdef __cplusplus
} //extern "C"
#endif

#endif


/* ALL values in [m] */

// vertex region
//#define ZVERTEX 0.1   // total length of vertex region
#define ZVERTEX 0.2   // total length of vertex region

// thickness of layer in radiation lengths
#define THICKNESS 0.015
//#define THICKNESS 0.01

// Pixel detector layers

// Geometry -> used for long track search
#define ZACCEPT 3.0   // geometric acceptance of tracker (in [m]]


/******** Experimental Setup ***********/

#include "geo_setup.h"

//#define SETUP_TSUNO_LOI

//#define SETUP_CENTRAL_36

#define BEAMPIPE
#define RBEAMPIPE 0.025
#define XBEAMPIPE (0.008/0.353)


const static double XBARR0={THICKNESS};   // p1
const static double XBARR1={THICKNESS};
const static double XBARR2={THICKNESS};



#if defined(ITK_EXT_SIMPLIFIED)
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 4
 const static double RLAYER2[NBARR2]   ={0.405  ,0.562  ,0.762  ,1.00};
 const static double ZMAXLAYER2[NBARR2]={1.4    ,1.4    ,1.4    ,1.4};
 const static double ZSMLAYER2[NBARR2] ={1e-3   ,1e-3   ,1e-3   ,1e-3};
 const static double TSMLAYER2[NBARR2] ={75.5e-6,75.5e-6,75.5e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif


#if defined(ITK_EXT_SIMPLIFIED_80TRI)
#define FIRSTTRIGGERLAYER 7
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.405  ,0.562  ,0.74 ,0.76 ,0.78 ,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4    ,1.4    ,1.4  ,1.4  ,1.4  ,1.4};
 const static double ZSMLAYER2[NBARR2] ={1e-3   ,1e-3   ,50e-6,50e-6,50e-6,1e-3};
 const static double TSMLAYER2[NBARR2] ={75.5e-6,75.5e-6,50e-6,50e-6,50e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif


// long pixel
#if defined(ITK_EXT_SIMPLIFIED_80TRI_AL)
#define FIRSTTRIGGERLAYER 7
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.405  ,0.562  ,0.74 ,0.76 ,0.78 ,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4    ,1.4    ,1.4  ,1.4  ,1.4  ,1.4};
 const static double ZSMLAYER2[NBARR2] ={1e-3   ,1e-3   ,100e-6,100e-6,100e-6,1e-3};
 const static double TSMLAYER2[NBARR2] ={75.5e-6,75.5e-6,50e-6,50e-6,50e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif




#if defined(ITK_EXT_SIMPLIFIED_80TRI_B)
#define FIRSTTRIGGERLAYER 7
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};

// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.405  ,0.562  ,0.744,0.760,0.776,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4    ,1.4    ,1.4  ,1.4  ,1.4  ,1.4};
 const static double ZSMLAYER2[NBARR2] ={1e-3   ,1e-3   ,50e-6,50e-6,50e-6,1e-3};
 const static double TSMLAYER2[NBARR2] ={75.5e-6,75.5e-6,50e-6,50e-6,50e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif



#if defined(ITK_EXT_SIMPLIFIED_80TRI_Z)
#define FIRSTTRIGGERLAYER 7
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.405  ,0.562  ,0.735,0.76 ,0.785,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4    ,1.4    ,1.4  ,1.4  ,1.4  ,1.4};
 const static double ZSMLAYER2[NBARR2] ={1e-3   ,1e-3   ,50e-6,50e-6,50e-6,1e-3};
 const static double TSMLAYER2[NBARR2] ={75.5e-6,75.5e-6,50e-6,50e-6,50e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif



#if defined(ITK_EXT_SIMPLIFIED_80TRI_Y)
#define FIRSTTRIGGERLAYER 7
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.405  ,0.562  ,0.730,0.760,0.790,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4    ,1.4    ,1.4  ,1.4  ,1.4  ,1.4};
 const static double ZSMLAYER2[NBARR2] ={1e-3   ,1e-3   ,50e-6,50e-6,50e-6,1e-3};
 const static double TSMLAYER2[NBARR2] ={75.5e-6,75.5e-6,50e-6,50e-6,50e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif


#if defined(ITK_EXT_SIMPLIFIED_80TRI_X)
#define FIRSTTRIGGERLAYER 7
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.405  ,0.562  ,0.725,0.760,0.795,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4    ,1.4    ,1.4  ,1.4  ,1.4  ,1.4};
 const static double ZSMLAYER2[NBARR2] ={1e-3   ,1e-3   ,50e-6,50e-6,50e-6,1e-3};
 const static double TSMLAYER2[NBARR2] ={75.5e-6,75.5e-6,50e-6,50e-6,50e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif


#if defined(ITK_EXT_SIMPLIFIED_80TRI_W)
#define FIRSTTRIGGERLAYER 7
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.405  ,0.562  ,0.720,0.760,0.800,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4    ,1.4    ,1.4  ,1.4  ,1.4  ,1.4};
 const static double ZSMLAYER2[NBARR2] ={1e-3   ,1e-3   ,50e-6,50e-6,50e-6,1e-3};
 const static double TSMLAYER2[NBARR2] ={75.5e-6,75.5e-6,50e-6,50e-6,50e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif


#if defined(ITK_EXT_SIMPLIFIED_80TRI_V)
#define FIRSTTRIGGERLAYER 7
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.405  ,0.562  ,0.715,0.760,0.805,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4    ,1.4    ,1.4  ,1.4  ,1.4  ,1.4};
 const static double ZSMLAYER2[NBARR2] ={1e-3   ,1e-3   ,50e-6,50e-6,50e-6,1e-3};
 const static double TSMLAYER2[NBARR2] ={75.5e-6,75.5e-6,50e-6,50e-6,50e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif


#if defined(ITK_EXT_SIMPLIFIED_80TRI_U)
#define FIRSTTRIGGERLAYER 7
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.405  ,0.562  ,0.710,0.760,0.810,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4    ,1.4    ,1.4  ,1.4  ,1.4  ,1.4};
 const static double ZSMLAYER2[NBARR2] ={1e-3   ,1e-3   ,50e-6,50e-6,50e-6,1e-3};
 const static double TSMLAYER2[NBARR2] ={75.5e-6,75.5e-6,50e-6,50e-6,50e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif






#if defined(ITK_EXT_SIMPLIFIED_60TRI_U)
#define FIRSTTRIGGERLAYER 6
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.405  ,0.528 ,0.578 ,0.628 ,0.762  ,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4    ,1.4   ,1.4   ,1.4   ,1.4    ,1.4};
 const static double ZSMLAYER2[NBARR2] ={1e-3   ,50e-6 ,50e-6 ,50e-6 ,1e-3   ,1e-3};
 const static double TSMLAYER2[NBARR2] ={75.5e-6,50e-6 ,50e-6 ,50e-6 ,75.5e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif

#if defined(ITK_EXT_SIMPLIFIED_60TRI_V)
#define FIRSTTRIGGERLAYER 6
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.405  ,0.533 ,0.578 ,0.623 ,0.762  ,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4    ,1.4   ,1.4   ,1.4   ,1.4    ,1.4};
 const static double ZSMLAYER2[NBARR2] ={1e-3   ,50e-6 ,50e-6 ,50e-6 ,1e-3   ,1e-3};
 const static double TSMLAYER2[NBARR2] ={75.5e-6,50e-6 ,50e-6 ,50e-6 ,75.5e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif

#if defined(ITK_EXT_SIMPLIFIED_60TRI_W)
#define FIRSTTRIGGERLAYER 6
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.405  ,0.538 ,0.578 ,0.618 ,0.762  ,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4    ,1.4   ,1.4   ,1.4   ,1.4    ,1.4};
 const static double ZSMLAYER2[NBARR2] ={1e-3   ,50e-6 ,50e-6 ,50e-6 ,1e-3   ,1e-3};
 const static double TSMLAYER2[NBARR2] ={75.5e-6,50e-6 ,50e-6 ,50e-6 ,75.5e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif

#if defined(ITK_EXT_SIMPLIFIED_60TRI_X)
#define FIRSTTRIGGERLAYER 6
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.405  ,0.543 ,0.578 ,0.613 ,0.762  ,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4    ,1.4   ,1.4   ,1.4   ,1.4    ,1.4};
 const static double ZSMLAYER2[NBARR2] ={1e-3   ,50e-6 ,50e-6 ,50e-6 ,1e-3   ,1e-3};
 const static double TSMLAYER2[NBARR2] ={75.5e-6,50e-6 ,50e-6 ,50e-6 ,75.5e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif

#if defined(ITK_EXT_SIMPLIFIED_60TRI_Y)
#define FIRSTTRIGGERLAYER 6
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.405  ,0.548 ,0.578 ,0.608 ,0.762  ,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4    ,1.4   ,1.4   ,1.4   ,1.4    ,1.4};
 const static double ZSMLAYER2[NBARR2] ={1e-3   ,50e-6 ,50e-6 ,50e-6 ,1e-3   ,1e-3};
 const static double TSMLAYER2[NBARR2] ={75.5e-6,50e-6 ,50e-6 ,50e-6 ,75.5e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif

#if defined(ITK_EXT_SIMPLIFIED_60TRI_Z)
#define FIRSTTRIGGERLAYER 6
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.405  ,0.553 ,0.578 ,0.603 ,0.762  ,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4    ,1.4   ,1.4   ,1.4   ,1.4    ,1.4};
 const static double ZSMLAYER2[NBARR2] ={1e-3   ,50e-6 ,50e-6 ,50e-6 ,1e-3   ,1e-3};
 const static double TSMLAYER2[NBARR2] ={75.5e-6,50e-6 ,50e-6 ,50e-6 ,75.5e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif



#if defined(ITK_EXT_SIMPLIFIED_60TRI_B)
#define FIRSTTRIGGERLAYER 6
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.405  ,0.562 ,0.578 ,0.594 ,0.762  ,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4    ,1.4   ,1.4   ,1.4   ,1.4    ,1.4};
 const static double ZSMLAYER2[NBARR2] ={1e-3   ,50e-6 ,50e-6 ,50e-6 ,1e-3   ,1e-3};
 const static double TSMLAYER2[NBARR2] ={75.5e-6,50e-6 ,50e-6 ,50e-6 ,75.5e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif

#if defined(ITK_EXT_SIMPLIFIED_60TRI_A)
#define FIRSTTRIGGERLAYER 6
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.405  ,0.558 ,0.578 ,0.598 ,0.762  ,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4    ,1.4   ,1.4   ,1.4   ,1.4    ,1.4};
 const static double ZSMLAYER2[NBARR2] ={1e-3   ,50e-6 ,50e-6 ,50e-6 ,1e-3   ,1e-3};
 const static double TSMLAYER2[NBARR2] ={75.5e-6,50e-6 ,50e-6 ,50e-6 ,75.5e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif



#if defined(ITK_EXT_SIMPLIFIED_40TRI_A)
#define FIRSTTRIGGERLAYER 5
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.385 ,0.405 ,0.425 ,0.562 ,0.762  ,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4   ,1.4   ,1.4   ,1.4   ,1.4    ,1.4};
 const static double ZSMLAYER2[NBARR2] ={50e-6 ,50e-6 ,50e-6 ,1e-3   ,1e-3   ,1e-3};
 const static double TSMLAYER2[NBARR2] ={50e-6 ,50e-6 ,50e-6 ,75.5e-6,75.5e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif


#if defined(ITK_EXT_SIMPLIFIED_40TRI_B)
#define FIRSTTRIGGERLAYER 5
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.393 ,0.409 ,0.425 ,0.562 ,0.762  ,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4   ,1.4   ,1.4   ,1.4   ,1.4    ,1.4};
 const static double ZSMLAYER2[NBARR2] ={50e-6 ,50e-6 ,50e-6 ,1e-3   ,1e-3   ,1e-3};
 const static double TSMLAYER2[NBARR2] ={50e-6 ,50e-6 ,50e-6 ,75.5e-6,75.5e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif



#if defined(ITK_EXT_SIMPLIFIED_40TRI_C)
#define FIRSTTRIGGERLAYER 5
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.405 ,0.415 ,0.425 ,0.562 ,0.762  ,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4   ,1.4   ,1.4   ,1.4   ,1.4    ,1.4};
 const static double ZSMLAYER2[NBARR2] ={50e-6 ,50e-6 ,50e-6 ,1e-3   ,1e-3   ,1e-3};
 const static double TSMLAYER2[NBARR2] ={50e-6 ,50e-6 ,50e-6 ,75.5e-6,75.5e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif


#if defined(ITK_EXT_SIMPLIFIED_40TRI_Z)
#define FIRSTTRIGGERLAYER 5
// pixel 0+1+2
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.039,0.075,0.155};  
 const static double ZMAXLAYER0[NBARR0]={1.20,1.2,0.74};
 const static double ZSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
 const static double TSMLAYER0[NBARR0]={50e-6,50e-6,50e-6};
// pixel 3+4
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.213,0.271};  
 const static double ZMAXLAYER1[NBARR1]={0.74,0.74};
 const static double ZSMLAYER1[NBARR1]={50e-6,50e-6};
 const static double TSMLAYER1[NBARR1]={50e-6,50e-6};
// strips
 #define NBARR2 6
 const static double RLAYER2[NBARR2]   ={0.390 ,0.415 ,0.440 ,0.562 ,0.762  ,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4   ,1.4   ,1.4   ,1.4   ,1.4    ,1.4};
 const static double ZSMLAYER2[NBARR2] ={50e-6 ,50e-6 ,50e-6 ,1e-3   ,1e-3   ,1e-3};
 const static double TSMLAYER2[NBARR2] ={50e-6 ,50e-6 ,50e-6 ,75.5e-6,75.5e-6,75.5e-6};

// pixel inner
 #define NRING0 5
 const static double RRING0[2]={0.15,0.32};   // minimum and maximum radius
 const static double ZRING0[NRING0]={0.823,0.90,1.00,1.10,1.20};
 const static double XRING0={THICKNESS};

// pixel outer
 #define NRING1 8
 const static double RRING1[2]={0.08,0.32};   // minimum and maximum radius
 const static double ZRING1[NRING1]={1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00};
 const static double XRING1={THICKNESS};
// strip
 #define NRING2 8
 const static double RRING2[2]={0.389,0.970};   // minimum and maximum radius
 const static double ZRING2[NRING2]={1.50,1.66,1.70,1.74,1.94,2.25,2.60,3.00};
 const static double XRING2={THICKNESS};
// const static double XRING2={0.001};

#endif





// trigger layers at 82cm
#if defined(SETUP_5PIXEL_TR82_50)
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.04,0.08,0.14};  
 const static double ZMAXLAYER0[NBARR0]={1.4,1.4,1.4};
 const static double ZSMLAYER0[NBARR0]={40e-6,40e-6,40e-6};
 const static double TSMLAYER0[NBARR0]={40e-6,40e-6,40e-6};

 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.220,0.300}; 
 const static double ZMAXLAYER1[NBARR1]={1.4,1.4};
 const static double ZSMLAYER1[NBARR1]={40e-6,40e-6};
 const static double TSMLAYER1[NBARR1]={40e-6,40e-6};

 #define NBARR2 6
const static double RLAYER2[NBARR2]={0.46,0.64,0.80,0.82,0.84,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4,1.4,1.4,1.4,1.4,1.4};
 const static double ZSMLAYER2[NBARR2]={40e-6,40e-6,50e-6,50e-6,50e-6,40e-6};
 const static double TSMLAYER2[NBARR2]={40e-6,40e-6,50e-6,50e-6,50e-6,40e-6};
#endif


// trigger layers at 82cm
#if defined(SETUP_5PIXEL_TR82_40)
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.04,0.08,0.14};  
 const static double ZMAXLAYER0[NBARR0]={1.4,1.4,1.4};
 const static double ZSMLAYER0[NBARR0]={40e-6,40e-6,40e-6};
 const static double TSMLAYER0[NBARR0]={40e-6,40e-6,40e-6};
 const static double XMLAYER0[NBARR0]={0.01,0.01,0.01};

 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.220,0.300}; 
 const static double ZMAXLAYER1[NBARR1]={1.4,1.4};
 const static double ZSMLAYER1[NBARR1]={40e-6,40e-6};
 const static double TSMLAYER1[NBARR1]={40e-6,40e-6};
 const static double XMLAYER1[NBARR1]={0.01,0.01};

 #define NBARR2 6
 const static double RLAYER2[NBARR2]={0.46,0.64,0.80,0.82,0.84,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4,1.4,1.4,1.4,1.4,1.4};
 const static double ZSMLAYER2[NBARR2]={40e-6,40e-6,40e-6,40e-6,40e-6,40e-6};
 const static double TSMLAYER2[NBARR2]={40e-6,40e-6,40e-6,40e-6,40e-6,40e-6};
 const static double XMLAYER2[NBARR2]={0.01,0.01,0.01,0.01,0.01,0.01};
#endif

// trigger layers at 82cm   gap 2.5cm
#if defined(SETUP_5PIXEL_TR82B_40)
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.04,0.08,0.14};  
 const static double ZMAXLAYER0[NBARR0]={1.4,1.4,1.4};
 const static double ZSMLAYER0[NBARR0]={40e-6,40e-6,40e-6};
 const static double TSMLAYER0[NBARR0]={40e-6,40e-6,40e-6};

 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.220,0.300}; 
 const static double ZMAXLAYER1[NBARR1]={1.4,1.4};
 const static double ZSMLAYER1[NBARR1]={40e-6,40e-6};
 const static double TSMLAYER1[NBARR1]={40e-6,40e-6};

 #define NBARR2 6
const static double RLAYER2[NBARR2]={0.46,0.64,0.795,0.82,0.845,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4,1.4,1.4,1.4,1.4,1.4};
 const static double ZSMLAYER2[NBARR2]={40e-6,40e-6,40e-6,40e-6,40e-6,40e-6};
 const static double TSMLAYER2[NBARR2]={40e-6,40e-6,40e-6,40e-6,40e-6,40e-6};
#endif

// trigger layers at 82cm   gap 2.5cm
#if defined(SETUP_5PIXEL_TR82A_40)
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.04,0.08,0.14};  
 const static double ZMAXLAYER0[NBARR0]={1.4,1.4,1.4};
 const static double ZSMLAYER0[NBARR0]={40e-6,40e-6,40e-6};
 const static double TSMLAYER0[NBARR0]={40e-6,40e-6,40e-6};

 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.220,0.300}; 
 const static double ZMAXLAYER1[NBARR1]={1.4,1.4};
 const static double ZSMLAYER1[NBARR1]={40e-6,40e-6};
 const static double TSMLAYER1[NBARR1]={40e-6,40e-6};

 #define NBARR2 6
const static double RLAYER2[NBARR2]={0.46,0.64,0.805,0.82,0.835,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.4,1.4,1.4,1.4,1.4,1.4};
 const static double ZSMLAYER2[NBARR2]={40e-6,40e-6,40e-6,40e-6,40e-6,40e-6};
 const static double TSMLAYER2[NBARR2]={40e-6,40e-6,40e-6,40e-6,40e-6,40e-6};
#endif







// geometry setting: CENTRAL6_EQUI
#if  defined(SETUP_CENTRAL_331)
 #define SETUP_BARR1
 #define SETUP_BARR2_3
 #define SETUP_BARR3_1
#elif  defined(SETUP_CENTRAL_332)
 #define SETUP_BARR1
 #define SETUP_BARR2_3
 #define SETUP_BARR3_2
#elif  defined(SETUP_CENTRAL_333)
 #define SETUP_BARR1
 #define SETUP_BARR2_3
 #define SETUP_BARR3_3
#elif  defined(SETUP_CENTRAL_323)
 #define SETUP_BARR1
 #define SETUP_BARR2_2
 #define SETUP_BARR3_3
#elif  defined(SETUP_CENTRAL_322)
 #define SETUP_BARR1
 #define SETUP_BARR2_2
 #define SETUP_BARR3_2
#elif  defined(SETUP_CENTRAL_34)
 #define SETUP_BARR1
#elif  defined(SETUP_CENTRAL_35)
 #define SETUP_BARR1
#elif  defined(SETUP_CENTRAL_36)
 #define SETUP_BARR1
#endif


// Barrel0
#if defined(SETUP_BARR1)
 #define NBARR0 3
 const static double RLAYER0[NBARR0]={0.04,0.07,0.10};  
 const static double ZMAXLAYER0[NBARR0]={1.0,1.0,1.0};
#endif

// Barrel1
#if defined(SETUP_BARR2_3)
 #define NBARR1 3
 const static double RLAYER1[NBARR1]={0.49,0.52,0.55}; 
 const static double ZMAXLAYER1[NBARR1]={1.25,1.25,1.25};
#elif defined(SETUP_BARR2_2)
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.505,0.535}; 
 const static double ZMAXLAYER1[NBARR1]={1.25,1.25};
#elif defined(SETUP_BARR2_1)
 #define NBARR1 1
 const static double RLAYER1[NBARR1]={0.52}; 
 const static double ZMAXLAYER1[NBARR1]={1.25};
#endif

// Barrel2
#if defined(SETUP_BARR3_3)
 #define NBARR2 3
 const static double RLAYER2[NBARR2]={0.94,0.97,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.5,1.5,1.5};

#elif defined(SETUP_BARR3_2)
 #define NBARR2 2
 const static double RLAYER2[NBARR2]={0.97,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.5,1.5};

#elif defined(SETUP_BARR3_1)
 #define NBARR2 1
 const static double RLAYER2[NBARR2]={1.0};
 const static double ZMAXLAYER2[NBARR2]={1.5};

#elif defined(SETUP_CENTRAL_36)
 #define NBARR1 3
 const static double RLAYER1[NBARR1]={0.25,0.40,0.55}; 
// const static double ZMAXLAYER1[NBARR1]={1.3,1.3,1.3};
 const static double ZMAXLAYER1[NBARR1]={3.0,3.0,3.0};
 #define NBARR2 3
 const static double RLAYER2[NBARR2]={0.70,0.85,1.0};
// const static double ZMAXLAYER2[NBARR2]={1.065,1.2825,1.5};
 const static double ZMAXLAYER2[NBARR2]={3.0,3.0,3.0};

#elif defined(SETUP_CENTRAL_35)
 #define NBARR1 3
 const static double RLAYER1[NBARR1]={0.28,0.46,0.64}; 
 const static double ZMAXLAYER1[NBARR1]={1.3,1.3,1.4};
 #define NBARR2 2
 const static double RLAYER2[NBARR2]={0.82,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.239,1.5};

#elif defined(SETUP_CENTRAL_34)
 #define NBARR1 2
 const static double RLAYER1[NBARR1]={0.325,0.55}; 
 const static double ZMAXLAYER1[NBARR1]={1.3,1.4};
 #define NBARR2 2
 const static double RLAYER2[NBARR2]={0.775,1.0};
 const static double ZMAXLAYER2[NBARR2]={1.5,1.5};

#endif


#ifndef NRING0
 // Ring0
 #define NRING0 0
 const static double RRING0[2]={0.0,0.0};   // minimum and maximum radius
 const static double ZRING0[NRING0]={};
 const static double XRING0={THICKNESS};
#endif

#ifndef NRING1
 // Ring1
 #define NRING1 0
 const static double RRING1[2]={0.0,0.0};   // minimum and maximum radius
 const static double ZRING1[NRING1]={};
 const static double XRING1={THICKNESS};
#endif

#ifndef NRING2
 // Ring2
 #define NRING2 0
 const static double RRING2[2]={0.0,0.0};   // minimum and maximum radius
 const static double ZRING2[NRING2]={};
 const static double XRING2={THICKNESS};
#endif


//Endcap1
#define NENDC1 0
const static double ZENDC1[NENDC1]={};
const static double R1ENDC1[NENDC1]={};
const static double R2ENDC1[NENDC1]={};
const static double XENDC1={THICKNESS};
//Endcap2
#define NENDC2 0
const static double ZENDC2[NENDC2]={};
const static double R1ENDC2[NENDC2]={};
const static double R2ENDC2[NENDC2]={};
const static double XENDC2={THICKNESS};

#define NBARRTOT (NBARR0+NBARR1+NBARR2)
#define NRINGTOT (NRING0+NRING1+NRING2)
#define NENDCTOT (NENDC1+NENDC2)





// geometry setting: ALL9
#ifdef  SETUP_ALL9

#define BEAMPIPE
#define RBEAMPIPE 0.025
#define XBEAMPIPE (0.008/0.353)

// Barrel0
#define NBARR0 3
const static double RLAYER0[NBARR0]={0.02856,0.04957,0.07035};
const static double ZMAXLAYER0[NBARR0]={0.4,0.55,0.7};
const static double XBARR0={THICKNESS};   // p2

// Barrel1
#define NBARR1 3
const static double RLAYER1[NBARR1]={0.49469,0.51530,0.53591};
const static double ZMAXLAYER1[NBARR1]={1.01,1.01,1.01};
const static double XBARR1={THICKNESS};

// Barrel2
#define NBARR2 3
const static double RLAYER2[NBARR2]={0.95992,0.98053,1.00115};
const static double ZMAXLAYER2[NBARR2]={1.51,1.51,1.51};
const static double XBARR2={THICKNESS};

//#define NBARRTOT (NBARR0+NBARR1+NBARR2)

// Ring1
#define NRING1 7
const static double RRING1[2]={0.47,0.55};   // minimum and maximum radius
const static double ZRING1[NRING1]={1.040,1.085,1.138,1.187,1.238,1.293,1.350};
const static double XRING1={THICKNESS};

// Ring2
#define NRING2 17
const static double RRING2[2]={0.90,1.01};   // minimum and maximum radius
const static double ZRING2[NRING2]={1.540,1.594,1.649,1.706,1.765,1.827,1.890,1.956,2.024,2.094,2.167,2.242,2.320,2.401,2.484,2.571,2.660};
const static double XRING2={THICKNESS};

#define NRINGTOT (NRING0+NRING1+NRING2)

//Endcap1
#define NENDC1 3
const static double ZENDC1[NENDC1]={1.40,1.48,1.56};
const static double R1ENDC1[NENDC1]={0.14,0.15,0.16};
const static double R2ENDC1[NENDC1]={0.55,0.55,0.55};
const static double XENDC1={THICKNESS};

//Endcap2
#define NENDC2 3
const static double ZENDC2[NENDC2]={2.76,2.88,3.00};
const static double R1ENDC2[NENDC2]={0.275,0.285,0.295};
const static double R2ENDC2[NENDC2]={1.01,1.01,1.01};
const static double XENDC2={THICKNESS};

#define NENDCTOT (NENDC1+NENDC2)

#endif


// hit smearing

//#define ZSMEAR  40E-6  // Pixel resolution in longitudinal direction
//#define TSMEAR  40E-6  // Pixel resolution in transverse direction
//#define ZSMEAR  80E-6  // Pixel resolution in longitudinal direction
//#define TSMEAR  80E-6  // Pixel resolution in transverse direction


#define ZSMEAR  40E-6  // Pixel resolution in longitudinal direction
//#define ZSMEAR  160E-6  // Pixel resolution in longitudinal direction
#define TSMEAR  40E-6  // Pixel resolution in transverse direction
//#define TSMEAR  20E-6  // Pixel resolution in transverse direction

#define XYSMEAR 40E-6  // Pixel resolution for discs in x-y



/* obsolete
#ifdef ZSMEAR
#define ZPIXBARR0 ZSMEAR
#define ZPIXBARR1 ZSMEAR
#define ZPIXBARR2 ZSMEAR
#else
#define ZPIXBARR0 0.0
#define ZPIXBARR1 0.0
#define ZPIXBARR2 0.0
#endif
*/
/* obsolete
#ifdef TSMEAR
#define TPIXBARR0 TSMEAR
#define TPIXBARR1 TSMEAR
#define TPIXBARR2 TSMEAR
#else
#define TPIXBARR0 0.0
#define TPIXBARR1 0.0
#define TPIXBARR2 0.0
#endif
*/

#ifdef XYSMEAR
#define PIXRING0 XYSMEAR
#define PIXRING1 XYSMEAR
#define PIXRING2 XYSMEAR

#define PIXENDC1 XYSMEAR
#define PIXENDC2 XYSMEAR

#else
#define PIXRING0 0.0
#define PIXRING1 0.0
#define PIXRING2 0.0

#define PIXENDC1 0.0
#define PIXENDC2 0.0
#endif

#include "field.h"

// hit efficiency
//#define HITEFF 0.9    // single hit efficiency

//#define NOISE_BARR   1e-6 // noise level


#ifdef MAIN
double rbarr[NBARRTOT],zbarr[NBARRTOT];
double r1ring[NRINGTOT],r2ring[NRINGTOT],zring[NRINGTOT];
double r1endc[NENDCTOT],r2endc[NENDCTOT],zendc[NENDCTOT];

double bms_barr[NBARRTOT];
double rres_barr[NBARRTOT];
double tres_barr[NBARRTOT];
double zres_barr[NBARRTOT];

double bms_ring[NRINGTOT];
double rres_ring[NRINGTOT];
double tres_ring[NRINGTOT];
double zres_ring[NRINGTOT];

double bms_endc[NENDCTOT];
double rres_endc[NENDCTOT];
double tres_endc[NENDCTOT];
double zres_endc[NENDCTOT];
#else
extern double rbarr[NBARRTOT],zbarr[NBARRTOT];
extern double r1ring[NRINGTOT],r2ring[NRINGTOT],zring[NRINGTOT];
extern double r1endc[NENDCTOT],r2endc[NENDCTOT],zendc[NENDCTOT];

extern double bms_barr[NBARRTOT];
extern double rres_barr[NBARRTOT];
extern double tres_barr[NBARRTOT];
extern double zres_barr[NBARRTOT];

extern double bms_ring[NRINGTOT];
extern double rres_ring[NRINGTOT];
extern double tres_ring[NRINGTOT];
extern double zres_ring[NRINGTOT];

extern double bms_endc[NENDCTOT];
extern double rres_endc[NENDCTOT];
extern double tres_endc[NENDCTOT];
extern double zres_endc[NENDCTOT];
#endif

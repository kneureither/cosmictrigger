#ifndef COSMICTRIGGER_KARIMAKIHELIXFIT_H
#define COSMICTRIGGER_KARIMAKIHELIXFIT_H

#define MAXLAYER 16
#include "basicDefines.h"


typedef struct KariFit
{
  float r3d;
  float rad;
  float dca;
  float phi;
  float tchi2n;
  float z0;
  float theta;
  float zchi2n;
} KARITRACK;

int karimakiHelixfit(KariFit &karires, int npoints, double xp[MAXLAYER], double yp[MAXLAYER], double zp[MAXLAYER],
        double phip[MAXLAYER], double thetas[MAXLAYER], double tres[MAXLAYER], double zres[MAXLAYER], double rres[MAXLAYER]);

template <typename T>
static int sign(const T val) {
    return (T(0) < val) - (val < T(0));
}

static void swapKariMomentum(KariFit &kari, float rec_r, int &count) {
    if(sign(kari.rad) != sign(rec_r)) {
        kari.rad = -kari.rad;
        kari.dca = -kari.dca;
        kari.r3d = -kari.r3d;
        count++;
    }
}

static void swapKariBField(KariFit &kari) {
    kari.rad = -kari.rad;
    kari.r3d = -kari.r3d;
    kari.dca = -kari.dca;
}

static void correctKariDirection(KariFit &kari) {

    if(kari.phi > 0) {
        //printf("\tKARI CORRECT OLD: \t rad %f \tr3d %f phi %f \t theta %f\n", kari.rad, kari.r3d, kari.phi, kari.theta);
        kari.rad = -kari.rad;
        kari.r3d = -kari.r3d;
        kari.phi = kari.phi - PI; //this is save as phi>0
//        kari.theta = PI - kari.theta;
        //printf("\tKARI CORRECT NEW: \t rad %f \tr3d %f phi %f \t theta %f\n", kari.rad, kari.r3d, kari.phi, kari.theta);
    }
}

#endif //COSMICTRIGGER_KARIMAKIHELIXFIT_H

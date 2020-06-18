//
// Created by Konstantin Neureither on 29.05.20.
//
#ifndef COSMICTRIGGER_RECONSTRUCTION_ACCURACY_H
#define COSMICTRIGGER_RECONSTRUCTION_ACCURACY_H

#ifndef PI
#define PI 3.1415926535
#endif //PI

#include <math.h>
#include "../3rdparty/karimaki/karimaki.h"


void correctKariDirection(KariFit &kari) {

    if(kari.phi > 0) {
        printf("\tKARI CORRECT OLD: \t rad %f \tr3d %f phi %f \t theta %f\n", kari.rad, kari.r3d, kari.phi, kari.theta);
        kari.rad = -kari.rad;
        kari.r3d = -kari.r3d;
        kari.phi = kari.phi - PI; //this is save as phi>0
//        kari.theta = PI - kari.theta;
        printf("\tKARI CORRECT NEW: \t rad %f \tr3d %f phi %f \t theta %f\n", kari.rad, kari.r3d, kari.phi, kari.theta);
    }
 }

 void swapKariMomentum(KariFit &kari, float rec_r, int &count) {
    if(sgn(kari.rad) != sgn(rec_r)) {
        kari.rad = -kari.rad;
        kari.dca = -kari.dca;
        kari.r3d = -kari.r3d;
        count++;
    }
}

//Declaration for function defined in "../3rdparty/karimaki/karimaki_hit.c"
int karimaki_hit(KariFit&, int , double *, double *, double *, double *, double*, double *, double *, double *);

#endif //COSMICTRIGGER_RECONSTRUCTION_ACCURACY_H
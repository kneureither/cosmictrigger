//
// Created by Konstantin Neureither on 29.05.20.
//
#ifndef COSMICTRIGGER_RECONSTRUCTION_ACCURACY_H
#define COSMICTRIGGER_RECONSTRUCTION_ACCURACY_H
#ifndef TRIPLET_HIT_ARRAY_LENGTH
#define TRIPLET_HIT_ARRAY_LENGTH 1024
#endif //TRIPLET_HIT_ARRAY_LENGTH
#endif //COSMICTRIGGER_RECONSTRUCTION_ACCURACY_H

#ifndef PI
#define PI 3.1415926535
#endif //PI

#include <math.h>
#include "../3rdparty/karimaki/karimaki.h"


int combineHits(std::vector<double> &xp, std::vector<double> &yp, std::vector<double> &zp,
                std::vector<int> &sids, std::vector<double> &phi_tracks, std::vector<double> &thetas, int nhits, int ntriplets,
                float (&x00)[TRIPLET_HIT_ARRAY_LENGTH], float (&x10)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&x01)[TRIPLET_HIT_ARRAY_LENGTH], float (&x20)[TRIPLET_HIT_ARRAY_LENGTH], float (&x21)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&y00)[TRIPLET_HIT_ARRAY_LENGTH], float (&y10)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&y01)[TRIPLET_HIT_ARRAY_LENGTH], float (&y20)[TRIPLET_HIT_ARRAY_LENGTH], float (&y21)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&z00)[TRIPLET_HIT_ARRAY_LENGTH], float (&z10)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&z01)[TRIPLET_HIT_ARRAY_LENGTH], float (&z20)[TRIPLET_HIT_ARRAY_LENGTH], float (&z21)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&sid00)[TRIPLET_HIT_ARRAY_LENGTH], float (&sid10)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&sid01)[TRIPLET_HIT_ARRAY_LENGTH], float (&sid20)[TRIPLET_HIT_ARRAY_LENGTH], float (&sid21)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&tan01)[TRIPLET_HIT_ARRAY_LENGTH], float (&tan12)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&lam01)[TRIPLET_HIT_ARRAY_LENGTH], float (&lam12)[TRIPLET_HIT_ARRAY_LENGTH]) {


    //These two loops add the basic ntriplets+2 hits
    for(int i=0; i<ntriplets; i++) {
    xp.push_back(x00[i]);
    yp.push_back(y00[i]);
    zp.push_back(z00[i]);
    sids.push_back((int) sid00[i]);
    phi_tracks.push_back(tan01[(i == 0 ? 0 : i-1)]);
    thetas.push_back(PI / 2 - lam01[(i == 0 ? 0 : i-1)]);
    }
    for(int i=ntriplets-2; i<ntriplets; i++) {
    xp.push_back(x20[i]);
    yp.push_back(y20[i]);
    zp.push_back(z20[i]);
    sids.push_back((int) sid20[i]);
    phi_tracks.push_back(tan12[i]);
    thetas.push_back(PI / 2 - lam12[i]);
    }
    int ncomhits = ntriplets+2;

    //when this is true, all hits in the segments are searched for additional duplicate hits
    bool ADD_REDUNDANT_HITS = false;
    //TODO extract tan and theta for these hits also

    if((ADD_REDUNDANT_HITS ? nhits > ntriplets+2 : false)) {
        printf("ntriplets = %d, nhits = %d\n", ntriplets, nhits);
        std::vector<double>::iterator it;
        for(int i=0; i<ntriplets; i++) {
            it = std::find(xp.begin(), xp.end(), x01[i]);
            if(it == xp.end()) {
                xp.push_back(x01[i]);
                yp.push_back(y01[i]);
                zp.push_back(z01[i]);
                sids.push_back((int) sid01[i]);
                ncomhits++;
                printf("Hit found in x01\n");
            }

            it = std::find(xp.begin(), xp.end(), x10[i]);
            if(it == xp.end()) {
                xp.push_back(x10[i]);
                yp.push_back(y10[i]);
                zp.push_back(z10[i]);
                sids.push_back((int) sid10[i]);
                ncomhits++;
                printf("Hit found in x10\n");
            }

            it = std::find(xp.begin(), xp.end(), x21[i]);
            if(it == xp.end()) {
                xp.push_back(x21[i]);
                yp.push_back(y21[i]);
                zp.push_back(z21[i]);
                sids.push_back((int) sid21[i]);
                ncomhits++;
                printf("Hit found in x21\n");
            }
        }
    }
    return ncomhits;
}

template <typename T>
unsigned int get_layer(T sid) {
    if(0 <= sid && sid < 1024) {
        return 0;
    } else if (1024 <= sid && sid < 2048) {
        return 1;
    } else if ((2000 <= sid && sid < 3000) || (10000 <= sid && sid < 11500) || (14000 <= sid && sid < 15200)) {
        return 2;
    } else if ((3000 <= sid && sid < 4000) || (11500 <= sid && sid < 12500) || (15200 <= sid && sid < 16500)) {
        return 3;
    } else {
        return 100;
    }
}

int counthitlayers(int nhits, float (&sid00)[TRIPLET_HIT_ARRAY_LENGTH], float (&sid11)[TRIPLET_HIT_ARRAY_LENGTH]) {
    int layerhits[4] = {0,0,0,0};

    for(int i=0; i<nhits-2; i++) {
        int layer=get_layer(sid00[i]);
        printf("layer %d sid %f", layer, sid00[i]);
        layerhits[layer] += (layerhits[layer] < 2 ? 1 : 0);
    }
    for(int i=nhits-4; i<nhits-2; i++) {
        int layer=get_layer(sid11[i]);
        printf("layer %d sid %f", layer, sid11[i]);
        layerhits[layer] += (layerhits[layer] < 2 ? 1 : 0);
    }
    return layerhits[0] + layerhits[1] + layerhits[2] + layerhits[3];
}


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


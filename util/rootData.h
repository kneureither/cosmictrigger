//
// Created by Konstantin Neureither on 18.06.20.
//

#ifndef COSMICTRIGGER_ROOTDATA_H
#define COSMICTRIGGER_ROOTDATA_H

#include <math.h>
#include "../karimaki/karimakiHelixfit.h"
#include "basicDefines.h"
#include "SegsTreeRead.h"

template <typename T>
static unsigned int get_layer(T sid) {
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

static unsigned int get_layer(float x, float y) {
    float r = sqrt(pow(x,2) + pow(y,2));
    const float layerBoundaries[5] = {0.0, 26.00, 51.0, 78.5, 100.00};
    if (layerBoundaries[0] <= r && r <= layerBoundaries[1]) {
        return 0;
    } else if (layerBoundaries[1] < r && r <= layerBoundaries[2]) {
        return 1;
    } else if (layerBoundaries[2] < r && r <= layerBoundaries[3]) {
        return 2;
    } else if (layerBoundaries[3] < r && r <= layerBoundaries[4]) {
        return 3;
    } else {
        return -1;
    }
}

static int combineHits(std::vector<double> &xp, std::vector<double> &yp, std::vector<double> &zp,
                std::vector<int> &sids, std::vector<double> &phi_tracks, std::vector<double> &thetas, int nhits, int ntriplets,
                float (&x00)[TRIPLET_HIT_ARRAY_LENGTH], float (&x10)[TRIPLET_HIT_ARRAY_LENGTH], float (&x11)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&x01)[TRIPLET_HIT_ARRAY_LENGTH], float (&x20)[TRIPLET_HIT_ARRAY_LENGTH], float (&x21)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&y00)[TRIPLET_HIT_ARRAY_LENGTH], float (&y10)[TRIPLET_HIT_ARRAY_LENGTH], float (&y11)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&y01)[TRIPLET_HIT_ARRAY_LENGTH], float (&y20)[TRIPLET_HIT_ARRAY_LENGTH], float (&y21)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&z00)[TRIPLET_HIT_ARRAY_LENGTH], float (&z10)[TRIPLET_HIT_ARRAY_LENGTH], float (&z11)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&z01)[TRIPLET_HIT_ARRAY_LENGTH], float (&z20)[TRIPLET_HIT_ARRAY_LENGTH], float (&z21)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&sid00)[TRIPLET_HIT_ARRAY_LENGTH], float (&sid10)[TRIPLET_HIT_ARRAY_LENGTH], float (&sid11)[TRIPLET_HIT_ARRAY_LENGTH],
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
                phi_tracks.push_back(tan01[(i == 0 ? 0 : i-1)]);
                thetas.push_back(PI / 2 - lam01[(i == 0 ? 0 : i-1)]);
                ncomhits++;
                printf("Hit found in x01\n");
            }

            it = std::find(xp.begin(), xp.end(), x11[i]);
            if(it == xp.end()) {
                xp.push_back(x11[i]);
                yp.push_back(y11[i]);
                zp.push_back(z11[i]);
                sids.push_back((int) sid11[i]);
                phi_tracks.push_back(tan01[i]);
                thetas.push_back(PI / 2 - lam01[i]);
                ncomhits++;
                printf("Hit found in x11\n");
            }

            it = std::find(xp.begin(), xp.end(), x21[i]);
            if(it == xp.end()) {
                xp.push_back(x21[i]);
                yp.push_back(y21[i]);
                zp.push_back(z21[i]);
                sids.push_back((int) sid21[i]);
                phi_tracks.push_back(tan12[i]);
                thetas.push_back(PI / 2 - lam12[i]);
                ncomhits++;
                printf("Hit found in x21\n");
            }
        }
    }
    return ncomhits;
}

static int combineHits(std::vector<double> &xp, std::vector<double> &yp, std::vector<double> &zp,
                       std::vector<int> &sids, std::vector<double> &phi_tracks, std::vector<double> &thetas,
                       SegsTreeRead &Segs, bool ADD_REDUNDANT_HITS) {

    //These two loops add the basic ntriplets+2 hits
    for(int i=0; i<Segs.rec_ntriplet; i++) {
        xp.push_back(Segs.x00[i]);
        yp.push_back(Segs.y00[i]);
        zp.push_back(Segs.z00[i]);
        sids.push_back((int) Segs.sid00[i]);
        phi_tracks.push_back(Segs.rec_tan01[(i == 0 ? 0 : i-1)]);
        thetas.push_back(PI / 2 - Segs.rec_lam01[(i == 0 ? 0 : i-1)]);
    }
    for(int i=Segs.rec_ntriplet-2; i<Segs.rec_ntriplet; i++) {
        xp.push_back(Segs.x20[i]);
        yp.push_back(Segs.y20[i]);
        zp.push_back(Segs.z20[i]);
        sids.push_back((int) Segs.sid20[i]);
        phi_tracks.push_back(Segs.rec_tan12[i]);
        thetas.push_back(PI / 2 - Segs.rec_lam12[i]);
    }
    int ncomhits = Segs.rec_ntriplet+2;

    //when this is true, all hits in the segments are searched for additional duplicate hits
//    bool ADD_REDUNDANT_HITS = false;
    //TODO extract tan and theta for these hits also

    if((ADD_REDUNDANT_HITS ? Segs.rec_nhit > Segs.rec_ntriplet+2 : false)) {
        printf("ntriplets = %d, nhits = %d\n", Segs.rec_ntriplet, Segs.rec_nhit);
        std::vector<double>::iterator it;
        for(int i=0; i<Segs.rec_ntriplet; i++) {
            it = std::find(xp.begin(), xp.end(), Segs.x01[i]);
            if(it == xp.end()) {
                xp.push_back(Segs.x01[i]);
                yp.push_back(Segs.y01[i]);
                zp.push_back(Segs.z01[i]);
                sids.push_back((int) Segs.sid01[i]);
                phi_tracks.push_back(Segs.rec_tan01[(i == 0 ? 0 : i-1)]);
                thetas.push_back(PI / 2 - Segs.rec_lam01[(i == 0 ? 0 : i-1)]);
                ncomhits++;
                printf("Hit found in x01\n");
            }

            it = std::find(xp.begin(), xp.end(), Segs.x11[i]);
            if(it == xp.end()) {
                xp.push_back(Segs.x11[i]);
                yp.push_back(Segs.y11[i]);
                zp.push_back(Segs.z11[i]);
                sids.push_back((int) Segs.sid11[i]);
                phi_tracks.push_back(Segs.rec_tan01[i]);
                thetas.push_back(PI / 2 - Segs.rec_lam01[i]);
                ncomhits++;
                printf("Hit found in x11\n");
            }

            it = std::find(xp.begin(), xp.end(), Segs.x21[i]);
            if(it == xp.end()) {
                xp.push_back(Segs.x21[i]);
                yp.push_back(Segs.y21[i]);
                zp.push_back(Segs.z21[i]);
                sids.push_back((int) Segs.sid21[i]);
                phi_tracks.push_back(Segs.rec_tan12[i]);
                thetas.push_back(PI / 2 - Segs.rec_lam12[i]);
                ncomhits++;
                printf("Hit found in x21\n");
            }
        }
    }
    return ncomhits;
}

static int combineBasicHits(std::vector<double> &xp, std::vector<double> &yp, std::vector<double> &zp,
                std::vector<double> &phi_tracks, std::vector<double> &thetas, std::vector<int> &layers, int nhits, int ntriplets,
                float (&x00)[TRIPLET_HIT_ARRAY_LENGTH], float (&x20)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&y00)[TRIPLET_HIT_ARRAY_LENGTH], float (&y20)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&z00)[TRIPLET_HIT_ARRAY_LENGTH], float (&z20)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&sid00)[TRIPLET_HIT_ARRAY_LENGTH], float (&sid20)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&tan01)[TRIPLET_HIT_ARRAY_LENGTH], float (&tan12)[TRIPLET_HIT_ARRAY_LENGTH],
                float (&lam01)[TRIPLET_HIT_ARRAY_LENGTH], float (&lam12)[TRIPLET_HIT_ARRAY_LENGTH]) {


    //These two loops add the basic ntriplets+2 hits
    for(int i=0; i<ntriplets; i++) {
        xp.push_back(x00[i]);
        yp.push_back(y00[i]);
        zp.push_back(z00[i]);
        layers.push_back((int) get_layer(xp[xp.size()-1], yp[yp.size()-1]));
        phi_tracks.push_back(tan01[(i == 0 ? 0 : i-1)]);
        thetas.push_back(PI / 2 - lam01[(i == 0 ? 0 : i-1)]);
    }
    for(int i=ntriplets-2; i<ntriplets; i++) {
        xp.push_back(x20[i]);
        yp.push_back(y20[i]);
        zp.push_back(z20[i]);
        layers.push_back((int) get_layer(xp[xp.size()-1], yp[yp.size()-1]));
        phi_tracks.push_back(tan12[i]);
        thetas.push_back(PI / 2 - lam12[i]);
    }
    int ncomhits = ntriplets+2;
    return ncomhits;
}

static int counthitlayers(int nhits, float (&sid00)[TRIPLET_HIT_ARRAY_LENGTH], float (&sid11)[TRIPLET_HIT_ARRAY_LENGTH]) {
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

static int get_charge_from_type(int type) {
    if (type == 4) {
        //mu-
        return 1;
    } else if (type == 3) {
        //mu+
        return -1;
    } else {
        //not a muon
        return 1;
    }
}

#endif //COSMICTRIGGER_ROOTDATA_H

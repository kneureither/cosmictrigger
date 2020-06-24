//
// Created by Konstantin Neureither on 07.06.20.
//

#ifndef COSMICTRIGGER_TRIGONOMETRY_H
#define COSMICTRIGGER_TRIGONOMETRY_H
#ifndef PI
#define PI 3.1415926535
#endif //PI

#include <assert.h>
#include <type_traits>




typedef struct HitWeights {
    float azimutal;
    float transversal;
}HITWE;

HitWeights getHitWeigths(float RMS, float xp, float yp, float phi_track, float theta) {
    float phi_hit;
    float dphi;
    HitWeights Angles;
    phi_hit = atan2(yp, xp);
    dphi= phi_track - phi_hit;

    Angles.azimutal = std::sin(dphi) * RMS;
    Angles.transversal = std::cos(theta) * RMS;

    return Angles;
}

void getResolutionArrs(int nhits, double RMS, std::vector<double> &tres, std::vector<double> &zres, std::vector<double> &rres,
                       std::vector<double> &xp, std::vector<double> &yp, std::vector<double> &phi_hits, double theta) {

    for(int i=0; i < nhits; i++) {
        double phi_track = atan2(yp[i], xp[i]);
        double dphi = phi_track - phi_hits[i];

        tres.push_back(abs(sin(dphi) * RMS));
        zres.push_back(abs(cos(theta) * RMS));
        rres.push_back(0.0);
    }
}

#endif //COSMICTRIGGER_TRIGONOMETRY_H
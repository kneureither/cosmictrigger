//
// Created by Konstantin Neureither on 09.06.20.
//
#include <assert.h>
#include <cstdlib>
#include <ctime>
#include <stdio.h>
#include <vector>
#include "../scripts/reconstruction_accuracy.h"

#ifndef TRIPLET_HIT_ARRAY_LENGTH
#define TRIPLET_HIT_ARRAY_LENGTH 1024
#endif //TRIPLET_HIT_ARRAY_LENGTH

void TESTncombinedhits();

int main (int argc, char *argv[]) {
    srand (static_cast <unsigned> (time(0)));

    TESTncombinedhits();


    printf("SUCCESS! Finished all tests.\n");
    return 0;
}

void TESTncombinedhits() {
    printf("STATUS : Running TESTncombinedhits()...\n");
    //float r2 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/X));
    std::vector<double> xp;
    std::vector<double> yp;
    std::vector<double> zp;
    std::vector<int> sids;
    std::vector<double> phi_tracks;
    std::vector<double> thetas;

    int ntriplets = 4;
    int nhits = 6;

    float x00[TRIPLET_HIT_ARRAY_LENGTH] = {1.0, 2.0, 3.0, 4.0};
    float x20[TRIPLET_HIT_ARRAY_LENGTH] = {3.0, 4.0, 5.0, 6.0};
    float y00[TRIPLET_HIT_ARRAY_LENGTH] = {11.0, 12.0, 13.0, 14.0};
    float y20[TRIPLET_HIT_ARRAY_LENGTH] = {13.0, 14.0, 15.0, 16.0};
    float z00[TRIPLET_HIT_ARRAY_LENGTH] = {21.0, 22.0, 23.0, 24.0};
    float z20[TRIPLET_HIT_ARRAY_LENGTH] = {23.0, 24.0, 25.0, 26.0};
    float sid00[TRIPLET_HIT_ARRAY_LENGTH];
    float sid20[TRIPLET_HIT_ARRAY_LENGTH];
    float tan01[TRIPLET_HIT_ARRAY_LENGTH] = {11.0, 12.0, 13.0};
    float tan12[TRIPLET_HIT_ARRAY_LENGTH] = {12.0, 13.0, 14.0};
    float lam01[TRIPLET_HIT_ARRAY_LENGTH];
    float lam12[TRIPLET_HIT_ARRAY_LENGTH];


    float x10[TRIPLET_HIT_ARRAY_LENGTH];
    float x01[TRIPLET_HIT_ARRAY_LENGTH];
    float x21[TRIPLET_HIT_ARRAY_LENGTH];
    float y10[TRIPLET_HIT_ARRAY_LENGTH];
    float y01[TRIPLET_HIT_ARRAY_LENGTH];
    float y21[TRIPLET_HIT_ARRAY_LENGTH];
    float z10[TRIPLET_HIT_ARRAY_LENGTH];
    float z01[TRIPLET_HIT_ARRAY_LENGTH];
    float z21[TRIPLET_HIT_ARRAY_LENGTH];
    float sid10[TRIPLET_HIT_ARRAY_LENGTH];
    float sid01[TRIPLET_HIT_ARRAY_LENGTH];
    float sid21[TRIPLET_HIT_ARRAY_LENGTH];

    int ncombinedhits = combineHits(xp, yp, zp, sids, phi_tracks, thetas, nhits, ntriplets,
                                    x00, x10, x01, x20, x21, y00, y10, y01, y20, y21, z00, z10, z01, z20, z21,
                                    sid00, sid10, sid01, sid20, sid21, tan01, tan12, lam01, lam12);


    assert(ncombinedhits == nhits);
    assert(phi_tracks[0] == tan01[0]);
    assert(phi_tracks[1] == tan01[0]);
    assert(phi_tracks[2] == tan01[1]);
    printf("\tasserting: phi_tracks[nhits - 1]=%f == tan12[ntriplets-2]=%f\n", phi_tracks[nhits - 1], tan12[ntriplets-2]);
    assert(phi_tracks[nhits - 1] == tan12[ntriplets-2]);
    printf("\tasserting: phi_tracks[nhits - 2]=%f == "
           "tan12[ntriplets-3]=%f\n", phi_tracks[nhits - 2], tan12[ntriplets-3]);
    assert(phi_tracks[nhits - 2] == tan12[ntriplets-3]);

}


void TESTkarimaki() {
    printf("STATUS : Running TESTkarimaki()...\n");

}


//
// Created by Konstantin Neureither on 29.05.20.
//

#ifndef COSMICTRIGGER_RECONSTRUCTION_ACCURACY_H
#define COSMICTRIGGER_RECONSTRUCTION_ACCURACY_H

#endif //COSMICTRIGGER_RECONSTRUCTION_ACCURACY_H


void combinehits(std::vector<float> &xp, std::vector<float> &yp, std::vector<float> &zp, int nhits,
        float (&x00)[TRIPLET_HIT_ARRAY_LENGTH], float (&x11)[TRIPLET_HIT_ARRAY_LENGTH],
        float (&y00)[TRIPLET_HIT_ARRAY_LENGTH], float (&y11)[TRIPLET_HIT_ARRAY_LENGTH],
        float (&z00)[TRIPLET_HIT_ARRAY_LENGTH], float (&z11)[TRIPLET_HIT_ARRAY_LENGTH]) {

    for(int i=0; i<nhits-2; i++) {
        xp.push_back(x00[i]);
        yp.push_back(y00[i]);
        zp.push_back(z00[i]);
//        xp.push_back(1);
//        yp.push_back(2);
//        zp.push_back(3);
    }
    for(int i=nhits-4; i<nhits-2; i++) {
        xp.push_back(x11[i]);
        yp.push_back(y11[i]);
        zp.push_back(z11[i]);
//        xp.push_back(4);
//        yp.push_back(5);
//        zp.push_back(6);
    }
}


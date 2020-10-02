//
// Created by Konstantin Neureither on 09.06.20.
//
#include <assert.h>
#include <cstdlib>
#include <ctime>
#include <stdio.h>
#include <vector>
#include "rootData.h"
#include "basicDefines.h"
#include "SegsTreeRead.h"
#include "SlimSegsTree.h"
#include <iostream>
#include "testTemplateData.h"
#include "testTemplateID.h"
#include "../Mu3eCosPat/include/TemplateBank.h"
#include "testTemplateBank.h"


#ifndef TRIPLET_HIT_ARRAY_LENGTH
#define TRIPLET_HIT_ARRAY_LENGTH 1024
#endif //TRIPLET_HIT_ARRAY_LENGTH

void TESTncombinedhits();
void TESTslimsegs();



int main (int argc, char *argv[]) {
    srand (static_cast <unsigned> (time(0)));

//    TESTncombinedhits();
//    TESTslimsegs();
//    TESTtemplateData();
//    testTemplateID();

    testTemplateBank();


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
    float tan01[TRIPLET_HIT_ARRAY_LENGTH] = {11.0, 12.0, 13.0, 14.0};
    float tan12[TRIPLET_HIT_ARRAY_LENGTH] = {12.0, 13.0, 14.0, 15.0};
    float lam01[TRIPLET_HIT_ARRAY_LENGTH] = {0.0,0.1,0.2,0.3};
    float lam12[TRIPLET_HIT_ARRAY_LENGTH] = {0.1,0.2,0.3,0.4};


    float x10[TRIPLET_HIT_ARRAY_LENGTH];
    float x01[TRIPLET_HIT_ARRAY_LENGTH];
    float x11[TRIPLET_HIT_ARRAY_LENGTH];
    float x21[TRIPLET_HIT_ARRAY_LENGTH];
    float y10[TRIPLET_HIT_ARRAY_LENGTH];
    float y01[TRIPLET_HIT_ARRAY_LENGTH];
    float y11[TRIPLET_HIT_ARRAY_LENGTH];
    float y21[TRIPLET_HIT_ARRAY_LENGTH];
    float z10[TRIPLET_HIT_ARRAY_LENGTH];
    float z01[TRIPLET_HIT_ARRAY_LENGTH];
    float z11[TRIPLET_HIT_ARRAY_LENGTH];
    float z21[TRIPLET_HIT_ARRAY_LENGTH];
    float sid10[TRIPLET_HIT_ARRAY_LENGTH];
    float sid01[TRIPLET_HIT_ARRAY_LENGTH];
    float sid21[TRIPLET_HIT_ARRAY_LENGTH];
    float sid11[TRIPLET_HIT_ARRAY_LENGTH];

    int ncombinedhits = combineHits(xp, yp, zp, sids, phi_tracks, thetas, nhits, ntriplets,
                                    x00, x10, x11, x01, x20, x21, y00, y10,y11, y01, y20, y21, z00, z10,z11, z01, z20, z21,
                                    sid00, sid10, sid11, sid01, sid20, sid21, tan01, tan12, lam01, lam12);


    assert(ncombinedhits == nhits);
    assert(phi_tracks[0] == tan01[0]);
    assert(phi_tracks[1] == tan01[0]);
    assert(phi_tracks[2] == tan01[1]);
    printf("\tasserting: phi_tracks[nhits - 1]=%f == tan12[ntriplets-1]=%f\n", phi_tracks[nhits - 1], tan12[ntriplets-1]);
    assert(phi_tracks[nhits - 1] == tan12[ntriplets-1]);
    printf("\tasserting: phi_tracks[nhits - 2]=%f == "
           "tan12[ntriplets-2]=%f\n", phi_tracks[nhits - 2], tan12[ntriplets-2]);
    assert(phi_tracks[nhits - 2] == tan12[ntriplets-2]);

}

void TESTkarimaki() {
    printf("STATUS : Running TESTkarimaki()...\n");

}

void TESTslimsegs() {
    printf("STATUS : Running TESTslimsegs()...\n");
    std::string rootfile = "data/SlimmedData/mu3e_slimmed_segs_test.root";

    // FILE FOR WRITING
    TFile toutF(rootfile.c_str(), "recreate");
    if (!toutF.IsOpen()) {
        std::cout << "[ERROR] File " << toutF.GetName() << " is not open!" << std::endl;
        exit(0);
    }
    TTree t_slim("SlimSegs","Tree with slimmed simulation data");

    //class representation for slimmed down 'slimSegs' Tree and Write functionality
    SlimSegsTreeWrite SlimSegsW = SlimSegsTreeWrite(&t_slim);
    std::vector<SlimSegsTreeWrite> SlimSegsVec;

    for(int i=0; i<10; i++) {
        SlimSegsW.reInitializeData();
        SlimSegsW.eventID = i;
        SlimSegsW.runID = rand();
        SlimSegsW.segsIndex = (unsigned int) rand();

        SlimSegsW.rec_nhit = (rand()%8);
        SlimSegsW.rec_p = static_cast <float> (rand());
        SlimSegsW.kari_p = static_cast <float> (rand());
        SlimSegsW.ncombinedhits = 2;
        for(int i=0; i<4; i++) {
            SlimSegsW.xp.push_back(static_cast <float> (rand()));
            SlimSegsW.yp.push_back(static_cast <float> (rand()));
            SlimSegsW.zp.push_back(static_cast <float> (rand()));
            SlimSegsW.layerp.push_back(rand());
        }
        SlimSegsW.fillTree();
        SlimSegsVec.push_back(SlimSegsW);
    }

    toutF.Write();
    toutF.Close();


    // FILE FOR READING
    TFile tinF(rootfile.c_str());
    if (!tinF.IsOpen()) {
        std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
        exit(0);
    }

    TTree *t_slimsegs;
    tinF.GetObject("SlimSegs", t_slimsegs);
    //class representation for segs tree and read functionality
    SlimSegsTreeRead SlimSegsR = SlimSegsTreeRead(t_slimsegs);

    for(int i = 0; i<10; i++) {

        SlimSegsR.getEntry(i);
        printf("asserting SlimSegsR.eventID == SlimSegsW.eventID : %d == %d\n", SlimSegsR.eventID, SlimSegsVec[i].eventID);
        assert(SlimSegsR.eventID == SlimSegsVec[i].eventID);
        printf("asserting SlimSegsR.runID == SlimSegsR[%d].runID : %d == %d\n", i, SlimSegsR.runID, SlimSegsVec[i].runID);
        assert(SlimSegsR.runID == SlimSegsVec[i].runID);
        printf("asserting SlimSegsR.segsIndex == SlimSegsR[%d].segsIndex : %d == %d\n", i, SlimSegsR.segsIndex, SlimSegsVec[i].segsIndex);
        assert(SlimSegsR.segsIndex == SlimSegsVec[i].segsIndex);
        printf("asserting SlimSegsR.rec_p == SlimSegsR[%d].rec_p : %f == %f\n", i, SlimSegsR.rec_p, SlimSegsVec[i].rec_p);
        assert(SlimSegsR.rec_p == SlimSegsVec[i].rec_p);
        printf("asserting SlimSegsR.kari_p == SlimSegsR[%d].kari_p : %f == %f\n\n", i, SlimSegsR.kari_p, SlimSegsVec[i].kari_p);
        assert(SlimSegsR.kari_p == SlimSegsVec[i].kari_p);
    }
}


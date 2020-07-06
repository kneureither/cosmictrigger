//
// Created by Konstantin Neureither on 30.06.20.
//

#ifndef COSMICTRIGGER_SEGSTREEREAD_H
#define COSMICTRIGGER_SEGSTREEREAD_H

#include "TFile.h"
#include "TTree.h"
#include "basicDefines.h"


class SegsTreeRead {
private:
public:
    TTree *tr_segs;

    void setBranches();
    explicit SegsTreeRead(TTree *t_segs);
    void getEntry(const int &index);

    //meta data
    unsigned int segs_entries;
    int rec_event;
    int rec_nhit;
    int rec_ntriplet;

    //monte carlo data
    int mc_tid;
    float mc_p;
    float mc_pt;
    float mc_theta;
    float mc_phi;
    float mc_lam;
    int mc_type;
    int mc_pid;

    //reconstruction data
    float rec_p;
    float rec_r;
    float rec_rt;
    float rec_zpca_z;
    float rec_zpca_r;

    //actual hit data of reconstructed segs
    float x00[TRIPLET_HIT_ARRAY_LENGTH];
    float x10[TRIPLET_HIT_ARRAY_LENGTH];
    float x20[TRIPLET_HIT_ARRAY_LENGTH];
    float x01[TRIPLET_HIT_ARRAY_LENGTH];
    float x11[TRIPLET_HIT_ARRAY_LENGTH];
    float x21[TRIPLET_HIT_ARRAY_LENGTH];

    float y00[TRIPLET_HIT_ARRAY_LENGTH];
    float y10[TRIPLET_HIT_ARRAY_LENGTH];
    float y20[TRIPLET_HIT_ARRAY_LENGTH];
    float y01[TRIPLET_HIT_ARRAY_LENGTH];
    float y11[TRIPLET_HIT_ARRAY_LENGTH];
    float y21[TRIPLET_HIT_ARRAY_LENGTH];

    float z00[TRIPLET_HIT_ARRAY_LENGTH];
    float z10[TRIPLET_HIT_ARRAY_LENGTH];
    float z20[TRIPLET_HIT_ARRAY_LENGTH];
    float z01[TRIPLET_HIT_ARRAY_LENGTH];
    float z11[TRIPLET_HIT_ARRAY_LENGTH];
    float z21[TRIPLET_HIT_ARRAY_LENGTH];

    float sid00[TRIPLET_HIT_ARRAY_LENGTH];
    float sid10[TRIPLET_HIT_ARRAY_LENGTH];
    float sid20[TRIPLET_HIT_ARRAY_LENGTH];
    float sid01[TRIPLET_HIT_ARRAY_LENGTH];
    float sid11[TRIPLET_HIT_ARRAY_LENGTH];
    float sid21[TRIPLET_HIT_ARRAY_LENGTH];

    float rec_tan01[TRIPLET_HIT_ARRAY_LENGTH];
    float rec_tan12[TRIPLET_HIT_ARRAY_LENGTH];
    float rec_lam01[TRIPLET_HIT_ARRAY_LENGTH];
    float rec_lam12[TRIPLET_HIT_ARRAY_LENGTH];

};

class SegsTreeReadPlus : public SegsTreeRead {
    //this class also calculates some further data from the segs data
public:
    //calc
    float rec_pt;
    float rec_phi;
    float rec_theta;

    //additional data
    float mc_p_corr;
    float mc_pt_corr;

    float mc_p_inv_corr;
    float mc_pt_inv_corr;

    float rec_p_inv;
    float rec_pt_inv;

    float p_inv_abs_error;
    float pt_inv_abs_error;

    SegsTreeReadPlus(TTree *segs)
    : SegsTreeRead(segs) {}

    void calcAdditionalData();
};


#endif //COSMICTRIGGER_SEGSTREEREAD_H

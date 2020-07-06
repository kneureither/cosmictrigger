//
// Created by Konstantin Neureither on 30.06.20.
//

#ifndef COSMICTRIGGER_SLIMSEGSTREE_H
#define COSMICTRIGGER_SLIMSEGSTREE_H

#define GET_ORIGINAL_HITS false

#include <vector>
#include "karimakiHelixfit.h"
#include "TFile.h"
#include "TTree.h"
#include "SegsTreeRead.h"

//for more compact and clearer function call when filling the TTree of the slimmed down segs data
struct SlimSegsMeta {
    unsigned int uEventID;
    unsigned int segsIndex;
    unsigned int runID;
};

//some further data from the kari fit (for calculation)
struct KariFitCalc : KariFit {
    float p = -9999.9;
    float pt = -9999.9;
};


/*This class represents the data structure of the slimmed down simulation data*/
/*From this class there are Write and Read classes derived which handle the write to or read from a ROOT::TTree object*/
class SlimSegsTree {
public:
    TTree *t_slimSegs;
    unsigned int entries;

    //meta data
    unsigned int eventID;
    unsigned int runID;
    unsigned int uEventID;
    unsigned int segsIndex;

    int rec_nhit;
    int rec_ntriplet;

    //monte carlo data
    int mc_tid;
    float mc_p;
    float mc_p_corr;
    float mc_pt;
    int mc_type;
    int mc_pid;

#if GET_ORIGINAL_HITS
    //hit data
    std::vector<float> x00;
    std::vector<float> x01;
    std::vector<float> x10;
    std::vector<float> x11;
    std::vector<float> x20;
    std::vector<float> x21;

    std::vector<float> y00;
    std::vector<float> y01;
    std::vector<float> y10;
    std::vector<float> y11;
    std::vector<float> y20;
    std::vector<float> y21;

    std::vector<float> z00;
    std::vector<float> z01;
    std::vector<float> z10;
    std::vector<float> z11;
    std::vector<float> z20;
    std::vector<float> z21;

    std::vector<int> sid00;
    std::vector<int> sid01;
    std::vector<int> sid10;
    std::vector<int> sid11;
    std::vector<int> sid20;
    std::vector<int> sid21;

    std::vector<float> rec_tan01;
    std::vector<float> rec_tan12;
    std::vector<float> rec_lam01;
    std::vector<float> rec_lam12;
#endif

    //ms fit data
    float rec_p;
    float rec_perr; //vlt
    float rec_r;
    float rec_rerr; //vlt
    float rec_rt;
    float rec_phi;
    float rec_theta;
    float rec_pt;
    float rec_dca_z;
    float rec_dca_r;

    //kari fit data
    float kari_rad;
    float kari_r3d;
    float kari_dca;
    float kari_z0;
    float kari_theta;
    float kari_phi;
    float kari_p;
    float kari_pt;
    float kari_tchi2;
    float kari_lchi2;

    //combined data
    std::vector<double> xp;
    std::vector<double> yp;
    std::vector<double> zp;
    std::vector<int> layerp;
    unsigned int ncombinedhits;

    void reInitializeData();

};

/* Derived class, that handles writing the SlimSegs data to a ROOT::TTree */
class SlimSegsTreeWrite : public SlimSegsTree {
public:
    explicit SlimSegsTreeWrite(TTree * slimSegs);
    void fillData(const SegsTreeReadPlus &Segs,
                  const SlimSegsMeta &Meta,
                  const KariFitCalc &Karires,
                  const unsigned int &ncombinedhits,
                  const std::vector<double> &xps, const std::vector<double> &zps,
                  const std::vector<double> &yps, const std::vector<int> &layerps);
};

/* Derived class, that handles reading the SlimSegs data from a ROOT::TTree */
class SlimSegsRead : public SlimSegsTree {
public:
    explicit SlimSegsRead(TTree * slimSegs);
};


#endif //COSMICTRIGGER_SLIMSEGSTREE_H

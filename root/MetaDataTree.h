//
// Created by Konstantin Neureither on 14.09.20.
//

#ifndef COSMICTRIGGER_METADATATREE_H
#define COSMICTRIGGER_METADATATREE_H

#include "TTree.h"

class MetaDataTreeFile {
public:
    TTree *tT_meta;
    TTree *tT_bgeff;

    //pattern meta data -- once per tree
    int dataset;
    int zBins[3];
    int wBins[3];
    char areaDescript[3][8];
    int mode;
    float efficiency;
    int eventcount;
    std::string mode_description;
    std::string *mode_description_ptr = nullptr;

    int bg_run;
    int max_muon_hits;
    int max_frame_nhits;
    int processed_frames;
    float tb_stopping_eff;
    unsigned int sp_count;
    float sp_target_ratio;

    //bg ana data
    float bg_discr_eff;
    float tb_training_eff;
    std::vector<float> frame_effs;
    std::vector<float> *frame_effs_ptr = nullptr;
    std::vector<int>frame_hits;
    std::vector<int> *frame_hits_ptr = nullptr;

    void reinitializeData();
    MetaDataTreeFile() {}
};

class MetaDataTreeRead : public MetaDataTreeFile {
public:
    MetaDataTreeRead(TTree* tT_meta);
    void setBranches();
    void getEntry(const int &index);
};

class MetaDataTreeWrite : public MetaDataTreeFile {
public:
    MetaDataTreeWrite(TTree *tT_meta, const int dataset, const int* zBins,const int* wBins,
                          char areaDescript[3][8], const int mode, const float efficiency, const int eventcount, std::string mode_description,
                          int bg_run,
                          int max_muon_hits,
                          int max_frame_nhits,
                          int processed_frames,
                          float tb_stopping_eff,
                          unsigned int sp_count,
                          float sp_target_ratio);
};


class BGAnaResTreeWrite : public MetaDataTreeWrite {
//    BGAnaResTreeWrite(TTree *tT_meta, TTree *tT_bgeff,
//                      const int dataset,
//                      const int* zBins, const int* wBins,
//                      char areaDescript[3][8],
//                      const int mode,
//                      const float efficiency,
//                      const int eventcount,
//                      std::string mode_description,
//                      int bg_run,
//                      int max_muon_hits,
//                      int max_frame_nhits,
//                      int processed_frames,
//                      float tb_stopping_eff,
//                      unsigned int sp_count,
//                      float sp_target_ratio) {
//
//        MetaDataTreeWrite(tT_meta, dataset, zBins, wBins, areaDescript[3][8],
//                mode, efficiency, eventcount, mode_description, bg_run,max_muon_hits,
//                int max_frame_nhits,
//                int processed_frames,
//                float tb_stopping_eff,
//                unsigned int sp_count,
//                float sp_target_ratio);
//    }
};

class BGAnaResTreeRead : public MetaDataTreeRead {
public:
    BGAnaResTreeRead(TTree *t_meta, TTree *t_bgeff);
    void setBranches();
    void getEntry(const int &index);
};


#endif //COSMICTRIGGER_METADATATREE_H
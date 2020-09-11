//
// Created by Konstantin Neureither on 11.09.20.
//

#ifndef COSMICTRIGGER_BACKGROUNDDATAFILE_H
#define COSMICTRIGGER_BACKGROUNDDATAFILE_H

#include "basicDefines.h"
#include "../Mu3eCosPat/include/TemplateData.h"


class BackgroundDataFile {
public:
    TTree *tT_meta;
    TTree *tT_frames;

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


    //actual templates -- every entry is a frame
    //hit data
    //TID data
    int rid_len;
    std::vector<unsigned short tid[TID_LEN]> rids;
    int max_cosmic_hits;
    unsigned short tid[TID_LEN] cosmic_hit;
    int nhits;

    void reinitializeData();
};


class BackgroundDataRead : public BackgroundDataFile {

};

class BackgroundDataWrite : public BackgroundDataFile {
    BackgroundDataWrite(TTree *tT_meta, TTree *tT_frames, const int dataset, const int *zBins, const int *wBins, char areaDescript[3][8],
                                                 const int mode, const float efficiency, const int eventcount, std::string mode_description);


};

#endif //COSMICTRIGGER_BACKGROUNDDATAFILE_H

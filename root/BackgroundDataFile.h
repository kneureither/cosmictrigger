//
// Created by Konstantin Neureither on 11.09.20.
//

#ifndef COSMICTRIGGER_BACKGROUNDDATAFILE_H
#define COSMICTRIGGER_BACKGROUNDDATAFILE_H

#include "basicDefines.h"
#include "../Mu3eCosPat/include/TemplateData.h"
#include "TTree.h"

using temidarr = std::vector<unsigned short>;
//typedef std::vector<unsigned short>  temidarr[TID_LEN];


class BackgroundDataFile {
public:
    TTree *tT_meta;
    TTree *tT_frames;

    //consolidated background meta data -- once per tree
    int zBins[3];
    int wBins[3];
    char areaDescript[3][8];
    int mode;
    float efficiency;
    int eventcount;
    std::string mode_description;
    std::string *mode_description_ptr = nullptr;


    //actual background templates -- every entry is a frame

    int bgtid_len;
    std::vector<temidarr> bgtids;
    std::vector<int> tidtypes;
    temidarr cosmic_track;
    int max_cosmic_hits;
    int nhit;

    void reinitializeData();
};


class BackgroundDataRead : public BackgroundDataFile {
private:
    std::vector<temidarr> *bgtids_ptr = nullptr;
    std::vector<int> *tidtype_ptr = nullptr;
public:
    BackgroundDataRead(TTree* tT_meta, TTree* tT_frames);
    void setBranches();
    void getEntry(const int &index);
};

class BackgroundDataWrite : public BackgroundDataFile {
public:
    void fillBGTIDData(std::vector<temidarr> &bgtids, std::vector<int> types, temidarr cosmic_track, int nhits, int max_cosmics_hits);

    BackgroundDataWrite(TTree *tT_meta, TTree *tT_frames, const int *zBins, const int *wBins, char areaDescript[3][8],
                                                 const int mode, const int eventcount, std::string mode_description);
};

#endif //COSMICTRIGGER_BACKGROUNDDATAFILE_H

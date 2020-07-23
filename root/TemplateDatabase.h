//
// Created by Konstantin Neureither on 15.07.20.
//

#ifndef COSMICTRIGGER_TEMPLATEDATABASE_H
#define COSMICTRIGGER_TEMPLATEDATABASE_H
#include "../Mu3eCosPat/TemplateData.h"
#include "TTree.h"

class TemplateDatabase {
public:
    TTree *tT_meta;
    TTree *tT_tid;

    //pattern meta data -- once per tree
    int dataset;
    int zBins[3];
    int wBins[3];
    char areaDescript[3][8];
    int mode;
    float efficiency;
    std::string mode_description;


    //actual templates -- every entry is one template
    //TID data
    int tid_len;
    short tid[TID_LEN];
    std::string tid_repr;
    int frequency;

    //maybe store these as mean and sigma values
    std::vector<int> nhit;
    std::vector<float> pt;
    std::vector<float> phi;
    std::vector<float> theta;
    std::vector<float> dca;

    void reinitializeData();
};

class TemplateDatabaseRead : public TemplateDatabase {
    TemplateDatabaseRead(TTree* tT_meta, TTree* tT_tid);
    void setBranches();
    void getEntry(const int &index);
};

class TemplateDatabaseWrite : private TemplateDatabase {
    TemplateDatabaseWrite(TTree *tT_meta, TTree *tT_tid, const int dataset, const int* zBins,const int* wBins, const char **areaDescript, const int mode,const float efficiency, std::string mode_description);
    void fillTIDData(short *tid, const int tid_len, const int &freq, std::vector<int> &nhit,
            std::vector<float> &p, std::vector<float> &phi, std::vector<float> &theta, std::vector<float> dca);

    void fillTIDData(short *tid, const int tid_len, const int &freq);
};


#endif //COSMICTRIGGER_TEMPLATEDATABASE_H

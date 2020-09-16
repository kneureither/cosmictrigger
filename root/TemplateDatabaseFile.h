//
// Created by Konstantin Neureither on 15.07.20.
//

#ifndef COSMICTRIGGER_TEMPLATEDATABASEFILE_H
#define COSMICTRIGGER_TEMPLATEDATABASEFILE_H
#include "../Mu3eCosPat/include/TemplateData.h"
#include "TTree.h"

class TemplateDatabaseFile {
public:
    TTree *tT_meta;
    TTree *tT_tid;

    //TODO add some important plots, such as eff, tcount, etc.

    //pattern meta data -- once per tree
    int dataset;
    int zBins[3];
    int wBins[3];
    char areaDescript[3][8];
    int mode;
    float training_efficiency;
    float stopping_efficiency;
    int training_events;
    unsigned int template_count;
    std::string mode_description;
    std::string *mode_description_ptr = nullptr;


    //actual templates -- every entry is one template
    //TID data
    int tid_len;
    unsigned short tid[TID_LEN];
    std::string tid_repr;
    char tid_repr_char[TID_LEN*4];
    int frequency;

    //maybe store these as mean and sigma values
    std::vector<int> nhit;
    std::vector<float> pt;
    std::vector<float> phi;
    std::vector<float> theta;
    std::vector<float> dca;
    std::vector<unsigned int> uEventIDs;

    void reinitializeData();
};

class TemplateDatabaseRead : public TemplateDatabaseFile {
public:
    TemplateDatabaseRead(TTree* tT_meta, TTree* tT_tid);
    void setBranches();
    void getEntry(const int &index);
};

class TemplateDatabaseWrite : private TemplateDatabaseFile {
public:
    TemplateDatabaseWrite(TTree *tT_meta, TTree *tT_tid, const int dataset, const int *zBins, const int *wBins,
                          char areaDescript[3][8], const int mode, const float efficiency, const int eventcount,
                          std::string mode_description, unsigned int template_count, const float stopping_efficiency);

    void fillTIDData(unsigned short *tid, const int tid_len, std::string tid_repr, const int &freq, std::vector<int> &nhit,
            std::vector<float> &p, std::vector<float> &phi, std::vector<float> &theta, std::vector<float> dca, std::vector<unsigned int> &uEventIDs);

    void fillTIDData(unsigned short *tid, const int &tid_len, std::string tid_repr, const int &freq);
};


#endif //COSMICTRIGGER_TEMPLATEDATABASEFILE_H

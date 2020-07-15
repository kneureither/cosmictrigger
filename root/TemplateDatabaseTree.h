//
// Created by Konstantin Neureither on 15.07.20.
//

#ifndef COSMICTRIGGER_TEMPLATEDATABASETREE_H
#define COSMICTRIGGER_TEMPLATEDATABASETREE_H
#include "../Mu3eCosPat/TemplateData.h"

class TemplateDatabaseTree {

    //pattern meta data -- once per tree
    int zBins;
    int wBins;
    int mode;
    float efficiency;
    std::vector<int> used_runs;
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
};


#endif //COSMICTRIGGER_TEMPLATEDATABASETREE_H

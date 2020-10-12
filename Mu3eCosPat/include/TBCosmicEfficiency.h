//
// Created by Konstantin Neureither on 12.10.20.
//

#ifndef COSMICTRIGGER_TBCOSMICEFFICIENCY_H
#define COSMICTRIGGER_TBCOSMICEFFICIENCY_H

#include "PatternEngine.h"
#include "TemplateBank.h"


class TBCosmicEfficiency {
private:
    int wBins;
    int zBins;
    int mode;
    float tb_cosmic_eff;


    std::vector<std::vector<unsigned int>> cosmic_spid_tracks;
    PatternEngine PE; //-> initializer list
public:
    TBCosmicEfficiency(int wBins, int zBins, int mode, int cosmic_dataset, std::string pathtoplots);
    float CalcTBCosmicEff(TemplateBank &TB, TIDLoadingFilter);
    float GetTBEff();
};


#endif //COSMICTRIGGER_TBCOSMICEFFICIENCY_H

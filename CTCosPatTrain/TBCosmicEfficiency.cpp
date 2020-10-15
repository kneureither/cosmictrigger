//
// Created by Konstantin Neureither on 12.10.20.
//

#include "TBCosmicEfficiency.h"

TBCosmicEfficiency::TBCosmicEfficiency(int wBins, int zBins, int mode, int cosmic_dataset, std::string pathtoplots)
    : PE {wBins, zBins, pathtoplots} {

}

float TBCosmicEfficiency::CalcTBCosmicEff(TemplateBank &TB, TIDLoadingFilter) {
    return 0;
}

float TBCosmicEfficiency::GetTBEff() {
    return 0;
}



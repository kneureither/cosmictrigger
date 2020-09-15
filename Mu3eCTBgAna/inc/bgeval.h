//
// Created by Konstantin Neureither on 14.09.20.
//

#ifndef COSMICTRIGGER_BGEVAL_H
#define COSMICTRIGGER_BGEVAL_H

#include "utilityFunctions.h"

std::string get_bgevalfile(int bgevents, int run, int dataset, int mode, int wbins, int zbins) {
    return "CosmicBackgroundEval_bgevents_" + get_padded_string(bgevents, 6, '0') + "_run_" + get_padded_string(run, 6, '0') +
        "_dataset" + get_string(dataset) + "_mode" + get_string(mode) +
        "zBins" + get_string(zbins) + "wBins" + get_string(wbins) + "_plots.root";
}

#endif //COSMICTRIGGER_BGEVAL_H

//
// Created by Konstantin Neureither on 14.09.20.
//

#ifndef COSMICTRIGGER_BGEVAL_H
#define COSMICTRIGGER_BGEVAL_H

#include "utilityFunctions.h"

std::string get_bgevalfile(int bgevents, int run, int dataset, int mode, int wbins, int zbins) {
    return "CosmicBackgroundEval_bgevents_" + get_padded_string(bgevents, 6, '0') + "_run_" + get_padded_string(run, 6, '0') +
        "_dataset" + get_string(dataset) + "_" + getfileidtag(mode, wbins, zbins) + "_plots.root";
}

std::string get_bgevalfile(int bg_events, int cosmic_eff_events, int bg_run, int cosmic_eff_dataset, int dataset, int mode, int wbins, int zbins) {
    return "CosmicTBGAna_bkgEv" + get_padded_string(bg_events, 6, '0') +
           "_cosEv" + get_padded_string(cosmic_eff_events, 6, '0') +
           "_bkgrun_" + get_padded_string(bg_run, 6, '0') +
           "_cosdst_" + get_padded_string(cosmic_eff_dataset, 6, '0') +
           "_dataset" + get_string(dataset) + "_" + getfileidtag(mode, wbins, zbins) + "_plots.root";
}

#endif //COSMICTRIGGER_BGEVAL_H

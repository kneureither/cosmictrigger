//
// Created by Konstantin Neureither on 30.10.20.
//

#include <vector>
#include "Configuration.h"
#include "cosmicTemplatesBgAna.h"

int main(int argc, char *argv[]) {

    Configuration CONFIG;
    CONFIG.BGANA_DATAPOINTS();

    int pattern_dataset = CONFIG.dataset;
    int background_dataset = CONFIG.background_run;
    bool append_to_outfile = true;
    std::vector<TIDLoadingFilter> filters = CONFIG.TmplBankFilter.filters;


    for (auto &curves : CONFIG.DBconfigDatapoints) {
        for (auto &config : curves) {
            std::cout << "(STATUS) : Running BG Eval for dataset " << pattern_dataset
                      << " CONFIG: WxZ bins " << config.wbins << " x " << config.zbins
                      << " | SPR " << config.spr << " | SPC " << config.spc << " | EFF " << config.stopp_eff
                      << std::endl;

            cosmicTemplatesBgAna(background_dataset, pattern_dataset, config.spc, config.spr, config.stopp_eff, append_to_outfile,
                                 filters, CONFIG.max_bg_frames, CONFIG.max_cosmic_events, CONFIG.max_bg_frame_nhits, CONFIG.cosmic_testing_dataset);

            append_to_outfile = true;
        }
    }
}
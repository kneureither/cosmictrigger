//
// Created by Konstantin Neureither on 25.08.20.
//
#include "inc/cosmicTemplatesBgAna.h"
#include "../CTCoreModules/inc/TemplateBank.h"
#include <iostream>
#include "Configuration.h"

int main() {
    Configuration CONFIG;
    CONFIG.BGANA();
    //only muon background
    int max_muon_hits = CONFIG.max_bg_frame_nhits;
    int background_dataset = CONFIG.background_run;
    int pattern_dataset = CONFIG.dataset;
    std::vector<float> spratios = CONFIG.sp_ratios;
    std::vector<int> SPcounts = CONFIG.sp_cnt;
    std::vector<float> tb_stopping_efficiencies = CONFIG.stopping_effs;
    std::vector<TIDLoadingFilter> filters = CONFIG.TmplBankFilter.filters;

    int append = 1;
    for(auto &spratio : spratios) {
        for(auto &spcout : SPcounts) {
            for (auto &tb_eff : tb_stopping_efficiencies) {
                std::cout << "STATUS : Running BG Eval for tb_eff: " << tb_eff << std::endl;
                cosmicTemplatesBgAna(background_dataset, pattern_dataset, spcout, spratio, tb_eff, append,
                                     filters, CONFIG.max_bg_frames, CONFIG.max_cosmic_events, CONFIG.max_bg_frame_nhits, CONFIG.cosmic_testing_dataset);
                append = true;
            }
        }
    }


    //each frame gets one cosmic
//    max_muon_hits = 4;
//    cosmicTemplatesBgEval(102, 400, 1, max_muon_hits);

    return 0;
}
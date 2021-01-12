//
// Created by Konstantin Neureither on 27.10.20.
//
#include "inc/cosmicTemplatesBgAna.h"
#include "../CTCoreModules/inc/TemplateBank.h"
#include <iostream>
#include "../CTCoreModules/Configuration.h"

int main() {
    Configuration CONFIG;
    CONFIG.BGANA_MULTI();
    //only muon background

    int bkg_run;
    int max_bg_frame_nhits;
    int max_bkg_frames;

    std::vector<int> bkg_runs = CONFIG.BkgFiles.bg_runs;
    std::vector<int> make_max_nhit_cut = CONFIG.BkgFiles.make_max_nhit_cut;

    int pattern_dataset = CONFIG.dataset;
    std::vector<float> spratios = CONFIG.sp_ratios;
    std::vector<int> SPcounts = CONFIG.sp_cnt;
    std::vector<float> tb_stopping_efficiencies = CONFIG.stopping_effs;
    std::vector<TIDLoadingFilter> filters = CONFIG.TmplBankFilter.filters;


    int append = 1;

    for(int datapoint=0; datapoint < bkg_runs.size(); datapoint++) {
        for(auto &make_cut : make_max_nhit_cut) {
            if(make_cut == 1) {
                max_bg_frame_nhits = CONFIG.BkgFiles.max_nhits[datapoint];
            } else {
                max_bg_frame_nhits = 0;
            }
            bkg_run = bkg_runs[datapoint];
            max_bkg_frames = CONFIG.BkgFiles.max_bkg_frames[datapoint];

            for(auto &spratio : spratios) {
                for(auto &spcout : SPcounts) {
                    for (auto &tb_eff : tb_stopping_efficiencies) {
                        std::cout << "\n--------\n(STATUS) : Running BG Eval for tb_eff " << tb_eff << " | bkg_run " << bkg_run
                        << " | make_cut " << make_cut << " | spr " << spratio << " | spc " << spcout << std::endl;
                        cosmicTemplatesBgAna(bkg_run, pattern_dataset, spcout, spratio, tb_eff, append,
                                             filters, max_bkg_frames, CONFIG.max_cosmic_events, max_bg_frame_nhits, CONFIG.cosmic_testing_dataset);
                        append = true;
                    }
                }
            }
        }
    }

    return 0;
}


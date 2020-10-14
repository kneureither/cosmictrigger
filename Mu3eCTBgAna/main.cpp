//
// Created by Konstantin Neureither on 25.08.20.
//
#include "inc/cosmicTemplatesBgAna.h"
#include "../Mu3eCosPat/include/TemplateBank.h"
#include <iostream>

int main() {

    //only muon background
    int max_muon_hits = 0;
    int background_dataset = 107;
    int pattern_dataset = 11;
    std::vector<float> spratios = {1};
    std::vector<int> SPcounts = {784};
    std::vector<float> tb_stopping_efficiencies = {0.6};
    std::vector<TIDLoadingFilter> filters = {ALL, NO_CENTER, CUT_ON_FREQ};

    int append = 1;
    for(auto &spratio : spratios) {
        for(auto &spcout : SPcounts) {
            for (auto &tb_eff : tb_stopping_efficiencies) {
                std::cout << "STATUS : Running BG Eval for tb_eff: " << tb_eff << std::endl;
                cosmicTemplatesBgAna(background_dataset, pattern_dataset, spcout, spratio, tb_eff, append,
                                     filters);
                append = true;
            }
        }
    }


    //each frame gets one cosmic
//    max_muon_hits = 4;
//    cosmicTemplatesBgEval(102, 400, 1, max_muon_hits);

    return 0;
}
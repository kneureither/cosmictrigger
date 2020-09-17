//
// Created by Konstantin Neureither on 25.08.20.
//
#include "inc/cosmicTemplatesBgEval.h"
#include <iostream>

int main() {

    //only muon background
    int max_muon_hits = 0;
    int background_dataset = 107;
    int pattern_dataset = 11;
    std::vector<float> spratios = {32};
    std::vector<int> SPcounts = {1568};
    std::vector<float> tb_stopping_efficiencies = {0.7};

    int append = 1;
    for(auto &spratio : spratios) {
        for(auto &spcout : SPcounts) {
            for (auto &tb_eff : tb_stopping_efficiencies) {
                std::cout << "STATUS : Running BG Eval for tb_eff: " << tb_eff << std::endl;
                cosmicTemplatesBgEval(background_dataset, pattern_dataset, spcout, spratio, tb_eff, append);
                append = true;
            }
        }
    }


    //each frame gets one cosmic
//    max_muon_hits = 4;
//    cosmicTemplatesBgEval(102, 400, 1, max_muon_hits);

    return 0;
}
//
// Created by Konstantin Neureither on 25.08.20.
//
#include "inc/cosmicTemplatesBgEval.h"
#include <iostream>

int main() {

    //only muon background
    int max_muon_hits = 0;
    int background_dataset = 106;
    float SPratio = 1;
    std::vector<int> SPcounts = {200,400,600,800};
    std::vector<float> tb_stopping_efficiencies = {0.7, 0.75, 0.8, 0.85, 0.9};

    for(auto &spcout : SPcounts) {
        for (auto &tb_eff : tb_stopping_efficiencies) {
            std::cout << "STATUS : Running BG Eval for tb_eff: " << tb_eff << std::endl;
            cosmicTemplatesBgEval(background_dataset, spcout, SPratio, max_muon_hits, tb_eff);
        }
    }

    //each frame gets one cosmic
//    max_muon_hits = 4;
//    cosmicTemplatesBgEval(102, 400, 1, max_muon_hits);

    return 0;
}
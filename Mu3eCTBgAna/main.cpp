//
// Created by Konstantin Neureither on 25.08.20.
//
#include "inc/cosmicTemplatesBgEval.h"

int main() {
    int max_muon_hits;

    //only muon background
    max_muon_hits = 0;
    cosmicTemplatesBgEval(102, 400, 1, max_muon_hits);

    //each frame gets one cosmic
//    max_muon_hits = 4;
//    cosmicTemplatesBgEval(102, 400, 1, max_muon_hits);

    return 0;
}
//
// Created by Konstantin Neureither on 27.08.20.
//

#include "inc/cosmicTemplatesBgPlots.h"
#include <iostream>

int main(int argc, char *argv[]) {
    int dataset;
    std::string filename;

    if(argc < 2) {
        std::cout << "ERROR: Error in argument! Usage: "
                     "Mu3eCosmicTemplatesEvalBGPlots <dataset number> <filename>" << std::endl;
        exit(0);
    } else {
        dataset = atoi(argv[1]);
        filename = std::string(argv[2]);
    }

    std::cout << "Producing combined plots for dataset " << dataset << " filename " << filename << "..." << std::endl;
    makeBgEvalPlots(dataset, 102, filename);
    return 0;
}
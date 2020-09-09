//
// Created by Konstantin Neureither on 27.08.20.
//

#include "inc/cosmicTemplatesBgPlots.h"
#include <iostream>
#include <string>

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

    int runidx = filename.find("run_");
    int run = stoi(filename.substr(runidx + 4, 6));

    std::cout << "Producing combined plots for dataset " << dataset << " filename " << filename << "..." << std::endl;
    makeBgEvalPlots(dataset, run, filename);
    return 0;
}
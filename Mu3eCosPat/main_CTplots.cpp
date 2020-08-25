//
// Created by Konstantin Neureither on 11.08.20.
//
#include "include/makePlots.h"
#include <iostream>

int main(int argc, char *argv[]) {
    int dataset;

    if(argc < 1) {
        std::cout << "ERROR: Error in argument! Usage: "
                     "Mu3eCosPatPlots <dataset number> " << std::endl;
        exit(0);
    } else {
        dataset = atoi(argv[1]);
    }

    std::cout << "Producing combined plots for dataset " << dataset << "..." << std::endl;
    makeCosPatPlots(dataset);
    return 0;
}
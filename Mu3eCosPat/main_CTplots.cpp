//
// Created by Konstantin Neureither on 11.08.20.
//
#include "include/makePlots.h"
#include <iostream>

int main(int argc, char *argv[]) {
    int dataset;
    int combination_id;

    if(argc < 2) {
        std::cout << "ERROR: Error in argument! Usage: "
                     "Mu3eCosPatPlots <dataset number> <sp ratio combination id>" << std::endl;
        exit(0);
    } else {
        dataset = atoi(argv[1]);
        combination_id = atoi(argv[2]);


    }

    std::cout << "Producing combined plots for dataset " << dataset << " with id " << combination_id << "..." << std::endl;
    makeCosPatPlots(dataset, combination_id);
    return 0;
}
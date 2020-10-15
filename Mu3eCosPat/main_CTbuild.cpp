//
// Created by Konstantin Neureither on 11.08.20.
//

#include "include/cosmicTemplatesBuild.h"
#include <iostream>

int main(int argc, char *argv[]) {
    int dataset;
    int spcount;
    float spratio;
    int combination_id;

    if(argc < 4) {
        std::cout << "ERROR: Error in argument! Usage: "
                     "Mu3eCosPatBuild <dataset number> <super pixel count (single station)> <SP ration width(phi):length(z)> <sp ratio id tag>" << std::endl;
        exit(0);
    } else {
        dataset = atoi(argv[1]);
        spcount = atoi(argv[2]);
        spratio = atof(argv[3]);
        combination_id = atoi(argv[4]);

    }

    std::cout << "Building Template Database for dataset " << dataset << "..." << std::endl;
    cosmicTemplatesBuild(dataset, spcount, spratio, combination_id, 0, false);
    return 0;
}
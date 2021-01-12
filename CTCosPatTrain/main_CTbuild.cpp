//
// Created by Konstantin Neureither on 11.08.20.
//

#include "inc/cosmicTemplatesBuild.h"
#include <iostream>

int main(int argc, char *argv[]) {
    /**
     * Builds exactly one template database file, parameters are given as commandline parameters.
     */

    int dataset;
    int spcount;
    float spratio;
    int combination_id;
    float max_eff;

    if(argc < 5) {
        std::cout << "ERROR: Error in argument! \n Usage: "
                     "BuildDBSingleConfig <dataset number> <SPC> <SPR> <stopping efficiency> <id tag outfile>" << std::endl;
        exit(0);
    } else {
        dataset = atoi(argv[1]);
        spcount = atoi(argv[2]);
        spratio = atof(argv[3]);
        max_eff = atof(argv[4]);
        combination_id = atoi(argv[5]);
    }

    std::cout << "Building Template Database for dataset " << dataset << "..." << std::endl;
    cosmicTemplatesBuild(dataset, spcount, spratio, combination_id, max_eff, false);
    return 0;
}
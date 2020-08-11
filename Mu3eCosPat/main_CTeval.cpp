//
// Created by Konstantin Neureither on 11.08.20.
//

#include "include/cosmicTemplatesEval.h"
#include <iostream>


int main(int argc, char *argv[]) {
    int dataset;
    int spcount;
    float spratio;

    if(argc < 3) {
        std::cout << "ERROR: Error in argument! Usage: "
                     "Mu3eCosPatEval <dataset number> <super pixel count (single station)> <SP ration width(phi):length(z)" << std::endl;
        exit(0);
    } else {
        dataset = atoi(argv[1]);
        spcount = atoi(argv[2]);
        spratio = atof(argv[3]);

    }

    std::cout << "Evaluating Template Database for dataset " << dataset << "..." << std::endl;
    cosmicTemplatesEval(dataset, spcount, spratio);
    return 0;
}
//
// Created by Konstantin Neureither on 11.08.20.
//

#include "cosmicEffTemplFilter.h"
#include <iostream>


int main(int argc, char *argv[]) {
    /**
    * This routine is just a basic a test for the method
    * getCosmicSIDtracks(...)
    *
    * TODO Move to test/ folder
    */

    int dataset = 13;
    int spcount = 1024;
    float spratio = 64;
    float stopping_eff = 0.8;

    std::cout << "(INFO)   : Evaluating Template Database for dataset " << dataset << "..." << std::endl;
    cosmicTemplatesEfficiency(dataset, spcount, spratio, stopping_eff);
    return 0;
}
//
// Created by Konstantin Neureither on 11.08.20.
//

#include "cosmicTemplatesEfficiency.h"
#include <iostream>


int main(int argc, char *argv[]) {
    int dataset = 12;
    int spcount = 3200;
    float spratio = 128;
    float stopping_eff = 0.5;

    std::cout << "(INFO)   : Evaluating Template Database for dataset " << dataset << "..." << std::endl;
    cosmicTemplatesEfficiency(dataset, spcount, spratio, stopping_eff);
    return 0;
}
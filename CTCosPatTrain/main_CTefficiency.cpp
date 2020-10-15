//
// Created by Konstantin Neureither on 11.08.20.
//

#include "cosmicTemplatesEfficiency.h"
#include <iostream>


int main(int argc, char *argv[]) {
    int dataset = 9;
    int spcount = 500;
    float spratio = 1;
    float stopping_eff = 0.6;

    std::cout << "(INFO)   : Evaluating Template Database for dataset " << dataset << "..." << std::endl;
    cosmicTemplatesEfficiency(dataset, spcount, spratio, stopping_eff);
    return 0;
}
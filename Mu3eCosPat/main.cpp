//
// Created by Konstantin Neureither on 16.06.20.
//

#include "TemplateBank.h"
#include "PatternEngineSingle.h"
#include "PatternEngine.h"
#include "cosmicTemplatesBuild.h"
#include "cosmicTemplatesEval.h"
#include "makePlots.h"


void processSegsPrototype(int, int);

const bool TEST_PES = false;
const bool TEST_PE = false;
const bool TEST_TB = false;
const bool TEST_SEGS_PROCESS = false;
const bool TEST_buildDB = false;
const bool TEST_readDB = true;


int main(int argc, char *argv[]) {

    std::vector<float> SPratios;

    //this must be set -->
    SPratios.push_back(0.1);
    SPratios.push_back(0.25);
    SPratios.push_back(0.5);
    SPratios.push_back(1);
    SPratios.push_back(2);
    SPratios.push_back(4);
    SPratios.push_back(10);
    int combination_id = 1; //will produce a separate file
    //   <-- up till here.

    int dataset;
    int spcount;

    if(argc < 2) {
        std::cout << "ERROR: Error in argument! Usage: "
                     "Mu3eCosPat <dataset number> <super pixel count (single station)>" << std::endl;
        exit(0);
    } else {
        dataset = atoi(argv[1]);
        spcount = atoi(argv[2]);
    }

    for(int i=0; i<SPratios.size(); i++) {
        std::cout << "Building Template Database for dataset " << dataset << " SPratio? " << SPratios[i] <<"..." << std::endl;
        cosmicTemplatesBuild(dataset, spcount, SPratios[i], combination_id);
    }

    std::cout << "Producing combined plots for dataset " << dataset << "..." << std::endl;
    makeCosPatPlots(dataset, combination_id);



    if(TEST_SEGS_PROCESS) {
        processSegsPrototype(14, 0);
    }

    if(TEST_PE) {
        PatternEngine PE(20, 100, "plots/Mu3eCosPat/");
        PE.displayBinBoundaries();
        PE.closePlot();
    }

    if(TEST_PES) {
        PatternEngineSingle PE(40, 10, 0);
        PE.displayBinBoundaries();
//    PE.testbinSearch();
//    PE.testCoordImpl();

        PE.getSuperPixel(1, 0, -100);
        PE.getSuperPixel(0, 0, -200);
        PE.getSuperPixel(20, 0, 0);
        PE.getSuperPixel(45, 0, 100);
        PE.getSuperPixel(55, 0, 200);
        PE.getSuperPixel(65, 0, 200);

        PE.testSPID();
        PE.displayBinWeightDistribution();

    }

    if(TEST_TB) {
        TemplateBank TB("plots/Mu3eCosPat");
    TB.testTemplateID();
//    TB.testFill();
//        TB.testGetMostPopTemplates();
//        TB.displayTemplatePopulationHistogram("test");

//        TB.testCheck();
    }

    return 0;
}
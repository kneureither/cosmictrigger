//
// Created by Konstantin Neureither on 16.06.20.
//

#include "TemplateBank.h"
#include "PatternEngineSingle.h"
#include "PatternEngine.h"
#include "cosmicTemplatesBuild.h"
#include "cosmicTemplatesEval.h"
#include "makePlots.h"
#include <vector>


void processSegsPrototype(int, int);

const bool TEST_PES = false;
const bool TEST_PE = false;
const bool TEST_TB = false;
const bool TEST_SEGS_PROCESS = false;
const bool TEST_buildDB = false;
const bool TEST_readDB = true;


int main(int argc, char *argv[]) {

    //this must be set -->
//    std::vector<float> sp_ratios = {0.1, 0.25, 0.5, 1, 2, 4, 10};
//    std::vector<float> sp_ratios = {0.25, 1, 4};
    std::vector<float> sp_ratios = {1,2,4,8};
//    std::vector<int> sp_count = {200,400,600,800, 1024};
    std::vector<int> sp_count = {300};
    std::vector<float> stopping_effs = {0.6};
    int combination_id = 4; //will produce a separate file
    int dataset = 9;
    //   <-- up till here.

    std::vector<int> cycle_plotting_order;


    if(argc < 0) {
        std::cout << "ERROR: Error in argument! Usage: "
                     "Mu3eCosPat <dataset number> <super pixel count (single station)>" << std::endl;
        exit(0);
    } else {
//        dataset = atoi(argv[1]);
//        spcount = atoi(argv[2]);
    }

    for(int i=0; i < sp_ratios.size(); i++) {
        for(int j=0; j < sp_count.size(); j++) {
            for(int n=0; n < stopping_effs.size(); n++) {
                std::cout << "Building Template Database for dataset " << dataset << " SPratio=" << sp_ratios[i]
                          << " SPcout=" << sp_count[j] << ".." << std::endl;
                cycle_plotting_order.push_back(cycle_plotting_order.size()+1);
                cosmicTemplatesBuild(dataset, sp_count[j], sp_ratios[i], combination_id, stopping_effs[n]);
            }
        }
    }

    std::cout << "Producing combined plots for dataset " << dataset << "..." << std::endl;
    makeCosPatPlots(dataset, combination_id, cycle_plotting_order);



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
        TemplateBank TB("plots/Mu3eCosPat", 0, 0, 0, 0);
    TB.testTemplateID();
//    TB.testFill();
//        TB.testGetMostPopTemplates();
//        TB.displayTemplatePopulationHistogram("test");

//        TB.testCheck();
    }

    return 0;
}
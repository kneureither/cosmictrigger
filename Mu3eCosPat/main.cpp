//
// Created by Konstantin Neureither on 16.06.20.
//

#include "TemplateBank.h"
#include "PatternEngineSingle.h"
#include "PatternEngine.h"
#include "buildCosmicTrajTemplatesScript.h"


void processSegsPrototype(int, int);

const bool TEST_PES = false;
const bool TEST_PE = false;
const bool TEST_TB = true;
const bool TEST_SEGS_PROCESS = false;
const bool TEST_buildDB = false;


int main(int argc, char *argv[]) {

    if(TEST_SEGS_PROCESS) {
        processSegsPrototype(10, 0);
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
        TemplateBank TB;
//    TB.testTemplateID();
//    TB.testFill();
        TB.testGetMostPopTemplates();
        TB.displayTemplatePopulationHistogram();
    }

    if(TEST_buildDB) {
        buildCosmicTemplatesScript(2);
    }

    return 0;
}
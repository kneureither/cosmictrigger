//
// Created by Konstantin Neureither on 16.06.20.
//

#include <vector>
#include <TFile.h>
#include <TROOT.h>
#include <algorithm>
#include <TTree.h>
#include <string>
#include <TH1.h>
#include <TGraph.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <cmath>
#include <assert.h>
#include "PatternEngine.h"

void processSegsPrototype(int, int);

int main(int argc, char *argv[]) {

    processSegsPrototype(10, 0);

//    PatternEngine PE(40, 10, 0);
//    PE.displayBinBoundaries();
////    PE.testbinSearch();
////    PE.testCoordImpl();
//
//    PE.getSuperPixel(1, 0, -100);
//    PE.getSuperPixel(0, 0, -200);
//    PE.getSuperPixel(20, 0, 0);
//    PE.getSuperPixel(45, 0, 100);
//    PE.getSuperPixel(55, 0, 200);
//    PE.getSuperPixel(65, 0, 200);
//
//    PE.displayBinWeightDistribution();


    return 0;
}
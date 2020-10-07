//
// Created by Konstantin Neureither on 06.08.20.
//

#include "cosmicTemplatesEval.h"

//basic stuff
#include <string>
#include <cassert>
#include <iostream>

//root
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include "cosmicTemplatesBuild.h"
#include "PatternEngine.h"
#include "TemplateBank.h"
#include "utilityFunctions.h"


void cosmicTemplatesEval(const int dataset, unsigned int centralTPcount, float spWZratio) {
    const std::string pathtodata = "data/TemplateData/";
    const std::string pathtoplots = "plots/Mu3eCosPatEval/";
    std::string pathtorunplots = pathtoplots + "dataset_" + get_padded_string(dataset, 3, '0') + "/";
    std::string pathtorundata = pathtodata + "dataset_" + get_padded_string(dataset, 3, '0') + "/";

    const bool MAKE_PLOT = true;
    const int MAX_ENTRIES = 0;
    const bool PRINTS = false;
    const int mode = 0;

    check_create_directory(pathtodata);
    check_create_directory(pathtoplots);
    check_create_directory(pathtorunplots);
    check_create_directory(pathtorundata);

    //Get the Pattern Engine and Template Manager
    const int spWbins = (int) sqrt((float) spWZratio * (float) centralTPcount);
    const int spZbins = (int) sqrt((float) centralTPcount / (float) spWZratio);
    std::cout << "\n -- PE config data:" << std::endl << "  wbins=" << spWbins << std::endl << "  zbins=" << spZbins << std::endl << std::endl;

    TemplateBank TB(pathtorunplots, 0, 0, 0, 0);
    TB.PRINTS = PRINTS;

    TB.readAMfromFile(pathtorundata, 0, ALL);
    TB.getMostPopulatedTemplates(50);

    TB.PlotTemplatePopulationHistogram();
}
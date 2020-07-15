//
// Created by Konstantin Neureither on 09.07.20.
//

//basic stuff
#include <string>
#include <cassert>
#include <iostream>


//root
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"



#include "buildCosmicTrajTemplatesScript.h"
#include "PatternEngine.h"
#include "TemplateBank.h"
#include "../root/SlimSegsTree.h"
#include "utilityFunctions.h"

void getReferenceHits(unsigned int *pInt, int nhit, unsigned int ncombinedhits);

void buildCosmicTemplatesScript(const int dataset) {
    const std::string pathtodata = "data/SlimmedData/";
    const std::string pathtoplots = "plots/Mu3eCosPat/";

    const bool MAKE_PLOT = true;
    const int MAX_ENTRIES = 10;
    const bool PRINTS = true;

    std::string runpadded = get_padded_string(dataset, 6, '0');
    std::string pathtorunplots = pathtoplots + "dataset_" + get_padded_string(dataset, 3, '0') + "/";
    std::string infile = pathtodata + "mu3e_slimmed_segs_" + get_padded_string(dataset, 6, '0') + ".root";

    check_create_directory(pathtodata);
    check_create_directory(pathtoplots);
    check_create_directory(pathtorunplots);

    //Get the Pattern Engine and Template Manager
    PatternEngine PE(20, 100, pathtorunplots);
    TemplateBank TB;
    PE.displayBinBoundaries(); //check if it worked and was initialized correctly.

    // FILE FOR READING
    TFile tinF(infile.c_str());
    if (!tinF.IsOpen()) {
        std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
        exit(0);
    }

    //TODO loop over all segs trees
    TTree *t_slimsegs;
    tinF.GetObject("SlimSegs", t_slimsegs);

    SlimSegsTreeRead SlimSegs = SlimSegsTreeRead(t_slimsegs);


    //stats data
    int p_fail_count = 0;
    int processed_entries = 0;
    int used_entries = 0;
    int too_many = 0;
    int twohittracks = 0;
    int failed_count = 0;


    for (unsigned int entryno = 0; entryno < (MAX_ENTRIES == 0 ? SlimSegs.entries : MAX_ENTRIES); entryno++) {
        SlimSegs.getEntry(entryno);

//        if(SlimSegs.rec_ntriplet < 7) {
//            continue;
//        }

        if(PRINTS) printf("\nBUILD COSMIC TRAJ ENTRY %d\n------------------------------\n\n", entryno);
        std::vector<unsigned int> SPIDs;
        std::vector<float> xpr;
        std::vector<float> ypr;
        std::vector<float> zpr;
        std::vector<int> layerpr;

        unsigned int SPID;
        bool enoughhits = false;

        // TODO
        // -- get hits from SlimSegs
        // -- only use outer layer hits (function needed)
        // -- get TIDs for these hits


        for(int i = 0; i<SlimSegs.layerp.size(); i++) {
            if(PRINTS) printf("  --layer=%d\t x=%f, y=%f, z=%f \n", SlimSegs.layerp[i], SlimSegs.xp[i], SlimSegs.yp[i], SlimSegs.zp[i]);
        }

        enoughhits = getSymmetricRefHits(xpr, ypr, zpr, layerpr, SlimSegs, 0);
        if(!enoughhits) {failed_count++; continue;}

        for(int i = 0; i < xpr.size(); i++) {
            SPID = PE.getSuperPixel(xpr[i], ypr[i], zpr[i]);
            SPIDs.push_back(SPID);
            if(PRINTS) printf(" -- SPIDS =%#X \n", SPIDs[i]);
        }

    }
}

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
    const int MAX_ENTRIES = 0;
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


    for (unsigned int entryno = 0; entryno < (MAX_ENTRIES == 0 ? SlimSegs.entries : MAX_ENTRIES); entryno++) {
        SlimSegs.getEntry(entryno);

//        if(SlimSegs.rec_ntriplet < 7) {
//            continue;
//        }

        if(PRINTS) printf("\nBUILD COSMIC TRAJ ENTRY %d\n------------------------------\n\n", entryno);
        std::vector<unsigned int> SPIDs;
        unsigned int SPID;

        // TODO
        // -- get hits from SlimSegs
        // -- only use outer layer hits (function needed)
        // -- get TIDs for these hits


        if( 2 <= SlimSegs.rec_ntriplet && SlimSegs.rec_ntriplet <= 10) {
            for(int i = 0; i<SlimSegs.layerp.size(); i++) {
                printf("  --layer=%d\t x=%f, y=%f, z=%f \n", SlimSegs.layerp[i], SlimSegs.xp[i], SlimSegs.yp[i], SlimSegs.zp[i]);
                if(SlimSegs.layerp[i] == 2 || SlimSegs.layerp[i] == 3) {
//                    SPID = PE.getSuperPixel(SlimSegs.xp[i], SlimSegs.yp[i], SlimSegs.zp[i]);
//                    SPIDs.push_back(SPID);
//                    printf("SID=%d, SIDhex=%#X, SIDs.size()=%d\n", SPID, SPID, SPIDs.size());
                }
            }
        } else {
            std::cout << "WARNING : triplet count out of bounds - is: " << SlimSegs.rec_ntriplet << std::endl;
        }

        getSymmetricRefHits(SPIDs, SlimSegs, PE);

        for(int i=0; i<SPIDs.size(); i++) {
            printf(" -- SPIDS =%#X \n", SPIDs[i]);
        }

//        assert(SPIDs.size() == 4);
//        assert(SlimSegs.rec_ntriplet < 7);

//        if(SlimSegs.rec_ntriplet == 2) { // 4-hit tracks
//
//        } else if (SlimSegs.rec_ntriplet == 3) { // 6-hit tracks
//
//        } else if (SlimSegs.rec_ntriplet == 4) { // 8-hit tracks
//
//        } else {
//            std::cout << "WARNING : More than 4 triplets" << std::endl;
//        }

    }
}

//
// Created by Konstantin Neureither on 11.09.20.
//

#include "TFile.h"
#include "TH1F.h"

#include <cassert>
#include <map>
#include <stdlib.h>

#include "../inc/cosmicTemplatesBgEval.h"
#include "utilityFunctions.h"
#include "Mu3eTree.h"
#include "../Mu3eCosPat/include/PatternEngine.h"
#include "../Mu3eCosPat/include/TemplateBank.h"
#include "../Mu3eCosPat/include/TemplateData.h"
#include "bgDataCombinatorics.h"
#include "../root/BackgroundDataFile.h"


void backgroundCombinatorics(const int run, unsigned int centralTPcount, float spWZratio) {
    /*
     * Read and analyse the mu3e mc hits
     * get the hits in xyz
     * assign sps to these hits
     * initialize template bank
     *
     */

    int MAX_ENTRIES = 0;
    int MAX_MUON_HITS = 0;
    int MAX_COSMIC_HITS = 0; //not implemented
    const bool RECREATE_FILE = false;
    const int MUONTYPE = 1;
    const int PRINTS = false;
    const int mode = 0;

    const std::string pathtoBGdata = "data/SimulationData/";
    const std::string pathtoTemplateData = "data/TemplateData/";
    const std::string pathtoProcessedBGdata = "data/BackgroundData/run_" + get_padded_string(run, 6, '0') +"/";
    const std::string pathtoplots = "plots/BackgroundCombinatorics/";
    const std::string infile = pathtoBGdata + "mu3e_run_" + get_padded_string(run, 6, '0') + ".root";
    const std::string pathtooutfile = pathtoplots + "bgrun_" + get_padded_string(run, 3, '0') + "/"; //this is where the root file is stored
    const std::string pathtorunplots = pathtooutfile +"/PDF/"; //this is where the pdf files are stored
//    const std::string pathtodatasettemplatedata = pathtoTemplateData + "dataset_" + get_padded_string(dataset, 3, '0') + "/";    const std::string pathtodatasettemplatedata = pathtoTemplateData + "dataset_" + get_padded_string(dataset, 3, '0') + "/";

    check_create_directory(pathtoBGdata);
    check_create_directory(pathtoTemplateData);
    check_create_directory(pathtoplots);
    check_create_directory(pathtorunplots);
    check_create_directory(pathtoProcessedBGdata);

    // FILE FOR READING
    TFile tinF(infile.c_str());
    if (!tinF.IsOpen()) {
        std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
        exit(0);
    }

    TTree *t_mu3e;
    tinF.GetObject("mu3e", t_mu3e);
    TTree *t_mu3e_mchits;
    tinF.GetObject("mu3e_mchits", t_mu3e_mchits);

    //Get the Pattern Engine and Template Manager
    const int spWbins = (int) sqrt((float) spWZratio * (float) centralTPcount);
    const int spZbins = (int) sqrt((float) centralTPcount / (float) spWZratio);
//    std::cout << "\n -- PE config data:" << std::endl << "  wbins=" << spWbins << std::endl << "  zbins=" << spZbins << std::endl << std::endl;

    Mu3eTree Mu3e = Mu3eTree(t_mu3e);
    Mu3eMChitsTree Mu3eMChits = Mu3eMChitsTree(t_mu3e_mchits);

    int bg_events = (MAX_ENTRIES == 0 ? Mu3e.my_entries : MAX_ENTRIES);
    int processed_frames = 0;

    const std::string bgdatatoutfile = pathtoProcessedBGdata + "BackgroundCombinatoricsData_run_" + get_padded_string(run, 6, '0') + "mode" + get_string(mode) + "zBins" + get_string(spZbins) + "wBins" + get_string(spWbins) + ".root";

    // FILE FOR WRITING BACKGROUND DATA
    TFile toutF(bgdatatoutfile.c_str());
    if (!toutF.IsOpen()) {
        std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
        exit(0);
    }

    TTree t_meta("t_meta", "Meta data of Pattern Engine and Background Combinatorics Computation");
    TTree t_frames("t_frames", "Full combinatorics of Background simulation frames as tids");

    PatternEngine PE(spWbins, spZbins, pathtorunplots);
    BackgroundDataWrite BGFile = BackgroundDataWrite(&t_meta, &t_frames, PE.ZBins, PE.WBins, PE.areaDescript, mode, bg_events, "normal_testing_mode");



    for(int frame=0; frame <= bg_events; frame++) {
        Mu3e.getEntry(frame);
        if (PRINTS) Mu3e.Print(frame);

        if(bg_events % 1000 == 0) std::cout << "STATUS : --> processing frame " << frame << std::endl;

        processed_frames++;

        std::vector<BGhit> bgframehits;
        BGSortedHits hits;
        std::vector<temidarr> TIDS;
        std::vector<int> types;
        temidarr cosmic_track = {0,0,0,0};
        std::map<unsigned short, int> SIDMem;
        BGhit BGHIT;
        unsigned short SID;
        int innerhits = 0;
        int outerhits = 0;
        int doublesid = 0;

        //get all hits from one frame
        assert(Mu3e.Nhit == Mu3e.hit_mc_i->size());
        for (int hitno = 0; hitno < Mu3e.Nhit; hitno++) {
            Mu3eMChits.getEntry((*Mu3e.hit_mc_i)[hitno]);
            if (PRINTS) Mu3eMChits.Print((*Mu3e.hit_mc_i)[hitno]);
            BGHIT.fill(Mu3eMChits.pos_g_x, Mu3eMChits.pos_g_y, Mu3eMChits.pos_g_z, 0);
            bgframehits.push_back(BGHIT);

            SID = (unsigned short) PE.getSuperPixel(BGHIT.x, BGHIT.y, BGHIT.z);
            if(SIDMem.count(SID) == 0) {
                SIDMem[SID] = 1;
            } else {
                doublesid++;
                continue;
            }

            int layer = PE.getLayerFromSPID(SID);

            if (layer == 3 && BGHIT.y >= 0) {
                hits.h0.push_back(SIDtype(SID, BGHIT.type));
            } else if (layer == 2 && BGHIT.y >= 0) {
                hits.h1.push_back(SIDtype(SID, BGHIT.type));
            } else if (layer == 2 && BGHIT.y < 0) {
                hits.h2.push_back(SIDtype(SID, BGHIT.type));
            } else if (layer == 3 && BGHIT.y < 0) {
                hits.h3.push_back(SIDtype(SID, BGHIT.type));
            } else {
                innerhits++;
            }
        }

        int cosmiccount[4] = {0, 0, 0, 0};
        int toomanycosmicscount = 0;
        outerhits = bgframehits.size() - innerhits;

//        int randindex = (rand() % cosmicpool);
//        temid COSMICTID = COSMICTIDs[randindex];

//        if(PRINTS) std::cout << "   -> got cosmic at index " << randindex << " TID=" << COSMICTID.toString() << std::endl;

//        if(SIDMem.count(COSMICTID.HIDS[0]) == 0) hits.h0.push_back(SIDtype(COSMICTID.HIDS[0], MUONTYPE));
//        if(SIDMem.count(COSMICTID.HIDS[1]) == 0) hits.h1.push_back(SIDtype(COSMICTID.HIDS[1], MUONTYPE));
//        if(SIDMem.count(COSMICTID.HIDS[2]) == 0) hits.h2.push_back(SIDtype(COSMICTID.HIDS[2], MUONTYPE));
//        if(SIDMem.count(COSMICTID.HIDS[3]) == 0) hits.h3.push_back(SIDtype(COSMICTID.HIDS[3], MUONTYPE));

//        hits.h0.push_back(SIDtype(COSMICTID.HIDS[0], MUONTYPE));
//        hits.h1.push_back(SIDtype(COSMICTID.HIDS[1], MUONTYPE));
//        hits.h2.push_back(SIDtype(COSMICTID.HIDS[2], MUONTYPE));
//        hits.h3.push_back(SIDtype(COSMICTID.HIDS[3], MUONTYPE));

        //Go through all possible combinations of hits and create corresponding Template IDs
        //for all hits in upper layer 3
        for (const auto &h0 : hits.h0) {
            cosmiccount[0] = (h0.type == MUONTYPE ? 1 : 0);

            //for all hits in upper layer 2
            for (const auto &h1 : hits.h1) {
                cosmiccount[1] = (h1.type == MUONTYPE ? 1 : 0);

                //for all hits in lower layer 2
                for (const auto &h2 : hits.h2) {
                    cosmiccount[2] = (h2.type == MUONTYPE ? 1 : 0);

                    //for all hits in lower layer 3
                    for (const auto &h3 : hits.h3) {
                        cosmiccount[3] = (h3.type == MUONTYPE ? 1 : 0);

                        if (cosmiccount[0] + cosmiccount[1] + cosmiccount[2] + cosmiccount[3] == MAX_MUON_HITS) {
                            temidarr TID1 = {h0.SID, h1.SID, h2.SID, h3.SID};
                            TIDS.push_back(TID1);
                            types.push_back(0);
                        } else {
                            toomanycosmicscount++;
                        }

                    }
                }
            }
        }
        BGFile.fillBGTIDData(TIDS, types, cosmic_track, bgframehits.size(), MAX_COSMIC_HITS);
    }

    std::cout << "STATUS : Background run " << run << " processed." << std::endl;
    std::cout << "INFO   : Processed " << processed_frames << " of " << Mu3e.my_entries << " Events.";

    PE.displayBinBoundaries();
    PE.displayBinWeightDistribution();
    PE.closePlot();

    toutF.Write();
    toutF.Close();

    std::cout << "STATUS : background templates were written to " + bgdatatoutfile << std::endl;
}

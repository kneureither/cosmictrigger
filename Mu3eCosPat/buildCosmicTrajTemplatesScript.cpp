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
    const bool PRINTS = false;

    std::string runpadded = get_padded_string(dataset, 6, '0');
    std::string pathtorunplots = pathtoplots + "dataset_" + get_padded_string(dataset, 3, '0') + "/";
    std::string infile = pathtodata + "mu3e_slimmed_segs_" + get_padded_string(dataset, 6, '0') + ".root";
//    std::string infile = pathtodata + "mu3e_test_slimmed_file_000000.root";

    check_create_directory(pathtodata);
    check_create_directory(pathtoplots);
    check_create_directory(pathtorunplots);

    //Get the Pattern Engine and Template Manager
    PatternEngine PE(10, 10, pathtorunplots);
    PE.PRINTS = PRINTS;
    TemplateBank TB(pathtorunplots);
    TB.PRINTS = PRINTS;

    // FILE FOR READING
    TFile tinF(infile.c_str());
    if (!tinF.IsOpen()) {
        std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
        exit(0);
    }

//    tinF.ls();
    int treecount = 0;
    TList *list = tinF.GetListOfKeys();
    TIter iter(list->MakeIterator());
    while(TObject* obj = iter()){
        treecount++;
//        TKey* theKey = (TKey*)obj;
//        cout<< "name/type of Key is: "<<theKey->GetName() << " / " << theKey->GetClassName() << endl;
//        theKey->Class()->Dump();
    }

    //stats data
    int p_fail_count = 0;
    int processed_entries = 0;
    int used_entries = 0;
    int too_many = 0;
    int twohittracks = 0;
    int failed_count = 0;

    for(int treeid = 1; treeid <= treecount; treeid++) {

        TTree *t_slimsegs;
        std::string treename = "SlimSegs;" + get_string(treeid);
        tinF.GetObject(treename.c_str(), t_slimsegs);

        SlimSegsTreeRead SlimSegs = SlimSegsTreeRead(t_slimsegs);

        for (unsigned int entryno = 0; entryno < (MAX_ENTRIES == 0 ? SlimSegs.entries : MAX_ENTRIES); entryno++) {
            SlimSegs.getEntry(entryno);
            processed_entries++;

            std::vector<unsigned int> SPIDs;
            std::vector<float> xpr;
            std::vector<float> ypr;
            std::vector<float> zpr;
            std::vector<int> layerpr;

//        if(SlimSegs.rec_ntriplet < 7) {
//            continue;
//        }

            if(PRINTS) printf("\nSLIM SEGS;%d ENTRY %d\n------------------------------\n\n",treeid, entryno);

            unsigned int SPID;
            bool enoughhits = 0;


            enoughhits = getSymmetricRefHits(xpr, ypr, zpr, layerpr, SlimSegs, 4 - (TID_LEN / 2));
            if(!enoughhits) {
                failed_count++;
                continue;
            }

            for(int i = 0; i<SlimSegs.layerp.size(); i++) {
//                SPID = PE.getSuperPixel(SlimSegs.xp[i], SlimSegs.yp[i], SlimSegs.zp[i]);
                if(PRINTS) printf("  --layer=%d\t x=%f, y=%f, z=%f \n", SlimSegs.layerp[i], SlimSegs.xp[i], SlimSegs.yp[i], SlimSegs.zp[i]);
            }
//            printf("\n");

            used_entries++;

            for(int i = 0; i < xpr.size(); i++) {
                SPID = PE.getSuperPixel(xpr[i], ypr[i], zpr[i]);
                SPIDs.push_back(SPID);
                if(PRINTS) printf(" -- SPIDS =%#X \n", SPIDs[i]);
            }

            TB.fillTemplate(&SPIDs[0], SPIDs.size(), SlimSegs.kari_p, SlimSegs.kari_dca, SlimSegs.kari_phi, SlimSegs.kari_theta);

        }
    }

    tinF.Close();


    //open new TFile for plots
    TFile * tF = new TFile((pathtorunplots +"TemplateBank_dataset_" +
            get_padded_string(dataset, 3, '0') + "_plots.root").c_str(), "recreate");
    if (!tF->IsOpen()) {
        std::cout << "[ERROR] File " << tF->GetName() << " is not open!" << std::endl;
    }

    //add some meta data for the
    int datast = dataset;
    TTree tT_met("MetadataTree","Metadata associated with these plots (SID config and dataset)");
    tT_met.Branch("dataset", &datast, "dataset/I");
    tT_met.Branch("area0Description", &PE.areaDescript[0], "area0Description/C");
    tT_met.Branch("area1Description", &PE.areaDescript[1], "area1Description/C");
    tT_met.Branch("area2Description", &PE.areaDescript[2], "area2Description/C");
    tT_met.Branch("wBins0", &PE.WBins[0], "wBins0/I");
    tT_met.Branch("wBins1", &PE.WBins[1], "wBins1/I");
    tT_met.Branch("wBins2", &PE.WBins[2], "wBins2/I");
    tT_met.Branch("zBins0", &PE.ZBins[0], "zBins0/I");
    tT_met.Branch("zBins1", &PE.ZBins[1], "zBins1/I");
    tT_met.Branch("zBins2", &PE.ZBins[2], "zBins2/I");
    tT_met.Branch("mode", &PE.mode, "mode/I");
    tT_met.Fill();
    tT_met.Write();

    //Make plots and add them to the root file
    PE.displayBinBoundaries(); //check if it worked and was initialized correctly.
    PE.displayBinWeightDistribution();
    PE.closePlot();
    TB.displayTemplatePopulationHistogram(PE.getModeTag());
    TB.displayEfficiency(PE.getModeTag());

    //close plot root file
    tF->Close();

    TB.getMostPopulatedTemplates(50);
    TB.writeAMtoFile(pathtorunplots, PE.ZBins, PE.WBins, PE.areaDescript, datast, PE.mode, "testing_mode_descr");

    std::cout << "\n\n>>>>> GENERAL STATS <<<<<\n\n";

    std::cout << " - total entries processed: " << processed_entries << endl;
    std::cout << "   of which " << used_entries << " (" << used_entries / (float) processed_entries *100 << "%) were used" << std::endl;
    std::cout << " - total trees processed: " << treecount << endl;
}

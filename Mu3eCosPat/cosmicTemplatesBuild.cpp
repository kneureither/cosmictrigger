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

#include "cosmicTemplatesBuild.h"
#include "PatternEngine.h"
#include "TemplateBank.h"
#include "utilityFunctions.h"
#include "MetaDataTree.h"

void getReferenceHits(unsigned int *pInt, int nhit, unsigned int ncombinedhits);

void cosmicTemplatesBuild(const int dataset, unsigned int centralTPcount, float spWZratio, int combination_id,
                          float max_efficiency) {
    const std::string pathtodata = "data/SlimmedData/";
    const std::string pathtoplots = "output/Mu3eCosPat/";

    const bool MAKE_PLOT = true;
    const int MAX_ENTRIES = 0;
    const bool PRINTS = false;
    float MAX_EFFICIENCY = max_efficiency;
    const bool WRITE_DB_FILE = true;

    std::string runpadded = get_padded_string(dataset, 6, '0');
    std::string pathtorunplots = pathtoplots + "dataset_" + get_padded_string(dataset, 3, '0') + "/";
    std::string pathtotemplatedb = "data/TemplateData/dataset_" + get_padded_string(dataset, 3, '0') + "/";
    std::string infile = pathtodata + "mu3e_slimmed_segs_" + get_padded_string(dataset, 6, '0') + ".root";
//    std::string infile = pathtodata + "mu3e_test_slimmed_file_000000.root";

    check_create_directory(pathtodata);
    check_create_directory(pathtoplots);
    check_create_directory(pathtorunplots);
    check_create_directory(pathtotemplatedb);

    //Get the Pattern Engine and Template Manager
    const int spWbins = (int) sqrt((float) spWZratio * (float) centralTPcount);
    const int spZbins = (int) sqrt((float) centralTPcount / (float) spWZratio);

    std::cout << "\n -- PE config data:" << std::endl << "  wbins=" << spWbins << std::endl << "  zbins=" << spZbins << std::endl << std::endl;

    PatternEngine PE(spWbins, spZbins, pathtorunplots);
    PE.PRINTS = PRINTS;
    TemplateBank TB(pathtorunplots, MAX_EFFICIENCY);
    TB.PRINTS = PRINTS;

    // FILE FOR READING
    TFile tinF(infile.c_str());
    if (!tinF.IsOpen()) {
        std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
        exit(0);
    }

    //stats data
    int p_fail_count = 0;
    int processed_entries = 0;
    int used_entries = 0;
    int too_many = 0;
    int twohittracks = 0;
    int failed_count = 0;
    unsigned int runID;
    bool stopping_point_reached = false;


    //get different cycles of trees in file
    int treecount = 0;
    int cycle = 0;
    std::string tree;
    std::vector<int> processed_run_ids;

    TList *list = tinF.GetListOfKeys();
    TIter iter(list->MakeIterator());

    //iterate over all cycles of trees separately
    while(TObject* obj = iter()){
        treecount++;
        TKey* theKey = (TKey*)obj;
        tree = theKey->GetName();
        cycle = theKey->GetCycle();
//        cout<< "name/type of Key is: "<<theKey->GetName() << " / " << theKey->GetClassName() << endl;
//        theKey->Class()->Dump();

        TTree * t_slimsegs;
        std::string treename = tree + ";" + get_string(cycle);
        tinF.GetObject(treename.c_str(), t_slimsegs);
        t_slimsegs->SetBranchAddress("runID", &runID);
        t_slimsegs->GetEntry(0);
        std::cout << "(INFO)   : Processing tree " << tree << " cycle " << cycle << " | run " << runID << " | TB eff " << TB.getEfficiency() << std::endl;
        std::cout << "(INFO)   : sp count " << centralTPcount << " | wbins " << spWbins << " | zbins " << spZbins << " | target eff " << max_efficiency << std::endl;

        SlimSegsTreeRead SlimSegs = SlimSegsTreeRead(t_slimsegs);

        SlimSegs.getEntry(0);
        if (std::find(processed_run_ids.begin(), processed_run_ids.end(),SlimSegs.runID)!=processed_run_ids.end()) {
            std::cout << "(WARNING): runID " << SlimSegs.runID << " has already been processed -- skipping cycle "<< cycle << std::endl;
            continue;
        } else {
            processed_run_ids.push_back(SlimSegs.runID);
        }
        int max_entry = (MAX_ENTRIES == 0 ? SlimSegs.entries : MAX_ENTRIES);

        for (unsigned int entryno = 0; entryno < max_entry; entryno++) {
            SlimSegs.getEntry(entryno);
            processed_entries++;

            if(entryno % 1000 == 0){
                int MAX_LEN = 50;
                float prog_perc = entryno /  (float) max_entry;
                std::string prog_bar_fill((int) (MAX_LEN * prog_perc), '=');
                std::string prog_bar_empty((int) (MAX_LEN * (1-prog_perc)), ' ');
                std::cout << "\r(STATUS) : " << "[" << prog_bar_fill << ">" << prog_bar_empty << "] ";
                std::cout << entryno /  (float) max_entry * 100 << "% | TB eff " << TB.getEfficiency() * 100 << "%" << std::flush;
            }

            std::vector<unsigned int> SPIDs;
            std::vector<float> xpr;
            std::vector<float> ypr;
            std::vector<float> zpr;
            std::vector<int> layerpr;

//        if(SlimSegs.rec_ntriplet < 7) {
//            continue;
//        }

            if(PRINTS) printf("\nSLIM SEGS;%d ENTRY %d\n------------------------------\n\n",treecount, entryno);

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
//            if(PRINTS) printf("\n");

            used_entries++;

            for(int i = 0; i < xpr.size(); i++) {
                SPID = PE.getSuperPixel(xpr[i], ypr[i], zpr[i]);
                SPIDs.push_back(SPID);
                if(PRINTS) printf(" -- SPIDS =%#X \n", SPIDs[i]);
            }

            stopping_point_reached = TB.fillTemplate(&SPIDs[0], SPIDs.size(), SlimSegs.kari_p, SlimSegs.kari_dca, SlimSegs.kari_phi, SlimSegs.kari_theta);

            if(stopping_point_reached){
                std::cout << std::endl << "(STATUS) : Reached efficiency stopping point! MAX_EFF=" << MAX_EFFICIENCY <<  " TB EFF=" << TB.getEfficiency() << std::endl;
                break;
            }
        }
        std::cout << std::endl;
        if(stopping_point_reached) break;
    }


    tinF.Close();


    //open new TFile for plots
    TFile * tF = new TFile((pathtorunplots +"TemplateBank_dataset_" +
            get_padded_string(dataset, 3, '0') + "_id" + get_padded_string(combination_id, 3, '0') + "_plots.root").c_str(), "update");
    if (!tF->IsOpen()) {
        std::cout << "[ERROR] File " << tF->GetName() << " is not open!" << std::endl;
    }

    //add some meta data for the
    int datast = dataset;
    float eff = TB.getEfficiency();
    int templatecount = TB.getTemplateCount();
    TTree tT_met("MetadataTree","Metadata associated with these plots (PE, TB config and dataset)");
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
    tT_met.Branch("training_efficiency", &eff, "training_efficiency/F");
    tT_met.Branch("template_count", &templatecount, "template_count/I");
    tT_met.Branch("processed_events", &processed_entries, "processed_events/I");
    tT_met.Branch("sp_count", &centralTPcount, "spcount/i");
    tT_met.Branch("sp_target_ratio", &spWZratio, "sp_target_ratio/F");
    tT_met.Fill();
    tT_met.Write();

    //Make plots and add them to the root file
    PE.displayBinBoundaries(); //check if PE worked and was initialized correctly.
    PE.displayBinWeightDistribution();
    PE.closePlot();
    TB.displayTemplatePopulationHistogram(PE.getModeTag());
    TB.displayEfficiency(PE.getModeTag());

    //close plot root file
    tF->Close();

    //do some template bank stuff, eg. write the TemplateBank to a root database file
//    TB.getMostPopulatedTemplates(50);

    //FIXME: ZBins and wbins are in wrong order!! TDB filenames are wrong --> now corrected, but problem in db
    if(WRITE_DB_FILE) TB.writeAMtoFile(pathtotemplatedb, PE.ZBins, PE.WBins, PE.areaDescript, datast, PE.mode, "testing_mode_descr");

    TB.plotFreqTimesTemplatecount(PE.getModeTag());

    std::cout << "\n\n[============ GENERAL STATS ============]\n";

    std::cout << " - total entries processed: " << processed_entries << endl;
    std::cout << "   of which " << used_entries << " (" << used_entries / (float) processed_entries *100 << "%) were used" << std::endl;
    std::cout << " - total trees processed: " << treecount << endl << endl << endl;
}

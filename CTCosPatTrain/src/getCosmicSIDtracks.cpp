//
// Created by Konstantin Neureither on 09.01.21.
//

#include "getCosmicSIDtracks.h"

//basic stuff
#include <string>
#include <cassert>
#include <iostream>

//root
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TLatex.h"

#include "cosmicTemplatesBuild.h"
#include "PatternEngine.h"
#include "TemplateBank.h"
#include "utilityFunctions.h"
#include "plots.h"


std::vector<std::vector<unsigned int>> getCosmicSIDtracks(int cosmic_testing_dataset, int centralTPcount, float spWZratio, int mode) {
    /**
     * reads a file from the Mu3e Sim in Cosmic Mode (Trirec Output), computes the cosmic TIDs and saves the result to a file, which can later be used for
     * efficiency evalutation. If the file already exists, then it just reads the file and returns the vector.
     *
     * This function can be used, if the TIDs from a cosmic simulation are needed, e.g. during Template Bank Training.
     * Especially if the same cosmic data is needed multiple times, it offers a speed advantage, as in only reads and
     * comines the actual cosmic sim file and then automatically refers to its own small combined file
     * "CosmicSIDtracks_[fileidtag].root" stored in data/CosmicSIDtrackData
     *
     */

    const bool PRINTS = false;

    /// Maybe these paths must be modified, depending on your environment.
    const std::string path_to_cosmic_inputdata = "data/SlimmedData/";
    const std::string path_to_cosmic_tids = "data/CosmicSIDtrackData/";
    const std::string path_to_plots = "output/5_CosmicSIDs/";


    const std::string path_to_plots_dataset =
            path_to_plots + "dataset_" + get_padded_string(cosmic_testing_dataset, 3, '0') + "/";
    const std::string path_to_cosmic_tids_dataset =
            path_to_cosmic_tids + "dataset_" + get_padded_string(cosmic_testing_dataset, 3, '0') + "/";

    check_create_directory(path_to_cosmic_inputdata);
    check_create_directory(path_to_cosmic_tids);
//    check_create_directory(path_to_plots);
//    check_create_directory(path_to_plots_dataset);
    check_create_directory(path_to_cosmic_tids_dataset);

    //Get the Pattern Engine and Template Manager
    const int spWbins = (int) sqrt((float) spWZratio * (float) centralTPcount);
    const int spZbins = (int) sqrt((float) centralTPcount / (float) spWZratio);
    PatternEngine PE(spWbins, spZbins, path_to_plots_dataset);

    std::string infile =
            path_to_cosmic_inputdata + "mu3e_slimmed_segs_" + get_padded_string(cosmic_testing_dataset, 6, '0') +
            ".root";
    std::string cosmic_sid_file = path_to_cosmic_tids_dataset + "CosmicSIDtracks_" +
                                  getfileidtag(cosmic_testing_dataset, mode, spWbins, spZbins) + ".root";
    std::vector <std::vector<unsigned int>> cosmic_spid_tracks;

    //First check if a file exists already
    TFile tcosTIDF(cosmic_sid_file.c_str());
    if (!tcosTIDF.IsOpen()) {
        // If there is no precomputed file, a Mu3eTriRec (cosmic mode) file is used to compute the TIDs.

        std::cout << "(INFO)   : No TID File " << cosmic_sid_file << " exists. Starting to create data!" << std::endl;

        //open file again with writing access
        tcosTIDF.Close();
        TFile tcosTIDF(cosmic_sid_file.c_str(), "recreate");

        // FILE FOR READING
        TFile tinF(infile.c_str());
        if (!tinF.IsOpen()) {
            std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
            exit(0);
        } else {
            std::cout << "(INFO)   : cosmic testing data from file " << infile << std::endl;
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
        while (TObject * obj = iter()) {
            treecount++;
            TKey *theKey = (TKey *) obj;
            tree = theKey->GetName();
            cycle = theKey->GetCycle();
            TTree *t_slimsegs;
            std::string treename = tree + ";" + get_string(cycle);
            tinF.GetObject(treename.c_str(), t_slimsegs);
            t_slimsegs->SetBranchAddress("runID", &runID);
            t_slimsegs->GetEntry(0);
            std::cout << "(INFO)   : Processing tree " << tree << " cycle " << cycle << " | run " << runID << std::endl;
            std::cout << "(INFO)   : sp count " << centralTPcount << " | wbins " << spWbins << " | zbins " << spZbins
                      << std::endl;

            SlimSegsTreeRead SlimSegs = SlimSegsTreeRead(t_slimsegs);

            SlimSegs.getEntry(0);
            if (std::find(processed_run_ids.begin(), processed_run_ids.end(), SlimSegs.runID) !=
                processed_run_ids.end()) {
                std::cout << "(WARNING): runID " << SlimSegs.runID << " has already been processed -- skipping cycle "
                          << cycle << std::endl;
                continue;
            } else {
                processed_run_ids.push_back(SlimSegs.runID);
            }


            int max_entry = SlimSegs.entries;
            for (unsigned int entryno = 0; entryno < max_entry; entryno++) {
                SlimSegs.getEntry(entryno);
                processed_entries++;

                //progress bar
                print_status_bar(entryno, max_entry, "building cosmic sid tracks",
                                 "processed entries " + get_string(processed_entries));

                std::vector<unsigned int> SPIDs;
                std::vector<float> xpr;
                std::vector<float> ypr;
                std::vector<float> zpr;
                std::vector<int> layerpr;

                unsigned int SPID;
                bool enoughhits = 0;

                // combine hits into array that suits template format
                enoughhits = getSymmetricRefHits(xpr, ypr, zpr, layerpr, SlimSegs, 4 - (TID_LEN / 2));
                if (!enoughhits) {
                    failed_count++;
                    continue;
                }
                used_entries++;

                //combining hits to super pixel hit array
                for (int i = 0; i < xpr.size(); i++) {
                    SPID = PE.getSuperPixel(xpr[i], ypr[i], zpr[i]);
                    SPIDs.push_back(SPID);
                    if (PRINTS) printf(" -- SPIDS =%#X \n", SPIDs[i]);
                }

                //store preprocessed tracks of cosmics
                cosmic_spid_tracks.push_back(SPIDs);

            }
            std::cout << std::endl;
        }
        tinF.Close();
        std::cout << "(STATUS) : Got cosmic test data. cosmic count: " << used_entries << " processed entries: "
                  << processed_entries << std::endl;

        tcosTIDF.cd();

        //add some meta data for the
        int wbins = spWbins;
        int zbins = spZbins;
        TTree tT_met("CosmicSIDMeta", "Metadata associated with these plots (PE config and dataset)");
        tT_met.Branch("cosmic_dataset", &cosmic_testing_dataset, "cosmic_dataset/I");
        tT_met.Branch("cosmic_testing_processed_entries", &processed_entries, "cosmic_processed_entries/I");
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
        tT_met.Branch("sp_count", &centralTPcount, "spcount/i");
        tT_met.Branch("sp_target_ratio", &spWZratio, "sp_target_ratio/F");
        tT_met.Fill();
        tT_met.Write();

        TTree tT_sids("CosmicSIDtracks", "SIDs of cosmic muon tracks (one cosmic track per entry)");
        std::vector<unsigned int> SPIDs;
        tT_sids.Branch("cosmic_track_sids", &SPIDs);

        for (int track = 0; track < cosmic_spid_tracks.size(); track++) {
            SPIDs = cosmic_spid_tracks[track];
            tT_sids.Fill();
        }

        tT_sids.Write();
        tcosTIDF.Close();

        return cosmic_spid_tracks;

    } else {
        //if there is a precomputed TID file, open it an read and return the TID data.

        std::cout << "(INFO)   : Opened file " << cosmic_sid_file << std::endl;

        TTree *tT_met;
        TTree *tT_sids;

        tcosTIDF.GetObject("CosmicSIDMeta", tT_met);
        tcosTIDF.GetObject("CosmicSIDtracks", tT_sids);

        int tf_cosmic_dataset;
        int tf_cosmic_testing_processed_entries;
        int tf_mode;
        int wBins;
        int zBins;
        unsigned int tf_sp_count;
        float tf_sp_target_ratio;

        std::vector<unsigned int> *SPIDs = nullptr;
        cosmic_spid_tracks.clear();


        tT_met->SetBranchAddress("cosmic_dataset", &tf_cosmic_dataset);
        tT_met->SetBranchAddress("cosmic_testing_processed_entries", &tf_cosmic_testing_processed_entries);
        tT_met->SetBranchAddress("wBins0", &wBins);
        tT_met->SetBranchAddress("zBins0", &zBins);
        tT_met->SetBranchAddress("sp_count", &tf_sp_count);
        tT_met->SetBranchAddress("sp_target_ratio", &tf_sp_target_ratio);
        tT_met->GetEntry(0);

        tT_sids->SetBranchAddress("cosmic_track_sids", &SPIDs);
        int entries = tT_sids->GetEntries();

        std::cout << "(INFO)   : Got data from configuration: wbins " << wBins << " | zbins " << zBins << " | mode "
                  << mode << std::endl;
        std::cout << "(INFO)   : Processed entries " << tf_cosmic_testing_processed_entries << " entries in file "
                  << entries << std::endl;

        for (int i = 0; i < entries; i++) {
            tT_sids->GetEntry(i);
            cosmic_spid_tracks.push_back(*SPIDs);
        }


        assert(cosmic_spid_tracks.size() == tf_cosmic_testing_processed_entries);

        return cosmic_spid_tracks;
    }

}
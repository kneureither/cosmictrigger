//
// Created by Konstantin Neureither on 06.08.20.
//

#include "cosmicTemplatesEfficiency.h"

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


void
cosmicTemplatesEfficiency(const int dataset, unsigned int centralTPcount, float spWZratio, float stopping_efficiency) {
    /*
     * Calculates the efficiency of a template database by adding some new cosmic data.
     * train_eff = #templ_matched / #templ_tested
     *
     * This is done for different template filters.
     *
     * input: template database configuration
     * output: training efficiencies fof different template filters (area, freq cut)
     */

    //result data
    std::vector<TIDLoadingFilter> filters = {ALL, CENTER_ONLY, RECURL_ONLY, MIXED_ONLY, NO_CENTER, CUT_ON_FREQ};
    std::vector<float> train_effs_total;
    std::vector<float> train_effs_relative;
    std::vector<std::vector<unsigned int>> cosmic_spid_tracks;

    int cosmic_testing_dataset = 30;
    const bool MAKE_PLOT = true;
    const int MAX_ENTRIES = 0;
    const bool PRINTS = false;
    const int mode = 0;

    //// INITIALIZE PATHS, TEMPLATE BANK

    const std::string pathtocosmicdata = "data/SlimmedData/";
    const std::string pathtotemplatedata = "data/TemplateData/";
    const std::string pathtooutput = "output/comsicTemplatesEfficiency/";
    const std::string pathtooutfile = pathtooutput + "dataset_" + get_padded_string(dataset, 3, '0') + "/";
    const std::string pathtoplots = pathtooutfile + "PDF/";
    const std::string pathtodatasettempldata = pathtotemplatedata + "dataset_" + get_padded_string(dataset, 3, '0') + "/";

    std::string infile = pathtocosmicdata + "mu3e_slimmed_segs_" + get_padded_string(cosmic_testing_dataset, 6, '0') + ".root";

    check_create_directory(pathtocosmicdata);
    check_create_directory(pathtotemplatedata);
    check_create_directory(pathtooutput);
    check_create_directory(pathtooutfile);
    check_create_directory(pathtoplots);
    check_create_directory(pathtodatasettempldata);

    //Get the Pattern Engine and Template Manager
    const int spWbins = (int) sqrt((float) spWZratio * (float) centralTPcount);
    const int spZbins = (int) sqrt((float) centralTPcount / (float) spWZratio);
//    std::cout << "\n -- PE config data:" << std::endl << "  wbins=" << spWbins << std::endl << "  zbins=" << spZbins << std::endl << std::endl;

    PatternEngine PE(spWbins, spZbins, pathtoplots);

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
    while(TObject* obj = iter()) {
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
        std::cout << "(INFO)   : sp count " << centralTPcount << " | wbins " << spWbins << " | zbins " << spZbins << std::endl;

        SlimSegsTreeRead SlimSegs = SlimSegsTreeRead(t_slimsegs);

        SlimSegs.getEntry(0);
        if (std::find(processed_run_ids.begin(), processed_run_ids.end(), SlimSegs.runID) != processed_run_ids.end()) {
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
            print_status_bar(entryno, max_entry, "building cosmic tids", "processed entries " + get_string(processed_entries));

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
        break;
    }
    tinF.Close();

    std::cout << "(STATUS) : Got cosmic test data. cosmic count: " << used_entries << " processed entries: " << processed_entries << std::endl;

    //open new TFile for plots
    TFile * tF = new TFile((pathtooutfile +"TemplateBankEfficiency_filter_" +
            getfileidtag(dataset, mode, spWbins, spZbins, stopping_efficiency)+ "_plots.root").c_str(), "recreate");
    if (!tF->IsOpen()) {
        std::cout << "[ERROR] File " << tF->GetName() << " is not open!" << std::endl;
    }

    float training_efficiency = 0;
    int template_count = 0;

    for(auto &filter : filters) {
        std::cout << "(STATUS) : Filter set to " << enum_to_string(filter) << std::endl;

        std::string mydirectory = "Plots_" + enum_to_string(filter);
        tF->mkdir(mydirectory.c_str());

        //init template bank
        TemplateBank TB(pathtoplots, dataset, mode, spWbins, spZbins);
        TB.readAMfromFile(pathtodatasettempldata, stopping_efficiency, filter);

        tF->cd(mydirectory.c_str());
        TB.SetPrints(PRINTS);
        TB.PlotTemplatePopulationHistogram();
        TB.PlotTemplateTypeDistribution();

        //get info to store in file later (maybe just do when filter==ALL)
        training_efficiency = TB.getEfficiency();
        template_count = TB.getInitialTemplateCount();


        //check the templates collected from file above
        for(auto &SPIDs : cosmic_spid_tracks) {
            TB.checkCosmicTemplate(&SPIDs[0], SPIDs.size(), filter);
        }

        //store data
        train_effs_total.push_back(TB.GetTrainEffTotal());
        train_effs_relative.push_back(TB.GetTrainEffRelative());
        std::cout << "(INFO)   : total eff " << TB.GetTrainEffTotal() << " | relative eff " << TB.GetTrainEffRelative() << std::endl;

        tF->cd();
    }

    tF->cd();

    //// Make train eff histogram

    int bins = filters.size() * 2;
    std::vector<float> binsX = {0};
    float binwidth=1/(float) bins;
    for(int i=0; i<bins*2; i++) {
        binsX.push_back(binsX[binsX.size()-1] + binwidth/10.0);
        binsX.push_back(binsX[binsX.size()-1] + 4*binwidth/10.0);
        binsX.push_back(binsX[binsX.size()-1] + 4*binwidth/10.0);
        binsX.push_back(binsX[binsX.size()-1] + binwidth/10.0);
    }
    auto *canvas = new TCanvas("cosmic efficiencies for different templ filters", "template type efficiencies", 900, 600);
    TH1F *h_templtypetotaleff = new TH1F("h_templtypetotaleff", "template type efficiencies", bins*2, &binsX[0]);
    TH1F *h_templtypereleff = new TH1F("h_templtypereleff", "template type efficiencies", bins*2, &binsX[0]);
    h_templtypetotaleff->SetStats(false);
    h_templtypereleff->SetStats(false);
    h_templtypetotaleff->GetXaxis()->SetNdivisions(bins, false);

    //fill with data
    int count=0;
    for(int i=2; i<bins*2; i+=4) {
        h_templtypetotaleff->SetBinContent(i, train_effs_total[count]);
        h_templtypereleff->SetBinContent(i+1, train_effs_relative[count]);
//        h_templtypetotaleff->GetXaxis()->SetBinLabel(i+1, enum_to_string(static_cast<TIDLoadingFilter>(i)).c_str());
//        h_templtypetotaleff->GetXaxis()->ChangeLabel(1+2*count,-1, -1, -1, -1,-1, enum_to_string(filters[count]));
        count++;
        h_templtypetotaleff->GetXaxis()->ChangeLabel(count*2,-1, -1, -1, -1,-1, enum_to_string(filters[count-1]));
        h_templtypetotaleff->GetXaxis()->ChangeLabel(count*2-1,-1, -1, -1, -1,-1, " ");
    }
    h_templtypetotaleff->GetXaxis()->ChangeLabel(count*2+1,-1, -1, -1, -1,-1, " ");
    h_templtypetotaleff->GetXaxis()->SetTickSize(0);


    //make it nice
    h_templtypetotaleff->GetYaxis()->SetTitle("cosmic efficiency");
    h_templtypetotaleff->GetXaxis()->SetTitle("template filter");
    h_templtypetotaleff->SetFillColor(kGreen);
    h_templtypetotaleff->SetLineColor(kGreen);
    h_templtypereleff->SetFillColor(kBlue);
//    h_templtypereleff->SetLineColor(kBlue);
    h_templtypetotaleff->SetFillStyle(3003);
    h_templtypereleff->SetFillStyle(3003);

    setPlottingStyle(h_templtypetotaleff); //general plotting style defined in plots.h

    h_templtypetotaleff->SetMaximum(1);
    h_templtypereleff->SetMaximum(1);
    h_templtypetotaleff->Draw();
    h_templtypereleff->Draw("SAME");

    //legend
    TLegend *legend = new TLegend(0.5, 0.8, 0.9, 0.9);
    legend->AddEntry(h_templtypetotaleff, "#epsilon_{training, total} for all cosmics");
    legend->AddEntry(h_templtypereleff, "#epsilon_{training, relative} for filter conditions");
    legend->SetTextSize(0.03);
    legend->Draw();

    TLatex tline1(.15,.81,("#it{#bf{TEMPLATE BANK} @ " + get_string(training_efficiency*100) + "% TRAINING EFF" +
                          " | TEMPL CNT " + get_string(template_count) + "}").c_str());
    tline1.SetTextFont(43);
    tline1.SetTextSize(10);
    tline1.SetNDC(kTRUE);
    tline1.Draw();

    TLatex tline2(.15,.78,("#it{#bf{CONFIG} WBINS " + get_string(spWbins) + " | ZBINS " + get_string(spZbins) +
                          " | DATASET " + get_string(dataset) + "}").c_str());
    tline2.SetTextFont(43);
    tline2.SetTextSize(10);
    tline2.SetNDC(kTRUE);
    tline2.Draw();

    TLatex tline3(.15,.75,("#it{COSMIC TESTING DATASET " + get_string(cosmic_testing_dataset) + " | COSMIC CNT " + get_string(processed_entries) + "}").c_str());
    tline3.SetTextFont(43);
    tline3.SetTextSize(10);
    tline3.SetNDC(kTRUE);
    tline3.Draw();
//
    canvas->Write();
    h_templtypetotaleff->Write();
    h_templtypereleff->Write();

    saveCanvas(canvas, ("TemplateTypeEffs_" + getfileidtag(dataset, mode, spWbins, spZbins, stopping_efficiency)).c_str(), pathtoplots);

    //add some meta data for the
    int datast = dataset;
    TTree tT_met("MetadataTree","Metadata associated with these plots (PE, TB config and dataset)");
    tT_met.Branch("dataset", &datast, "dataset/I");
    tT_met.Branch("cosmic_testing_dataset", &cosmic_testing_dataset, "cosmic_testing_dataset/I");
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
    tT_met.Branch("training_efficiency", &training_efficiency, "training_efficiency/F");
    tT_met.Branch("stopping_efficiency", &stopping_efficiency, "stopping_efficiency/F");
    tT_met.Branch("template_count", &template_count, "template_count/I");
    tT_met.Branch("sp_count", &centralTPcount, "spcount/i");
    tT_met.Branch("sp_target_ratio", &spWZratio, "sp_target_ratio/F");
    tT_met.Fill();
    tT_met.Write();

    tF->Close();


}


std::vector<std::vector<unsigned int>> produceCosmicSIDtracks(int cosmic_testing_dataset, int centralTPcount, float spWZratio, int mode) {
    /**
     * reads a cosmic file, computes the TIDs and saves the result to a file, which can later be used for
     * efficiency evalutation. If the file already exists, then it just reads the file and returns the vector.
     */

    const bool PRINTS = false;

    const std::string path_to_cosmic_inputdata = "data/SlimmedData/";
    const std::string path_to_cosmic_tids = "data/CosmicSIDtrackData/";
    const std::string path_to_plots = "output/produceCosmicSIDtracks/";
    const std::string path_to_plots_dataset = path_to_plots + "dataset_" + get_padded_string(cosmic_testing_dataset, 3, '0') + "/";
    const std::string path_to_cosmic_tids_dataset = path_to_cosmic_tids + "dataset_" + get_padded_string(cosmic_testing_dataset, 3, '0') + "/";

    check_create_directory(path_to_cosmic_inputdata);
    check_create_directory(path_to_cosmic_tids);
    check_create_directory(path_to_plots);
    check_create_directory(path_to_plots_dataset);
    check_create_directory(path_to_cosmic_tids_dataset);

    //Get the Pattern Engine and Template Manager
    const int spWbins = (int) sqrt((float) spWZratio * (float) centralTPcount);
    const int spZbins = (int) sqrt((float) centralTPcount / (float) spWZratio);
    PatternEngine PE(spWbins, spZbins, path_to_plots_dataset);

    std::string infile = path_to_cosmic_inputdata + "mu3e_slimmed_segs_" + get_padded_string(cosmic_testing_dataset, 6, '0') + ".root";
    std::string cosmic_sid_file = path_to_cosmic_tids_dataset + "CosmicSIDtracks_" + getfileidtag(cosmic_testing_dataset, mode, spWbins, spZbins)  + ".root";
    std::vector<std::vector<unsigned int>> cosmic_spid_tracks;

    //First check if a file exists already
    TFile tcosTIDF(cosmic_sid_file.c_str());
    if (!tcosTIDF.IsOpen()) {
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
        while(TObject* obj = iter()) {
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
            std::cout << "(INFO)   : sp count " << centralTPcount << " | wbins " << spWbins << " | zbins " << spZbins << std::endl;

            SlimSegsTreeRead SlimSegs = SlimSegsTreeRead(t_slimsegs);

            SlimSegs.getEntry(0);
            if (std::find(processed_run_ids.begin(), processed_run_ids.end(), SlimSegs.runID) != processed_run_ids.end()) {
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
                print_status_bar(entryno, max_entry, "building cosmic sid tracks", "processed entries " + get_string(processed_entries));

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
        std::cout << "(STATUS) : Got cosmic test data. cosmic count: " << used_entries << " processed entries: " << processed_entries << std::endl;

        tcosTIDF.cd();

        //add some meta data for the
        int wbins = spWbins;
        int zbins = spZbins;
        TTree tT_met("CosmicSIDMeta","Metadata associated with these plots (PE config and dataset)");
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

        TTree tT_sids("CosmicSIDtracks","SIDs of cosmic muon tracks (one cosmic track per entry)");
        std::vector<unsigned int> SPIDs;
        tT_sids.Branch("cosmic_track_sids", &SPIDs);

        for(int track=0; track < cosmic_spid_tracks.size(); track++) {
            SPIDs = cosmic_spid_tracks[track];
            tT_sids.Fill();
        }

        tT_sids.Write();
        tcosTIDF.Close();

        return cosmic_spid_tracks;

    } else {

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

        tT_sids->SetBranchAddress("cosmic_track_sids",&SPIDs);
        int entries = tT_sids->GetEntries();

        std::cout << "(INFO)   : Got data from configuration: wbins " << wBins << " | zbins " << zBins << " | mode " << mode << std::endl;
        std::cout << "(INFO)   : Processed entries " << tf_cosmic_testing_processed_entries << " entries in file " <<  entries << std::endl;

        for(int i=0; i<entries; i++) {
            tT_sids->GetEntry(i);
            cosmic_spid_tracks.push_back(*SPIDs);
        }


        assert(cosmic_spid_tracks.size() == tf_cosmic_testing_processed_entries);

        return cosmic_spid_tracks;
    }

}
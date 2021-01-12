//
// Created by Konstantin Neureither on 06.08.20.
//

#include "cosmicEffTemplFilter.h"

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
#include "getCosmicSIDtracks.h"


void cosmicTemplatesEfficiency(const int dataset, unsigned int centralTPcount, float spWZratio, float stopping_efficiency) {
    /**
     * Calculates the efficiency of a template database by adding some new cosmic data.
     * train_eff = #templ_matched / #templ_tested
     *
     * This is done for different template filters and a histogram showing template filter acceptance and efficiency
     * is produced and stored at the path specified in variable std::string pathtoplotsfilter.
     *
     * input: template database configuration
     * output: plot of training efficiencies of different template filters (area, freq cut)
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

    gStyle->SetLegendBorderSize(0);

    //// INITIALIZE PATHS, TEMPLATE BANK

    const std::string pathtocosmicdata = "data/SlimmedData/";
    const std::string pathtotemplatedata = "data/TemplateData/";
    const std::string pathtooutput = "output/4_cosmicEffCheck/";
    const std::string pathtooutfile = pathtooutput + "dataset_" + get_padded_string(dataset, 3, '0') + "/";
    const std::string pathtoplots = pathtooutfile + "PDF/";
    const std::string pathtoplotsfilter = pathtoplots /*+ "FilterEffs/"*/;
    const std::string pathtodatasettempldata = pathtotemplatedata + "dataset_" + get_padded_string(dataset, 3, '0') + "/";

    std::string infile = pathtocosmicdata + "mu3e_slimmed_segs_" + get_padded_string(cosmic_testing_dataset, 6, '0') + ".root";

    check_create_directory(pathtocosmicdata);
    check_create_directory(pathtotemplatedata);
    check_create_directory(pathtooutput);
    check_create_directory(pathtooutfile);
    check_create_directory(pathtoplots);
    check_create_directory(pathtodatasettempldata);
    check_create_directory(pathtoplotsfilter);

    //Get the Pattern Engine and Template Manager
    const int spWbins = (int) sqrt((float) spWZratio * (float) centralTPcount);
    const int spZbins = (int) sqrt((float) centralTPcount / (float) spWZratio);
//    std::cout << "\n -- PE config data:" << std::endl << "  wbins=" << spWbins << std::endl << "  zbins=" << spZbins << std::endl << std::endl;

    PatternEngine PE(spWbins, spZbins, pathtoplots);

    //stats data
    int processed_entries = 0;
    int used_entries = 0;

    cosmic_spid_tracks= getCosmicSIDtracks(cosmic_testing_dataset, centralTPcount, spWZratio, mode);
    processed_entries = cosmic_spid_tracks.size();


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

        std::string label = enum_to_string(filters[count-1]);
        std::replace( label.begin(), label.end(), '_', ' ');
        h_templtypetotaleff->GetXaxis()->ChangeLabel(count*2,-1, -1, -1, -1,-1, label);
        h_templtypetotaleff->GetXaxis()->ChangeLabel(count*2-1,-1, -1, -1, -1,-1, " ");
    }
    h_templtypetotaleff->GetXaxis()->ChangeLabel(count*2+1,-1, -1, -1, -1,-1, " ");
    h_templtypetotaleff->GetXaxis()->SetTickSize(0);


    //make it nice
    h_templtypetotaleff->GetYaxis()->SetTitle("cosmic efficiency");
    h_templtypetotaleff->GetXaxis()->SetTitle("template filter");
    h_templtypetotaleff->SetFillColor(kRed-4);
    h_templtypetotaleff->SetLineColor(kRed-4);
    h_templtypereleff->SetLineColor(kCyan+4);
    h_templtypereleff->SetFillColor(kCyan+4);
    h_templtypetotaleff->SetFillStyle(3005);
    h_templtypereleff->SetFillStyle(3004);

    setPlottingStyle(h_templtypetotaleff); //general plotting style defined in plots.h

    h_templtypetotaleff->SetMaximum(1);
    h_templtypereleff->SetMaximum(1);
    h_templtypetotaleff->Draw();
    h_templtypereleff->Draw("SAME");

    float leg_x = 0.14;
    float leg_y = 0.85;
    float spacing = 0.04;

    //legend
    TLegend *legend = new TLegend(leg_x + 0.39, 0.8, 0.89, leg_y+0.03);
    legend->AddEntry(h_templtypetotaleff, "cosmic detection efficiency #epsilon^{cosmic}_{detection}");
    legend->AddEntry(h_templtypereleff, "acceptance of filter #Alpha^{cosmic}_{filter}");
    legend->SetTextSize(0.03);
    legend->Draw();

    TLatex tline1(leg_x,leg_y,("#it{#bf{TEMPLATE BANK} @ " + get_string(training_efficiency*100).substr(0,2) + "% TRAINING EFF" +
                          " | T-CNT " + get_string(template_count) + "}").c_str());
    tline1.SetTextFont(43);
    tline1.SetTextSize(14);
    tline1.SetNDC(kTRUE);
    tline1.Draw();

    TLatex tline2(leg_x, leg_y-spacing,("#it{#bf{CONFIG} BINS " + get_string(spWbins) + "#times" + get_string(spZbins) +
                          " | DATASET " + get_string(dataset) + "}").c_str());
    tline2.SetTextFont(43);
    tline2.SetTextSize(14);
    tline2.SetNDC(kTRUE);
    tline2.Draw();

//    TLatex tline3(.15,.75,("#it{COSMIC TESTING DATASET " + get_string(cosmic_testing_dataset) + " | COSMIC CNT " + get_string(processed_entries) + "}").c_str());
//    tline3.SetTextFont(43);
//    tline3.SetTextSize(14);
//    tline3.SetNDC(kTRUE);
//    tline3.Draw();
//
    canvas->Write();
    h_templtypetotaleff->Write();
    h_templtypereleff->Write();

    saveCanvas(canvas, ("FilterTypeEffs_" + getfileidtag(dataset, mode, spWbins, spZbins, stopping_efficiency)).c_str(), pathtoplotsfilter);

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
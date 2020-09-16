//
// Created by Konstantin Neureither on 27.08.20.
//

#include <string>
#include <cassert>
#include <iostream>
#include <vector>

//root
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "THStack.h"

//custom
#include "../inc/cosmicTemplatesBgPlots.h"
#include "utilityFunctions.h"
#include "plots.h"

void makeBgEvalPlots(const int dataset, const int bgrun, std::string filename) {
    const std::string pathtodata = "plots/Mu3eCosPatBgEval/";
    const std::string pathtoplots = "plots/Mu3eCosPatBgEvalPlots/";
    std::string pathtorunplots = pathtoplots + "dataset_" + get_padded_string(dataset, 3, '0')+ "_run_" + get_padded_string(bgrun, 3, '0') + "/";
    std::string pathtorundata = pathtodata + "bgrun_" + get_padded_string(bgrun, 3, '0') + "/";

    const bool MAKE_PLOT = true;
    const bool PRINTS = false;
    const int delete_cycle=0;
    std::vector<int> cycles = {1,3};

    gStyle->SetTitleFontSize(0.06);

    check_create_directory(pathtodata);
    check_create_directory(pathtoplots);
    check_create_directory(pathtorunplots);
    check_create_directory(pathtorundata);

    //open file with plots from Pattern Engine and Template Bank
    TFile tinF((pathtorundata + filename).c_str());
    if (!tinF.IsOpen()) {
        std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
        exit(0);
    }

    //root file meta data
    TTree *tree_meta;
    TH1F *h_bgeff;
    std::vector<TH1F*> h_bgeffs;
    std::string ltext;

    THStack *hs = new THStack("hs","");

    int mode;
    int bg_events;
    int bg_run;
    int max_muon_hits;
    int max_frame_nhits;
    int processed_frames;
    int data_ds;
    int wBins[3];
    int zBins[3];
    int SPcount;
    float spWZratio;

    int colpalette[10] = {632,600,427,420,410,414,601,603,861,854};

    gStyle->SetPalette(kBlueGreenYellow);

    TCanvas *canvas = new TCanvas("canvas", "Background match efficiency result", 900, 600);
//    canvas->SetLeftMargin(0.15);
//    canvas->SetRightMargin(0.15);
    canvas->SetGrid(1, 1);
    canvas->SetTicks(1, 1);
    canvas->SetLogy(1);

    auto *pad0 = new TPad("background match efficiency", "background match efficiency", 0, 0, 1,1);
//    pad0->SetLogy(1);
    pad0->Draw();
    pad0->cd();


    //get different cycles of trees in file
    int cycle;


    for(const auto cycle : cycles) {

        std::cout << "STATUS : Processing cycle " << cycle << std::endl;

        std::string treename = "MetadataTree;" + get_string(cycle);
        tinF.GetObject(treename.c_str(), tree_meta);
        tree_meta->SetBranchAddress("bg_events", &bg_events);
        tree_meta->SetBranchAddress("bg_run", &bg_run);
        tree_meta->SetBranchAddress("max_muon_hits", &max_muon_hits);
        tree_meta->SetBranchAddress("max_frame_nhits", &max_frame_nhits);
        tree_meta->SetBranchAddress("processed_frames", &processed_frames);
        tree_meta->SetBranchAddress("mode", &mode);
        tree_meta->SetBranchAddress("wBins0", &wBins[0]);
        tree_meta->SetBranchAddress("zBins0", &zBins[0]);
        tree_meta->GetEntry(0);

        SPcount = wBins[0] * zBins[0];
        spWZratio = (float) wBins[0] / (float) zBins[0];

        std::cout << "  -- SPcount=" << SPcount << " SPratio=" << spWZratio << " wBins=" << wBins[0] << " zBins=" << zBins[0] << std::endl;
        std::cout << "     max muon hits=" << max_muon_hits << std::endl;
        std::cout << "     max frame nhits=" << max_frame_nhits << std::endl;

        tinF.GetObject(("h_bgeff;" + get_string(cycle)).c_str(), h_bgeff);

        TH1F *h_bgeffclone = (TH1F*) h_bgeff->Clone();
        h_bgeffs.push_back(h_bgeffclone);

        hs->Add(h_bgeffclone);

        //make the plots
        h_bgeffclone->SetMarkerStyle(5);

        float stddev = h_bgeffclone->GetStdDev();
        float mean = h_bgeffclone->GetMean();

        if(max_muon_hits != 0) {
            ltext = "#splitline{#bf{FRAMES WITH MUON} (max hits included: " + get_string(max_muon_hits) + ")}"+
                    "{#it{mean:}  " + get_string(mean) + "  #it{std dev:}  " +get_string(stddev) + "}";
        } else {
//            ltext = "#splitline{#bf{FRAMES WITHOUT MOUN}}"
//                    "{#it{mean:}  " + get_string(mean) + "  #it{std dev:}  " + get_string(stddev) + "}";
            ltext = "#splitline{#bf{FALSE POSITIVE} #it{max nhit:} " + get_string(max_frame_nhits) + " #it{# frames: }" + get_string(processed_frames) +"}"
                    "{#it{mean:}  " + get_string(mean) + "  #it{std dev:}  " + get_string(stddev) + "}";
        }
        h_bgeffclone->SetTitle(ltext.c_str());
    }




    std::cout << hs->GetNhists() << std::endl;

    std::string lline1 = "Processed frames of background data: " + get_string(bg_events);
    std::string lline2 = "#it{run:}    #bf{" + get_string(bg_run) + "}   #it{sp aspect ratio:}  #bf{" + (spWZratio < 1 ? "1:" + get_string(1/spWZratio) : get_string(spWZratio) + ":1") + "}";
    std::string lline3 = "#it{bins_{z}:}   #bf{" + get_string(zBins[0]) + "}   #it{bins_{w}:} #bf{" + get_string(wBins[0]) + "}  #it{dataset:} #bf{" + get_string(dataset) + "}";

    auto legend = new TLegend(0.55,0.8 - 0.05*hs->GetNhists(),0.9,0.9);
    legend->SetTextFont(43);
    legend->SetTextSize(15);
    legend->SetHeader(lline1.c_str(), "c"); //"C" centers header)

    legend->SetTextSize(10);
    legend->AddEntry((TObject*)0, ("#splitline{" + lline2 + "}{" + lline3 + "}").c_str(), "");

    for(int cycle=0; cycle < hs->GetNhists(); cycle++) {
        std::cout << hs->GetNhists() << std::endl;
        TH1F* h_current = (TH1F*) (hs->GetHists()->At(cycle));
        legend->AddEntry(h_current,h_current->GetTitle(),"p");
    }

    hs->Draw("nostack PLC PMC");

    hs->GetXaxis()->SetTitle("training_efficiency #epsilon = #frac{#tau_{matched}}{#tau_{combinatorics}}");
    hs->GetXaxis()->SetLabelFont(43);
    hs->GetXaxis()->SetLabelSize(14);
    hs->GetXaxis()->SetTitleFont(53);
    hs->GetXaxis()->SetTitleSize(14);
    hs->GetXaxis()->SetTitleOffset(1.4);
    hs->GetXaxis()->CenterTitle (false);

    hs->GetYaxis()->SetTitle("normalized distribution");
    hs->GetYaxis()->SetLabelFont(43);
    hs->GetYaxis()->SetLabelSize(14);
    hs->GetYaxis()->SetTitleFont(53);
    hs->GetYaxis()->SetTitleSize(14);
    hs->GetYaxis()->SetTitleOffset(1.6);
    hs->GetYaxis()->CenterTitle(false);

    canvas->Modified();

    legend->SetBorderSize(1);
    legend->Draw();

//    canvas->SaveAs((pathtorunplots + "CosPatBGeff.pdf").c_str());
    saveCanvas(canvas, "CosPatPlots_dataset_" + get_string(dataset) + "_bgrun_" + get_padded_string(bgrun, 6, '0') + "_bgevents_" + get_string(bg_events), pathtorunplots);

    tinF.Close();

}
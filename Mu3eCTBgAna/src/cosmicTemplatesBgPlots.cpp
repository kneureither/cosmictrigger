//
// Created by Konstantin Neureither on 27.08.20.
//

#include <string>
#include <cassert>
#include <iostream>

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
    const int MAX_ENTRIES = 0;
    const bool PRINTS = false;
    const int delete_cycle=0;

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
    TH1F h_disceff[4];

    THStack *hs = new THStack("hs","");

    int mode;
    int bg_events;
    int bg_run;
    int max_muon_hits;
    int data_ds;
    int wBins[3];
    int zBins[3];
    int SPcount;
    float spWZratio;

    int colpalette[10] = {632,600,427,420,410,414,601,603,861,854};

    TCanvas *canvas = new TCanvas("canvas", "Template Bank Pattern Result", 900, 600);
//    canvas->SetLeftMargin(0.15);
//    canvas->SetRightMargin(0.15);
    canvas->SetGrid(1, 1);
    canvas->SetTicks(1, 1);
    canvas->SetLogx(1);

    auto *pad0 = new TPad("background match efficiency", "background match efficiency", 0, 0, 1,1);
    pad0->Draw();

    auto legend = new TLegend(0.3,0.6,0.7,0.9);

    //get different cycles of trees in file
    int treecount = 0;
    int cycle = 0;
    std::string tree;

//    TList *list = tinF.GetListOfKeys();
//    TIter iter(list->MakeIterator());
//    //iterate over all cycles of trees separately = different SP configs
//    while (TObject *obj = iter()) {
//        TKey *theKey = (TKey *) obj;
//        tree = theKey->GetName();
//        cycle = theKey->GetCycle();
//
//        if(delete_cycle != 0) {
//            if(cycle == delete_cycle) {
//                std::string object_to_remove = tree + ";" + get_string(delete_cycle);
//                std::cout << " now deleting object " << tree << " cycle " << cycle << " complete name " << object_to_remove << std::endl;
//                gDirectory->Delete(object_to_remove.c_str());
//                continue;
//            }
//        }
//
//        std::cout << tree << std::endl;
//
//        if(tree != "MetadataTree") continue;
//        treecount++;
//        std::cout << treecount << std::endl;

    for(cycle=1; cycle <= 3; cycle++) {

        std::cout << "STATUS : Processing tree " << tree << " cycle " << cycle << std::endl;

        std::string treename = "MetadataTree;" + get_string(cycle);
        tinF.GetObject(treename.c_str(), tree_meta);
        tree_meta->SetBranchAddress("bg_events", &bg_events);
        tree_meta->SetBranchAddress("bg_run", &bg_run);
        tree_meta->SetBranchAddress("max_muon_hits", &max_muon_hits);
        tree_meta->SetBranchAddress("mode", &mode);
        tree_meta->SetBranchAddress("wBins0", &wBins[0]);
        tree_meta->SetBranchAddress("zBins0", &zBins[0]);
        tree_meta->GetEntry(0);

        SPcount = wBins[0] * zBins[0];
        spWZratio = (float) wBins[0] / (float) zBins[0];

        std::cout << "  -- SPcount=" << SPcount << " SPratio=" << spWZratio << " wBins=" << wBins[0] << " zBins=" << zBins[0] << std::endl;
        std::cout << "     max muon hits=" << max_muon_hits << std::endl;

        tinF.GetObject(("h_bgeff;" + get_string(cycle)).c_str(), h_bgeff);

        h_disceff[cycle] = *h_bgeff;


        //make the plots
        pad0->cd();
        h_disceff[cycle].SetLineColor(colpalette[treecount-1]);
        h_disceff[cycle].SetMarkerStyle(5);
        h_disceff[cycle].SetMarkerColor(colpalette[treecount-1]);

        h_disceff[cycle].GetXaxis()->SetTitle("efficiency #epsilon_{templates matched/generated}");
        h_disceff[cycle].GetXaxis()->SetLabelFont(43);
        h_disceff[cycle].GetXaxis()->SetLabelSize(14);
        h_disceff[cycle].GetXaxis()->SetTitleFont(63);
        h_disceff[cycle].GetXaxis()->SetTitleSize(14);
        h_disceff[cycle].GetXaxis()->SetTitleOffset(1.4);
        h_disceff[cycle].GetXaxis()->CenterTitle(false);

        h_disceff[cycle].GetYaxis()->SetTitle("normalized distribution");
        h_disceff[cycle].GetYaxis()->SetLabelFont(43);
        h_disceff[cycle].GetYaxis()->SetLabelSize(14);
        h_disceff[cycle].GetYaxis()->SetTitleFont(63);
        h_disceff[cycle].GetYaxis()->SetTitleSize(11);
        h_disceff[cycle].GetYaxis()->SetTitleOffset(1.6);
        h_disceff[cycle].GetYaxis()->CenterTitle(false);


        hs->Add(&h_disceff[cycle]);

        std::string ltext="muon hits=" + get_string(max_muon_hits);

        legend->AddEntry(&h_disceff[cycle],ltext.c_str(),"l");
    }

    std::string lline1 = "Analysed frames: " + get_string(bg_events);
    std::string lline2 = "run: " + get_string(bg_run) + " SP ratio:" + (spWZratio < 1 ? "1:" + get_string(1/spWZratio) : get_string(spWZratio) + ":1");
    std::string lline3 = "bins_{z}=" + get_string(zBins[0]) + " bins_{w}=" + get_string(wBins[0]);

//    legend->SetHeader(("#splitline{" + lline1 + "}{" + lline2 + "}{" + lline3 + "}").c_str(), "C"); //"C" centers header)
    legend->SetHeader(lline1.c_str(), "C"); //"C" centers header)

    hs->Draw();

    legend->SetTextFont(43);
    legend->SetTextSize(10);
    legend->Draw();

    std::cout << "check" << std::endl;

    canvas->SaveAs((pathtorunplots + "CosPatBGeff.pdf").c_str());
//    saveCanvas(canvas, "CosPatPlots_dataset_" + get_string(dataset) + "_bgrun_" + get_padded_string(bgrun, 6, '0'), pathtorunplots);

    tinF.Close();

}
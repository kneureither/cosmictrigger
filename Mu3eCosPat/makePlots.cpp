//
// Created by Konstantin Neureither on 22.07.20.
//
//basic stuff
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

//custom
#include "include/makePlots.h"
#include "utilityFunctions.h"
#include "plots.h"

void makeCosPatPlots(const int dataset, const int combination_id) {
    const std::string pathtodata = "plots/Mu3eCosPat/";
    const std::string pathtoplots = "plots/Mu3eCosPatPlots/";
    std::string pathtorunplots = pathtoplots + "dataset_" + get_padded_string(dataset, 3, '0') + "/";
    std::string pathtorundata = pathtodata + "dataset_" + get_padded_string(dataset, 3, '0') + "/";

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
    TFile tinF((pathtorundata + "TemplateBank_dataset_" + get_padded_string(dataset, 3, '0') + "_id" + get_padded_string(combination_id, 3, '0') + "_plots.root").c_str(), "update");
    if (!tinF.IsOpen()) {
        std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
        exit(0);
    }

    //root file meta data
    TTree *tree_meta;

    float efficiency;
    int templ_count;
    int mode;
    int processed_events;
    int data_ds;
    int wBins[3];
    int zBins[3];
    int SPcount;
    float spWZratio;

    TH1F *h_templfreq;
    TGraph *g_efficiency;
    TGraph *g_tnumber;

    int colpalette[10] = {433,435,427,420,410,414,601,603,861,854};

    TCanvas *canvas = new TCanvas("canvas", "Template Bank Pattern Result", 900, 600);
//    canvas->SetLeftMargin(0.15);
//    canvas->SetRightMargin(0.15);
    canvas->SetGrid(1, 1);
    canvas->SetTicks(1, 1);
    canvas->SetLogx(1);

//    auto *pad0 = new TPad("title", "title", 0, 0.3, 1, 0.99);
//    pad0->SetLogx(0);
//    pad0->Draw();

    auto *pad1 = new TPad("template efficiency", "template efficiency", 0, 0.3, 1, 0.99);
    pad1->SetLogx(0);
    pad1->Draw();

    auto legend = new TLegend(0.5,0.1,0.9,0.5);

    auto *pad2 = new TPad("template count", "template count", 0, 0, 1, 0.3);
    pad2->SetLogx(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGrid(1,5);
    pad2->Draw();


    //get different cycles of trees in file
    int treecount = 0;
    int cycle = 0;
    std::string tree;

    TList *list = tinF.GetListOfKeys();
    TIter iter(list->MakeIterator());
    //iterate over all cycles of trees separately = different SP configs
    while (TObject *obj = iter()) {
        treecount++;
        TKey *theKey = (TKey *) obj;
        tree = theKey->GetName();
        cycle = theKey->GetCycle();

        if(delete_cycle != 0) {
            if(cycle == delete_cycle) {
                std::string object_to_remove = tree + ";" + get_string(delete_cycle);
                std::cout << " now deleting object " << tree << " cycle " << cycle << " complete name " << object_to_remove << std::endl;
                gDirectory->Delete(object_to_remove.c_str());
                continue;
            }
        }

        if(tree != "MetadataTree") continue;

        std::cout << "STATUS : Processing tree " << tree << " cycle " << cycle << std::endl;

        std::string treename = "MetadataTree;" + get_string(cycle);
        tinF.GetObject(treename.c_str(), tree_meta);
        tree_meta->SetBranchAddress("dataset", &data_ds);
        tree_meta->SetBranchAddress("efficiency", &efficiency);
        tree_meta->SetBranchAddress("templ_count", &templ_count);
        tree_meta->SetBranchAddress("processed_events", &processed_events);
        tree_meta->SetBranchAddress("mode", &mode);
        tree_meta->SetBranchAddress("wBins0", &wBins[0]);
        tree_meta->SetBranchAddress("zBins0", &zBins[0]);
        tree_meta->GetEntry(0);

        SPcount = wBins[0] * zBins[0];
        spWZratio = (float) wBins[0] / (float) zBins[0];

        std::cout << "  -- SPcount=" << SPcount << " SPratio=" << spWZratio << " wBins=" << wBins[0] << " zBins=" << zBins[0] << std::endl;

        tinF.GetObject(("g_efficiency;" + get_string(cycle)).c_str(), g_efficiency);
        tinF.GetObject(("g_tnumber;" + get_string(cycle)).c_str(), g_tnumber);
        tinF.GetObject(("h_templfreq;" + get_string(cycle)).c_str(), h_templfreq);


        //make the plots
        pad1->cd();
        g_efficiency->SetLineColor(colpalette[treecount-1]);
        g_efficiency->SetMarkerStyle(5);
        g_efficiency->SetMarkerColor(colpalette[treecount-1]);

        g_efficiency->GetXaxis()->SetTitle("Number of training events");
        g_efficiency->GetXaxis()->SetLabelFont(43);
        g_efficiency->GetXaxis()->SetLabelSize(14);
        g_efficiency->GetXaxis()->SetTitleFont(63);
        g_efficiency->GetXaxis()->SetTitleSize(14);
        g_efficiency->GetXaxis()->SetTitleOffset(1.4);
        g_efficiency->GetXaxis()->CenterTitle(true);

        g_efficiency->GetYaxis()->SetTitle("efficiency #epsilon = #frac{tmpl_{matched}}{tmpl_{total}}");
        g_efficiency->GetYaxis()->SetLabelFont(43);
        g_efficiency->GetYaxis()->SetLabelSize(14);
        g_efficiency->GetYaxis()->SetTitleFont(63);
        g_efficiency->GetYaxis()->SetTitleSize(11);
        g_efficiency->GetYaxis()->SetTitleOffset(1.6);
        g_efficiency->GetYaxis()->CenterTitle(false);

//        TPaveText *pt = (TPaveText*)(g_efficiency->GetTitle()); pt->SetTextSize(0.1);


        g_efficiency->Draw((treecount == 1 ? "AL" : "SAME"));
        std::string ltext="bins_{w}=" + get_string(wBins[0]) +
                "   bins_{z}=" + get_string(zBins[0]) +
                "   #epsilon =" + get_string(efficiency * 100).substr(0, 4) + "%" +
                "   cnt=" + get_string(SPcount) +
                "   ratio W/Z=" + (spWZratio < 1 ? "1:" + get_string(1/spWZratio) : get_string(spWZratio) + ":1");

        legend->AddEntry(g_efficiency,ltext.c_str(),"l");

        pad2->cd();
        g_tnumber->SetLineColor(colpalette[treecount-1]);
        g_tnumber->GetYaxis()->SetRangeUser(0, 750000);

        g_tnumber->SetTitle("");

        g_tnumber->GetXaxis()->SetTitle("Number of training events");
        g_tnumber->GetXaxis()->SetLabelFont(43);
        g_tnumber->GetXaxis()->SetLabelSize(14);
        g_tnumber->GetXaxis()->SetTitleOffset(3.3);
        g_tnumber->GetXaxis()->SetTitleFont(63);
        g_tnumber->GetXaxis()->SetTitleSize(14);
        g_tnumber->GetXaxis()->CenterTitle(true);
        g_tnumber->GetXaxis()->SetTickLength(0.05);

        g_tnumber->GetYaxis()->SetTitle("number of templates");
        g_tnumber->GetYaxis()->SetLabelFont(43);
        g_tnumber->GetYaxis()->SetLabelSize(14);
        g_tnumber->GetYaxis()->SetTitleFont(63);
        g_tnumber->GetYaxis()->SetTitleSize(11);
        g_tnumber->GetYaxis()->SetTitleOffset(1.6);
        g_tnumber->GetYaxis()->SetNdivisions(3, 5, 0, false);
        g_tnumber->GetYaxis()->CenterTitle(false);


        g_tnumber->Draw((treecount == 1 ? "AL" : "SAME"));
    }

    pad1->cd();
//    legend->SetHeader(("Super Pixel count (central area): " + get_string(SPcount)).c_str(),"C"); // option "C" allows to center the header
    legend->SetTextFont(43);
    legend->SetTextSize(10);
    legend->Draw();
    saveCanvas(canvas, "CosPatPlots_dataset_" +get_string(dataset) + "_id" + get_padded_string(combination_id, 3, '0'), pathtorunplots);

    tinF.Close();

}



//read the data set analysis root file and combine eff plots with SP config data as TLatex, make one single nice pdf.
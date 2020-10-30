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
#include "TMultiGraph.h"
#include "TF1.h"
#include "TMath.h"

//custom
#include "../inc/cosmicTemplatesTrainingPlots.h"
#include "utilityFunctions.h"
#include "plots.h"

#define LEGEND_ON_TOP false

void makeCosPatPlots(const int dataset, const int combination_id, std::vector<int> cycle_plotting_order,
                     std::string filelabel="") {
    const std::string pathtodata = "output/Mu3eCosPat/";
    const std::string pathtoplots = "output/Mu3eCosPatPlots/";
    std::string pathtorunplots = pathtoplots + "dataset_" + get_padded_string(dataset, 3, '0') + "/";
    std::string pathtorundata = pathtodata + "dataset_" + get_padded_string(dataset, 3, '0') + "/";

    const bool MAKE_PLOT = true;
    const int MAX_ENTRIES = 0;
    const bool PRINTS = false;
    const int delete_cycle=0;
    const bool FILTER = false;
//    const bool LEGEND_ON_TOP = false;

    const bool SET_LOGX = 0;
    const bool FIT = true;
    const int LINE_W = 2;

    gStyle->SetLegendBorderSize(0);


    gStyle->SetLegendBorderSize(0);
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetPalette(kThermometer);
    gStyle->SetLegendTextSize(0.027);
    gStyle->SetLineWidth(1);
    gStyle->SetHistLineWidth(3);
    gStyle->SetMarkerStyle(5);

    check_create_directory(pathtodata);
    check_create_directory(pathtoplots);
    check_create_directory(pathtorunplots);
    check_create_directory(pathtorundata);

    //open file with plots from Pattern Engine and Template Bank
    TFile tinF((pathtorundata + get_tbtrainingoutfile(dataset, combination_id)).c_str(), "update");
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
    TGraph *g_eff_tnumber;
    auto g_efficiencies  = new TMultiGraph();
    auto g_tnumbers  = new TMultiGraph();
    auto g_eff_tnumbers  = new TMultiGraph();



    TCanvas *canvas = new TCanvas("canvas", "Template Bank Pattern Result", 900, 750);
    canvas->SetGrid(0, 0);
    canvas->SetTicks(1, 1);

    auto *pad1 = new TPad("template efficiency", "template efficiency", 0, 0.3, 1, 0.99);
    pad1->SetLogx(SET_LOGX);
//    pad1->SetGrid(1,5);
    pad1->Draw();

#if LEGEND_ON_TOP
        auto legend = new TLegend(0.6,0.55,0.9,0.9);
#else
        auto legend = new TLegend(0.5,0.15,0.88,0.25+cycle_plotting_order.size()*0.05);
        legend->SetTextAlign(13);
#endif

    auto *pad2 = new TPad("template count", "template count", 0, 0, 1, 0.3);
    pad2->SetLogx(SET_LOGX);
    pad2->SetBottomMargin(0.2);
//    pad2->SetGrid(1,5);
    pad2->Draw();


    //get different cycles of trees in file
    int treecount = 0;
    std::string tree;

    if(cycle_plotting_order.size() == 0) {
        std::cout << "STATUS : Processing complete file chronologically " << std::endl;

        TList *list = tinF.GetListOfKeys();
        TIter iter(list->MakeIterator());
        //iterate over all cycles of trees separately = different SP configs
        while (TObject *obj = iter()) {
            treecount++;
            TKey *theKey = (TKey *) obj;
            tree = theKey->GetName();
            int cycle = theKey->GetCycle();
            if(tree == "MetadataTree") cycle_plotting_order.push_back(cycle);

        }
    } else {
        std::cout << "STATUS : Processing file in custom order - it is: " << std::endl;
        std::cout << "     -> ";
        for( auto &cycle : cycle_plotting_order) std::cout << cycle << " - ";
        std::cout << std::endl << std::endl;
    }


    for(auto &cycle : cycle_plotting_order) {

        if (delete_cycle != 0) {
            if (cycle == delete_cycle) {
                std::string object_to_remove = tree + ";" + get_string(delete_cycle);
                std::cout << " now deleting object " << tree << " cycle " << cycle << " complete name "
                          << object_to_remove << std::endl;
                gDirectory->Delete(object_to_remove.c_str());
                continue;
            }
        }

        std::cout << "STATUS : Processing tree " << tree << " cycle " << cycle << std::endl;

        std::string treename = "MetadataTree;" + get_string(cycle);
        tinF.GetObject(treename.c_str(), tree_meta);
        tree_meta->SetBranchAddress("dataset", &data_ds);
        tree_meta->SetBranchAddress("training_efficiency", &efficiency);
        tree_meta->SetBranchAddress("template_count", &templ_count);
        tree_meta->SetBranchAddress("processed_events", &processed_events);
        tree_meta->SetBranchAddress("mode", &mode);
        tree_meta->SetBranchAddress("wBins0", &wBins[0]);
        tree_meta->SetBranchAddress("zBins0", &zBins[0]);
        tree_meta->GetEntry(0);

        SPcount = wBins[0] * zBins[0];
        spWZratio = (float) wBins[0] / (float) zBins[0];

        std::cout << "  -- SPcount=" << SPcount << " SPratio=" << spWZratio << " wBins=" << wBins[0] << " zBins="
                  << zBins[0] << std::endl;


        if (FILTER) {
            if (1040 < SPcount && SPcount < 1500) {

            } else {
                continue;
            }
        }


        tinF.GetObject(("g_efficiency;" + get_string(cycle)).c_str(), g_efficiency);
        tinF.GetObject(("g_tnumber;" + get_string(cycle)).c_str(), g_tnumber);
        tinF.GetObject(("h_templfreq;" + get_string(cycle)).c_str(), h_templfreq);
        tinF.GetObject(("g_eff_templ;" + get_string(cycle)).c_str(), g_eff_tnumber);


        //make the plots
        pad1->cd();
        g_efficiency->SetMarkerStyle(5);
        g_efficiency->SetLineWidth(LINE_W);

        g_efficiency->GetXaxis()->SetTitle("Number of training events");
        g_efficiency->GetXaxis()->SetLabelFont(43);
        g_efficiency->GetXaxis()->SetLabelSize(14);
        g_efficiency->GetXaxis()->SetTitleFont(63);
        g_efficiency->GetXaxis()->SetTitleSize(14);
        g_efficiency->GetXaxis()->SetTitleOffset(1.4);

        g_efficiency->GetYaxis()->SetTitle("efficiency #epsilon^{cosmic}_{train}");
        g_efficiency->GetYaxis()->SetLabelFont(43);
        g_efficiency->GetYaxis()->SetLabelSize(14);
        g_efficiency->GetYaxis()->SetTitleFont(63);
        g_efficiency->GetYaxis()->SetTitleSize(11);
        g_efficiency->GetYaxis()->SetTitleOffset(1.6);
        g_efficiency->GetYaxis()->CenterTitle(false);

//        if(FIT) {
//            TF1 *f1 = new TF1("f1", "1-[1]*exp(((-1)*x*[2]/100000.))", 6000000, 10000000);
////        TF1 *f1 = new TF1("f1", "1-[1]*exp(((-1)*x*[2]/100000.))", 500000, 15000000);
//            g_efficiency->Fit("f1");
//        }


        //add to multi graph
        g_efficiencies->Add(g_efficiency, "L");

        std::string ltext="#it{" + get_string(wBins[0]) + "x" + get_string(zBins[0]) + "}  |  " +
                          "#it{" + get_string(SPcount) + "}  |  " +
                          "#it{" + (spWZratio < 1 ? "1:" + get_string(1/spWZratio) : get_string(spWZratio) + ":1") + "}  |  "+
                          "#bf{#epsilon_{train}} #it{" + get_string(efficiency * 100).substr(0, 4) + "%}";


        legend->AddEntry(g_efficiency,ltext.c_str(),"L");

        pad2->cd();
        g_tnumber->SetTitle("");
        g_tnumber->SetLineWidth(LINE_W);

        g_tnumber->GetXaxis()->SetTitle("Number of training events");
        g_tnumber->GetXaxis()->SetLabelFont(43);
        g_tnumber->GetXaxis()->SetLabelSize(16);
        g_tnumber->GetXaxis()->SetTitleOffset(3.3);
        g_tnumber->GetXaxis()->SetTitleFont(63);
        g_tnumber->GetXaxis()->SetTitleSize(14);
        g_tnumber->GetXaxis()->SetTickLength(0.05);

        g_tnumber->GetYaxis()->SetTitle("number of templates");
        g_tnumber->GetYaxis()->SetLabelFont(43);
        g_tnumber->GetYaxis()->SetLabelSize(16);
        g_tnumber->GetYaxis()->SetTitleFont(63);
        g_tnumber->GetYaxis()->SetTitleSize(13);

        g_tnumber->GetYaxis()->CenterTitle(false);

        //add to multi graph
        g_tnumbers->Add(g_tnumber, "L");

        /// graph eff - tnumber
        g_eff_tnumber->SetTitle("");
        g_eff_tnumber->SetLineWidth(LINE_W);
        g_eff_tnumber->SetMarkerStyle(5);
        g_eff_tnumber->SetMaximum(1);

        g_eff_tnumber->GetYaxis()->CenterTitle(false);

        if(FIT) {
            TF1 *f1 = new TF1("f1", "1-[0]*exp(((-1)*x*[1]/100000.))", 0, 5000000);
//        TF1 *f1 = new TF1("f1", "1-[1]*exp(((-1)*x*[2]/100000.))", 500000, 15000000);
            g_eff_tnumber->Fit("f1");
        }

        //add to multi graph
        g_eff_tnumbers->Add(g_eff_tnumber, "L");
    }

    TF1 *f2= new TF1("f2", "1-[0]*exp(((-1)*x*[1]/100000.))", 6000000, 10000000);
    f2->SetLineWidth(LINE_W);
    f2->SetLineColor(kRed);

    legend->AddEntry(f2, "Fit function f(x) = a #upoint * exp( #minus b #upoint N_{tmpl})", "L");


    // plot the overview multigraphs

    pad2->cd(); // template count
    g_tnumbers->GetYaxis()->SetRangeUser(0, g_tnumbers->GetYaxis()->GetXmax());
    g_tnumbers->GetXaxis()->SetRangeUser(0, g_efficiencies->GetXaxis()->GetXmax());
    g_tnumbers->GetYaxis()->SetNdivisions(4, 5, 0, true);
    g_tnumbers->GetXaxis()->SetLabelSize(0.08);
    g_tnumbers->GetYaxis()->SetLabelSize(0.08);
    g_tnumbers->GetXaxis()->SetTitle("# training events");
    g_tnumbers->GetYaxis()->SetTitle("# templates");
    g_tnumbers->GetXaxis()->SetTitleSize(0.08);
    g_tnumbers->GetYaxis()->SetTitleSize(0.08);
    g_tnumbers->GetXaxis()->SetTitleFont(52);
    g_tnumbers->GetYaxis()->SetTitleFont(52);
    g_tnumbers->GetYaxis()->SetMaxDigits(2);
    g_tnumbers->GetYaxis()->SetTitleOffset(0.4);

    g_tnumbers->Draw("A PLC PMC");

    pad1->cd(); // efficiency
    g_efficiencies->GetYaxis()->SetRangeUser(0, 1);
    g_efficiencies->GetXaxis()->SetTitle("# training events");
    g_efficiencies->GetYaxis()->SetTitle("cosmic training efficiency #epsilon^{cosmic}_{training}");
    g_efficiencies->GetXaxis()->SetTitleFont(52);
    g_efficiencies->GetYaxis()->SetTitleFont(52);
    g_efficiencies->Draw("A PLC PMC");

    std::string lheadtext="#bf{SP BINS} #it{W#timesZ} | #bf{SP CNT} | #bf{SP RATIO} #it{W:Z} | #bf{TRAIN EFF}";

    legend->SetHeader(lheadtext.c_str(),""); // option "C" allows to center the header
    legend->SetTextAlign(13);
    legend->Draw("");
    saveCanvas(canvas, "CosPatPlots_dataset_" + get_string(dataset) +"_" + filelabel +
    "_id" + get_padded_string(combination_id, 3, '0') + "_overview", pathtorunplots);


    // only efficiency plot

    TCanvas *canvas2 = new TCanvas("canvas2", "Efficiency", 900, 600);
    canvas2->SetGrid(1, 1);
    canvas2->SetTicks(1, 1);

    auto *pad3 = new TPad("template efficiency", "template efficiency", 0, 0, 1, 0.99);
    pad3->SetLogx(0);
    pad3->Draw();
    pad3->cd();

    legend->Draw("");
    g_efficiencies->Draw("A PLC PMC");
    saveCanvas(canvas2, "CosPatPlots_dataset_" + get_string(dataset) +"_" + filelabel + "_id" + get_padded_string(combination_id, 3, '0') + "_Efficiency", pathtorunplots);

    // only template number plot

    TCanvas *canvas3 = new TCanvas("canvas3", "Template Count", 900, 500);
    canvas3->SetGrid(1, 1);
    canvas3->SetTicks(1, 1);

    auto *pad4 = new TPad("template count", "template count", 0, 0, 1, 0.99);
    pad4->SetLogx(0);
    pad4->Draw();
    pad4->cd();

    g_tnumbers->GetYaxis()->SetNdivisions(10, 5, 0, true);
    g_tnumbers->GetXaxis()->SetLabelSize(0.03);
    g_tnumbers->GetYaxis()->SetLabelSize(0.03);
    g_tnumbers->GetXaxis()->SetTitleSize(0.03);
    g_tnumbers->GetYaxis()->SetTitleSize(0.03);
    g_tnumbers->GetXaxis()->SetTitleFont(52);
    g_tnumbers->GetYaxis()->SetTitleFont(52);
    g_tnumbers->GetYaxis()->SetMaxDigits(1);
    g_tnumbers->GetYaxis()->SetTitleOffset(1.6);

    g_tnumbers->Draw("A PLC PMC");
    legend->Draw("C");
    saveCanvas(canvas3, "CosPatPlots_dataset_" + get_string(dataset) +"_" + filelabel + "_id" + get_padded_string(combination_id, 3, '0') + "_TemplateCount", pathtorunplots);


    // eff vs. templ count plot

    TCanvas *c_eff_templ = new TCanvas("c_eff_templ", "Template Bank Pattern Result", 900, 500);
    c_eff_templ->SetGrid(1, 1);
    c_eff_templ->SetTicks(1, 1);
    c_eff_templ->SetLogx(1);

    auto *pad5 = new TPad("eff-template_count", "eff over template count", 0, 0, 1, 1);
    pad5->SetLogx(0);
    pad5->Draw();
    pad5->cd();

    g_eff_tnumbers->GetXaxis()->SetTitle("# templates");
    g_eff_tnumbers->GetYaxis()->SetTitle("efficiency #epsilon^{cosmic}_{train}");
    g_eff_tnumbers->SetMaximum(1);

    g_eff_tnumbers->GetXaxis()->SetTitleFont(52);
    g_eff_tnumbers->GetXaxis()->SetLabelFont(42);
    g_eff_tnumbers->GetXaxis()->SetTitleOffset(1);
    g_eff_tnumbers->GetXaxis()->SetTitleSize(0.04);
    g_eff_tnumbers->GetXaxis()->SetLabelSize(0.035);

    g_eff_tnumbers->GetYaxis()->SetTitleFont(52);
    g_eff_tnumbers->GetYaxis()->SetLabelFont(42);
    g_eff_tnumbers->GetYaxis()->SetTitleOffset(1);

    g_eff_tnumbers->GetYaxis()->SetTitleSize(0.04);
    g_eff_tnumbers->GetYaxis()->SetLabelSize(0.035);

    g_eff_tnumbers->Draw("A PMC");
    legend->Draw("");
    saveCanvas(c_eff_templ, "CosPatPlots_dataset_" + get_string(dataset) +"_" + filelabel + "_id" + get_padded_string(combination_id, 3, '0') + "_EffvsTemplateCount", pathtorunplots);

    tinF.Close();

}



//read the data set analysis root file and combine eff plots with SP config data as TLatex, make one single nice pdf.
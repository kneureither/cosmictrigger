//
// Created by Konstantin Neureither on 08.11.20.
//


//
// Created by Konstantin Neureither on 07.11.20.
//

#include "utilityFunctions.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TPad.h"
#include "TCanvas.h"
#include "plots.h"

void SetStyle(TH1F* h, std::string name) {
    h->SetTitle(("cosmic muon track " + name).c_str());
    h->SetName(name.c_str());
    h->GetYaxis()->SetTitleOffset(1.6);
    h->GetYaxis()->SetTitleFont(52);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);

    h->GetXaxis()->SetTitleOffset(1);
    h->GetXaxis()->SetTitleFont(52);
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.045);

    h->SetLineWidth(2);
}


void PlotPEconfigs() {
    int dataset = 13;
    int combination_id = 0;
    int cycle = 20;
    const std::string pathtoplots = "output/Mu3eCosPat/";
    const std::string pathtodatasetplots = pathtoplots + "dataset_" + get_padded_string(dataset, 3, '0') + "/";

    check_create_directory(pathtoplots);
    check_create_directory(pathtodatasetplots);

//    gStyle->SetPalette(kThermometer);
    gStyle->SetMarkerStyle(23);
    gStyle->SetHistLineWidth(4);
//    gStyle->SetStatH(0.3);
//    gStyle->SetStatW(0.3);
//    gStyle->SetStatFontSize(0.2);


    // FILE FOR READING
    TFile tinF((pathtodatasetplots + "CTCosPatBuild_dataset_013_id" + get_padded_string(combination_id, 3, '0') + "_plots.root").c_str());
    if (!tinF.IsOpen()) {
        std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
        exit(0);
    }

    std::cout << "opened file" << std::endl;


    std::string hist_title = "h_templfreq";
    std::string h_title_cycle = hist_title + ";" + get_string(cycle);

    TH1F *h;
    tinF.GetObject(h_title_cycle.c_str(), h);


    TCanvas *c1 = new TCanvas("c1", "c1", 900, 650);
    c1->SetTicks(1, 1);
    auto *pad1 = new TPad("binweights", "binweights", 0, 0, 1, 1);
    pad1->Draw();

    pad1->cd();
    pad1->SetLogy(1);
    pad1->SetLeftMargin(1.5);
    pad1->SetRightMargin(0.1);
    h->SetStats(0);

    h->SetTitle("");
//    h->GetYaxis()->SetTitleOffset();
    h->GetYaxis()->SetTitle("# templates");
    h->GetYaxis()->SetTitleFont(52);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);

    h->GetXaxis()->SetTitleOffset(1);
    h->GetXaxis()->SetTitleFont(52);
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetRangeUser(0, 50);


    h->SetLineWidth(2);

    h->Draw("COLZ");


    std::string l1 = "#it{#bf{SP BINS }}#it{32#times32 @ 60% TRAIN EFF}";
    std::string l2 = "#it{#bf{TEMPLATES} " + get_string(h->GetEntries()) + "}";
    drawAdditionalInfoBlock(pad1,0.37, 0.8, l1, l2);

    float leg_x=0.5;
    float leg_y=0.8;
    float spacing = 0.05;

    TLatex tline1(leg_x,leg_y,l1.c_str());
    tline1.SetTextFont(43);
    tline1.SetTextSize(20);
    tline1.SetNDC(kTRUE);
    tline1.Draw();
    TLatex tline2(leg_x,leg_y-spacing,l2.c_str());
    tline2.SetTextFont(43);
    tline2.SetTextSize(20);
    tline2.SetNDC(kTRUE);
    tline2.Draw();

    saveCanvas(c1, "plots/PE_Plots_" + h_title_cycle, pathtodatasetplots);

    tinF.Close();
}

int main(int argc, char *argv[]) {
    PlotPEconfigs();
}


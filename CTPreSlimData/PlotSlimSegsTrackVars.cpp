//
// Created by Konstantin Neureither on 07.11.20.
//

#include "utilityFunctions.h"
#include "TFile.h"
#include "TH1F.h"
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


void PlotSlimSegsTrackVars() {
    int dataset = 15;
    const std::string pathtoplots = "output/Mu3eSlimSegs/";
    const std::string pathtodatasetplots = pathtoplots + "dataset_" + get_padded_string(dataset, 3, '0') + "/";

    check_create_directory(pathtoplots);
    check_create_directory(pathtodatasetplots);


    // FILE FOR READING
    TFile tinF((pathtodatasetplots + "SlimSegsOutput.root").c_str());
    if (!tinF.IsOpen()) {
        std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
        exit(0);
    }

    std::cout << "opened file" << std::endl;

    TH1F *h_theta;
    TH1F *h_phi;
    TH1F *h_z0;
    TH1F *h_dca;
    TH1F *h_pt;

    tinF.GetObject("h_thetakari", h_theta);
    tinF.GetObject("h_phikari", h_phi);
    tinF.GetObject("h_z0kari", h_z0);
    tinF.GetObject("h_dcakari", h_dca);
    tinF.GetObject("h_ptkari", h_pt);



    gStyle->SetPalette(kThermometer);
    gStyle->SetMarkerStyle(23);
    gStyle->SetHistLineWidth(4);
    gStyle->SetStatH(0.3);
    gStyle->SetStatW(0.3);
    gStyle->SetStatFontSize(0.2);


    TCanvas *c1 = new TCanvas("c1", "canvasphitheta", 900, 450);
    c1->SetTicks(1, 1);
    auto *pad1 = new TPad("phi", "phi", 0, 0, 0.5, 1);
    pad1->Draw();

    auto *pad2 = new TPad("theta", "theta", 0.5, 0, 1, 1);
    pad2->Draw();

    pad1->cd();
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.15);
    pad1->SetBottomMargin(0.1);
    SetStyle(h_phi, "#phi");
    h_phi->Draw();

    pad2->cd();
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.15);
    pad2->SetBottomMargin(0.1);
    SetStyle(h_theta, "#Theta");
    h_theta->Draw();

    saveCanvas(c1, "Angles_Reconstruction", pathtodatasetplots);


    TCanvas *c2 = new TCanvas("c2", "canvasz0dca", 600, 450);
    c1->SetTicks(1, 1);
    auto *pad3 = new TPad("phi", "phi", 0, 0, 0.33, 1);
    pad3->Draw();

    auto *pad4 = new TPad("theta", "theta", 0.33, 0, 0.66, 1);
    pad4->Draw();

    auto *pad5 = new TPad("pt", "pt", 0.66, 0, 1, 1);
    pad5->Draw();

    pad3->cd();
//    pad3->SetRightMargin(0.15);
    pad3->SetLeftMargin(0.15);
    pad3->SetBottomMargin(0.1);
    SetStyle(h_dca, "DCA");
    h_dca->SetStats(0);
    h_dca->Draw();

    pad4->cd();
//    pad4->SetRightMargin(0.15);
    pad4->SetLeftMargin(0.15);
    pad4->SetBottomMargin(0.1);
    SetStyle(h_z0, "z0");
    h_z0->SetStats(0);
    h_z0->Draw();

    pad5->cd();
    pad5->SetRightMargin(0.15);
    pad5->SetBottomMargin(0.1);
    h_pt->GetYaxis()->SetTitle("count");
    SetStyle(h_pt, "p_{t}");
    h_pt->SetStats(0);
    h_pt->GetXaxis()->SetMaxDigits(2);
    h_pt->Draw();

    saveCanvas(c2, "lengths_Reconstruction", pathtodatasetplots);

    tinF.Close();
}

int main(int argc, char *argv[]) {
    PlotSlimSegsTrackVars();
}
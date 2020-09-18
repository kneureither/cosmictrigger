//
// Created by Konstantin Neureither on 14.09.20.
//

#include "cosmicTemplatesBgROC.h"
#include "TFile.h"
#include "TMultiGraph.h"
#include "TH1F.h"
#include "TStyle.h"

#include <map>
#include <stdlib.h>

#include "utilityFunctions.h"
#include "Mu3eTree.h"
#include "MetaDataTree.h"
#include "bgeval.h"
#include "plots.h"
#include "../Mu3eCosPat/include/PatternEngine.h"
#include "../Mu3eCosPat/include/TemplateData.h"


void cosmicTemplatesBgROC() {
    /*
     * Read the BGEval output files from multiple root files
     * put the data into  two vectors (bg eff and train eff)
     * Create a roc curve for these plots
     */

    const bool RECREATE_FILE = true;
    const int PRINTS = false;
    const int mode = 0;
    const int run = 107;
    const int dataset = 12;
    const int bgevents = 99900;

    std::vector<int> SPcounts = { 3200};
    std::vector<float> SPratios = { 128};
    std::vector<float> tb_stopping_effs = {0.5, 0.6 };

    const std::string pathtoplots = "output/Mu3eCosPatBgEval/dataset_" + get_padded_string(dataset, 3, '0') + "/";
    const std::string pathtooutfile =
            pathtoplots + "bgrun_" + get_padded_string(run, 3, '0') + "/"; //this is where the root file is stored
    const std::string pathtorunplots = pathtooutfile + "PDF/"; //this is where the pdf files are stored

    check_create_directory(pathtoplots);
    check_create_directory(pathtorunplots);
    check_create_directory(pathtooutfile);

    std::vector<float> training_effs;
    std::vector<float> bg_discr_effs;
    std::vector<std::string> ltexts;
    auto g_effROCs = new TMultiGraph();

    gStyle->SetPalette(kBlueGreenYellow);

    TCanvas *canvas = new TCanvas("canvas", "Background Evaluation ROC curve", 900, 600);
    canvas->SetTicks(1, 1);
    auto *pad1 = new TPad("bg eff vs training eff", "bg eff vs training eff", 0, 0, 1, 0.99);
    pad1->SetLogx(0);
//    pad1->SetGrid(1,5);
    pad1->Draw();
    auto legend = new TLegend(0.55, 0.3, 0.9, 0.55);

    for (int spcount : SPcounts) {
        for (auto &spratio : SPratios) {

//            if(spratio == 32 && spcount == 1500){
//                spcount = 1568;
//            }


            // get the filename
            // Get the Pattern Engine and Template Manager
            const int spWbins = (int) sqrt((float) spratio * (float) spcount);
            const int spZbins = (int) sqrt((float) spcount / (float) spratio);
//            std::cout << "(STATUS) : Trying to open config | bgevents " << bgevents << " | run " << run << " | wbins" << spWbins << " | zbins " << spZbins << std::endl;

            // FILE FOR READING BG EVAL DATA
            std::string infile = get_bgevalfile(bgevents, run, dataset, mode, spWbins, spZbins);
            std::cout << "(STATUS) : reading file " << infile << std::endl;
            TFile tinF((pathtooutfile + infile).c_str());
            if (!tinF.IsOpen()) {
                std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
                exit(0);
            }

            TTree *t_meta;
            TTree *t_eff;

            training_effs.clear();
            bg_discr_effs.clear();

            BGAnaResTreeRead *BGAna;

            for (auto &stopp_eff : tb_stopping_effs) {
                std::string pathinfile = "trainingEff" + get_string(stopp_eff) + "/";
                tinF.GetObject((pathinfile + "BackgroundEfficiency").c_str(), t_eff);
                tinF.GetObject((pathinfile + "MetadataTree").c_str(), t_meta);

                BGAna = new BGAnaResTreeRead(t_meta, t_eff);

                training_effs.push_back(BGAna->tb_training_eff);
                bg_discr_effs.push_back(BGAna->bg_discr_eff);

                std::cout << "STATUS : sp count " << spcount << " | sp ratio " << BGAna->sp_target_ratio << " | tb_eff " << BGAna->tb_training_eff << " | bg_eff " << BGAna->bg_discr_eff << std::endl;
            }

            std::string ltext="#it{" + get_string(BGAna->wBins[0]) + "x" + get_string(BGAna->zBins[0]) + "}  |  " +
                              "#it{" + get_string(BGAna->wBins[0] * BGAna->zBins[0]) + "}  |  " +
                              "#it{" + (spratio < 1 ? "1:" + get_string(1/spratio) : get_string(spratio) + ":1") + "}";

            TGraph *g_effroc = new TGraph(training_effs.size(), &training_effs[0], &bg_discr_effs[0]);
            legend->AddEntry(g_effroc, ltext.c_str(), "lp");
            g_effroc->SetMarkerStyle(23);
            g_effROCs->Add(g_effroc, "PL");

        }
    }

    g_effROCs->GetXaxis()->SetTitle("training #epsilon");
    g_effROCs->GetXaxis()->SetTitleFont(53);
    g_effROCs->GetXaxis()->SetTitleSize(14);
    g_effROCs->GetXaxis()->SetTitleOffset(1.6);

    g_effROCs->GetYaxis()->SetTitle("background discr. #epsilon");
    g_effROCs->GetYaxis()->SetTitleFont(53);
    g_effROCs->GetYaxis()->SetTitleSize(14);
    g_effROCs->GetYaxis()->SetTitleOffset(1.6);
    g_effROCs->Draw("A PLC PMC");

    std::string lheadtext="#bf{SPBINS} #it{WxZ} | #bf{SPCOUNT} | #bf{SPRATIO} #it{W:Z}";
    legend->SetTextSize(0.02);
    legend->SetHeader(lheadtext.c_str(),"C"); // option "C" allows to center the header
    legend->Draw("C");


    saveCanvas(canvas, "BgEvalROC_dataset_" + get_string(dataset) + "_spc_" + get_string(SPcounts[0]), pathtorunplots);
}

//
// Created by Konstantin Neureither on 15.09.20.
//

#include "BgEvalEff.h"

#include "cosmicTemplatesBgROC.h"
#include "TFile.h"
#include "TMultiGraph.h"
#include "TH1F.h"

#include <map>
#include <stdlib.h>

#include "utilityFunctions.h"
#include "Mu3eTree.h"
#include "MetaDataTree.h"
#include "bgeval.h"
#include "plots.h"
#include "../Mu3eCosPat/include/PatternEngine.h"
#include "../Mu3eCosPat/include/TemplateData.h"


void BgEvalEff() {
    /*
     * Read the BGEval output files from multiple root files
     * put the data into  two vectors (bg eff and train eff)
     * Create a roc curve for these plots
     */

    const bool RECREATE_FILE = true;
    const int PRINTS = false;
    const int mode = 0;
    const int run = 107;
    const int dataset = 9;
    const int bgevents = 99900;

    std::vector<int> SPcounts = {300, 512};
    std::vector<float> SPratios = {1,2,4,8};
    std::vector<float> tb_stopping_effs = {0.6};
    float tb_stopping_eff = tb_stopping_effs[0];

    const std::string pathtoplots = "output/Mu3eCosPatBgEval/";
    const std::string pathtooutfile =
            pathtoplots + "bgrun_" + get_padded_string(run, 3, '0') + "/"; //this is where the root file is stored
    const std::string pathtorunplots = pathtooutfile + "PDF/"; //this is where the pdf files are stored

    check_create_directory(pathtoplots);
    check_create_directory(pathtorunplots);
    check_create_directory(pathtooutfile);

    std::vector<float> training_effs;
    std::vector<float> bg_discr_effs;
    std::vector<float> spcounts; //must be float for graph
    std::vector<float> templ_count;
    std::vector<float> training_eventcount;
    std::vector<std::string> ltexts;
    auto g_effROCs = new TMultiGraph();

    TCanvas *canvas = new TCanvas("canvas", "Background Evaluation Eff", 900, 600);
    canvas->SetTicks(1, 1);
    auto *pad1 = new TPad("bg eff vs sp count", "bg eff vs sp count", 0, 0, 1, 0.99);
    pad1->SetLogx(0);
//    pad1->SetGrid(1,5);
    pad1->Draw();
    auto legend = new TLegend(0.1, 0.1, 0.35, 0.25);


    for (auto &spratio : SPratios) {
        training_effs.clear();
        bg_discr_effs.clear();
        spcounts.clear();
        templ_count.clear();
        training_eventcount.clear();

        for (auto &spcount : SPcounts) {

            // get the filename
            // Get the Pattern Engine and Template Manager
            const int spWbins = (int) sqrt((float) spratio * (float) spcount);
            const int spZbins = (int) sqrt((float) spcount / (float) spratio);

            // FILE FOR READING BG EVAL DATA
            std::string infile = get_bgevalfile(bgevents, run, dataset, mode, spWbins, spZbins);
            TFile tinF((pathtooutfile + infile).c_str());
            if (!tinF.IsOpen()) {
                std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
                exit(0);
            }

            TTree *t_meta;
            TTree *t_eff;

            BGAnaResTreeRead *BGAna;

            std::string pathinfile = "trainingEff" + get_string(tb_stopping_eff) + "/";
            tinF.GetObject((pathinfile + "BackgroundEfficiency").c_str(), t_eff);
            tinF.GetObject((pathinfile + "MetadataTree").c_str(), t_meta);

            BGAna = new BGAnaResTreeRead(t_meta, t_eff);

            training_effs.push_back(BGAna->tb_training_eff);
            bg_discr_effs.push_back(BGAna->bg_discr_eff);
            spcounts.push_back((float) BGAna->sp_count);
            templ_count.push_back((float) BGAna->template_count);
            training_eventcount.push_back((float) BGAna->tb_training_eventcount);

            std::cout << "STATUS : sp count " << spcount << " | train events " << BGAna->tb_training_eventcount << " | bg events " << BGAna->bg_events << " | templ count " << BGAna->template_count;
            std::cout << " | tb_eff " << BGAna->tb_training_eff << " | bg_eff " << BGAna->bg_discr_eff << std::endl;
        }


//        std::string ltext="#it{" + get_string(BGAna->wBins[0]) + "x" + get_string(BGAna->zBins[0]) + "}  |  " +
//                          "#it{" + get_string(BGAna->wBins[0] * BGAna->zBins[0]) + "}  |  " +
//                          "#it{" + (spratio < 1 ? "1:" + get_string(1/spratio) : get_string(spratio) + ":1") + "}";
        std::string ltext = "#bf{TB EFF} " + get_string(tb_stopping_eff) + + " | #bf{SP RATIO} #it{" "#it{" + (spratio < 1 ? "1:" + get_string(1/spratio) : get_string(spratio) + ":1") + "}}";

        TGraph *g_effroc = new TGraph(spcounts.size(), &spcounts[0], &bg_discr_effs[0]);
        legend->AddEntry(g_effroc, ltext.c_str());
        g_effroc->SetMarkerColor(5);
        g_effROCs->Add(g_effroc, "PL"); //no markers shown
    }
    g_effROCs->GetXaxis()->SetTitle("# sp central area");
    g_effROCs->GetXaxis()->SetTitleFont(53);
    g_effROCs->GetXaxis()->SetTitleSize(14);
    g_effROCs->GetXaxis()->SetTitleOffset(1.6);

    g_effROCs->GetYaxis()->SetTitle("background discr. #epsilon");
    g_effROCs->GetYaxis()->SetTitleFont(53);
    g_effROCs->GetYaxis()->SetTitleSize(14);
    g_effROCs->GetYaxis()->SetTitleOffset(1.6);
    g_effROCs->Draw("A PLC PMC");

    std::string lheadtext="#bf{SP CONFIGURATION} #it{DATASET " + get_string(dataset) + "}";
    legend->SetTextSize(0.02);
    legend->SetHeader(lheadtext.c_str(),"C"); // option "C" allows to center the header
    legend->Draw("C");


    saveCanvas(canvas, "BgEvalEffSPC_" + get_string(dataset), pathtorunplots);
}

//
// Created by Konstantin Neureither on 14.09.20.
//

#include "cosmicTemplatesBgROC.h"
#include "bgeval.h"
#include "TFile.h"
#include "TMultiGraph.h"
#include "plots.h"


//
// Created by Konstantin Neureither on 25.08.20.
//

#include "TFile.h"
#include "TH1F.h"

#include <cassert>
#include <map>
#include <stdlib.h>

#include "../inc/cosmicTemplatesBgEval.h"
#include "utilityFunctions.h"
#include "Mu3eTree.h"
#include "MetaDataTree.h"
#include "../Mu3eCosPat/include/PatternEngine.h"
#include "../Mu3eCosPat/include/TemplateBank.h"
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
    const int dataset = 9;
    const int bgevents = 40000;

    std::vector<int> SPcounts = {200, 400, 600, 800};
    std::vector<float> SPratios = {1.0};
    std::vector<float> tb_stopping_effs = {0.7, 0.75, 0.8, 0.85};

    const std::string pathtoplots = "output/Mu3eCosPatBgEval/";
    const std::string pathtooutfile =
            pathtoplots + "bgrun_" + get_padded_string(run, 3, '0') + "/"; //this is where the root file is stored
    const std::string pathtorunplots = pathtooutfile + "/PDF/"; //this is where the pdf files are stored

    check_create_directory(pathtoplots);
    check_create_directory(pathtorunplots);
    check_create_directory(pathtooutfile);

    std::vector<float> training_effs;
    std::vector<float> bg_discr_effs;
    std::vector<std::string> ltexts;
    auto g_effROCs = new TMultiGraph();

    TCanvas *canvas = new TCanvas("canvas", "Background Evaluation ROC curve", 900, 600);
    canvas->SetTicks(1, 1);
    auto *pad1 = new TPad("bg eff vs training eff", "bg eff vs training eff", 0, 0, 1, 0.99);
    pad1->SetLogx(0);
//    pad1->SetGrid(1,5);
    pad1->Draw();
    auto legend = new TLegend(0.6, 0.1, 0.9, 0.7);

    for (auto &spcount : SPcounts) {
        for (auto &spratio : SPratios) {

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

            training_effs.clear();
            bg_discr_effs.clear();

            for (auto &stopp_eff : tb_stopping_effs) {
                std::string pathinfile = "trainingEff" + get_string(stopp_eff) + "/";
                tinF.GetObject((pathinfile + "BackgroundEfficiency").c_str(), t_eff);
                tinF.GetObject((pathinfile + "MetadataTree").c_str(), t_meta);

                BGAnaResTreeRead BGAna = BGAnaResTreeRead(t_meta, t_eff);

                training_effs.push_back(BGAna.tb_training_eff);
                bg_discr_effs.push_back(BGAna.bg_discr_eff);

                std::cout << "STATUS : sp count " << spcount << " | tb_eff " << BGAna.tb_training_eff << " | bg_eff " << BGAna.bg_discr_eff << std::endl;
            }

            ltexts.push_back("some legend text");
            TGraph *g_effroc = new TGraph(training_effs.size(), &training_effs[0], &bg_discr_effs[0]);
            g_effroc->SetMarkerColor(5);
            g_effROCs->Add(g_effroc, "PL");
            
        }
    }

    g_effROCs->Draw("A PLC PMC");

    saveCanvas(canvas, "BgEvalROC_" + get_string(dataset), pathtorunplots);
}

//
// Created by Konstantin Neureither on 15.09.20.
//

#include "../inc/bgAnaPlots_EffSPcount.h"

#include "../inc/bgAnaPlots_ROC.h"
#include "TFile.h"
#include "TMultiGraph.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TStyle.h"

#include <map>
#include <stdlib.h>

#include "utilityFunctions.h"
#include "Mu3eTree.h"
#include "MetaDataTree.h"
#include "bgeval.h"
#include "plots.h"
#include "PatternEngine.h"
#include "TemplateData.h"
#include "Configuration.h"
#include "TemplateBank.h"


void BgAnaPlots_EffSPcount() {
    /**
     * Read the BGEval output files from multiple root files
     * put the data into  two vectors (bg eff and sp count)
     * Create Plot for train eff / sp count and templ count / sp count
     */

    Configuration CONFIG;
    CONFIG.set2_spcount().BGANA_PLOT_SPC();

    const bool RECREATE_FILE = true;
    const int PRINTS = CONFIG.prints;
    const int mode = CONFIG.mode;
    const int run = CONFIG.background_run;
    const int dataset = CONFIG.dataset;
    const int bgevents = CONFIG.max_bg_frames;

    std::vector<int> SPcounts = CONFIG.sp_res;
    std::vector<float> SPratios = CONFIG.sp_ratios;
    std::vector<float> tb_stopping_effs = CONFIG.stopping_effs;
    float tb_stopping_eff = tb_stopping_effs[0];
    float max_frame_nhits = CONFIG.max_bg_frame_nhits;
    std::vector<TIDLoadingFilter> filters = CONFIG.TmplBankFilter.filters;

    const std::string pathtoplots = "output/Mu3eCosPatBgEval/dataset_" + get_padded_string(dataset, 3, '0') + "/";
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
    auto g_eff_spress = new TMultiGraph();
    auto g_tcounts_spress = new TMultiGraph();

    gStyle->SetPalette(kRust);

    TCanvas *canvas = new TCanvas("canvas", "Background Evaluation Eff", 900, 750);
    canvas->SetTicks(1, 1);
    auto *pad1 = new TPad("bg eff vs sp count", "bg eff vs sp count", 0, 0.3, 1, 0.99);
    pad1->SetLogx(0);
//    pad1->SetGrid(1,5);
    pad1->Draw();
    gStyle->SetLegendBorderSize(0);
    auto legend1 = new TLegend(0.5, 0.12, 0.89, 0.35);

    auto *pad2 = new TPad("template count vs sp count", "template count vs sp count", 0, 0.01, 1, 0.3);
    pad2->Draw();
    auto legend2 = new TLegend(0.1, 0., 0.35, 0.9);

    int temp_spcount;
    float tb_train_eff_total;
    float tb_train_eff_relative;

    for(auto &filter : filters) {

        for (auto &tb_stopping_eff : tb_stopping_effs) {
            for (auto &spratio : SPratios) {

                training_effs.clear();
                bg_discr_effs.clear();
                spcounts.clear();
                templ_count.clear();
                training_eventcount.clear();

                float tb_train_eff_total=101;
                float tb_train_eff_relative=101;

                for (int spcount : SPcounts) {
                    temp_spcount = spcount;


                    if (spratio == 32 && spcount == 500) {
                        spcount = 1568;
                    }

                    // get the filename
                    // Get the Pattern Engine and Template Manager
                    const int spWbins = (int) sqrt((float) spratio * (float) spcount);
                    const int spZbins = (int) sqrt((float) spcount / (float) spratio);

                    // FILE FOR READING BG EVAL DATA
                    std::string infile = get_bgevalfile(bgevents, CONFIG.max_cosmic_events, run, CONFIG.cosmic_testing_dataset, dataset, mode, spWbins, spZbins );
                    TFile tinF((pathtooutfile + infile).c_str());
                    if (!tinF.IsOpen()) {
                        std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
                        exit(0);
                    }

                    TTree *t_meta;
                    TTree *t_eff;

                    BGAnaResTreeRead *BGAna;

                    std::string pathinfile =
                            enum_to_string(filter) + "/trainingEff" + get_string(tb_stopping_eff) + "/";
                    tinF.GetObject((pathinfile + "BackgroundEfficiency").c_str(), t_eff);
                    tinF.GetObject((pathinfile + "MetadataTree").c_str(), t_meta);

                    BGAna = new BGAnaResTreeRead(t_meta, t_eff);

                    tb_train_eff_relative = BGAna->train_eff_rel;
                    tb_train_eff_total = BGAna->train_eff_total;
                    training_effs.push_back(BGAna->train_eff_total);
                    bg_discr_effs.push_back(BGAna->bg_discr_eff);
                    spcounts.push_back((float) BGAna->wBins[0] * BGAna->zBins[0]);
                    templ_count.push_back((float) BGAna->template_count);
                    training_eventcount.push_back((float) BGAna->tb_training_eventcount);

                    //not implemented
//                    std::cout << filter << " BG " << BGAna->filter;
//                    assert(filter == BGAna->filter);

                    max_frame_nhits = BGAna->max_frame_nhits;

                    std::cout << "STATUS : filter " << enum_to_string(filter) << " | sp count " << BGAna->wBins[0] * BGAna->zBins[0] << " | train events "
                              << BGAna->tb_training_eventcount
                              << " | bg events " << BGAna->bg_events << " | templ count " << BGAna->template_count;
                    std::cout << " | tb_target_eff " << BGAna->tb_stopping_eff << " | tb_eff " << tb_train_eff_total << " | bg_eff " << BGAna->bg_discr_eff
                              << std::endl;
                }


                //        std::string ltext="#it{" + get_string(BGAna->wBins[0]) + "x" + get_string(BGAna->zBins[0]) + "}  |  " +
                //                          "#it{" + get_string(BGAna->wBins[0] * BGAna->zBins[0]) + "}  |  " +
                //                          "#it{" + (spratio < 1 ? "1:" + get_string(1/spratio) : get_string(spratio) + ":1") + "}";

                std::string ltext = "#bf{TB EFF} #it{" + get_string(tb_train_eff_total).substr(0, 4) + " [" + get_string(tb_train_eff_relative).substr(0,4) + " ]} | #bf{SP RATIO} W:Z #it{" +
                                    (spratio < 1 ? "1:" + get_string(1 / spratio) : get_string(spratio) + ":1")
                                    +"} | #bf{FLTR} #it{" + enum_to_string(filter) + "}";

                TGraph *g_eff_spres = new TGraph(spcounts.size(), &spcounts[0], &bg_discr_effs[0]);
                legend1->AddEntry(g_eff_spres, ltext.c_str());
                g_eff_spres->SetMarkerStyle(23);
                g_eff_spress->Add(g_eff_spres, "PL"); //no markers shown

                TGraph *g_tcount_spres = new TGraph(spcounts.size(), &spcounts[0], &templ_count[0]);
                //        legend1->AddEntry(g_tcount_spres, ltext.c_str());
                g_tcount_spres->SetMarkerStyle(23);
                g_tcounts_spress->Add(g_tcount_spres, "PL"); //no markers shown
            }
        }
    }

    //set eff sp res mutliplot style
    pad1->cd();
    g_eff_spress->GetXaxis()->SetTitle("# sp central area");
    g_eff_spress->GetXaxis()->SetTitleFont(53);
    g_eff_spress->GetXaxis()->SetTitleSize(14);
    g_eff_spress->GetXaxis()->SetTitleOffset(1.6);

    g_eff_spress->GetYaxis()->SetTitle("background discr. #epsilon");
    g_eff_spress->GetYaxis()->SetTitleFont(53);
    g_eff_spress->GetYaxis()->SetTitleSize(14);
    g_eff_spress->GetYaxis()->SetTitleOffset(1.9);

    expandYaxisRange(g_eff_spress);
    g_eff_spress->Draw("A PLC PMC");

    std::string lheadtext="#bf{SP CONFIGURATION} #it{DATASET " + get_string(dataset) + "}";
    legend1->SetTextSize(0.02);
    legend1->SetHeader(lheadtext.c_str(), "C"); // option "C" allows to center the header
    legend1->Draw("C");

    std::string effs = "";
    for(auto &eff : tb_stopping_effs)  effs = effs + get_string(100*eff) + "% ";

    std::string l1 = "#it{#bf{COSMIC TRIGGER SIM @ }}#it{" + effs + "TRAINING EFF}";
    std::string l2 = "#it{#bf{MU3E SIM} RUN " + get_string(run) + " | EVENTS " + get_string(bgevents) + " | MAX HITS " + get_string(max_frame_nhits) + "}";
    drawAdditionalInfoBlock(pad1,0.37, 0.8, l1, l2);

    float leg_x=0.15;
    float leg_y=0.75;
    float spacing = 0.05;

    TLatex tline1(leg_x,leg_y,l1.c_str());
    tline1.SetTextFont(43);
    tline1.SetTextSize(14);
    tline1.SetNDC(kTRUE);
    tline1.Draw();
    TLatex tline2(leg_x,leg_y-spacing,l2.c_str());
    tline2.SetTextFont(43);
    tline2.SetTextSize(14);
    tline2.SetNDC(kTRUE);
    tline2.Draw();





    //set template count mutliplot style
    pad2->cd();
    pad2->SetBottomMargin(0.16);
    g_tcounts_spress->GetXaxis()->SetTitle("# sp central area");
    g_tcounts_spress->GetXaxis()->SetTitleSize(0.06);
    g_tcounts_spress->GetXaxis()->SetLabelSize(0.08);
    g_tcounts_spress->GetXaxis()->SetTitleFont(52);
    g_tcounts_spress->GetXaxis()->SetLabelFont(42);
    g_tcounts_spress->GetXaxis()->SetTitleOffset(1);
    g_tcounts_spress->GetYaxis()->SetTitleOffset(1.2);

    g_tcounts_spress->GetYaxis()->SetTitle("# templates in db");
    g_tcounts_spress->GetYaxis()->SetTitleSize(0.06);
    g_tcounts_spress->GetYaxis()->SetLabelSize(0.08);
    g_tcounts_spress->GetYaxis()->SetTitleFont(52);
    g_tcounts_spress->GetYaxis()->SetLabelFont(42);
    g_tcounts_spress->GetYaxis()->SetTitleOffset(0.7);
    expandYaxisRange(g_tcounts_spress);
    g_tcounts_spress->Draw("A PLC PMC");


    std::string ratios = CONFIG.RatiosToString();
    std::string stopp_effs = CONFIG.StoppingEffsToString();
    saveCanvas(canvas, "Plots_bgAnaSPcount_dst_" + get_string(dataset) + CONFIG.set_description +  "_spratio_" + ratios + "_tbeff_" + stopp_effs, pathtorunplots);
}

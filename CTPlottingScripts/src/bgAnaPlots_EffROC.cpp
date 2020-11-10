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

#define USE_RATE true



void BgAnaPlots_ROCdp() {
    /**
     * Read the BGEval output files from multiple root files
     * Using the DBDatapoints vector from Config file
     * put the data into  two vectors (bg eff and sp count)
     * Create Plot for train eff / sp count and templ count / sp count
     */

    /// Config data
    Configuration CONFIG;
    CONFIG.BGANA_PLOT_ROC_DATAPOINTS();

    const bool RECREATE_FILE = true;
    const int PRINTS = CONFIG.prints;
    const int mode = CONFIG.mode;
    const int run = CONFIG.background_run;
    const int dataset = CONFIG.dataset;
    const int bgevents = CONFIG.max_bg_frames;

#if USE_RATE
    const bool DRAW_LEGEND = true;
    const bool LOG_Y = true;
#else
    const bool DRAW_LEGEND = true;
    const bool LOG_Y = false;
#endif


    const std::string pathtoplots = "output/Mu3eCosPatBgEval/dataset_" + get_padded_string(dataset, 3, '0') + "/";
    const std::string pathtooutfile =
            pathtoplots + "bgrun_" + get_padded_string(run, 3, '0') + "/"; //this is where the root file is stored
    const std::string pathtorunplots = pathtooutfile + "PDF/"; //this is where the pdf files are stored

    check_create_directory(pathtoplots);
    check_create_directory(pathtorunplots);
    check_create_directory(pathtooutfile);


    ///Plotting data
    // len=datapoints per curve)
    std::vector<float> training_effs;
    std::vector<float> bg_discr_effs;
    std::vector<float> spcounts;
    std::vector<float> templ_count;
    std::vector<float> training_eventcount;
    std::vector<float> red_rates;

    // len=curves
    std::vector<float> tb_stopping_effs;

    // some values
    float tb_train_eff_total;
    float tb_train_eff_relative;
    int max_frame_nhits;

    //file label
    std::string sprs = "";
    std::string effs = "";

    // graphs
    std::vector<std::string> ltexts;
    auto g_eff_ROCs = new TMultiGraph();

    ///root section
    gStyle->SetPalette(kRust);
    gStyle->SetLegendBorderSize(0);

    TCanvas *canvas = new TCanvas("canvas", "Background Evaluation Eff ROC", 900, 600);
    canvas->SetTicks(1, 1);

    auto *pad1 = new TPad("bg eff vs train eff", "bg eff vs train eff", 0, 0.0, 1, 0.99);
    pad1->SetLogx(0);
    pad1->SetTicks(1,1);

    if(LOG_Y) pad1->SetLogy(1);
    pad1->Draw();

#if USE_RATE
    float leg_x=0.45;
    float leg_y=0.78;
    float spacing = 0.04;
    auto legend1 = new TLegend(leg_x, 0.58, 0.85, leg_y-0.08 );
#else
    float leg_x=0.15;
    float leg_y=0.33;
    float spacing = 0.04;
    auto legend1 = new TLegend(leg_x, 0.13, 0.5, leg_y-0.08 );
#endif


    int curve_idx = 0;
    for(auto &curve : CONFIG.DBconfigDatapoints) {
        for(auto &filter : CONFIG.TmplBankFilter.filters) {

            int run_idx=0;
            std::cout << "-- curve " << curve_idx  << std::endl;

            // prepare data for graphs
            training_effs.clear();
            bg_discr_effs.clear();
            spcounts.clear();
            templ_count.clear();
            training_eventcount.clear();
            red_rates.clear();

            //tmp values
            int spc;
            float spr;
            int zbins;
            int wbins;
            float tb_stopping_eff;
            float red_rate;

            float tb_train_eff_total=101;
            float tb_train_eff_relative=101;

            for (auto& config : curve) {
                spc = config.spc;
                spr = config.spr;
                wbins = config.wbins;
                zbins = config.zbins;
                tb_stopping_eff = config.stopp_eff;

                // get the filename

                // FILE FOR READING BG EVAL DATA
                std::string infile = get_bgevalfile(bgevents, CONFIG.max_cosmic_events, CONFIG.max_bg_frame_nhits,
                                                    run, CONFIG.cosmic_testing_dataset, dataset, mode, config.wbins, config.zbins );
                TFile tinF((pathtooutfile + infile).c_str());
                if (!tinF.IsOpen()) {
                    std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
                    exit(0);
                }

                TTree *t_meta;
                TTree *t_eff;

                BGAnaResTreeRead *BGAna;

                std::string pathinfile =
                        enum_to_string(filter) + "/trainingEff" + get_string(config.stopp_eff) + "/";
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
                max_frame_nhits = BGAna->max_frame_nhits;
                red_rate = 1/(1-BGAna->bg_discr_eff);
                red_rates.push_back(red_rate);


                std::cout << "STATUS : filter " << enum_to_string(filter)
                          << " | sp count " << config.spc
                          << " | train events " << BGAna->tb_training_eventcount
                          << " | bg events " << BGAna->bg_events
                          << " | templ count " << BGAna->template_count
                          << " | tb_target_eff " << BGAna->tb_stopping_eff
                          << " | tb_eff " << tb_train_eff_total
                          << " | bg_eff " << BGAna->bg_discr_eff
                          << " | tb_acc  " << tb_train_eff_relative
                          << " | red_rate " << red_rate
                          << std::endl;

                // checks
                assert(config.spc == BGAna->wBins[0] * BGAna->zBins[0]);

                if(curve_idx==0) {
                    sprs += "_" + get_string(config.spr);
                    tb_stopping_effs.push_back(tb_stopping_eff);
                }
                if(run_idx==0) effs += "_" + get_string(config.stopp_eff);
                run_idx++;

            }

            tb_stopping_effs.push_back(tb_stopping_eff);

            std::string ltext = "#bf{SPC} #it{" + get_string(spc)
                                + "} | #bf{BINS} #it{" + get_string(wbins) + "#times" + get_string(zbins)
                                //                                + "} | "
                                //                                + " [" + get_string(tb_train_eff_relative).substr(0,4) + " ]"
                                //                                + "} | #bf{SPR} W:Z #it{" + (spr < 1 ? "1:" + get_string(1 / spr) : get_string(spr) + ":1")
                                //                                + "} | #bf{Z SPs} fixed @ #it{" + get_string(zbins) + " bins [" + get_string(400 / zbins)+ "mm]"
                                //                                + "} | #bf{W SPs} fixed @ #it{" + get_string(wbins) + " bins"
                                + ( filter != ALL ? ("} | #bf{FLTR} #it{" + enum_to_string(filter)) : "") + "}";

            for(auto &entry : training_effs)  std::cout << entry << " ";
            std::cout << std::endl;
            for(auto &entry : bg_discr_effs)  std::cout << entry << " ";
            std::cout << std::endl;
            for(auto &entry : red_rates)  std::cout << entry << " ";
            std::cout << std::endl;

#if USE_RATE
            TGraph *g_effroc = new TGraph(training_effs.size(), &training_effs[0], &red_rates[0]);
#else
            TGraph *g_effroc = new TGraph(training_effs.size(), &training_effs[0], &bg_discr_effs[0]);
#endif
            g_effroc->SetMarkerStyle(23);
            g_effroc->SetMarkerSize(2);
            g_effroc->SetLineWidth(2);
            if(filter == ALL) g_effroc->SetLineStyle(0);
            if(filter == NO_CENTER) g_effroc->SetLineStyle(7);
            if(filter == CUT_ON_FREQ) g_effroc->SetLineStyle(2);

            legend1->AddEntry(g_effroc, ltext.c_str(), "PL");
            g_eff_ROCs->Add(g_effroc, "PL");

            curve_idx++;

        }
    }

    //set eff sp res mutliplot style
    pad1->cd();
    g_eff_ROCs->GetXaxis()->SetTitle("cosmic efficiency #epsilon_{training}^{cosmic} ");
    g_eff_ROCs->GetXaxis()->SetTitleFont(53);
    g_eff_ROCs->GetXaxis()->SetTitleSize(16);
    g_eff_ROCs->GetXaxis()->SetTitleOffset(1.4);
    g_eff_ROCs->GetXaxis()->SetLimits(0, 1);
    g_eff_ROCs->GetXaxis()->SetRangeUser(0, 1);

#if USE_RATE
    g_eff_ROCs->GetYaxis()->SetTitle("frame selectivity");
#else
    g_eff_ROCs->GetYaxis()->SetTitle("background rejection #epsilon");
#endif
    g_eff_ROCs->GetYaxis()->SetTitleFont(53);
    g_eff_ROCs->GetYaxis()->SetTitleSize(16);
    g_eff_ROCs->GetYaxis()->SetTitleOffset(1.6);

    if(!LOG_Y) expandYaxisRange(g_eff_ROCs);

#if USE_RATE
#else
    g_eff_ROCs->SetMaximum(1);
    g_eff_ROCs->SetMinimum(0.99);
#endif
    g_eff_ROCs->Draw("A PLC PMC");

    std::string lheadtext="#bf{SP CONFIGURATION} #it{DATASET " + get_string(dataset) + "}";
    legend1->SetTextSize(0.025);
//    legend1->SetHeader("", "C"); // option "C" allows to center the header
    if(DRAW_LEGEND) legend1->Draw("C");


    std::string effs2 = "";
    if(tb_stopping_effs.size() < 5) {
        for(auto &eff : tb_stopping_effs) {
            std::string tag=get_string(100*eff);
            if(effs2.find(tag) == std::string::npos) {
                effs2 = effs2 + get_string(100*eff) + "% ";
            }
        }
    } else {
        effs2 = get_string(100*tb_stopping_effs[0]) + "% TO " + get_string(100*tb_stopping_effs[tb_stopping_effs.size()-1]) +"% ";
    }


//    std::string l1 = "#it{#bf{COSMIC TRIGGER SIM} @  }#it{" + effs2 + "TRAINING EFF}";
    std::string l1 = "#it{#bf{COSMIC TRIGGER SIM} @  }#it{UP TO 80% TRAINING EFF}";
    std::string l2 = "#it{#bf{MU3E BKG SIM} RUN " + get_string(run) + " (MICHEL) | EVENTS " + get_string(bgevents)
            + (max_frame_nhits != 0 ? (" | MAX HITS " + get_string(max_frame_nhits)) : "") + "}";
    drawAdditionalInfoBlock(pad1,0.37, 0.8, l1, l2);

    TLatex tline1(leg_x,leg_y,l1.c_str());
    tline1.SetTextFont(43);
    tline1.SetTextSize(14);
    tline1.SetNDC(kTRUE);
    if(DRAW_LEGEND) tline1.Draw();
    TLatex tline2(leg_x,leg_y-spacing,l2.c_str());
    tline2.SetTextFont(43);
    tline2.SetTextSize(14);
    tline2.SetNDC(kTRUE);
    if(DRAW_LEGEND) tline2.Draw();


    std::string ratios = CONFIG.RatiosToString();
    std::string stopp_effs = CONFIG.StoppingEffsToString();

#if USE_RATE
    saveCanvas(canvas, "Plots_bgAnaEffROC_RATE_dst_" + get_string(dataset) + CONFIG.set_description + "_" + CONFIG.FiltersToString() + "_spratio_" + sprs + "_tbeff_" + effs, pathtorunplots);
#else
    saveCanvas(canvas, "Plots_bgAnaEffROC_dst_" + get_string(dataset) + CONFIG.set_description + "_" + CONFIG.FiltersToString() + "_spratio_" + sprs + "_tbeff_" + effs, pathtorunplots);
#endif
}

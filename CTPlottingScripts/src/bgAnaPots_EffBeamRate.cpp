//
// Created by Konstantin Neureither on 17.10.20.
//

#include "bgAnaPots_EffBeamRate.h"


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


void BgAnaPlots_EffBeamRate() {
    /**
     * Read the BGEval output files from multiple root files
     * put the data into  two vectors (bg eff and beam rate (i.e. run number)
     * Create Plot for train eff / sp count and templ count / sp count
     */

    Configuration CONFIG;
    CONFIG.bkg_rates();

    const bool RECREATE_FILE = true;
    const int PRINTS = CONFIG.prints;
    const int mode = CONFIG.mode;
    const int dataset = CONFIG.dataset;

    const std::string basepathtoplots = "output/Mu3eCosPatBgEval/dataset_" + get_padded_string(dataset, 3, '0') + "/";
    const std::string pathtooutfile = basepathtoplots + "BgBeamRateEff/";
    const std::string pathtorunplots = pathtooutfile + "PDF/"; //this is where the pdf files are stored

    check_create_directory(basepathtoplots);
    check_create_directory(pathtooutfile);
    check_create_directory(pathtorunplots);

    //// Plotting data
    //dim = (num curves x datapoints per curve)
    std::vector<std::vector<float>> bg_discr_effs;
    std::vector<std::vector<float>> beam_rates;

    //dim = (num curves)
    std::vector<float> training_effs_total;
    std::vector<float> training_effs_relative;
    std::vector<float> stopping_effs;
    std::vector<float> spcounts_real; //must be float for graph
    std::vector<float> tmpl_counts;
    std::vector<float> training_eventcounts;
    std::vector<std::string> ltexts;


    ////root section
    //set style attributes
    gStyle->SetPalette(kThermometer);
    gStyle->SetMarkerStyle(23);
    gStyle->SetHistLineWidth(4);

    //define data and plots
    auto g_eff_beamrates = new TMultiGraph();

    TCanvas *canvas = new TCanvas("canvas", "Background Evaluation Eff vs. Beamrates", 900, 550);
    canvas->SetTicks(1, 1);
    auto *pad1 = new TPad("bgeff_vs_beamrate", "bg eff vs beam rate", 0, 0, 1, 0.99);
    pad1->Draw();
    gStyle->SetLegendBorderSize(0);
    auto legend1 = new TLegend(0.15, 0.15, 0.5, 0.4);


    //read data section
    float tb_train_eff_total;
    float tb_train_eff_relative;

    //get data for axis labels
    assert(1 == CONFIG.sp_cnt.size());
    assert(1 == CONFIG.sp_ratios.size());
    float spWbins;
    float spZbins;
    float plot_spcount;
    float plot_spratio;

    int curve_idx=0;

    // spc and spr loops are designed to be flexibel later. Acutally this plot should be produced with only one SPconfig
    for(auto &spcount : CONFIG.sp_cnt) {
        for(auto &spratio : CONFIG.sp_ratios) {

            plot_spcount = spcount;
            plot_spratio = spratio;
            spWbins = (int) sqrt((float) spratio * (float) spcount);
            spZbins = (int) sqrt((float) spcount / (float) spratio);

            for(auto &filter : CONFIG.TmplBankFilter.filters) {
                for(auto &tb_stopping_eff : CONFIG.stopping_effs) {
                    for(auto &make_nhit_cut : CONFIG.BkgFiles.make_max_nhit_cut) {
                        int run_idx = 0;
                        float train_eff_total = 0;
                        float train_eff_rel = 0;
                        int sp_count_real = 0;
                        int tmpl_count = 0;
                        int training_eventcnt = 0;

                        std::vector<float> curve_bg_discr_eff;
                        std::vector<float> curve_beam_rates;

                        for(auto &bkg_run : CONFIG.BkgFiles.bg_runs) {

                            // Store parameters
                            int max_bkg_frames = CONFIG.BkgFiles.max_bkg_frames[run_idx];
                            int max_bkg_nhits = (make_nhit_cut ? CONFIG.BkgFiles.max_nhits[run_idx] : 0);
                            int beam_rate = CONFIG.BkgFiles.beam_rates[run_idx] / 2.466;

                            // FILE FOR READING BG EVAL DATA
                            std::string infile = get_bgevalfile(max_bkg_frames, CONFIG.max_cosmic_events, max_bkg_nhits, bkg_run,
                                                                CONFIG.cosmic_testing_dataset, dataset, mode, spWbins, spZbins);
                            TFile tinF((basepathtoplots + + "bgrun_" + get_padded_string(bkg_run, 3, '0') + "/" + infile).c_str());
                            if (!tinF.IsOpen()) {
                                std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
                                exit(0);
                            }

                            //get objects from file here
                            TTree *t_meta;
                            TTree *t_eff;

                            BGAnaResTreeRead *BGAna;

                            std::string pathinfile =
                                    enum_to_string(filter) + "/trainingEff" + get_string(tb_stopping_eff) + "/";
                            tinF.GetObject((pathinfile + "BackgroundEfficiency").c_str(), t_eff);
                            tinF.GetObject((pathinfile + "MetadataTree").c_str(), t_meta);

                            BGAna = new BGAnaResTreeRead(t_meta, t_eff);

                            //append to vectors here
                            curve_bg_discr_eff.push_back(BGAna->bg_discr_eff);
                            curve_beam_rates.push_back(beam_rate);


                            tb_train_eff_relative = BGAna->train_eff_rel;
                            tb_train_eff_total = BGAna->train_eff_total;
                            sp_count_real = BGAna->wBins[0] * BGAna->zBins[0];
                            tmpl_count = BGAna->template_count;
                            training_eventcnt = BGAna->tb_training_eventcount;

                            std::cout << "(STATUS) : Got datapoint for bg_run " << bkg_run
                                      << " | beam rate "  << beam_rate << std::endl;
                            std::cout << "(INFO)   : filter " << enum_to_string(filter)
                                      << " | sp count " << sp_count_real
                                      << " | train events " << training_eventcnt
                                      << " | max bg frames " << BGAna->bg_events
                                      << " | templ count " << BGAna->template_count
                                      << " | tb_stopping_eff " << BGAna->tb_stopping_eff
                                      << " | tb_eff " << tb_train_eff_total
                                      << " | bg_eff " << BGAna->bg_discr_eff
                                      << " | nhit_cut " << (make_nhit_cut == false ? "[none]" : get_string(max_bkg_nhits))
                                      << std::endl;

                            run_idx++;

                        }

                        //make space for new curve
                        bg_discr_effs.push_back(curve_bg_discr_eff);
                        beam_rates.push_back(curve_beam_rates);

                        //store meta data per curve
                        training_effs_total.push_back(tb_train_eff_total);
                        training_effs_relative.push_back(tb_train_eff_relative);
                        spcounts_real.push_back(sp_count_real);
                        tmpl_counts.push_back(tmpl_count);
                        training_eventcounts.push_back(training_eventcnt);

                        std::string ltext;

                        if(make_nhit_cut) {
                            ltext = "#bf{CUT} frames with upper 15% nhits";
                        } else {
                            ltext = "#bf{TB EFF [ACC]} #it{" + get_string(tb_train_eff_total).substr(0, 4) + " [" + get_string(tb_train_eff_relative).substr(0,4) + "]" +
                                                /*"} | #bf{SPRATIO} W:Z #it{" + (spratio < 1 ? "1:" + get_string(1 / spratio) : get_string(spratio) + ":1") +*/
                                                /*"} | #bf{BINS} W#timesZ #it{" + get_string(spWbins) + "#times" + get_string(spZbins) + */
                                                "} | #bf{FLTR} #it{" + enum_to_string(filter) + "}";
                            legend1->AddEntry((TObject*)0, "", "");
                        }


                        // make legend label

                        // create TGraph
                        TGraph *g_eff_beamrate = new TGraph(bg_discr_effs[curve_idx].size(), &beam_rates[curve_idx][0], &bg_discr_effs[curve_idx][0]);
                        g_eff_beamrate->SetLineWidth(2);
                        g_eff_beamrate->SetMarkerSize(2);
                        g_eff_beamrate->SetLineStyle(make_nhit_cut ? 2 : 1);
                        legend1->AddEntry(g_eff_beamrate, ltext.c_str());
                        g_eff_beamrates->Add(g_eff_beamrate, "PL");

                        curve_idx++;
                    }
                }
            }
        }
    }

    std::cout << curve_idx << std::endl;


    //// Make the Multigraph
    //set multiplot style
    pad1->cd();
    g_eff_beamrates->GetXaxis()->SetTitle("Mu3e simulation [normal mode] / target stopping rate");
    g_eff_beamrates->GetXaxis()->SetTitleFont(53);
    g_eff_beamrates->GetXaxis()->SetTitleSize(16);
    g_eff_beamrates->GetXaxis()->SetTitleOffset(1.6);
    g_eff_beamrates->GetXaxis()->SetMaxDigits(2);

    g_eff_beamrates->GetYaxis()->SetTitle(("#bf{CONFIG} bins w#timesz " + get_string(spWbins) + "#times" + get_string(spZbins) + "  /    background rejection #epsilon").c_str());
    g_eff_beamrates->GetYaxis()->SetTitleFont(53);
    g_eff_beamrates->GetYaxis()->SetTitleSize(16);
    g_eff_beamrates->GetYaxis()->SetTitleOffset(1.4);

    expandYaxisRange(g_eff_beamrates);
    g_eff_beamrates->GetXaxis()->SetLimits(2e7, 1e8);
    g_eff_beamrates->GetXaxis()->SetRangeUser(2e7, 1e8);
    g_eff_beamrates->Draw("A PLC PMC");

    std::string lheadtext="#bf{SP CONFIG} #it{DST " + get_string(dataset) + "}";
    legend1->SetTextSize(0.03);
//    legend1->SetHeader(lheadtext.c_str(), "C"); // option "C" allows to center the header
    legend1->Draw("C");




    // Add another lable into the plot

//    std::string effs = "";
//    for(auto &eff : CONFIG.stopping_effs)  effs = effs + get_string(100*eff) + "% ";
//
//    std::string l1 = "#it{#bf{COSMIC TRIGGER SIM @ }}#it{" + effs + "TRAINING EFF}";
////    std::string l2 = "#it{#bf{MU3E SIM} RUN " + get_string(run) + " | EVENTS " + get_string(bgevents) + " | MAX HITS " + get_string(max_frame_nhits) + "}";
//    std::string l2 = "#it{#bf{PE CONFIG} SP RES " + get_string(CONFIG.sp_cnt[0]) + " | WBINS " + get_string() + " | ZBINS " + get_string(max_frame_nhits) + "}";
//    drawAdditionalInfoBlock(pad1,0.37, 0.8, l1, l2);
//
//    float leg_x=0.15;
//    float leg_y=0.75;
//    float spacing = 0.05;
//
//    TLatex tline1(leg_x,leg_y,l1.c_str());
//    tline1.SetTextFont(43);
//    tline1.SetTextSize(14);
//    tline1.SetNDC(kTRUE);
//    tline1.Draw();
//    TLatex tline2(leg_x,leg_y-spacing,l2.c_str());
//    tline2.SetTextFont(43);
//    tline2.SetTextSize(14);
//    tline2.SetNDC(kTRUE);
//    tline2.Draw();


    std::string filters = CONFIG.FiltersToString();
    std::string spcounts = CONFIG.ResolutionsToString();
    std::string ratios = CONFIG.RatiosToString();
    std::string stopp_effs = CONFIG.StoppingEffsToString();
    saveCanvas(canvas, "Plots_bgAnaBeamRate_dst_" + get_string(dataset) + "_flt_" + filters + "_spc_" + spcounts +  "_spratio_" + ratios + "_tbeff_" + stopp_effs, pathtorunplots);
}

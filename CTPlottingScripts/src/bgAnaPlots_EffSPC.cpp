//
// Created by Konstantin Neureither on 15.09.20.
//

//root
#include "TFile.h"
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TStyle.h"

#include <map>
#include <stdlib.h>

//project files
#include "utilityFunctions.h"
#include "Mu3eTree.h"
#include "MetaDataTree.h"
#include "plots.h"
#include "PatternEngine.h"
#include "TemplateData.h"
#include "../../CTCoreModules/Configuration.h"
#include "TemplateBank.h"



void BgAnaPlots_SPC() {
    /**
     * Read the BGEval output files from multiple root files
     * Using the DBDatapoints vector from Config file
     * put the data into  two vectors (bg eff and sp count)
     * Create Plot for train eff / sp count and templ count / sp count
     */

    /// Config data
    Configuration CONFIG;
    CONFIG.PLOT_BGANA_SPC_DATAPOINTS();

    const bool RECREATE_FILE = true;
    const int PRINTS = CONFIG.prints;
    const int mode = CONFIG.mode;
    const int run = CONFIG.background_run;
    const int dataset = CONFIG.dataset;
    const int bgevents = CONFIG.max_bg_frames;

    const std::string pathtoplots = "output/3_BKGEvaluation/dataset_" + get_padded_string(dataset, 3, '0') + "/";
    const std::string pathtooutfile =
            pathtoplots + "bgrun_" + get_padded_string(run, 3, '0') + "/"; //this is where the root file is stored
    const std::string pathtorunplots = pathtooutfile + "PDF_SPCplots/"; //this is where the pdf files are stored

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
    auto g_eff_spress = new TMultiGraph();
    auto g_tcounts_spress = new TMultiGraph();

    float leg_x=0.45;
    float leg_y=0.4;
    float spacing = 0.05;

    ///root section
    gStyle->SetPalette(kRust);

    TCanvas *canvas = new TCanvas("canvas", "Background Evaluation Eff", 900, 750);
    canvas->SetTicks(1, 1);

    auto *pad1 = new TPad("bg eff vs sp count", "bg eff vs sp count", 0, 0.3, 1, 0.99);
    pad1->SetLogx(0);
    pad1->SetTicks(1, 1);
    pad1->Draw();
    gStyle->SetLegendBorderSize(0);
    auto legend1 = new TLegend(leg_x, 0.15, 0.85, leg_y-0.1);

    auto *pad2 = new TPad("template count vs sp count", "template count vs sp count", 0, 0.01, 1, 0.3);
    pad2->Draw();
    pad2->SetTicks(1, 1);
    auto legend2 = new TLegend(0.1, 0., 0.35, 0.9);

    int curve_idx = 0;
    for(auto &filter : CONFIG.TmplBankFilter.filters) {
        for(auto &curve : CONFIG.DBconfigCurveDatapoints) {

            int run_idx=0;

            // prepare data for graphs
            training_effs.clear();
            bg_discr_effs.clear();
            spcounts.clear();
            templ_count.clear();
            training_eventcount.clear();

            //tmp values
            int spc;
            float spr;
            int zbins;
            int wbins;
            float tb_stopping_eff;

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


                std::cout << "STATUS : filter " << enum_to_string(filter)
                          << " | sp count " << config.spc
                          << " | train events " << BGAna->tb_training_eventcount
                          << " | bg events " << BGAna->bg_events
                          << " | templ count " << BGAna->template_count
                          << " | tb_target_eff " << BGAna->tb_stopping_eff
                          << " | tb_eff " << tb_train_eff_total
                          << " | bg_eff " << BGAna->bg_discr_eff
                          << std::endl;

                // checks
                assert(config.spc == BGAna->wBins[0] * BGAna->zBins[0]);

                if(curve_idx==0) sprs += "_" + get_string(config.spr);
                if(run_idx==0) effs += "_" + get_string(config.stopp_eff);
                run_idx++;

            }

            tb_stopping_effs.push_back(tb_stopping_eff);

            std::string ltext = "#bf{TB EFF} #epsilon_{training}^{cosmic} #it{" + get_string(tb_stopping_eff).substr(0, 3)
//                                + " [" + get_string(tb_train_eff_relative).substr(0,4) + " ]"
//                                + "} | #bf{SPR} W:Z #it{" + (spr < 1 ? "1:" + get_string(1 / spr) : get_string(spr) + ":1")
//                                + "} | #bf{Z SPs} fixed @ #it{" + get_string(zbins) + " bins [" + get_string(400 / zbins)+ "mm]"
                                + "} | #bf{W SPs} fixed @ #it{" + get_string(wbins) + " bins"
                                + ( filter != ALL ? ("} | #bf{FLTR} #it{" + enum_to_string(filter)) : "") + "}";

            TGraph *g_eff_spres = new TGraph(spcounts.size(), &spcounts[0], &bg_discr_effs[0]);
            legend1->AddEntry(g_eff_spres, ltext.c_str());
            g_eff_spres->SetMarkerStyle(23);
            g_eff_spres->SetMarkerSize(2);
            g_eff_spres->SetLineWidth(2);
            std::cout << tb_stopping_eff << std::endl;
            if(tb_stopping_eff < 0.7){
//                g_eff_spres->SetLineStyle(7);
            }

//            TF1 *f1 = new TF1("f1", "1-[0]*exp([1] * x)");
//            g_eff_spres->Fit("f1");

            g_eff_spress->Add(g_eff_spres, "PL"); //no markers shown



            TGraph *g_tcount_spres = new TGraph(spcounts.size(), &spcounts[0], &templ_count[0]);
            g_tcount_spres->SetMarkerStyle(23);
            g_tcount_spres->SetMarkerSize(2);
            g_tcount_spres->SetLineWidth(2);
//            if(tb_stopping_eff < 0.7) g_tcount_spres->SetLineStyle(7);
            g_tcounts_spress->Add(g_tcount_spres, "PL"); //no markers shown

            curve_idx++;

        }
    }

    //set eff sp res mutliplot style
    pad1->cd();
    g_eff_spress->GetXaxis()->SetTitle("super pixel count [W bins * Z bins] ");
    g_eff_spress->GetXaxis()->SetTitleFont(53);
    g_eff_spress->GetXaxis()->SetTitleSize(16);
    g_eff_spress->GetXaxis()->SetTitleOffset(1.8);

    g_eff_spress->GetXaxis()->SetLimits(512, 4096);
    g_eff_spress->GetXaxis()->SetRange(512, 4096);

    g_eff_spress->GetYaxis()->SetTitle("background rejection #epsilon");
    g_eff_spress->GetYaxis()->SetTitleFont(53);
    g_eff_spress->GetYaxis()->SetTitleSize(16);
    g_eff_spress->GetYaxis()->SetTitleOffset(1.9);

    expandYaxisRange(g_eff_spress);
    g_eff_spress->Draw("A PLC PMC");

    std::string lheadtext="#bf{SP CONFIGURATION} #it{DATASET " + get_string(dataset) + "}";
    legend1->SetTextSize(0.025);
//    legend1->SetHeader("", "C"); // option "C" allows to center the header
    legend1->Draw("C");

    std::string effs2 = "";
    for(auto &eff : tb_stopping_effs) {
        std::string tag=get_string(100*eff);
        if(effs2.find(tag) == std::string::npos) {
            effs2 = effs2 + get_string(100*eff) + "% ";
        }
    }

    std::string l1 = "#it{#bf{COSMIC TRIGGER SIM @ }}#it{" + effs2 + "TRAINING EFF}";
    std::string l2 = "#it{#bf{MU3E BKG SIM} RUN " + get_string(run) + " | EVENTS " + get_string(bgevents)
            + (max_frame_nhits != 0 ? (" | MAX HITS " + get_string(max_frame_nhits)) : "") + "}";
    drawAdditionalInfoBlock(pad1,0.37, 0.8, l1, l2);

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
    pad2->SetBottomMargin(0.20);
    g_tcounts_spress->GetXaxis()->SetTitle("super pixel count [W bins * Z bins] ");
    g_tcounts_spress->GetXaxis()->SetTitleSize(0.075);
    g_tcounts_spress->GetXaxis()->SetLabelSize(0.08);
    g_tcounts_spress->GetXaxis()->SetTitleFont(52);
    g_tcounts_spress->GetXaxis()->SetLabelFont(42);
    g_tcounts_spress->GetXaxis()->SetTitleOffset(1.3);

    g_tcounts_spress->GetXaxis()->SetLimits(512, 4096);
    g_tcounts_spress->GetXaxis()->SetRange(512, 4096);


    g_tcounts_spress->GetYaxis()->SetTitle("# templates in db");
    g_tcounts_spress->GetYaxis()->SetTitleSize(0.08);
    g_tcounts_spress->GetYaxis()->SetLabelSize(0.08);
    g_tcounts_spress->GetYaxis()->SetTitleFont(52);
    g_tcounts_spress->GetYaxis()->SetLabelFont(42);
    g_tcounts_spress->GetYaxis()->SetTitleOffset(0.5);
    expandYaxisRange(g_tcounts_spress);
    g_tcounts_spress->Draw("A PLC PMC");


    std::string ratios = CONFIG.RatiosToString();
    std::string stopp_effs = CONFIG.StoppingEffsToString();
    saveCanvas(canvas, "Plots_bgAnaEffSpc_dst_" + get_string(dataset) + CONFIG.set_description +  "_spratio_" + sprs + "_tbeff_" + effs, pathtorunplots);
}

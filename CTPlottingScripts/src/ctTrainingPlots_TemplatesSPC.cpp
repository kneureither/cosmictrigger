//
// Created by Konstantin Neureither on 28.10.20.
//

#include "../inc/ctTrainingPlots_TemplatesSPC.h"
#include "Configuration.h"
#include "TStyle.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"

#include "plots.h"

void ctTrainingPlots_TemplatesSPC() {

    Configuration CONFIG;
    CONFIG.set13_spc_tmpl_plot();

    const int mode = CONFIG.mode;
    const int dataset = CONFIG.dataset;

    const std::string basepathtoplots = "output/Mu3eCosPatPlots/dataset_" + get_padded_string(dataset, 3, '0') + "/";
    const std::string pathtooutfile = basepathtoplots + "TempltCntSPC/";
    const std::string pathtorunplots = pathtooutfile + "PDF/"; //this is where the pdf files are stored
    const std::string pathtotempldb = "data/TemplateData/dataset_" + get_padded_string(dataset, 3, '0') + "/";

    check_create_directory(basepathtoplots);
    check_create_directory(pathtooutfile);
    check_create_directory(pathtorunplots);
    check_create_directory(pathtotempldb);

    //// Plotting data
    //dim = (num curves x datapoints per curve)
    std::vector<std::vector<double>> templ_counts;
    std::vector<std::vector<double>> sp_counts;
    std::vector<std::vector<double>> sp_wbins;
    std::vector<std::vector<double>> sp_zbins;

    //dim = (num curves)
    std::vector<double> stopping_effs;
    std::vector<std::string> ltexts;

    ////root section
    //set style attributes
    gStyle->SetPalette(kThermometer);
    gStyle->SetMarkerStyle(23);
    gStyle->SetHistLineWidth(4);

    int colors[2] = {kRed-2, kCyan+4};

    //define data and plots
    auto g_templ_spc_multi = new TMultiGraph();
    auto g_templ_spc_fits = new TMultiGraph();

    TCanvas *canvas = new TCanvas("canvas", "Template count vs. super pixel count", 900, 550);
    canvas->SetTicks(1, 1);
    auto *pad1 = new TPad("tmpl_cnt_vs_spc", "template count vs sp count", 0, 0, 1, 0.99);
    pad1->Draw();
    gStyle->SetLegendBorderSize(0);
    auto legend1 = new TLegend(0.12, 0.6, 0.4, 0.88);

    //file label
    std::string spcounts = "";
    std::string effs = "";

    //root file reading data
    unsigned int template_count;
    float training_eff;
    float stopping_eff;
    int training_events;

    int curve_idx=0;

    for(auto &curve : CONFIG.DBconfigDatapoints) {

        std::vector<double> tcounts;
        std::vector<double> spcs;
        std::vector<double> wbins;
        std::vector<double> zbins;

        int run_idx=0;

        for(auto &config : curve) {

            // FILE FOR READING BG EVAL DATA
            std::string infile = "CosmicPatternDatabase_" + getfileidtag(dataset, mode, config.wbins, config.zbins, config.stopp_eff) + ".root";
            TFile tinF((pathtotempldb + infile).c_str());
            if (!tinF.IsOpen()) {
                std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
                exit(0);
            }

            // get data
            TTree *t_config;
            tinF.GetObject("ConfigTree", t_config);
            t_config->SetBranchAddress("template_count", &template_count);
            t_config->SetBranchAddress("training_efficiency", &training_eff);
            t_config->SetBranchAddress("stopping_efficiency", &stopping_eff);
            t_config->SetBranchAddress("training_events", &training_events);
            t_config->GetEntry(0);

            // store data
            tcounts.push_back((double) template_count);
            spcs.push_back((double) config.spc);
            wbins.push_back((double) config.wbins);
            zbins.push_back((double) config.zbins);

            std::cout << "(STATUS) : Got datapoint for stopping eff " << config.stopp_eff
                      << " | spc "  << config.spc << std::endl;
            std::cout << "(INFO)   : template count " << template_count
                      << " | train events " << training_events
                      << " | wbins " << config.wbins
                      << " | zbins " << config.zbins
                      << " | tb_eff " << training_eff
                      << std::endl;

            if(curve_idx==0) spcounts += "_" + get_string(config.spc);
            if(run_idx==0) effs += "_" + get_string(config.stopp_eff);

            run_idx++;
        }

        //use last curves data to store
        stopping_effs.push_back(stopping_eff);

        templ_counts.push_back(tcounts);
        sp_counts.push_back(spcs);
        sp_wbins.push_back(wbins);
        sp_zbins.push_back(zbins);

        std::string ltext = "#bf{TB EFF} #epsilon^{cosmic}_{training} = #it{" + get_string(stopping_eff) +
                            "} | #bf{Z SPs} fixed @ #it{" + get_string(zbins[0]) + " bins [" + get_string(400 / zbins[0]) + " mm]" + "}";



        // fill graphs and mutligraphs here

        // create TGraph
        TGraph *g_tmpl_spc = new TGraph(sp_counts[curve_idx].size(), &sp_counts[curve_idx][0], &templ_counts[curve_idx][0]);
        g_tmpl_spc->SetLineWidth(2);
        g_tmpl_spc->SetLineStyle(1);
        g_tmpl_spc->SetMarkerSize(2);
        g_tmpl_spc->SetMarkerColor(colors[curve_idx]);
        if(curve_idx != 0) legend1->AddEntry((TObject*)0, "", "");
        legend1->AddEntry(g_tmpl_spc, ltext.c_str(), "P");
        g_templ_spc_multi->Add(g_tmpl_spc, "P");

        //fit the curve
        TF1  *f1 = new TF1("f1","[0]*x+[1]*x*x+[2]*x*x*x+[3]*x*x*x*x",0,4096);
        TF1  *f2 = new TF1("f2","[0]*x+[1]*x*x",0,4096);
//        TF1  *f3 = new TF1("f3","[0]*TMath::Exp(x/[1]) - 1",0,4096);  //Does not converge
        g_tmpl_spc->Fit("f1",0);
        g_tmpl_spc->Fit("f2",0);
        TF1* fit = g_tmpl_spc->GetFunction("f2");
        fit->SetLineColor(colors[curve_idx]);
//        std::string fitlegend = "Fit f(x) = " + get_string(fit->GetParameter(0))+  " x + " + get_string(fit->GetParameter(1)) + " x^{2}";
        std::string fitlegend = "Fit f(x) = ax + bx^{2}";
        legend1->AddEntry(fit, fitlegend.c_str(), "L");



        curve_idx++;
    }

    std::cout << "(STATUS) : Added curve " << curve_idx << std::endl;



    //// Make the Multigraph
    //set multiplot style
    pad1->cd();
    g_templ_spc_multi->GetXaxis()->SetTitle("[W#timesZ bins] / super pixel count ");
    g_templ_spc_multi->GetXaxis()->SetTitleFont(53);
    g_templ_spc_multi->GetXaxis()->SetTitleSize(16);
    g_templ_spc_multi->GetXaxis()->SetTitleOffset(1.6);
    g_templ_spc_multi->GetXaxis()->SetMaxDigits(2);
    g_templ_spc_multi->GetXaxis()->SetLimits(0, 4096);
    g_templ_spc_multi->GetXaxis()->SetRangeUser(0, 4096);
    g_templ_spc_multi->GetXaxis()->SetNdivisions(8, 5, 0, false);
    g_templ_spc_multi->GetXaxis()->SetMaxDigits(4);

    g_templ_spc_multi->GetYaxis()->SetTitle("# templates");
    g_templ_spc_multi->GetYaxis()->SetTitleFont(53);
    g_templ_spc_multi->GetYaxis()->SetTitleSize(16);
    g_templ_spc_multi->GetYaxis()->SetTitleOffset(1.1);

    expandYaxisRange(g_templ_spc_multi);
    g_templ_spc_multi->Draw("A");


//    std::string lheadtext="#bf{SP CONFIG} #it{DATA SET " + get_string(dataset) + "}";
    std::string lheadtext="";
    legend1->SetTextSize(0.03);
    legend1->SetHeader(lheadtext.c_str(), "C"); // option "C" allows to center the header
    legend1->Draw("C");


    saveCanvas(canvas, "Plots_ctTrainingTemplSPC_dst_" + get_string(dataset) + "_spc" + spcounts + "_effs" + effs , pathtorunplots);

}
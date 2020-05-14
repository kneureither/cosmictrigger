#include <TFile.h>
#include <TROOT.h>
#include <algorithm>
#include <TTree.h>
#include <map>
#include <string>
#include <TH1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <fstream>
#include <cmath>
#include "../util/utility_functions.h"
#include "../util/custom_types.h"

using std::cout;
using std::endl;

const std::string pathtodata = "../data/";
const std::string pathtoplots = "../plots/";
const int run = 11;
const bool DETAILED_PRINTS = true;
const bool MAKE_PLOT = true;

void reconstruction_accuracy() {

    gStyle->SetLegendBorderSize(1);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(42);
    gStyle->SetLegendTextSize(0.03);

    std::string runpadded = get_padded_string(run, 6, '0');
    std::string infile1 = pathtodata + "mu3e_run_" + runpadded + ".root";
    std::string infile2 = pathtodata + "mu3e_run_" + runpadded + "_trirec_cosmic.root";

    TFile f1(infile1.c_str());
    TFile f2(infile2.c_str());

    TTree *t_mu3e;
    f1.GetObject("mu3e", t_mu3e);
    TTree *t_segs;
    f2.GetObject("segs", t_segs);
    TTree *t_frames;
    f2.GetObject("frames", t_frames);

    //data in mu3e tree
    unsigned int mu3e_entries = t_mu3e->GetEntries();
    int mu3e_nhits;
    int header[4]; //event, run, type (empty), setup (emtpy)

    std::vector<int> *pixelid = nullptr;
    std::vector<int> *trajids = nullptr;
    std::vector<int> *timestamp = nullptr;
    std::vector<double> *traj_px = nullptr;
    std::vector<double> *traj_py = nullptr;
    std::vector<double> *traj_pz = nullptr;

    t_mu3e->SetBranchAddress("Header", &header);
    t_mu3e->SetBranchAddress("Nhit", &mu3e_nhits);

//    t_mu3e->SetBranchAddress("hit_pixelid",&pixelid);
//    t_mu3e->SetBranchAddress("hit_timestamp", &timestamp);
    t_mu3e->SetBranchAddress("traj_ID", &trajids);
    t_mu3e->SetBranchAddress("traj_px", &traj_px);
    t_mu3e->SetBranchAddress("traj_py", &traj_py);
    t_mu3e->SetBranchAddress("traj_pz", &traj_pz);

    cout << "Branches set for mu3e..." << endl;


    //data for trirec result tree segs
    unsigned int segs_entries = t_segs->GetEntries();
    int rec_event;
    int rec_nhit;
    int rec_trajid;

    float mc_p;
    float mc_pt;
    float mc_theta;
    float mc_phi;
    float mc_lam;
    float mc_vpca_offset;
    float mc_vpca_phi;

    float rec_p;
    float rec_r;
    float rec_tan01[2];
    float rec_lam01[2];

    //for calculation
    double rec_phi;
    double rec_theta;
    double rec_pt;

    t_segs->SetBranchAddress("eventId", &rec_event);
    t_segs->SetBranchAddress("nhit", &rec_nhit);

    t_segs->SetBranchAddress("p", &rec_p);
    t_segs->SetBranchAddress("r", &rec_r);
    t_segs->SetBranchAddress("tan01", &rec_tan01);
    t_segs->SetBranchAddress("lam01", &rec_lam01);

    t_segs->SetBranchAddress("mc_tid", &rec_trajid);
    t_segs->SetBranchAddress("mc_p", &mc_p);
    t_segs->SetBranchAddress("mc_pt", &mc_pt);
    t_segs->SetBranchAddress("mc_theta", &mc_theta);
    t_segs->SetBranchAddress("mc_phi", &mc_phi);
    t_segs->SetBranchAddress("mc_lam", &mc_lam);
    t_segs->SetBranchAddress("mc_vpca_offset", &mc_vpca_offset);
    t_segs->SetBranchAddress("mc_vpca_phi", &mc_vpca_phi);

    cout << "Branches set for segs..." << endl;


    //data for frames tree
    unsigned int frames_entries = t_frames->GetEntries();
    //TODO Maybe switch from segs tree to frames tree
    cout << "Branches set for frames..." << endl;

    //stats data definitions
    int p_fail_count = 0;
    std::vector<float> p_rel_errors;
    std::vector<float> p_inv_rel_errors;
    std::vector<float> p_over_pmcs;
    std::vector<float> rec_rs;
    std::vector<float> rec_p_corrs;
    std::vector<float> rec_ps;
    std::vector<float> rec_inv_ps;
    std::vector<float> rec_pts;
    std::vector<float> rec_phis;
    std::vector<float> rec_thetas;

    std::vector<float> mc_dcas;
    std::vector<float> mc_z_dcas;
    std::vector<float> mc_ps;
    std::vector<float> mc_pts;
    std::vector<int> sim_nhits;
    std::vector<int> rec_nhits;
    std::vector<float> p_inv_rel_errors_hits[6];
    int rec_hits_count[6] = {0};

    unsigned int mu3e_index = 1;
    t_mu3e->GetEntry(mu3e_index);
    mu3e_index++;

    for (unsigned int i = 0; i < segs_entries; i++) {
        t_segs->GetEntry(i);

        //find corresponding entry in segs tree
        while (rec_event > header[0]) {
            t_mu3e->GetEntry(mu3e_index);
            mu3e_index++;
        }

        for (unsigned int hit = 0; hit < (unsigned int) mu3e_nhits; hit++) {
//            PXID pixid = process_pixel_id((*pixelid)[hit]);
            unsigned int trajid = (*trajids)[hit];
            double trajpx = (*traj_px)[hit];
        }

        //get trajectory data from sim tree
        unsigned int trajid = (*trajids)[0];
        double trajpx = (*traj_px)[0];
        double trajpy = (*traj_py)[0];
        double trajpz = (*traj_pz)[0];
        double trajp = sqrt(pow(trajpx, 2) + pow(trajpy, 2) + pow(trajpz, 2));

        //Theta, phi and traverse p of reconstruction
        rec_pt = rec_p * std::cos((rec_lam01[0]));
        rec_phi = rec_tan01[0];
        rec_theta = 3.14159*0.5 - rec_lam01[0];

        //TODO Find good criteria to sort out, these numbers are kind of random
        //mc_p != 0 this should be one
//        if ((trajp / mc_p < 0.99 || trajp / mc_p > 1.01)) {
        if (mc_p == 0) {
            //check, if the mc_p corresponds to calculated impuls from px, py, pz from sim file
            p_fail_count++;
        } else {
            if (true) {
                //if the impulse is correct, compare reconstruction with mc truth info

                //TODO Use direction provided by phi to correct sign of impulse
                //compare impulses of reconstruction and mc trut

                float rec_p_corr = rec_p * sgn(rec_r);
                float mc_inv_p = 1. / mc_p;
                float rec_inv_p = 1. / rec_p_corr;
                float p_inv_rel_error = (mc_inv_p - rec_inv_p);
                float p_rel_error = (mc_p - rec_p_corr);
                float p_over_pmc = rec_p_corr / mc_p;
                float mc_z_dca = std::sin(mc_vpca_phi) * mc_vpca_offset;

                rec_ps.push_back(rec_p);
                rec_inv_ps.push_back(rec_inv_p);
                rec_p_corrs.push_back(rec_p_corr);
                mc_ps.push_back(mc_p);
                p_inv_rel_errors.push_back(p_inv_rel_error);
                p_rel_errors.push_back(p_rel_error);
                p_over_pmcs.push_back(p_over_pmc);

                rec_phis.push_back(rec_phi);
                rec_thetas.push_back(rec_theta);
                rec_pts.push_back(rec_pt);
                rec_rs.push_back(rec_r);
                mc_dcas.push_back(mc_vpca_offset);
                mc_z_dcas.push_back(mc_z_dca);
                mc_pts.push_back(mc_pt);

                sim_nhits.push_back(mu3e_nhits);
                rec_nhits.push_back(rec_nhit);


                switch(rec_nhit) {
                    case 4:
                        p_inv_rel_errors_hits[0].push_back(p_inv_rel_error);
                        break;
                    case 5:
                        p_inv_rel_errors_hits[1].push_back(p_inv_rel_error);
                        break;
                    case 6:
                        p_inv_rel_errors_hits[2].push_back(p_inv_rel_error);
                        break;
                    case 7:
                        p_inv_rel_errors_hits[3].push_back(p_inv_rel_error);
                        break;
                    case 8:
                        p_inv_rel_errors_hits[4].push_back(p_inv_rel_error);
                        break;
                    default:
                        p_inv_rel_errors_hits[5].push_back(p_inv_rel_error);
                }

                ////PRINT SECTION PER ENTRY IN TREE
                if(DETAILED_PRINTS) {
                    //data from reconstruction
                    cout << "rec_event: " << rec_event << " mu3e_event: " << header[0];
                    cout << " mc_p: " << mc_p << " rec_p: " << rec_p*sgn(rec_r) << "rec_nhit: "<< rec_nhit <<"\t\t";

                    //data from simulation
//                    cout << "\t" << "traj id; (px, py, pz) " <<  trajid << "; (" << trajpx << ", " << trajpy << ", " << trajpz << ") ";
//                    cout << "p_calc: " << trajp << "\t";
//                    cout << "\tp_calc / mc_p: " << double(trajp / mc_p);

                    //mc and rec deviations
                    cout << "pt_mc= " << mc_pt << " pt_rec= " << rec_pt << " phi_mc= " << mc_phi << " rec_phi= " << rec_phi;
                    cout << " mc_theta= " << mc_theta << " rec_theta= " << rec_theta << "z-dca=" << mc_z_dca <<"\t\t";

                    cout << "(1/mc_p-1/rec_p)/(1/mc_p): " << p_inv_rel_error*100 << " %" << endl;
                }
            }
        }
    }

    //TODO Print histogram of nhits and rel error of p reconstruction

    if (MAKE_PLOT) {
        //impulses
        const float LEFT_BOUNDARY = -3e4;
        const float RIGHT_BOUNDARY = 3e4 ;
        const int BIN_COUNT = 50;

        std::stringstream filename;
        std::stringstream run_info;
        run_info << "reconstruction-accuracy_run" << run;
        std::string run_info_str = (run_info.str()).c_str();

//        std::stringstream namecompletepdf;
//        namecompletepdf << run_info << ".pdf";
//        std::string namecompletepdf_str = (namecompletepdf.str()).c_str();


        ///FILLING THE HISTOGRAMS

        TH1F *h_p_relerror = new TH1F("h_p_relerror", "p_{rec} p_{truth} deviation", BIN_COUNT, -2.0, 1.0);
        labelAxis(h_p_relerror, "p - p_{mc}) [unit]", "count");
        fillHistWithVector(h_p_relerror, p_rel_errors);

        TH1F *h_pinv_relerror = new TH1F("h_pinv_relerror", "p_{truth}^{-1} #minus p_{rec}^{-1} deviation",
                BIN_COUNT, -0.0002, 0.0002);
        labelAxis(h_pinv_relerror, "p^{-1} - p_{mc}^{-1} [unit]", "count");
        fillHistWithVector(h_pinv_relerror, p_inv_rel_errors);

        TH1F *h_poverpmc = new TH1F("p / p_{mc} relation", "p / p_{mc}", 30, -0, 2);
        labelAxis(h_poverpmc, "#frac{p}{p_{mc}} [unit]", "count");
        fillHistWithVector(h_poverpmc, p_over_pmcs);

        TH1F *h_nhits = new TH1F("h_nhits", "hits per frame", 20, 0., 20.);
        labelAxis(h_nhits, "number of hits", "count");
        fillHistWithVector(h_nhits, sim_nhits);

        TH1F *h_rec_nhits = new TH1F("h_rec_nhits", "hits per reconstructed frame", 20, 0., 20.);
        labelAxis(h_rec_nhits, "number of hits", "count");
        fillHistWithVector(h_rec_nhits, rec_nhits);

        TH1F *h_p = new TH1F("h_p", "muon impulses", 30, LEFT_BOUNDARY, RIGHT_BOUNDARY);
        labelAxis(h_p, "p [unit]", "count");
        fillHistWithVector(h_p, rec_ps);

        TH1F *h_p_corr = new TH1F("h_p_corr", "muon impulses (sgn corrected)", 30, 0, RIGHT_BOUNDARY);
        labelAxis(h_p_corr, "p_{rec-corr} [unit]", "count");
        fillHistWithVector(h_p_corr, rec_p_corrs);

        TH1F *h_pt = new TH1F("h_pt", "muon traverse impulses", 30, LEFT_BOUNDARY, RIGHT_BOUNDARY);
        labelAxis(h_pt, "p_{t} [unit]", "count");
        fillHistWithVector(h_pt, rec_pts);

        TH1F *h_pmc = new TH1F("h_pmc", "muon impulses monte carlo", 30, LEFT_BOUNDARY, RIGHT_BOUNDARY);
        labelAxis(h_pmc, "p_{mc} [unit]", "count");
        fillHistWithVector(h_pmc, mc_ps);

        TH1F *h_ptmc = new TH1F("h_ptmc", "muon traverse impulses monte carlo", 30, LEFT_BOUNDARY, RIGHT_BOUNDARY);
        labelAxis(h_ptmc, "p_{t-mc} [unit]", "count");
        fillHistWithVector(h_ptmc, mc_pts);

        TH1F *h_dcamc = new TH1F("h_dcamc", "dca monte carlo [unit]", 30, 0, 100);
        labelAxis(h_dcamc, "distance of closest approach", "count");
        fillHistWithVector(h_dcamc, mc_dcas);

        TH1F *h_phi = new TH1F("h_phi", "reconstruction phi", 30, -3, 0);
        labelAxis(h_phi, "#Phi", "count");
        fillHistWithVector(h_phi, rec_phis);

        TH1F *h_theta = new TH1F("h_theta", "reconstruction theta", 30, 0, 3);
        labelAxis(h_theta, "#Theta", "count");
        fillHistWithVector(h_theta, rec_thetas);


        ///FILLING SCATTER PLOTS

        TGraph *g_pdev_phi = new TGraph(p_inv_rel_errors.size(),&p_inv_rel_errors[0], &rec_phis[0]);
        g_pdev_phi->SetTitle("#Phi over p_{rec}^{-1} #minus p_{mc}^{-1} correlation");
        labelAxis(g_pdev_phi, "p_{rec}^{-1} #minus p_{mc}^{-1} [unit]", "#Phi");
        setGraphRange(g_pdev_phi, -4e-4, 4e-3, -3, 0);

        TGraph *g_pdev_dca = new TGraph(p_inv_rel_errors.size(),&p_inv_rel_errors[0], &mc_dcas[0]);
        g_pdev_dca->SetTitle("DCA over p_{rec}^{-1} #minus p_{mc}^{-1} correlation");
        labelAxis(g_pdev_dca, "p_{rec}^{-1} #minus p_{mc}^{-1} [unit]", "dca [unit]");
        setGraphRange(g_pdev_dca, -4e-4, 4e-3, 0, 50);

        TGraph *g_pdev_z_dca = new TGraph(p_inv_rel_errors.size(),&p_inv_rel_errors[0], &mc_z_dcas[0]);
        g_pdev_z_dca->SetTitle("Z-DCA over p_{rec}^{-1} #minus p_{mc}^{-1} correlation");
        labelAxis(g_pdev_z_dca, "p_{rec}^{-1} #minus p_{mc}^{-1} [unit]", "z-dca [mm]");
        //setGraphRange(g_pdev_z_dca, -4e-4, 4e-3, 0, 50);

        TGraph *g_pdev_p = new TGraph(p_rel_errors.size(),&p_inv_rel_errors[0], &rec_ps[0]);
        g_pdev_p->SetTitle("p_{rec} over p_{rec}^{-1} #minus p_{mc}^{-1} correlation");
        labelAxis(g_pdev_p,"p_{rec}^{-1} #minus p_{mc}^{-1} [unit]", "p_{rec} [unit]");
        setGraphRange(g_pdev_p, -4e-4, 4e-3, -1.5e5, 1.5e5);

        TGraph *g_pdev_preccorr = new TGraph(p_rel_errors.size(),&p_inv_rel_errors[0], &rec_p_corrs[0]);
        g_pdev_preccorr->SetTitle("p_{rec-corr} over p_{rec}^{-1} #minus p_{mc}^{-1} correlation (sgn corrected)");
        labelAxis(g_pdev_preccorr,"p_{rec}^{-1} #minus p_{mc}^{-1} [unit]", "p_{rec-corr} [unit]");
        setGraphRange(g_pdev_preccorr, -4e-4, 4e-3, 0, 1.5e5);

        TGraph *g_prec_pmc = new TGraph(rec_ps.size(),&mc_ps[0],&rec_ps[0]);
        g_prec_pmc->SetTitle("p_{rec} over p_{mc} correlation");
        labelAxis(g_prec_pmc, "p_{mc} [unit]", "p_{rec} [unit]");
        setGraphRange(g_prec_pmc, 0, 150e3, -1.5e5, 1.5e5);

        TGraph *g_preccorr_pmc = new TGraph(rec_ps.size(),&mc_ps[0],&rec_p_corrs[0]);
        g_preccorr_pmc->SetTitle("p_{rec-corr} over p_{mc} correlation (sgn corrected)");
        labelAxis(g_preccorr_pmc, "p_{mc} [unit]", "p_{rec-corr} [unit]");
        setGraphRange(g_preccorr_pmc, 0, 150e3, 0, 1.5e5);

        TGraph *g_prec_rrec = new TGraph(rec_ps.size(),&rec_ps[0],&rec_rs[0]);
        g_prec_rrec->SetTitle("r_{rec} over p_{mc} correlation");
        labelAxis(g_prec_rrec, "p_{rec} [unit]", "r_{rec} [unit]");
        setGraphRange(g_prec_rrec, -4e6,4e6,-10e6, 10e6);

        // x-axis p_rec / p_mc plots
        TGraph *g_ppmc_phi = new TGraph(p_over_pmcs.size(),&p_over_pmcs[0],&rec_phis[0]);
        g_ppmc_phi->SetTitle("#phi over p_{rec-corr}/p_{mc} correlation");
        labelAxis(g_ppmc_phi, "#frac{p_{rec-corr}}{p_{mc}} [unit]", "#phi");
        setGraphRange(g_ppmc_phi, 0,2,-3, 0);

        TGraph *g_ppmc_dca = new TGraph(p_over_pmcs.size(),&p_over_pmcs[0],&mc_dcas[0]);
        g_ppmc_dca->SetTitle("DCA over p_{rec-corr}/p_{mc} correlation");
        labelAxis(g_ppmc_dca, "#frac{p_{rec-corr}}{p_{mc}} [unit]", "dca [mm]");
        setGraphRange(g_ppmc_dca, 0,2,0, 60);

        TGraph *g_ppmc_p = new TGraph(p_over_pmcs.size(),&p_over_pmcs[0],&rec_p_corrs[0]);
        g_ppmc_p->SetTitle("p_{rec-corr} over p_{rec-corr}/p_{mc} correlation");
        labelAxis(g_ppmc_p, "#frac{p_{rec-corr}}{p_{mc}} [unit]", "p_{rec-corr} [mm]");
        setGraphRange(g_ppmc_p, 0,2,0, 3e4);

        // p_rec - p_rec corr comparison plots





        //NHIT histograms
        std::vector<TH1F*> h_pinv_relerrorhits;
        for(int i = 0; i < 6; i++) {
            std::stringstream histname, histtitle;
            histname << "h_pinv_relerror_hits_" << i;
            histtitle << "p_{rec}^{-1} relative errors, nhit=" << i + 4;
            TH1F *h = new TH1F((histname.str()).c_str(), (histtitle.str()).c_str(), BIN_COUNT, LEFT_BOUNDARY, RIGHT_BOUNDARY);
            h_pinv_relerrorhits.push_back(h);

            //clear strings
            histname.str(std::string());
            histtitle.str(std::string());
            for(int j=0; j < p_inv_rel_errors_hits[i].size(); j++) {
                h_pinv_relerrorhits[i]->Fill(p_inv_rel_errors_hits[i][j]);
            }
        }

        ///PLOTTING THE HISTOGRAMS

        //overview canvas containing all plots
        auto *c_multi = new TCanvas("cmulti", "cmulti", 1200, 1200);
        c_multi->SetWindowPosition(0, 400);

        c_multi->Divide(4,4);

        {
            //first row
            c_multi->cd(1);
            gPad->SetLogy(1);
            gPad->SetLeftMargin(0.15);
            h_p->Draw();

            c_multi->cd(2);
            gPad->SetLogy(1);
            gPad->SetLeftMargin(0.15);
            h_pmc->Draw();

            c_multi->cd(3);
            gPad->SetLogy(1);
            gPad->SetLeftMargin(0.15);
            h_pt->Draw();

            c_multi->cd(4);
            gPad->SetLogy(1);
            gPad->SetLeftMargin(0.15);
            h_ptmc->Draw();

            //second row
            c_multi->cd(5);
            gPad->SetLogy(1);
            gPad->SetLeftMargin(0.15);
            h_dcamc->Draw();

            c_multi->cd(6);
            gPad->SetLogy(1);
            gPad->SetLeftMargin(0.15);
            h_phi->Draw();

            c_multi->cd(7);
            gPad->SetLogy(1);
            gPad->SetLeftMargin(0.15);
            h_nhits->Draw();

            c_multi->cd(8);
            gPad->SetLogy(1);
            gPad->SetLeftMargin(0.15);
            h_rec_nhits->Draw();

            //third row
            c_multi->cd(9);
            gPad->SetLogy(1);
            gPad->SetLeftMargin(0.15);
            h_theta->Draw();

            c_multi->cd(10);
            gPad->SetLogy(1);
            gPad->SetLeftMargin(0.15);
            h_pinv_relerror->Draw();

            c_multi->cd(11);
            gPad->SetLogy(1);
            gPad->SetLeftMargin(0.15);
            h_poverpmc->Draw();

            c_multi->cd(12);
            gPad->SetLeftMargin(0.15);
            g_prec_rrec->Draw("AP");


            //fourth row
            c_multi->cd(13);
            gPad->SetLeftMargin(0.15);
            g_pdev_phi->Draw("ap");

            c_multi->cd(14);
            gPad->SetLeftMargin(0.15);
            g_pdev_dca->Draw("ap");

            c_multi->cd(15);
            gPad->SetLeftMargin(0.15);
            g_pdev_p->Draw("ap");

            c_multi->cd(16);
            gPad->SetLeftMargin(0.15);
            g_prec_pmc->Draw("AP");


            filename << pathtoplots << run_info_str << "_overview.pdf";
            c_multi->SaveAs((filename.str()).c_str());
            filename.str(std::string());

            filename << pathtoplots << run_info_str << "_plots.pdf(";
            c_multi->Print((filename.str()).c_str(), "pdf");
            filename.str(std::string());
        }

        auto *c_multi2 = new TCanvas("cmulti2", "cmulti2", 1200, 600);
        c_multi2->SetWindowPosition(0, 400);

        c_multi2->Divide(3,1);
        {
            c_multi2->cd(1);
            gPad->SetLeftMargin(0.15);
            g_ppmc_phi->Draw("AP");

            gPad->Update();
            TLine *l=new TLine(1,gPad->GetUymin(),1,gPad->GetUymax());
            l->SetLineColor(kBlue);
            l->Draw();

            c_multi2->cd(2);
            gPad->SetLeftMargin(0.15);
            g_ppmc_dca->Draw("AP");

            gPad->Update();
            TLine *l1=new TLine(1,gPad->GetUymin(),1,gPad->GetUymax());
            l1->SetLineColor(kBlue);
            l1->Draw();

            c_multi2->cd(3);
            gPad->SetLeftMargin(0.15);
            g_ppmc_p->Draw("AP");

            gPad->Update();
            TLine *l2=new TLine(1,gPad->GetUymin(),1,gPad->GetUymax());
            l2->SetLineColor(kBlue);
            l2->Draw();

            filename << pathtoplots << run_info_str << "_poverpmc.pdf";
            c_multi2->SaveAs((filename.str()).c_str());
            filename.str(std::string());

            filename << pathtoplots << run_info_str << "_plots.pdf";
            c_multi2->Print((filename.str()).c_str(), "pdf");
            filename.str(std::string());
        }

        auto *c_multi3 = new TCanvas("cmulti3", "cmulti3", 900, 600);
        c_multi3->SetWindowPosition(0, 400);

        c_multi3->Divide(3,2);
        {
            c_multi3->cd(1);
            gPad->SetLogy(1);
            gPad->SetLeftMargin(0.15);
            h_p->Draw();

            c_multi3->cd(2);
            gPad->SetLeftMargin(0.15);
            g_prec_pmc->Draw("AP");

            c_multi3->cd(3);
            gPad->SetLeftMargin(0.15);
            g_pdev_p->Draw("AP");

            c_multi3->cd(4);
            gPad->SetLogy(1);
            gPad->SetLeftMargin(0.15);
            h_p_corr->Draw();

            c_multi3->cd(5);
            gPad->SetLeftMargin(0.15);
            g_preccorr_pmc->Draw("AP");

            c_multi3->cd(6);
            gPad->SetLeftMargin(0.15);
            g_pdev_preccorr->Draw("AP");

            filename << pathtoplots << run_info_str << "_prec_vs_corr.pdf";
            c_multi3->SaveAs((filename.str()).c_str());
            filename.str(std::string());

            filename << pathtoplots << run_info_str << "_plots.pdf";
            c_multi3->Print((filename.str()).c_str(), "pdf");
            filename.str(std::string());
        }

        auto *c_single1 = new TCanvas("csinlge1", "csinlge1", 600, 600);
        c_single1->SetWindowPosition(0, 400);

        c_single1->SetLeftMargin(0.15);
        g_pdev_z_dca->Draw("AP");

        filename << pathtoplots << run_info_str << "_plots.pdf)";
        c_single1->Print((filename.str()).c_str(), "pdf");
        filename.str(std::string());





        ////###########################################################################

        //group plot canvas p_err 1/p_err and nhit hist
        auto  *c1 = new TCanvas("c", "c", 1200, 600);
        c1->SetWindowPosition(0, 400);

        c1->Divide(3,1);
        c1->cd(1);
        gPad->SetLogy(1);
        gPad->SetLeftMargin(0.15);
        h_p_relerror->Draw();

        c1->cd(2);
        gPad->SetLogy(1);
        gPad->SetLeftMargin(0.15);
        h_pinv_relerror->Draw();

        c1->cd(3);
        gPad->SetLogy(1);
        gPad->SetLeftMargin(0.15);
        h_nhits->Draw();

//        auto legend = new TLegend(30,20);
//        legend->SetHeader("Legend and run data","C"); // option "C" allows to center the header
//        std::stringstream rundata;
//        rundata << "run=" << run << ", events=" << mu3e_entries;
//        legend->AddEntry((TObject*)0, (rundata.str()).c_str(), "");
//        rundata.str(std::string());
//        rundata  << "count p_fail/p_rec=" << p_fail_count << "/" << segs_entries;
//        legend->AddEntry((TObject*)0, (rundata.str()).c_str(), "");
//        legend->Draw();

        filename << pathtoplots << run_info_str << "_hist_3plots_err_nhits.pdf";
        c1->SaveAs((filename.str()).c_str());
        filename.str(std::string());


        //Single histogram for different nhits
        auto  *c2 = new TCanvas();
        c2->SetWindowPosition(0, 500 );
        c2->SetLogy(1);
        auto legend1 = new TLegend();

        h_pinv_relerror->SetFillStyle(3001);
        h_pinv_relerror->SetFillColor(16);
        h_pinv_relerror->SetLineColor(4);
        h_pinv_relerror->SetName("1/p reconstruction error");

        h_pinv_relerror->DrawClone("");

        for(int i=0; i < h_pinv_relerrorhits.size(); i++) {
            h_pinv_relerrorhits[i]->SetLineColor(i+1);
            h_pinv_relerrorhits[i]->SetFillStyle(3001);
            h_pinv_relerrorhits[i]->SetFillColor(i+1);
            h_pinv_relerrorhits[i]->Draw("same");
            std::stringstream histtitle;
            histtitle << "nhit=" << i + 4;
            legend1->AddEntry(h_pinv_relerrorhits[i], (histtitle.str()).c_str(), "l");
        }

        legend1->AddEntry(h_pinv_relerror, "all hits", "l");
        legend1->DrawClone("");

        filename << pathtoplots << run_info_str << "_hist_perr_nhits.pdf";
        c2->SaveAs((filename.str()).c_str());
        filename.str(std::string());


        //Single histogram for nhits
        auto  *c3 = new TCanvas();
        c3->SetWindowPosition(0, 500 );
        c3->SetLogy(1);
        auto legend2 = new TLegend();

        h_nhits->SetLineColor(4);
        h_nhits->SetFillStyle(16);
        h_nhits->SetName("nhits for reconstructed tracks");
        h_nhits->DrawClone("");

        filename << pathtoplots << run_info_str << "_hist_nhits.pdf";
        c3->SaveAs((filename.str()).c_str());
        filename.str(std::string());

        auto *c4 = new TCanvas();
        c4->SetWindowPosition(0, 500 );
        auto legend4 = new TLegend();
        g_pdev_phi->SetTitle("#phi - p_{err} correlation");
        //"#frac{1/p_{err} #minus 1/p_{mc}}{1/p_{mc}}"
        g_pdev_phi->GetXaxis()->SetTitle("p_{inv_err}");
        g_pdev_phi->GetXaxis()->SetRangeUser(-10.,10.);
        g_pdev_phi->GetYaxis()->SetTitle("#phi_{mc}");
        g_pdev_phi->Draw("ap");
        filename << pathtoplots << run_info_str << "_perr-phi.pdf";
        c4->SaveAs((filename.str()).c_str());
        filename.str(std::string());

        auto *c5 = new TCanvas();
        c5->SetWindowPosition(0, 500 );
        auto legend5 = new TLegend();
        g_pdev_dca->SetTitle("Dca - p_{err} correlation");
        g_pdev_dca->GetXaxis()->SetTitle("p_{inv_err}");
        g_pdev_dca->GetYaxis()->SetTitle("#phi_{mc}");
        g_pdev_dca->GetXaxis()->SetRangeUser(-7.,7.);
        g_pdev_dca->GetYaxis()->SetRangeUser(0.,60.);
        g_pdev_dca->Draw("ap");
        filename << pathtoplots << run_info_str << "_perr-dca.pdf";
        c5->SaveAs((filename.str()).c_str());
        filename.str(std::string());

        auto *c6 = new TCanvas();
        c6->SetWindowPosition(0, 500 );
        auto legend6 = new TLegend();
        g_pdev_p->SetTitle("p_{rec} - p_{err} correlation");
        g_pdev_p->GetYaxis()->SetTitle("p_{inv_err}");
        g_pdev_p->GetXaxis()->SetTitle("#p_{rec}");
        g_pdev_p->GetXaxis()->SetRangeUser(-1.e5,1.e5);
        g_pdev_p->GetYaxis()->SetRangeUser(-2.,2.);
        g_pdev_p->Draw("ap");
        filename << pathtoplots << run_info_str << "_p-perr.pdf";
        c6->SaveAs((filename.str()).c_str());
        filename.str(std::string());

        auto *c7 = new TCanvas("c7", "c7", 500,500);
        c7->SetWindowPosition(0, 500 );
//        c7->SetLogy(1);
//        c7->SetLogx(1);
        auto legend7 = new TLegend();
        g_prec_pmc->SetTitle("p_{rec} - p_{mc} correlation");
        g_prec_pmc->GetXaxis()->SetMaxDigits(2);
        g_prec_pmc->GetYaxis()->SetMaxDigits(2);
        g_prec_pmc->GetYaxis()->SetTitle("p_{mc}");
        g_prec_pmc->GetXaxis()->SetTitle("p_{rec}");
        g_prec_pmc->GetXaxis()->SetRangeUser(0,5.e4);
        g_prec_pmc->GetYaxis()->SetRangeUser(0,5.e4);
        g_prec_pmc->Draw("AP");
        filename << pathtoplots << run_info_str << "_p-pmc.pdf";
        c7->SaveAs((filename.str()).c_str());
        filename.str(std::string());

    }

    ////PRINT THE STATS
    float p_inv_rec_error_mean = vector_mean(p_inv_rel_errors);

    cout << endl << endl << "---General Stats---\n" << endl;
    cout << "trajp / p_mc != 1 fail count: " << p_fail_count << ", total: " << segs_entries << ", rate: "
         << (p_fail_count / (float) segs_entries) * 100 << " %" << endl;
    cout << "1 / p reconstruction error mean: "  << p_inv_rec_error_mean * 100 << "% " << endl;

    for(int i=0; i < 5; i++) {
        if(i < 4) {
            cout << i+4 << " hits count: " << p_inv_rel_errors_hits[i].size() << endl;
        } else {
            cout << "other: " << p_inv_rel_errors_hits[i].size() << endl;
        }
    }

}

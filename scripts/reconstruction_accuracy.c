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
//    std::vector<double> *rec_tan01 = nullptr;
    float mc_theta;
    float mc_phi;
    float mc_lam;
    float mc_vpca_offset;
    float mc_vpca_phi;


    float rec_p;
    float rec_r;
    float rec_tan01[2];
    float rec_lam01[2];
//    std::vector<double> *rec_lam01 = nullptr;

    //for calculation
    double rec_phi;
    double rec_theta;
    double rec_pt;

    t_segs->SetBranchAddress("eventId", &rec_event);
    t_segs->SetBranchAddress("nhit", &rec_nhit);

    t_segs->SetBranchAddress("p", &rec_p);
    t_segs->SetBranchAddress("r", &rec_r);

    t_segs->SetBranchAddress("mc_tid", &rec_trajid);
    t_segs->SetBranchAddress("mc_p", &mc_p);
    t_segs->SetBranchAddress("mc_pt", &mc_pt);
    t_segs->SetBranchAddress("mc_theta", &mc_theta);
    t_segs->SetBranchAddress("mc_phi", &mc_phi);
    t_segs->SetBranchAddress("mc_lam", &mc_lam);
    t_segs->SetBranchAddress("mc_vpca_offset", &mc_vpca_offset);
    t_segs->SetBranchAddress("mc_vpca_phi", &mc_vpca_phi);


    t_segs->SetBranchAddress("tan01", &rec_tan01);
    t_segs->SetBranchAddress("lam01", &rec_lam01);

    cout << "Branches set for segs..." << endl;


    //data for frames tree
    unsigned int frames_entries = t_frames->GetEntries();
    //TODO Maybe switch from segs tree to frames tree
    cout << "Branches set for frames..." << endl;

    //stats data definitions
    int p_fail_count = 0;
    std::vector<float> p_rel_errors;
    std::vector<float> p_inv_rel_errors;
    std::vector<float> rec_phis;
    std::vector<float> mc_dcas;
    std::vector<float> mc_ps;
    std::vector<float> rec_ps;
    std::vector<float> rec_inv_ps;
    std::vector<int> nhits;
    std::vector<int> rec_nhits;
    std::vector<float> p_inv_rel_errors_hits[6];
    int rec_hits_count[6] = {0};


    unsigned int mu3e_index = 1;
    t_mu3e->GetEntry(mu3e_index);
    mu3e_index++;

//    for (unsigned int i = 0; i < segs_entries; i++) {
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

        //TODO get theta and phi data for reconstruction
        //Theta, phi and traverse p of reconstruction
        rec_pt = rec_p * std::cos((rec_lam01[0]));
        rec_phi = rec_tan01[0];
        rec_theta = 3.14159*0.5 - rec_lam01[0];

        //TODO Find good criteria to sort out, these numbers are kind of random
        //mc_p != 0 this should be one
        if ((trajp / mc_p < 0.99 || trajp / mc_p > 1.01)) {
            //check, if the mc_p corresponds to calculated impuls from px, py, pz from sim file
            p_fail_count++;
        } else {
            if (true) {
                //if the impulse is correct, compare reconstruction with mc truth info

                //TODO Use direction provided by phi to correct sign of impulse
                //compare impulses of reconstruction and mc truth

                float p_rel_error = ((mc_p) - (rec_p*sgn(rec_r))) / (mc_p);
                p_rel_errors.push_back(p_rel_error);

                float mc_p_inv = 1. / mc_p;
                float rec_p_inv = 1. /rec_p * sgn(rec_r);
                float p_inv_rel_error = ((mc_p_inv) - (rec_p_inv)) / (mc_p_inv);
                p_inv_rel_errors.push_back(p_inv_rel_error);

                rec_phis.push_back(rec_phi);
                mc_dcas.push_back(mc_vpca_offset);

                nhits.push_back(mu3e_nhits);
                rec_nhits.push_back(rec_nhit);
                rec_ps.push_back(rec_p*sgn(rec_r));
                mc_ps.push_back(mc_p);
                rec_inv_ps.push_back(rec_p_inv);

                switch(rec_nhit) {
                    case 4:
                        rec_hits_count[0]++;
                        p_inv_rel_errors_hits[0].push_back(p_inv_rel_error);
                        break;
                    case 5:
                        rec_hits_count[1]++;
                        p_inv_rel_errors_hits[1].push_back(p_inv_rel_error);
                        break;
                    case 6:
                        rec_hits_count[2]++;
                        p_inv_rel_errors_hits[2].push_back(p_inv_rel_error);
                        break;
                    case 7:
                        rec_hits_count[3]++;
                        p_inv_rel_errors_hits[3].push_back(p_inv_rel_error);
                        break;
                    case 8:
                        rec_hits_count[4]++;
                        p_inv_rel_errors_hits[4].push_back(p_inv_rel_error);
                        break;
                    default:
                        rec_hits_count[5]++;
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
                    cout << " mc_theta= " << mc_theta << " rec_theta= " << rec_theta << "\t\t";

                    cout << "(1/mc_p-1/rec_p)/(1/mc_p): " << p_inv_rel_error*100 << " %" << endl;
                }
            }
        }
    }

    //TODO Print histogram of nhits and rel error of p reconstruction

    if (MAKE_PLOT) {
        const float LEFT_BOUNDARY = -2.0;
        const float RIGHT_BOUNDARY = 1.0 ;
        const int BIN_COUNT = 50;

        std::stringstream filename;
        std::stringstream run_info;
        run_info << "run" << run << "_frames" << mu3e_entries << "_recevents" << segs_entries;
        std::string run_info_str = (run_info.str()).c_str();


        ///FILLING THE HISTOGRAMS

        TH1F *h1 = new TH1F("h1", "p reconstruction errors", BIN_COUNT, LEFT_BOUNDARY, RIGHT_BOUNDARY);
        for (int i = 0; i < p_rel_errors.size(); i++) {
            h1->Fill(p_rel_errors[i]);
        }

        TH1F *h2 = new TH1F("h2", "1/p reconstruction errors", BIN_COUNT, LEFT_BOUNDARY, RIGHT_BOUNDARY);
        for (int i = 0; i < p_inv_rel_errors.size(); i++) {
            h2->Fill(p_inv_rel_errors[i]);
        }

        TH1F *h3 = new TH1F("h3", "number of hits per run", 20, 0., 20.);
        for (int i = 0; i < nhits.size(); i++) {
            h3->Fill(nhits[i]);
        }

        ///FILLING SCATTER PLOTS

        TGraph *g_phip = new TGraph(p_inv_rel_errors.size(),&p_inv_rel_errors[0],&rec_phis[0]);
        TGraph *g_dcap = new TGraph(p_inv_rel_errors.size(),&p_inv_rel_errors[0],&mc_dcas[0]);
        TGraph *g_perrp = new TGraph(p_rel_errors.size(),&rec_ps[0],&p_rel_errors[0]);
        TGraph *g_pmcp = new TGraph(rec_ps.size(),&rec_ps[0],&mc_ps[0]);


        //NHIT histograms
        std::vector<TH1F*> hist_hits;
        for(int i = 0; i < 6; i++) {
            std::stringstream histname, histtitle;
            histname << "h2" << i;
            histtitle << "1/p rec error, nhit=" << i + 4;
            TH1F *hh = new TH1F((histname.str()).c_str(), (histtitle.str()).c_str(), BIN_COUNT, LEFT_BOUNDARY, RIGHT_BOUNDARY);
            hist_hits.push_back(hh);
            histname.str(std::string());
            histtitle.str(std::string());
            for(int j=0; j < p_inv_rel_errors_hits[i].size(); j++) {
                hist_hits[i]->Fill(p_inv_rel_errors_hits[i][j]);
            }
        }

        ///PLOTTING THE HISTOGRAMS

        //group plot canvas p err 1/p err and nhit hist
        auto  *c1 = new TCanvas("c", "c", 1200, 600);
        c1->SetWindowPosition(0, 400 );
        c1->SetLogy(1);


        c1->Divide(3,1);
        c1->cd(1);
        c1->SetLogy(1);
        h1->Draw();
        c1->cd(2);
        c1->SetLogy(1);
        h2->Draw();
        c1->cd(3);
        c1->SetLogy(1);
        h3->Draw();

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

        h2->SetFillStyle(3001);
        h2->SetFillColor(16);
        h2->SetLineColor(4);
        h2->SetName("1/p reconstruction error");

        h2->DrawClone("");

        for(int i=0; i < hist_hits.size(); i++) {
            hist_hits[i]->SetLineColor(i+1);
            hist_hits[i]->SetFillStyle(3001);
            hist_hits[i]->SetFillColor(i+1);
            hist_hits[i]->Draw("same");
            std::stringstream histtitle;
            histtitle << "nhit=" << i + 4;
            legend1->AddEntry(hist_hits[i], (histtitle.str()).c_str(), "l");
        }

        legend1->AddEntry(h2, "all hits", "l");
        legend1->DrawClone("");

        filename << pathtoplots << run_info_str << "_hist_perr_nhits.pdf";
        c2->SaveAs((filename.str()).c_str());
        filename.str(std::string());


        //Single histogram for nhits
        auto  *c3 = new TCanvas();
        c3->SetWindowPosition(0, 500 );
        c3->SetLogy(1);
        auto legend2 = new TLegend();

        h3->SetLineColor(4);
        h3->SetFillStyle(16);
        h3->SetName("nhits for reconstructed tracks");
        h3->DrawClone("");

        filename << pathtoplots << run_info_str << "_hist_nhits.pdf";
        c3->SaveAs((filename.str()).c_str());
        filename.str(std::string());

        auto *c4 = new TCanvas();
        c4->SetWindowPosition(0, 500 );
        auto legend4 = new TLegend();
        g_phip->SetTitle("#phi - p_{err} correlation");
        //"#frac{1/p_{err} #minus 1/p_{mc}}{1/p_{mc}}"
        g_phip->GetXaxis()->SetTitle("p_{inv_err}");
        g_phip->GetXaxis()->SetRangeUser(-10.,10.);
        g_phip->GetYaxis()->SetTitle("#phi_{mc}");
        g_phip->Draw("ap");
        filename << pathtoplots << run_info_str << "_perr-phi.pdf";
        c4->SaveAs((filename.str()).c_str());
        filename.str(std::string());

        auto *c5 = new TCanvas();
        c5->SetWindowPosition(0, 500 );
        auto legend5 = new TLegend();
        g_dcap->SetTitle("Dca - p_{err} correlation");
        g_dcap->GetXaxis()->SetTitle("p_{inv_err}");
        g_dcap->GetYaxis()->SetTitle("#phi_{mc}");
        g_dcap->GetXaxis()->SetRangeUser(-7.,7.);
        g_dcap->GetYaxis()->SetRangeUser(0.,60.);
        g_dcap->Draw("ap");
        filename << pathtoplots << run_info_str << "_perr-dca.pdf";
        c5->SaveAs((filename.str()).c_str());
        filename.str(std::string());

        auto *c6 = new TCanvas();
        c6->SetWindowPosition(0, 500 );
        auto legend6 = new TLegend();
        g_perrp->SetTitle("p_{rec} - p_{err} correlation");
        g_perrp->GetYaxis()->SetTitle("p_{inv_err}");
        g_perrp->GetXaxis()->SetTitle("#p_{rec}");
        g_perrp->GetXaxis()->SetRangeUser(-1.e5,1.e5);
        g_perrp->GetYaxis()->SetRangeUser(-2.,2.);
        g_perrp->Draw("ap");
        filename << pathtoplots << run_info_str << "_p-perr.pdf";
        c6->SaveAs((filename.str()).c_str());
        filename.str(std::string());

        auto *c7 = new TCanvas("c7", "c7", 500,500);
        c7->SetWindowPosition(0, 500 );
//        c7->SetLogy(1);
//        c7->SetLogx(1);
        auto legend7 = new TLegend();
        g_pmcp->SetTitle("p_{rec} - p_{mc} correlation");
        g_pmcp->GetXaxis()->SetMaxDigits(2);
        g_pmcp->GetYaxis()->SetMaxDigits(2);
        g_pmcp->GetYaxis()->SetTitle("p_{mc}");
        g_pmcp->GetXaxis()->SetTitle("p_{rec}");
        g_pmcp->GetXaxis()->SetRangeUser(0,5.e4);
        g_pmcp->GetYaxis()->SetRangeUser(0,5.e4);
        g_pmcp->Draw("AP");
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
        cout << i+4 << " hits count: " << rec_hits_count[i] << "\t" << p_inv_rel_errors_hits[i].size() << endl;
    }
    cout << "other: " << rec_hits_count[5] << endl;
}

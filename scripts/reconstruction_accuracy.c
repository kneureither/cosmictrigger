#include <TFile.h>
#include <TROOT.h>
#include <algorithm>
#include <TTree.h>
#include <map>
#include <string>
#include <iostream>
#include <TH1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <fstream>
#include "../util/utility_functions.h"
#include "../util/custom_types.h"


using std::cout;
using std::endl;

const std::string pathtodata = "../data/";
const int run = 7;
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

    TFile
    f1(infile1.c_str());
    TFile
    f2(infile2.c_str());

    TTree *t_mu3e;
    f1.GetObject("mu3e", t_mu3e);
    TTree *t_segs;
    f2.GetObject("segs", t_segs);

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
    float rec_p;

    t_segs->SetBranchAddress("eventId", &rec_event);
    t_segs->SetBranchAddress("nhit", &rec_nhit);
    t_segs->SetBranchAddress("mc_tid", &rec_trajid);
    t_segs->SetBranchAddress("mc_p", &mc_p);
    t_segs->SetBranchAddress("p", &rec_p);

    cout << "Branches set for segs..." << endl;

    //stats data definitions
    int p_fail_count = 0;
    std::vector<float> p_rel_errors;
    std::vector<float> p_inv_rel_errors;
    std::vector<int> nhits;
    std::vector<int> rec_nhits;

    int rec_hits_count[6] = {0};
    std::vector<float> p_inv_rel_errors_hits[6];


    unsigned int mu3e_index = 1;
    t_mu3e->GetEntry(mu3e_index);
    mu3e_index++;

    for (unsigned int i = 0; i < segs_entries; i++) {
        t_segs->GetEntry(i);

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


        //TODO Find good criteria to sort out, these numbers are kind of random
        //mc_p != 0 this should be one
        if (trajp / mc_p < 0.99 || trajp / mc_p > 1.01) {
            //check, if the mc_p corresponds to calculated impuls from px, py, pz from sim file
            p_fail_count++;
        } else {
            if (true) {
                //if the impulse is correct, compare reconstruction with mc truth info

                //compare impulses of reconstruction and mc truth
                float p_rel_error = ((mc_p) - (rec_p)) / (mc_p);
                p_rel_errors.push_back(p_rel_error);

                nhits.push_back(mu3e_nhits);
                rec_nhits.push_back(rec_nhit);

                float mc_p_inv = 1. / mc_p;
                float rec_p_inv = 1. /rec_p;
                float p_inv_rel_error = (mc_p_inv - rec_p_inv) / mc_p_inv;
                p_inv_rel_errors.push_back(p_inv_rel_error);

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
                    cout << " mc_p: " << mc_p << " rec_p: " << rec_p << "rec_nhit: "<< rec_nhit <<"\t\t";

                    //data from simulation
                    cout << "\t" << "traj id; (px, py, pz) " <<  trajid << "; (" << trajpx << ", " << trajpy << ", " << trajpz << ") ";
                    cout << "p_calc: " << trajp << "\t";

                    cout << "\tp_calc / mc_p: " << double(trajp / mc_p) << endl;
//                cout << " relative error: "  << p_rel_error << endl;
                }
            }
        }
    }

    //TODO Print histogram of nhits and rel error of p reconstruction

    const float LEFT_BOUNDARY = -4;
    const float RIGHT_BOUNDARY = 4;
    const int BIN_COUNT = 20;

    if (MAKE_PLOT) {
        TH1F *h1 = new TH1F("h1", "p reconstruction errors", BIN_COUNT, LEFT_BOUNDARY, RIGHT_BOUNDARY);
        for (int i = 0; i < p_rel_errors.size(); i++) {
            h1->Fill(p_rel_errors[i]);
        }
        TH1F *h2 = new TH1F("h2", "1/p reconstruction errors", BIN_COUNT, LEFT_BOUNDARY, RIGHT_BOUNDARY);
        for (int i = 0; i < p_inv_rel_errors.size(); i++) {
            h2->Fill(p_inv_rel_errors[i]);
        }

        //NHIT histograms

        TH1F *h21 = new TH1F("h21", "1/p reconstruction errors", BIN_COUNT, LEFT_BOUNDARY, RIGHT_BOUNDARY);
        for (int i = 0; i < p_inv_rel_errors_hits[0].size(); i++) {
            h21->Fill(p_inv_rel_errors_hits[0][i]);
        }
        TH1F *h22 = new TH1F("h22", "1/p reconstruction errors", BIN_COUNT, LEFT_BOUNDARY, RIGHT_BOUNDARY);
        for (int i = 0; i < p_inv_rel_errors_hits[1].size(); i++) {
            h22->Fill(p_inv_rel_errors_hits[1][i]);
        }
        TH1F *h23 = new TH1F("h23", "1/p reconstruction errors", BIN_COUNT, LEFT_BOUNDARY, RIGHT_BOUNDARY);
        for (int i = 0; i < p_inv_rel_errors_hits[2].size(); i++) {
            h23->Fill(p_inv_rel_errors_hits[2][i]);
            cout << "data: " << p_inv_rel_errors_hits[2][i]<< endl;
        }
        TH1F *h24 = new TH1F("h24", "1/p reconstruction errors", BIN_COUNT, LEFT_BOUNDARY, RIGHT_BOUNDARY);
        for (int i = 0; i < p_inv_rel_errors_hits[3].size(); i++) {
            h24->Fill(p_inv_rel_errors_hits[3][i]);
        }
        TH1F *h25 = new TH1F("h25", "1/p reconstruction errors", BIN_COUNT, LEFT_BOUNDARY, RIGHT_BOUNDARY);
        for (int i = 0; i < p_inv_rel_errors_hits[4].size(); i++) {
            h25->Fill(p_inv_rel_errors_hits[4][i]);
        }
        TH1F *h26 = new TH1F("h26", "1/p reconstruction errors", BIN_COUNT, LEFT_BOUNDARY, RIGHT_BOUNDARY);
        for (int i = 0; i < p_inv_rel_errors_hits[5].size(); i++) {
            h26->Fill(p_inv_rel_errors_hits[5][i]);
        }


        std::vector<TH1F*> hist_hits;
        for(i = 0; i < p_inv_rel_errors_hits.size(); i++) {

            TH1F *hh = new THF1("h");
            hist_hits.push_back(hh);
        }

        TH1F *h3 = new TH1F("h3", "number of hits per run", 20, 0., 20.);
        for (int i = 0; i < nhits.size(); i++) {
            h3->Fill(nhits[i]);
        }

        //Single plots
        auto  *c2 = new TCanvas("c1", "c1");
        c2->SetWindowPosition(0, 500 );
        c2->SetLogy(1);
        h2->SetLineColor(4);
        h2->SetName("1/p reconstruction error");
        h2->Draw("same");

        h23->SetLineColor(5);
        h23->SetName("nhit = 6");
        h23->Draw("same");

        auto legend1 = new TLegend();
        legend1->AddEntry(h1, "p rec error", "");

        legend1->Draw();


        // Second plot
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

        auto legend = new TLegend(30,20);
        legend->SetHeader("Legend and run data","C"); // option "C" allows to center the header
        std::stringstream rundata;
        rundata << "run=" << run << ", events=" << mu3e_entries;
        legend->AddEntry((TObject*)0, (rundata.str()).c_str(), "");
        rundata.str(std::string());
//        rundata  << "count p_fail/p_rec=" << p_fail_count << "/" << segs_entries;
//        legend->AddEntry((TObject*)0, (rundata.str()).c_str(), "");
        legend->Draw();
    }

    float p_inv_rec_error_mean = vector_mean(p_inv_rel_errors);

    ////PRINT THE STATS
    cout << endl << endl << "---General Stats---\n" << endl;
    cout << "trajp / p_mc != 1 fail count: " << p_fail_count << ", total: " << segs_entries << ", rate: "
         << (p_fail_count / (float) segs_entries) * 100 << " %" << endl;
    cout << "1 / p reconstruction error mean: "  << p_inv_rec_error_mean * 100 << "% " << endl;

    for(int i=0; i < 5; i++) {
        cout << i+4 << " hits count: " << rec_hits_count[i] << "\t" << p_inv_rel_errors_hits[i].size() << endl;
    }
    cout << "other: " << rec_hits_count[5] << endl;
}

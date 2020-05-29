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

void read_mu3e() {

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

    cout << "event\tmc_p" << endl;

    int last_event = 0;

    unsigned int mu3e_index = 1;
    for (unsigned int i = 0; i < segs_entries; i++) {
        t_segs->GetEntry(i);

//        if(rec_event == last_event + 1) {
        if(rec_event < 1000 && rec_event > 900) {
            cout << rec_event << "\t" << mc_p << endl;
        }

        last_event = rec_event;

    }
}

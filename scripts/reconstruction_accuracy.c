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
const int run = 6;

void reconstruction_accuracy() {
    std::string runpadded = get_padded_string(run, 6, '0');

    std::string infile1 = pathtodata + "mu3e_run_" + runpadded + ".root";
    std::string infile2 = pathtodata + "mu3e_run_" + runpadded + "_trirec_cosmic2.root";

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

    //stats data
    int p_fail_count = 0;


    unsigned int mu3e_index = 1;
    for (unsigned int i = 0; i < segs_entries; i++) {
        t_segs->GetEntry(i);

        while (true) {
            t_mu3e->GetEntry(mu3e_index);
            mu3e_index++;

            if (rec_event > header[0]) {
                continue;
            } else {
                break;
            }
        }

        for (unsigned int hit = 0; hit < (unsigned int) mu3e_nhits; hit++) {
            //pixelid calc
            unsigned int trajid = (*trajids)[hit];
            double trajpx = (*traj_px)[hit];
        }

        //get traj data from sim tree
        unsigned int trajid = (*trajids)[0];
        double trajpx = (*traj_px)[0];
        double trajpy = (*traj_py)[0];
        double trajpz = (*traj_pz)[0];

        double trajp = sqrt(pow(trajpx, 2) + pow(trajpy, 2) + pow(trajpz, 2));

        if (trajp / mc_p < 0.99 || trajp / mc_p > 1.01) {
            p_fail_count++;
        }

        //data from reconstruction
        cout << "rec_event: " << rec_event << " mu3e_event: " << header[0];
        cout << " mc_p: " << mc_p << " rec_p: " << rec_p << " \t\t";

        //data from simulation
        cout << "\t" << "traj id; (px, py, pz) " <<  trajid << "; (" << trajpx << ", " << trajpy << ", " << trajpz << ") ";
        cout << "p: (calc)" << trajp << "\t";

        cout <<
        double(trajp / mc_p);
        cout << endl;


//        for(unsigned int j = 0; j < (*trajids).size(); j++) {
//            // for each entry in t raj ids
//            // not important, as always index 0 is the desired data.
//            unsigned int trajid = (*trajids)[j];
//            double trajpx = (*traj_px)[j];
//            double trajpy = (*traj_py)[j];
//            double trajpz = (*traj_pz)[j];
//
//            double trajp = sqrt(pow(trajpx, 2) + pow(trajpy, 2) + pow(trajpz, 2));
//
//        }

        //Entry / event / mc_p / rec_p / calc_p
        // 1/p vergleichen
    }

    cout << endl << endl << "---General Stats---\n" << endl;
    cout << "p_mc and trajp fail count: " << p_fail_count << " rate: " << p_fail_count / segs_entries * 100 << " %"
         << endl;

}

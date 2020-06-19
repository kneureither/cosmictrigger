//
// Created by Konstantin Neureither on 18.06.20.
//
#define TRIPLET_HIT_ARRAY_LENGTH 1024

#ifndef PI
#define PI 3.1415926535
#endif //PI

#ifndef DEBUG
#define DEBUG false
#endif //DEBUG

#include <TFile.h>
#include <TROOT.h>
#include <algorithm>
#include <TTree.h>
#include <map>
#include <string>
#include <TH1.h>
#include <TGraph.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <filesystem>
#include <cmath>
#include <assert.h>
#include "../util/utilityFunctions.h"
#include "../util/plots.h"
#include "../util/trigonometry.h"
#include "../util/rootData.h"
#include "PatternEngine.h"

using std::cout;
using std::endl;

void processSegsPrototype(int run, int FILTER) {

    const std::string pathtodata = "../data/";
    const std::string pathtoplots = "../plots/Mu3eCosPat/processSegsPrototype/";
    const bool MAKE_PLOT = true;
    const int MAX_ENTRIES = 0;

    check_create_directory(pathtodata);
    check_create_directory(pathtoplots);
    std::string pathtorunplots = pathtoplots + "run_" + get_padded_string(run, 3, '0') + "/";
    check_create_directory(pathtorunplots);

    PatternEngine PE(100, 20, 0, pathtorunplots);
    PE.displayBinBoundaries();

    std::string filtertag;

    std::string runpadded = get_padded_string(run, 6, '0');
    std::string infile2 = pathtodata + "mu3e_run_" + runpadded + "_trirec_cosmic.root";
    TFile f(infile2.c_str());
    TTree *t_segs;
    f.GetObject("segs", t_segs);

    //data for trirec result in "segs" tree
    unsigned int segs_entries = t_segs->GetEntries();
    int rec_event;
    int rec_nhit;
    int rec_ntriplet;
    int rec_trajid;

    int mc_type;
    int mc_pid;
    float mc_p;
    float mc_pt;
    float mc_theta;
    float mc_phi;
    float mc_lam;

    float rec_p;
    float rec_r;
    float rec_rt;
    float rec_tan01[TRIPLET_HIT_ARRAY_LENGTH];
    float rec_tan12[TRIPLET_HIT_ARRAY_LENGTH];
    float rec_lam01[TRIPLET_HIT_ARRAY_LENGTH];
    float rec_lam12[TRIPLET_HIT_ARRAY_LENGTH];
    float rec_zpca_x;
    float rec_zpca_y;
    float rec_zpca_z;
    float rec_zpca_r;

    float x00[TRIPLET_HIT_ARRAY_LENGTH];
    float x10[TRIPLET_HIT_ARRAY_LENGTH];
    float x20[TRIPLET_HIT_ARRAY_LENGTH];
    float x01[TRIPLET_HIT_ARRAY_LENGTH];
    float x11[TRIPLET_HIT_ARRAY_LENGTH];
    float x21[TRIPLET_HIT_ARRAY_LENGTH];

    float y00[TRIPLET_HIT_ARRAY_LENGTH];
    float y10[TRIPLET_HIT_ARRAY_LENGTH];
    float y20[TRIPLET_HIT_ARRAY_LENGTH];
    float y01[TRIPLET_HIT_ARRAY_LENGTH];
    float y11[TRIPLET_HIT_ARRAY_LENGTH];
    float y21[TRIPLET_HIT_ARRAY_LENGTH];

    float z00[TRIPLET_HIT_ARRAY_LENGTH];
    float z10[TRIPLET_HIT_ARRAY_LENGTH];
    float z20[TRIPLET_HIT_ARRAY_LENGTH];
    float z01[TRIPLET_HIT_ARRAY_LENGTH];
    float z11[TRIPLET_HIT_ARRAY_LENGTH];
    float z21[TRIPLET_HIT_ARRAY_LENGTH];

    float sid00[TRIPLET_HIT_ARRAY_LENGTH];
    float sid10[TRIPLET_HIT_ARRAY_LENGTH];
    float sid20[TRIPLET_HIT_ARRAY_LENGTH];
    float sid01[TRIPLET_HIT_ARRAY_LENGTH];
    float sid11[TRIPLET_HIT_ARRAY_LENGTH];
    float sid21[TRIPLET_HIT_ARRAY_LENGTH];

    t_segs->SetBranchAddress("eventId", &rec_event);
    t_segs->SetBranchAddress("nhit", &rec_nhit);
    t_segs->SetBranchAddress("n", &rec_ntriplet);

    t_segs->SetBranchAddress("p", &rec_p);
    t_segs->SetBranchAddress("r", &rec_r);
    t_segs->SetBranchAddress("rt", &rec_rt);
    t_segs->SetBranchAddress("tan01", &rec_tan01);
    t_segs->SetBranchAddress("tan12", &rec_tan12);
    t_segs->SetBranchAddress("lam01", &rec_lam01);
    t_segs->SetBranchAddress("lam12", &rec_lam12);
    t_segs->SetBranchAddress("zpca_z", &rec_zpca_z);
    t_segs->SetBranchAddress("zpca_x", &rec_zpca_x);
    t_segs->SetBranchAddress("zpca_y", &rec_zpca_y);
    t_segs->SetBranchAddress("zpca_r", &rec_zpca_r);

    t_segs->SetBranchAddress("mc_type", &mc_type);
    t_segs->SetBranchAddress("mc_pid", &mc_pid);
    t_segs->SetBranchAddress("mc_tid", &rec_trajid);
    t_segs->SetBranchAddress("mc_p", &mc_p);
    t_segs->SetBranchAddress("mc_pt", &mc_pt);
    t_segs->SetBranchAddress("mc_theta", &mc_theta);
    t_segs->SetBranchAddress("mc_phi", &mc_phi);
    t_segs->SetBranchAddress("mc_lam", &mc_lam);

    t_segs->SetBranchAddress("x00", &x00);
    t_segs->SetBranchAddress("x10", &x10);
    t_segs->SetBranchAddress("x20", &x20);
    t_segs->SetBranchAddress("x01", &x01);
    t_segs->SetBranchAddress("x11", &x11);
    t_segs->SetBranchAddress("x21", &x21);

    t_segs->SetBranchAddress("y00", &y00);
    t_segs->SetBranchAddress("y10", &y10);
    t_segs->SetBranchAddress("y20", &y20);
    t_segs->SetBranchAddress("y01", &y01);
    t_segs->SetBranchAddress("y11", &y11);
    t_segs->SetBranchAddress("y21", &y21);

    t_segs->SetBranchAddress("z00", &z00);
    t_segs->SetBranchAddress("z10", &z10);
    t_segs->SetBranchAddress("z20", &z20);
    t_segs->SetBranchAddress("z01", &z01);
    t_segs->SetBranchAddress("z11", &z11);
    t_segs->SetBranchAddress("z21", &z21);

    t_segs->SetBranchAddress("sid00", &sid00);
    t_segs->SetBranchAddress("sid10", &sid10);
    t_segs->SetBranchAddress("sid20", &sid20);
    t_segs->SetBranchAddress("sid01", &sid01);
    t_segs->SetBranchAddress("sid11", &sid11);
    t_segs->SetBranchAddress("sid21", &sid21);

    cout << "Branches set for segs..." << endl;

    //stats data
    int p_fail_count = 0;
    int processed_entries = 0;

    for (unsigned int i = 0; i < (MAX_ENTRIES == 0 ? segs_entries : MAX_ENTRIES); i++) {
        t_segs->GetEntry(i);

        //Theta, phi and traverse p of reconstruction
        float rec_pt = rec_p * std::cos((rec_lam01)[0]);
        float rec_phi = (rec_tan01)[0];
        float rec_theta = PI*0.5 - (rec_lam01)[0];

        //Consolidate data
        std::vector<double> xp;
        std::vector<double> yp;
        std::vector<double> zp;
        std::vector<int> sids;
        std::vector<double> phi_hits;
        std::vector<double> thetas;

        int ncombinedhits = combineHits(xp, yp, zp, sids, phi_hits, thetas, rec_nhit, rec_ntriplet,
                                        x00, x10, x01, x20, x21, y00, y10, y01, y20, y21, z00, z10, z01, z20, z21,
                                        sid00, sid10, sid01, sid20, sid21, rec_tan01, rec_tan12, rec_lam01, rec_lam12);

        if(true) {
            for(int i = 0; i < ncombinedhits; i++) {
                //cout << "\tHits sorted for kari:  sid=" << sids[i] << " x=" << xp[i] << "  y=" << yp[i] << "  z=" << zp[i] <<" phi=" << phi_hits[i] << " theta=" << thetas[i] << endl;
            }
        }

        if (mc_p == 0 || mc_pt == 0 || rec_p == 0 || rec_pt == 0) {
            p_fail_count++;
        } else {

            //Filtering entries //FIXME filtertag is defined over and over again
            bool choice;
            switch(FILTER) {
                case 0:
                    choice = choice = mc_type % 10 == 3 || mc_type % 10 == 4;
                    filtertag = "muonsonly";
                    break;
                case 1:
//                    choice = pt_kari_inv_err / mc_inv_pt < -1.8 && pt_kari_inv_err / mc_inv_pt > -2.2;
                    choice = false;
                    filtertag = "kari-outliers";
                    break;
                case 2:
                    choice = mc_type == 3;
                    filtertag = "mu+";
                    break;
                case 3:
                    choice = mc_type == 4;
                    filtertag = "mu-";
                    break;
                case 4:
//                    choice = sgn(karires.rad) != sgn(mc_p_corr);
                    choice = false;
                    filtertag = "p-kari-mc-sgn";
                    break;
                case 5:
//                    choice = karires.phi > 0 && karires.phi < PI/2 && (pt_kari_inv_err / mc_inv_pt < -1.8 && pt_kari_inv_err / mc_inv_pt > -2.2);
                    choice = false;
                    filtertag = "special1";
                    break;
                case 6:
                    choice = mc_type == 3 || mc_type == 4;
                    filtertag = "planemuonsonly";
                    break;
                case 7:
                    choice = true;
                    filtertag = "all";
                    break;
                default:
                    choice = true;
                    filtertag = "";
            }

            if (choice) {
                processed_entries++;

                for(int i = 0; i<ncombinedhits; i++) {
//                    printf("Getting SP Coordinates for Pixel %d: x=%f, y=%f, z=%f\n", i, xp.at(i),yp.at(i),zp.at(i));
                    int SID = PE.getSuperPixel(xp[i], yp[i], zp[i]);
//                    printf("SID=%d\n", SID);
                }
            }

        }
    }

    PE.displayBinWeightDistribution();
    PE.closePlot();

    printf("Entries processed: %d", processed_entries);

}

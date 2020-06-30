//
// Created by Konstantin Neureither on 25.06.20.
//

#ifndef DEBUG
#define DEBUG false
#endif //DEBUG

// Basic imports
#include <iostream>
#include <string>
#include <stdlib.h>

// ROOT-associated
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"

// own files
#include "slimSegsData.h"
#include "plots.h"
#include "rootData.h"
#include "utilityFunctions.h"
#include "SegsRepresentation.h"

using std::cout;
using std::endl;


//Declaration for function defined in "../3rdparty/karimaki/karimaki_hit.c"
int karimaki_hit(KariFit&, int , double *, double *, double *, double *, double*, double *, double *, double *);

void slimSegsData(/*std::string outputfile, int run, const bool appendToFile */) {

    int run = 10;
    const bool appendToFile = false;
    const std::string outputfile = "../data/SlimmedData/mu3e_test_slimmed_file_000000.root";

    const std::string pathtodata = "../data/SimulationData/";
    const std::string pathtoplots = "../plots/Mu3eSlimSegs/";
    const std::string pathtorunplots = pathtoplots + "run_" + get_padded_string(run, 3, '0') + "/";
    const std::string infile = pathtodata + "mu3e_run_" + get_padded_string(run, 6, '0'); + "_trirec_cosmic.root";

    const bool RECONSTRUCTION_PRINTS = false;
    const bool HIT_PRINTS = DEBUG;
    const bool MAKE_PLOT = true;
    const int MAX_ENTRIES = 0;

    std::string filtertag;

    check_create_directory(pathtodata);
    check_create_directory(pathtoplots);
    check_create_directory(pathtorunplots);

    TFile tF(outputfile.c_str(), "recreate");
    if (!tF.IsOpen()) {
        std::cout << "[ERROR] File " << tF.GetName() << " is not open!" << std::endl;
        exit(0);
    }

    TTree t_slim("SlimSegs","Tree with slimmed simulation data");

    //// FILE FOR READING

    TFile f2(infile.c_str());
    TTree *t_segs;
    f2.GetObject("segs", t_segs);
    TTree *t_frames;
    f2.GetObject("frames", t_frames);

    ////################# DATA FOR ROOT FILES ################################

    SegsRepresentationAndCalc Segs = SegsRepresentationAndCalc(t_segs);

//    //data for trirec result in "segs" tree
//    unsigned int segs_entries = t_segs->GetEntries();
//    int rec_event;
//    int rec_nhit;
//    int rec_ntriplet;
//    int rec_trajid;
//
//    float mc_p;
//    float mc_pt;
//    float mc_theta;
//    float mc_phi;
//    float mc_lam;
//    int mc_type;
//    int mc_pid;
//
//    float rec_p;
//    float rec_r;
//    float rec_rt;
//    float rec_tan01[TRIPLET_HIT_ARRAY_LENGTH];
//    float rec_tan12[TRIPLET_HIT_ARRAY_LENGTH];
//    float rec_lam01[TRIPLET_HIT_ARRAY_LENGTH];
//    float rec_lam12[TRIPLET_HIT_ARRAY_LENGTH];
//    float rec_zpca_x;
//    float rec_zpca_y;
//    float rec_zpca_z;
//    float rec_zpca_r;
//
//    float x00[TRIPLET_HIT_ARRAY_LENGTH];
//    float x10[TRIPLET_HIT_ARRAY_LENGTH];
//    float x20[TRIPLET_HIT_ARRAY_LENGTH];
//    float x01[TRIPLET_HIT_ARRAY_LENGTH];
//    float x11[TRIPLET_HIT_ARRAY_LENGTH];
//    float x21[TRIPLET_HIT_ARRAY_LENGTH];
//
//    float y00[TRIPLET_HIT_ARRAY_LENGTH];
//    float y10[TRIPLET_HIT_ARRAY_LENGTH];
//    float y20[TRIPLET_HIT_ARRAY_LENGTH];
//    float y01[TRIPLET_HIT_ARRAY_LENGTH];
//    float y11[TRIPLET_HIT_ARRAY_LENGTH];
//    float y21[TRIPLET_HIT_ARRAY_LENGTH];
//
//    float z00[TRIPLET_HIT_ARRAY_LENGTH];
//    float z10[TRIPLET_HIT_ARRAY_LENGTH];
//    float z20[TRIPLET_HIT_ARRAY_LENGTH];
//    float z01[TRIPLET_HIT_ARRAY_LENGTH];
//    float z11[TRIPLET_HIT_ARRAY_LENGTH];
//    float z21[TRIPLET_HIT_ARRAY_LENGTH];
//
//    float sid00[TRIPLET_HIT_ARRAY_LENGTH];
//    float sid10[TRIPLET_HIT_ARRAY_LENGTH];
//    float sid20[TRIPLET_HIT_ARRAY_LENGTH];
//    float sid01[TRIPLET_HIT_ARRAY_LENGTH];
//    float sid11[TRIPLET_HIT_ARRAY_LENGTH];
//    float sid21[TRIPLET_HIT_ARRAY_LENGTH];
//
//    t_segs->SetBranchAddress("eventId", &rec_event);
//    t_segs->SetBranchAddress("nhit", &rec_nhit);
//    t_segs->SetBranchAddress("n", &rec_ntriplet);
//
//    t_segs->SetBranchAddress("p", &rec_p);
//    t_segs->SetBranchAddress("r", &rec_r);
//    t_segs->SetBranchAddress("rt", &rec_rt);
//    t_segs->SetBranchAddress("tan01", &rec_tan01);
//    t_segs->SetBranchAddress("tan12", &rec_tan12);
//    t_segs->SetBranchAddress("lam01", &rec_lam01);
//    t_segs->SetBranchAddress("lam12", &rec_lam12);
//    t_segs->SetBranchAddress("zpca_z", &rec_zpca_z);
//    t_segs->SetBranchAddress("zpca_x", &rec_zpca_x);
//    t_segs->SetBranchAddress("zpca_y", &rec_zpca_y);
//    t_segs->SetBranchAddress("zpca_r", &rec_zpca_r);
//
//    t_segs->SetBranchAddress("mc_tid", &rec_trajid);
//    t_segs->SetBranchAddress("mc_p", &mc_p);
//    t_segs->SetBranchAddress("mc_pt", &mc_pt);
//    t_segs->SetBranchAddress("mc_theta", &mc_theta);
//    t_segs->SetBranchAddress("mc_phi", &mc_phi);
//    t_segs->SetBranchAddress("mc_lam", &mc_lam);
//    t_segs->SetBranchAddress("mc_type", &mc_type);
//    t_segs->SetBranchAddress("mc_pid", &mc_pid);
//
//    t_segs->SetBranchAddress("x00", &x00);
//    t_segs->SetBranchAddress("x10", &x10);
//    t_segs->SetBranchAddress("x20", &x20);
//    t_segs->SetBranchAddress("x01", &x01);
//    t_segs->SetBranchAddress("x11", &x11);
//    t_segs->SetBranchAddress("x21", &x21);
//
//    t_segs->SetBranchAddress("y00", &y00);
//    t_segs->SetBranchAddress("y10", &y10);
//    t_segs->SetBranchAddress("y20", &y20);
//    t_segs->SetBranchAddress("y01", &y01);
//    t_segs->SetBranchAddress("y11", &y11);
//    t_segs->SetBranchAddress("y21", &y21);
//
//    t_segs->SetBranchAddress("z00", &z00);
//    t_segs->SetBranchAddress("z10", &z10);
//    t_segs->SetBranchAddress("z20", &z20);
//    t_segs->SetBranchAddress("z01", &z01);
//    t_segs->SetBranchAddress("z11", &z11);
//    t_segs->SetBranchAddress("z21", &z21);
//
//    t_segs->SetBranchAddress("sid00", &sid00);
//    t_segs->SetBranchAddress("sid10", &sid10);
//    t_segs->SetBranchAddress("sid20", &sid20);
//    t_segs->SetBranchAddress("sid01", &sid01);
//    t_segs->SetBranchAddress("sid11", &sid11);
//    t_segs->SetBranchAddress("sid21", &sid21);
//
//    cout << "Branches set for segs..." << endl;



    //// further data for stats graphs and calculations
    int p_fail_count = 0;
    int processed_entries = 0;

    std::vector<float> p_rel_errors;
    std::vector<float> p_inv_rel_errors;
    std::vector<float> p_inv_kari_errors;
    std::vector<float> pt_inv_kvsms_errors;
    std::vector<float> pt_inv_errors;
    std::vector<float> p_over_pmcs;
    std::vector<float> rec_rs;
    std::vector<float> rec_rts;
    std::vector<float> rec_p_corrs;
    std::vector<float> rec_ps;
    std::vector<float> rec_inv_ps;
    std::vector<float> rec_inv_pts;
    std::vector<float> rec_pts;
    std::vector<float> rec_pt_corrs;
    std::vector<float> rec_phis;
    std::vector<float> rec_thetas;
    std::vector<float> rec_dca_rs;
    std::vector<float> rec_dca_xs;
    std::vector<float> rec_dca_ys;
    std::vector<float> rec_dca_zs;

    std::vector<float> kari_pts;
    std::vector<float> kari_inv_pts;
    std::vector<float> kari_r3ds;
    std::vector<float> kari_rads;
    std::vector<float> kari_inv_rads;
    std::vector<float> kari_dcas;
    std::vector<float> kari_phis;
    std::vector<float> kari_tchi2ns;
    std::vector<float> kari_z0s;
    std::vector<float> kari_thetas;
    std::vector<float> kari_zchi2ns;

    std::vector<float> mc_dcas;
    std::vector<float> mc_z_dcas;
    std::vector<float> mc_phi_dcas;
    std::vector<float> mc_ps;
    std::vector<float> mc_p_corrs;
    std::vector<float> mc_inv_ps;
    std::vector<float> mc_inv_pts;
    std::vector<float> mc_pt_corrs;
    std::vector<float> mc_pts;
    std::vector<float> mc_phis;
    std::vector<int> sim_nhits;
    std::vector<int> rec_nhits;
    std::vector<int> mc_types;
    std::vector<float> rec_nhits_float;

    //data for karifit
    KariFit karires;


    double RMS = 80 * 1e-3 / sqrt(12); //RMS in mm
    const float BFIELD = 1.0;

    double tres[TRIPLET_HIT_ARRAY_LENGTH] = {RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS};
    double zres[TRIPLET_HIT_ARRAY_LENGTH] = {RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS};
    double rres[TRIPLET_HIT_ARRAY_LENGTH] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};



    for (unsigned int i = 0; i < (MAX_ENTRIES == 0 ? Segs.segs_entries : MAX_ENTRIES); i++) {
        Segs.t_segs->GetEntry(i);
        Segs.calcAdditionalData();

        #if DEBUG
        printf("\nSEGS ENTRY\n------------------------------\n\n");
        #endif

        //Theta, phi and traverse p of reconstruction
        float rec_pt = Segs.rec_p * std::cos((Segs.rec_lam01)[0]);
        float rec_phi = (Segs.rec_tan01)[0];
        float rec_theta = PI*0.5 - (Segs.rec_lam01)[0];

        //Do Karimaki Helix fit
        std::vector<double> xp;
        std::vector<double> yp;
        std::vector<double> zp;
        std::vector<int> layerp;
        std::vector<double> phi_hitp;
        std::vector<double> thetap;

        int ncombinedhits = combineBasicHits(xp, yp, zp, phi_hitp, thetap, layerp, Segs.rec_nhit, Segs.rec_ntriplet,
                                             Segs.x00, Segs.x20, Segs.y00, Segs.y20, Segs.z00,Segs.z20,Segs.sid00,
                                             Segs.sid20,Segs.rec_tan01,Segs.rec_tan12,Segs.rec_lam01,Segs.rec_lam12);

//        getResolutionArrs(ncombinedhits, RMS, tres, zres, rres, xp, yp, phi_hitp, (double) rec_theta);
//        assert(tres.size() == ncombinedhits);

#if DEBUG
        //print triplet data from segs tree
        cout << "\tsegs tree data:" << endl;

        for(int j=0; j<rec_ntriplet; j++) {
            cout << "\ttriplet index" << j << endl;
            cout << "\ttriplet 00: " << "\tsid=" << sid00[j] << "\tx00=" << x00[j] << "\ty00=" << y00[j]<< endl;
            cout << "\ttriplet 10: " << "\tsid=" << sid10[j] << "\tx10=" << x10[j] << "\ty10=" << y10[j] << endl;
            cout << "\ttriplet 20: " << "\tsid=" << sid20[j] << "\tx20=" << x20[j] << "\ty20=" << y20[j] << endl;
            cout << "\ttriplet 01: " << "\tsid=" << sid01[j] << "\tx01=" << x01[j] << "\ty01=" << y01[j] << endl;
            cout << "\ttriplet 11: " << "\tsid=" << sid11[j] << "\tx11=" << x11[j] << "\ty11=" << y11[j] << endl;
            cout << "\ttriplet 21: " << "\tsid=" << sid21[j] << "\tx21=" << x21[j] << "\ty21=" << y21[j] << endl;
            cout << "\tangles    : " << "\ttan01=" << rec_tan01[j] << "\ttan12=" << rec_tan12[j] <<  endl;

        }
        cout << endl;

        for(int i = 0; i < ncombinedhits; i++) {
            cout << "\tHits sorted for kari: x=" << xp[i] << "  y=" << yp[i] << "  z=" << zp[i] <<" phi=" << phi_hits[i] << " theta=" << thetas[i] << endl;
        }
#endif

        karimaki_hit(karires, ncombinedhits, &xp[0], &yp[0], &zp[0], &phi_hitp[0], &thetap[0], &tres[0], &zres[0], &rres[0]);
        correctKariDirection(karires);

        if (Segs.mc_p == 0 || Segs.mc_pt == 0 || Segs.rec_p == 0 || rec_pt == 0) {
            p_fail_count++;
        } else {
            float mc_p_corr;
            float mc_pt_corr;

            if (FILTER == 0 || FILTER == 6) {
                mc_p_corr = (Segs.mc_p) * (Segs.mc_type % 10 == 3 ? -1 : 1);
                mc_pt_corr = (Segs.mc_pt) * (Segs.mc_type % 10 == 3 ? -1 : 1);
            } else {
                mc_p_corr = (Segs.mc_p) * sgn(Segs.mc_pid);
                mc_pt_corr = (Segs.mc_pt) * sgn(Segs.mc_pid);
            }

            float rec_p_corr = Segs.rec_p * sgn(Segs.rec_r);
            float rec_pt_corr = rec_pt;
            float kari_pt = 0.3 * karires.rad * BFIELD;

            float mc_inv_p = 1. / mc_p_corr;
            float mc_inv_pt = 1. / mc_pt_corr;

            float rec_inv_p = 1. / Segs.rec_p;
            float rec_inv_pt = 1. / rec_pt;

            float kari_inv_pt = 1. / kari_pt;

            float p_inv_abs_error = (rec_inv_p - mc_inv_p); //this will mostly be used as estimator for deviation
            float p_abs_error = (mc_p_corr - Segs.rec_p);
            float p_over_pmc = Segs.rec_p / mc_p_corr;
            float pt_inv_error = (1./rec_pt_corr) - (1./mc_pt_corr);
            float pt_kari_inv_err = kari_inv_pt - (1./mc_pt_corr);
            float pt_kvsms_inv_err = kari_inv_pt - rec_inv_pt;


            int choice = true;
            if (choice) {

//                printf("mc_type? %d \t mc_pid = %d\n", mc_type, mc_pid);
//                assert(!(mc_type == 4 && sgn(mc_pid) == -1));
//                assert(!(mc_type == 3 && sgn(mc_pid) == 1));

                //calculated data
                rec_inv_ps.push_back(rec_inv_p);
                rec_p_corrs.push_back(rec_p_corr);
                rec_inv_pts.push_back(rec_inv_p);
                p_inv_rel_errors.push_back(p_inv_abs_error);
                p_rel_errors.push_back(p_abs_error);
                pt_inv_errors.push_back(pt_inv_error);
                p_over_pmcs.push_back(p_over_pmc);
                rec_pts.push_back(rec_pt);
                rec_phis.push_back(rec_phi);
                rec_thetas.push_back(rec_theta);

                mc_p_corrs.push_back(mc_p_corr);
                mc_pt_corrs.push_back(mc_pt_corr);
                mc_inv_pts.push_back(mc_inv_pt);
                mc_inv_ps.push_back(mc_inv_p);

                //reconstruction data
                rec_ps.push_back(Segs.rec_p);
                rec_rs.push_back(Segs.rec_r);
                rec_rts.push_back(Segs.rec_rt);
                rec_dca_rs.push_back(Segs.rec_zpca_r);
                rec_dca_zs.push_back(Segs.rec_zpca_z);

                //kari data
                kari_r3ds.push_back(karires.r3d);
                kari_rads.push_back(karires.rad);
                kari_dcas.push_back(karires.dca);
                kari_phis.push_back(karires.phi);
                kari_tchi2ns.push_back(karires.tchi2n);
                kari_z0s.push_back(karires.z0);
                kari_thetas.push_back(karires.theta);
                kari_zchi2ns.push_back(karires.zchi2n);

                kari_pts.push_back(kari_pt);
                kari_inv_pts.push_back(kari_inv_pt);
                p_inv_kari_errors.push_back(pt_kari_inv_err);
                pt_inv_kvsms_errors.push_back(pt_kvsms_inv_err);

                //monte carlo data
                mc_ps.push_back(Segs.mc_p);
                mc_pts.push_back(Segs.mc_pt);
                mc_phis.push_back(Segs.mc_phi);
                mc_types.push_back(Segs.mc_type);

                //meta data
                rec_nhits.push_back(Segs.rec_nhit);

                processed_entries++;

                ////PRINT SECTION PER ENTRY IN TREE
                if(RECONSTRUCTION_PRINTS) {
                    //data from reconstruction
                    cout << "rec_event: " << Segs.rec_event;

                    cout << " mc_p: " << Segs.mc_p << " rec_p: " << Segs.rec_p*sgn(Segs.rec_r) << "rec_nhit: "<< Segs.rec_nhit <<"\t\t";
                    cout << "rec_r= " << Segs.rec_r << " mc_type=" << Segs.mc_type;

                    //mc and rec deviations
                    cout << "pt_mc= " << Segs.mc_pt << " pt_rec= " << rec_pt << " phi_mc= " << Segs.mc_phi << " rec_phi= " << rec_phi;
                    cout << " mc_theta= " << Segs.mc_theta << " rec_theta= " << rec_theta << "z-dca=" << Segs.mc_z_dca <<"\t\t";

                    cout << "(1/mc_p-1/rec_p)/(1/mc_p): " << p_inv_abs_error*100 << " %" << endl;
                }
            }
        }

    }


}


//
// Created by Konstantin Neureither on 25.06.20.
//

#define MAKE_PLOT false

// Basic imports
#include <iostream>
#include <string>
#include <stdlib.h>

// ROOT-associated
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"

// own files
#include "slimSegsDataScript.h"
#include "plots.h"
#include "rootData.h"
#include "utilityFunctions.h"
#include "SegsTreeRead.h"
#include "SlimSegsTree.h"
#include "basicDefines.h"

using std::cout;
using std::endl;

void slimSegsDataScript(std::string outputfilename, const int run, const bool appendToFile) {

    const int outnum = 0;
    const int MAX_ENTRIES = 0;

    const std::string pathtoplots = "plots/Mu3eSlimSegs/";
    const std::string pathtodata = "data/SimulationData/";
    const std::string pathtorunplots = pathtoplots + "run_" + get_padded_string(run, 3, '0') + "/";
    const std::string infile = pathtodata + "mu3e_run_" + get_padded_string(run, 6, '0') + "_trirec_cosmic.root";
//    const std::string outputfile = pathtodata + "mu3e_slimmed_segs_" + get_padded_string(outnum, 6, '0') + ".root";
    const std::string outputfile = "data/SlimmedData/" + outputfilename;

    check_create_directory(pathtodata);
    check_create_directory(pathtoplots);
    check_create_directory(pathtorunplots);


    // FILE FOR WRITING
    TFile toutF(outputfile.c_str(), (appendToFile ? "update" : "recreate"));
    if (!toutF.IsOpen()) {
        std::cout << "[ERROR] File " << toutF.GetName() << " is not open!" << std::endl;
        exit(0);
    }

    TTree t_slim("SlimSegs","Tree with slimmed simulation data");

    // FILE FOR READING
    TFile tinF(infile.c_str());
    if (!tinF.IsOpen()) {
        std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
        exit(0);
    }

    TTree *t_segs;
    tinF.GetObject("segs", t_segs);
    TTree *t_frames;
    tinF.GetObject("frames", t_frames);


    ////################# DATA FOR ROOT FILES ################################

    //class representation for segs tree and read functionality
    SegsTreeReadPlus Segs = SegsTreeReadPlus(t_segs);
    //class representation for slimmed down 'slimSegs' Tree and Write functionality
    SlimSegsTreeWrite SlimSegs = SlimSegsTreeWrite(&t_slim);

    //data for karifit
    KariFitCalc karires;

    //meta data
    SlimSegsMeta meta;

    //// further data for stats graphs and calculations
    int p_fail_count = 0;
    int processed_entries = 0;

#if MAKE_PLOT
        //vectors for making control plots
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
#endif

    //detector specific resolution data for helix fitting in karimakiHelixfit.cpp
    double RMS = 80 * 1e-3 / sqrt(12); //RMS in mm
    const float BFIELD = 1.0;
    double tres[TRIPLET_HIT_ARRAY_LENGTH] = {RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS};
    double zres[TRIPLET_HIT_ARRAY_LENGTH] = {RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS};
    double rres[TRIPLET_HIT_ARRAY_LENGTH] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

    //iterate through entries in segs tree
    for (unsigned int i = 0; i < (MAX_ENTRIES == 0 ? Segs.my_entries : MAX_ENTRIES); i++) {
        printf("SEGS ENTRY %d ------- \n", i);

        //update data in Segs class
        Segs.getEntry(i);

        //set meta data
        meta.uEventID = run + 1000*Segs.rec_event;
        meta.segsIndex = i;
        meta.runID = run;

        //only look at entries where the reconstruction did not fail
        if (Segs.mc_p == 0 || Segs.mc_pt == 0 || Segs.rec_p == 0 || Segs.rec_pt == 0) {
            p_fail_count++;
        } else {
            Segs.calcAdditionalData();

            ////Do Karimaki Helix fit
            //define variables that store sorted hits
            std::vector<double> xp;
            std::vector<double> yp;
            std::vector<double> zp;
            std::vector<int> layerp;
            std::vector<double> phi_hitp;
            std::vector<double> thetap;

            //gather the hits for particle from fit triplets, sorted as vectors
            int ncombinedhits = combineBasicHits(xp, yp, zp, phi_hitp, thetap, layerp, Segs.rec_nhit, Segs.rec_ntriplet,
                                                 Segs.x00, Segs.x20, Segs.y00, Segs.y20, Segs.z00,Segs.z20,Segs.sid00,
                                                 Segs.sid20,Segs.rec_tan01,Segs.rec_tan12,Segs.rec_lam01,Segs.rec_lam12);

            //calculate momentum from radius
            karires.p = 0.3 * karires.r3d * BFIELD;
            karires.pt = 0.3 * karires.rad * BFIELD;

/*some debug prints*/
#if DEBUG
            //print triplet data from segs tree
            cout << "\tsegs tree data:" << endl;

            for(int j=0; j<Segs.rec_ntriplet; j++) {
                cout << "\ttriplet index" << j << endl;
                cout << "\ttriplet 00: " << "\tsid=" << Segs.sid00[j] << "\tx00=" << Segs.x00[j] << "\ty00=" << Segs.y00[j]<< endl;
                cout << "\ttriplet 10: " << "\tsid=" << Segs.sid10[j] << "\tx10=" << Segs.x10[j] << "\ty10=" << Segs.y10[j] << endl;
                cout << "\ttriplet 20: " << "\tsid=" << Segs.sid20[j] << "\tx20=" << Segs.x20[j] << "\ty20=" << Segs.y20[j] << endl;
                cout << "\ttriplet 01: " << "\tsid=" << Segs.sid01[j] << "\tx01=" << Segs.x01[j] << "\ty01=" << Segs.y01[j] << endl;
                cout << "\ttriplet 11: " << "\tsid=" << Segs.sid11[j] << "\tx11=" << Segs.x11[j] << "\ty11=" << Segs.y11[j] << endl;
                cout << "\ttriplet 21: " << "\tsid=" << Segs.sid21[j] << "\tx21=" << Segs.x21[j] << "\ty21=" << Segs.y21[j] << endl;
                cout << "\tangles    : " << "\ttan01=" << Segs.rec_tan01[j] << "\ttan12=" << Segs.rec_tan12[j] <<  endl;

            }
            cout << endl;

            for(int i = 0; i < ncombinedhits; i++) {
                cout << "\tHits sorted for kari: x=" << xp[i] << "  y=" << yp[i] << "  z=" << zp[i] <<" phi=" << phi_hitp[i] << " theta=" << thetap[i] << " layer=" << layerp[i] << endl;
            }
#endif


            //do the helix fit
            karimakiHelixfit(karires, ncombinedhits, &xp[0], &yp[0], &zp[0], &phi_hitp[0], &thetap[0], &tres[0],
                             &zres[0], &rres[0]);
            //this corrects the direction according to phi, as the cosmic come from above
            correctKariDirection(karires);

            int choice = Segs.mc_type == 3 || Segs.mc_type == 4;
            if (choice) {
                processed_entries++;

                //write all gathered and calculated data to the SlimSegs tree
                SlimSegs.fillData(Segs, meta, karires, ncombinedhits,
                                  xp, zp, yp, layerp);


                float kari_inv_pt = 1. / karires.pt;
                float kari_inv_rad = 1. / karires.rad;

                float pt_kari_inv_err = kari_inv_pt - (1./Segs.mc_pt_inv_corr);
                float pt_kvsms_inv_err = kari_inv_pt - Segs.rec_p_inv;


                //Check some data if Segs->SlimSegs transfer worked
                //This is not fail safe, just a check of the principle
                assert(!(Segs.mc_type == 4 && sgn(Segs.mc_pid) == -1));
                assert(!(Segs.mc_type == 3 && sgn(Segs.mc_pid) == 1));

                assert(SlimSegs.eventID == Segs.rec_event);
                assert(SlimSegs.rec_p == Segs.rec_p);
                assert(SlimSegs.rec_ntriplet == Segs.rec_ntriplet);
//                for(int i=0; i<Segs.rec_ntriplet; i++) {
//                    assert(SlimSegs.x00[i] == Segs.x00[i]);
//                    assert(SlimSegs.x10[i] == Segs.x10[i]);
//                }

#if MAKE_PLOT
                    //calculated data
//                rec_inv_ps.push_back(Segs.rec_inv_p);
//                rec_p_corrs.push_back(rec_p_corr);
//                rec_inv_pts.push_back(rec_inv_p);
//                p_inv_rel_errors.push_back(p_inv_abs_error);
//                p_rel_errors.push_back(p_abs_error);
//                pt_inv_errors.push_back(pt_inv_error);
//                p_over_pmcs.push_back(p_over_pmc);
//                rec_pts.push_back(rec_pt);
//                rec_phis.push_back(rec_phi);
//                rec_thetas.push_back(rec_theta);
//
//                mc_p_corrs.push_back(mc_p_corr);
//                mc_pt_corrs.push_back(mc_pt_corr);
//                mc_inv_pts.push_back(mc_inv_pt);
//                mc_inv_ps.push_back(mc_inv_p);

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

                    kari_pts.push_back(karires.pt);
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
#endif
                }
            }
        }

    SlimSegs.t_slimSegs->Print();
    toutF.Write();
    toutF.Close();
    tinF.Close();

    std::cout << "\n\n>>>>> GENERAL STATS <<<<<\n\n";
    std::cout << " - P Fail count was " << p_fail_count << " of " << Segs.my_entries << " entries in total. (";
    std::cout << p_fail_count / (float)Segs.my_entries * 100 << " %)" << endl;
    std::cout << " - total entries processed: " << processed_entries << endl;


    //// Kontrollplots

}


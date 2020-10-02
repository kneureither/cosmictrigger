//
// Created by Konstantin Neureither on 25.06.20.
//

#define MAKE_PLOT true

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

void slimSegsDataScript(const int dataset, const int run, const bool appendToFile) {

    const int MAX_ENTRIES = 0;

    const std::string pathtoplots = "output/Mu3eSlimSegs/";
    const std::string pathtodata = "data/SimulationData/";
    const std::string pathtodatasetplots = pathtoplots + "dataset_" + get_padded_string(dataset, 3, '0') + "/";
    const std::string infile = pathtodata + "mu3e_run_" + get_padded_string(run, 6, '0') + "_trirec_cosmic.root";
//    const std::string outputfile = pathtodata + "mu3e_slimmed_segs_" + get_padded_string(dataset, 6, '0') + ".root";
    const std::string outputfile = "data/SlimmedData/mu3e_slimmed_segs_" + get_padded_string(dataset, 6, '0') + ".root";

    check_create_directory(pathtodata);
    check_create_directory(pathtoplots);
    check_create_directory(pathtodatasetplots);

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
        std::vector<float> p_inv_rel_errors;
        std::vector<float> p_inv_kari_errors;
        std::vector<float> pt_inv_kvsms_errors;
        std::vector<float> pt_inv_errors;
        std::vector<float> p_over_pmcs;
        std::vector<float> rec_p_corrs;
        std::vector<float> rec_ps;
        std::vector<float> rec_rs;
        std::vector<float> rec_inv_ps;
        std::vector<float> rec_inv_pts;
        std::vector<float> rec_pts;
        std::vector<float> rec_pt_corrs;
        std::vector<float> rec_phis;
        std::vector<float> rec_thetas;
        std::vector<float> rec_dca_rs;
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
        std::vector<int> rec_nhits;
        std::vector<int> mc_types;
#endif

    //detector specific resolution data for helix fitting in karimakiHelixfit.cpp
    double RMS = 80 * 1e-3 / sqrt(12); //RMS in mm
    const float BFIELD = 1.0;
    double tres[TRIPLET_HIT_ARRAY_LENGTH] = {RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS};
    double zres[TRIPLET_HIT_ARRAY_LENGTH] = {RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS};
    double rres[TRIPLET_HIT_ARRAY_LENGTH] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

    //iterate through entries in segs tree
    int max_entries = (MAX_ENTRIES == 0 ? Segs.my_entries : MAX_ENTRIES);

    for (unsigned int i = 0; i < max_entries; i++) {
//        printf("SEGS ENTRY %d ------- \n", i);

//        if(((i % (int) pow((float) 10, (float) std::floor(log10(i))) == 0) && (i >= 1000))
//           || (i >= 1000000 && i % 100000 == 0)) {
//
//        }

        if(i % 1000 == 0) {
            print_status_bar(i, max_entries, "processing run " + get_string(run), "");
//            std::cout << "\r(STATUS) : processing entry " << i << " of " << Segs.my_entries << std::flush;
        }

        //update data in Segs class
        Segs.getEntry(i);
        Segs.calcAdditionalData();

        //set meta data
        meta.uEventID = run + 1000*Segs.rec_event;
        meta.segsIndex = i;
        meta.runID = run;

        //only look at entries where the reconstruction did not fail
        if (Segs.mc_p == 0 || Segs.mc_pt == 0 || Segs.rec_p == 0 || Segs.rec_pt == 0) {
            p_fail_count++;
        } else {

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
            //according to bfield direction
            swapKariBField(karires);

            //calculate momentum from radius
            karires.p = 0.3 * karires.r3d * BFIELD;
            karires.pt = 0.3 * karires.rad * BFIELD;

            int choice = Segs.mc_type == 3 || Segs.mc_type == 4;
            if (choice) {
                processed_entries++;

                //write all gathered and calculated data to the SlimSegs tree
                SlimSegs.fillData(Segs, meta, karires, ncombinedhits,
                                  xp, zp, yp, layerp);


                float kari_inv_pt = 1. / karires.pt;
                float kari_inv_rad = 1. / karires.rad;

                float pt_kari_inv_err = kari_inv_pt - Segs.mc_pt_inv_corr;
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
                    rec_inv_ps.push_back(Segs.rec_p_inv);
                    rec_inv_pts.push_back(Segs.rec_pt_inv);
                    p_inv_rel_errors.push_back(Segs.p_inv_abs_error);
                    pt_inv_errors.push_back(Segs.pt_inv_abs_error);
                    rec_pts.push_back(Segs.rec_pt);
                    rec_phis.push_back(Segs.rec_phi);
                    rec_thetas.push_back(Segs.rec_theta);

                    mc_p_corrs.push_back(Segs.mc_p_corr);
                    mc_pt_corrs.push_back(Segs.mc_pt_corr);
                    mc_inv_pts.push_back(Segs.mc_pt_inv_corr);
                    mc_inv_ps.push_back(Segs.mc_p_inv_corr);

                    //reconstruction data
                    rec_ps.push_back(Segs.rec_p);
                    rec_rs.push_back(Segs.rec_r);
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
                    kari_inv_rads.push_back(kari_inv_rad);
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

#if MAKE_PLOT
    const float LEFT_BOUNDARY_P = -3e4;
    const float RIGHT_BOUNDARY_P = 3e4 ;
    const int BIN_COUNT = 50;

    const std::string plottingfile = pathtodatasetplots + "Mu3eSlimSegsScript_run_" + get_padded_string(run, 6, '0') +
                                     "_dataset_" + get_padded_string(dataset, 3, '0') + ".pdf";


    TH1F *h_rec_nhits = new TH1F("h_rec_nhits", "hits per frame of trirec output", 20, 0., 20.);
    labelAxis(h_rec_nhits, "number of hits", "count");
    fillHistWithVector(h_rec_nhits, rec_nhits);

    //Trajectory data Monte Carlo
    TH1F *h_pmc = new TH1F("h_pmc", "muon momenta monte carlo (charge corrrected)", 30, LEFT_BOUNDARY_P, RIGHT_BOUNDARY_P);
    labelAxis(h_pmc, "p_{mc} * charge_{mc} [MeV]", "count");
    fillHistWithVector(h_pmc, mc_p_corrs);

    TH1F *h_ptmc = new TH1F("h_ptmc", "muon traverse momenta monte carlo (charge corrrected)", 30, LEFT_BOUNDARY_P, RIGHT_BOUNDARY_P);
    labelAxis(h_ptmc, "p_{t-mc} * charge_{mc} [MeV]", "count");
    fillHistWithVector(h_ptmc, mc_pt_corrs);

    TH1F *h_phimc = new TH1F("h_phimc", "#Phi monte carlo (charge corrrected)", 30, -3.2, 1);
    labelAxis(h_phimc, "#Phi", "count");
    fillHistWithVector(h_phimc, mc_phis);

    TH1F *h_mctype = new TH1F("h_mctype", "Particle types of simulation", 55, 0, 55);
    labelAxis(h_mctype, "particle type", "count");
    fillHistWithVector(h_mctype, mc_types);

    //Trajectory data reconstructed
    TH1F *h_p = new TH1F("h_p", "muon momenta", 30, LEFT_BOUNDARY_P, RIGHT_BOUNDARY_P);
    labelAxis(h_p, "p_{rec} [MeV]", "count");
    fillHistWithVector(h_p, rec_ps);

    TH1F *h_p_corr = new TH1F("h_p_corr", "muon momenta (charge corrected)", 30, 0, RIGHT_BOUNDARY_P);
    labelAxis(h_p_corr, "p_{rec-corr} [MeV^{-1}]", "count");
    fillHistWithVector(h_p_corr, rec_p_corrs);

    TH1F *h_pt = new TH1F("h_pt", "muon transverse momentum", 30, LEFT_BOUNDARY_P, RIGHT_BOUNDARY_P);
    labelAxis(h_pt, "p_{t} [MeV]", "count");
    fillHistWithVector(h_pt, rec_pts);

    TH1F *h_phi = new TH1F("h_phi", "reconstruction #phi (var: tan01)", 30, -3.2, 3.2);
    labelAxis(h_phi, "#phi", "count");
    fillHistWithVector(h_phi, rec_phis);

    TH1F *h_theta = new TH1F("h_theta", "reconstruction #Theta_{msfit}", 30, -3.2, 3.2 );
    labelAxis(h_theta, "#Theta_{msfit}", "count");
    fillHistWithVector(h_theta, rec_thetas);

    //rec pt
    TH1F *h_pt_inv_err = new TH1F("h_pt_inv_err", "pt_{rec}^{-1} #minus pt_{mc}^{-1} transverse momentum error", 30, -5e-4,5e-4);
    labelAxis(h_pt_inv_err, "pt_{rec}^{-1} #minus pt_{mc}^{-1} [MeV^{-1}]", "count");
    fillHistWithVector(h_pt_inv_err, pt_inv_errors);

    TH1F *h_invptrec = new TH1F("h_invptrec", "reconstruction pt_{rec}^{-1}", 30, -5e-4,5e-4);
    labelAxis(h_invptrec, "pt_{rec}^{-1} [MeV^{-1}]", "count");
    fillHistWithVector(h_invptrec, rec_inv_pts);

    TH1F *h_invptmc = new TH1F("h_invptmc", "monte carlo pt_{mc}^{-1}", 30, -5e-4,5e-4);
    labelAxis(h_invptmc, "pt_{mc}^{-1} [MeV^{-1}]", "count");
    fillHistWithVector(h_invptmc, mc_inv_pts);


    //DCAs reconstruction
    TH1F *h_zdca = new TH1F("h_zdca", "dca_{reconstruction} along z-axis", 30, -600, 600);
    labelAxis(h_zdca, "dca_{z} [mm]", "count");
    fillHistWithVector(h_zdca, rec_dca_zs);

    TH1F *h_rdca = new TH1F("h_rdca", "dca_{reconstruction} r", 30, -100, 100);
    labelAxis(h_rdca, "dca_{r} [mm]", "count");
    fillHistWithVector(h_rdca, rec_dca_rs);

    //KARIMAKI Results
    TH1F *h_r3dkari = new TH1F("h_r3dkari", "rad_{kari} 3D ", 30, -1e5, 1e5);
    labelAxis(h_r3dkari, "3D radius [mm]", "count");
    fillHistWithVector(h_r3dkari, kari_r3ds);

    TH1F *h_rinvkari = new TH1F("h_rinvkari", "rad_{kari}^{-1} 2D ", 30, -0.0003, 0.0003);
    labelAxis(h_rinvkari, "rt^{-1} [mm^{-1}]", "count");
    fillHistWithVector(h_rinvkari, kari_inv_rads);

    TH1F *h_phikari = new TH1F("h_phikari", "#phi_{kari}", 30, -3.2, 3.2);
    labelAxis(h_phikari, "#phi", "count");
    fillHistWithVector(h_phikari, kari_phis);

    TH1F *h_dcakari = new TH1F("h_dcakari", "DCA_{kari}", 30, -100, 100);
    labelAxis(h_dcakari, "DCA [mm]", "count");
    fillHistWithVector(h_dcakari, kari_dcas);

    TH1F *h_tchi2nkari = new TH1F("h_tchi2nkari", "transverse #chi^{2}_{kari}", 30, 0, 5);
    labelAxis(h_tchi2nkari, "#chi^{2}", "count");
    fillHistWithVector(h_tchi2nkari, kari_tchi2ns);

    TH1F *h_z0kari = new TH1F("h_z0kari", "z0_{kari}", 30, -600, 600);
    labelAxis(h_z0kari, "z0 [mm]", "count");
    fillHistWithVector(h_z0kari, kari_z0s);

    TH1F *h_thetakari = new TH1F("h_thetakari", "#Theta_{kari}", 30, -3.2, 3.2);
    labelAxis(h_thetakari, "#Theta", "count");
    fillHistWithVector(h_thetakari, kari_thetas);

    TH1F *h_zchi2kari = new TH1F("h_zchi2kari", "longitudinal #chi^{2}_{kari}", 30,0, 5);
    labelAxis(h_zchi2kari, "#chi^{2}", "count");
    fillHistWithVector(h_zchi2kari, kari_zchi2ns);

    TH1F *h_ptkari = new TH1F("h_ptkari", "pt_{kari} transverse momentum", 30, -40000, 40000);
    labelAxis(h_ptkari, "pt [MeV]", "count");
    fillHistWithVector(h_ptkari, kari_pts);

    TH1F *h_ptkari_inv = new TH1F("h_ptkari_inv", "karimaki pt^{-1}_{kari}", 30, -1e-3, 1e-3);
    labelAxis(h_ptkari_inv, "pt^{-1} [MeV^{-1}]", "count");
    fillHistWithVector(h_ptkari_inv, kari_inv_pts);

    TH1F *h_pterr_kari = new TH1F("h_pterr_kari", "pt^{-1}_{kari} #minus pt^{-1}_{mc} transverse momentum error", 30, -0.0005, 0.0005);
    labelAxis(h_pterr_kari, "pt^{-1}_{kari} #minus pt^{-1}_{mc} [MeV^{-1}]", "count");
    fillHistWithVector(h_pterr_kari, p_inv_kari_errors);




    ///FILLING SCATTER PLOTS

    TGraph *g_pdev_phi = new TGraph(p_inv_rel_errors.size(), &rec_phis[0],&p_inv_rel_errors[0]);
    g_pdev_phi->SetTitle("p_{rec} #minus p_{mc}^{-1} over #Phi correlation");
    labelAxis(g_pdev_phi, "#Phi", "p_{rec}^{-1} #minus p_{mc}^{-1} [MeV^{-1}]");
    setGraphRange(g_pdev_phi, -3.2, 0, -4e-4, 4e-4);

    TGraph *g_pdev_dca = new TGraph(p_inv_rel_errors.size(), &rec_dca_rs[0],&p_inv_rel_errors[0]);
    g_pdev_dca->SetTitle("p_{rec}^{-1} #minus p_{mc}^{-1} over r-dca_{rec} correlation");
    labelAxis(g_pdev_dca, "r-dca_{rec} [mm]", "p_{rec}^{-1} #minus p_{mc}^{-1} [MeV^{-1}]");
    setGraphRange(g_pdev_dca, -80, 80, -4e-4, 4e-4);

    TGraph *g_pdev_z_dca = new TGraph(p_inv_rel_errors.size(), &rec_dca_zs[0],&p_inv_rel_errors[0]);
    g_pdev_z_dca->SetTitle(" p_{rec}^{-1} #minus p_{mc}^{-1} over z-dca_{rec} correlation");
    labelAxis(g_pdev_z_dca, "z-dca_{rec} [mm]", "p_{rec}^{-1} #minus p_{mc}^{-1} [MeV^{-1}]");
    setGraphRange(g_pdev_z_dca, -500, 500, -4e-4, 4e-4);

    TGraph *g_invpt_invptmc = new TGraph(rec_inv_pts.size(),&mc_inv_pts[0],&rec_inv_pts[0]);
    g_invpt_invptmc->SetTitle("pt_{rec}^{-1} over pt_{mc}^{-1} correlation");
    labelAxis(g_invpt_invptmc, "pt_{mc}^{-1} [MeV^{-1}]", "pt_{rec}^{-1} [MeV^{-1}]");
    setGraphRange(g_invpt_invptmc, -5e-4,5e-4,-1e-3, 1e-3);

    TGraph *g_ptdev_ptmc = new TGraph(pt_inv_errors.size(),&mc_inv_pts[0],&pt_inv_errors[0]);
    g_ptdev_ptmc->SetTitle("pt_{rec}^{-1} #minus pt_{mc}^{-1} over pt_{mc}^{-1} correlation");
    labelAxis(g_ptdev_ptmc, "pt_{mc}^{-1} [MeV^{-1}]", "pt_{rec}^{-1} #minus pt_{mc}^{-1} [MeV^{-1}]");
    setGraphRange(g_ptdev_ptmc,-5e-4,5e-4, -1e-3,1e-3);


    // KARI Scatter plots

    //pkari-pmc
    TGraph *g_ptkari_ptmc = new TGraph(kari_inv_pts.size(),&mc_inv_pts[0],&kari_inv_pts[0]);
    g_ptkari_ptmc->SetTitle("pt_{kari}^{-1} over pt_{mc}^{-1} correlation");
    labelAxis(g_ptkari_ptmc, "pt_{mc}^{-1} [MeV^{-1}]", "pt_{kari}^{-1} [MeV^{-1}]");
    setGraphRange(g_ptkari_ptmc,-5e-4,5e-4, -1e-3,1e-3);

    //perr pkari
    TGraph *g_ptkaridev_ptmc = new TGraph(p_inv_kari_errors.size(),&mc_inv_pts[0],&p_inv_kari_errors[0]);
    g_ptkaridev_ptmc->SetTitle("pt_{kari}^{-1} #minus pt_{mc}^{-1} over pt_{mc}^{-1} correlation");
    labelAxis(g_ptkaridev_ptmc, "pt_{mc}^{-1} [MeV^{-1}]", "pt_{kari}^{-1} #minus pt_{mc}^{-1} [MeV^{-1}]");
    setGraphRange(g_ptkaridev_ptmc,-5e-4,5e-4, -1e-3,1e-3);

    TGraph *g_karizdca_mczdca = new TGraph(kari_z0s.size(),&rec_dca_zs[0],&kari_z0s[0]);
    g_karizdca_mczdca->SetTitle("z-dca_{kari} over z-dca_{rec}  correlation");
    labelAxis(g_karizdca_mczdca, "z-dca_{rec}", "z-dca_{kari}");
    setGraphRange(g_karizdca_mczdca,-600, 600, -600, 600);

    //pterr kari - rec over pt rec
    TGraph *g_ptdevkvsms_ptms = new TGraph(pt_inv_kvsms_errors.size(),&rec_inv_pts[0],&pt_inv_kvsms_errors[0]);
    g_ptdevkvsms_ptms->SetTitle("pt_{kari}^{-1} #minus pt_{msfit}^{-1} over pt_{msfit}^{-1} correlation");
    labelAxis(g_ptdevkvsms_ptms, "#pt_{msfit}^{-1} [MeV^{-1}]", "pt_{kari}^{-1} #minus pt_{msfit}^{-1} [MeV^{-1}]");
    setGraphRange(g_ptdevkvsms_ptms,-8e-4,8e-4, -1e-3,1e-3);

    //perr rdca
    TGraph *g_ptkaridev_rdca = new TGraph(p_inv_kari_errors.size(),&kari_dcas[0],&p_inv_kari_errors[0]);
    g_ptkaridev_rdca->SetTitle("pt_{kari}^{-1} #minus pt_{mc}^{-1} over r-dca_{kari} correlation");
    labelAxis(g_ptkaridev_rdca, "r-dca_{kari} [mm]", "pt_{kari}^{-1} #minus pt_{mc}^{-1} [MeV^{-1}]");
    setGraphRange(g_ptkaridev_rdca,-80,80, -4e-4,4e-4);

    //perr zdca
    TGraph *g_ptkaridev_zdca = new TGraph(p_inv_kari_errors.size(),&kari_z0s[0],&p_inv_kari_errors[0]);
    g_ptkaridev_zdca->SetTitle("pt_{kari}^{-1} #minus pt_{mc}^{-1} over z-dca_{kari} correlation");
    labelAxis(g_ptkaridev_zdca, "z-dca_{kari} [mm]", "pt_{kari}^{-1} #minus pt_{mc}^{-1} [MeV^{-1}]");
    setGraphRange(g_ptkaridev_zdca,-500,500, -4e-4,4e-4);

    //perr phi_mc
    TGraph *g_ptkaridev_mcphi = new TGraph(p_inv_kari_errors.size(),&mc_phis[0],&p_inv_kari_errors[0]);
    g_ptkaridev_mcphi->SetTitle("pt_{kari}^{-1} #minus pt_{mc}^{-1} over #Phi_{mc} correlation");
    labelAxis(g_ptkaridev_mcphi, " #Phi_{mc}", "pt_{kari}^{-1} #minus pt_{mc}^{-1} [MeV^{-1}]");
    setGraphRange(g_ptkaridev_mcphi,-3.2, 3.2, -1e-3,1e-3);

    ////###### PLOTS ##################################################

    TGraph * graphs[4];
    TH1F * hists[4];

    graphs[0] = g_ptdev_ptmc; graphs[1]= g_ptkaridev_ptmc;
    makeSimpleMultiCanvas(1, 2, 2, graphs, plottingfile + "(");

    graphs[0] = g_invpt_invptmc; graphs[1]= g_ptkari_ptmc;
    makeSimpleMultiCanvas(1, 2, 2, graphs, plottingfile);

    hists[0] = h_pt_inv_err; hists[1]= h_pterr_kari;
    makeSimpleMultiCanvas(1, 2, 2, hists, plottingfile);

    hists[0] =h_invptmc; hists[1]= h_invptrec; hists[2]= h_ptkari_inv;
    makeSimpleMultiCanvas(1, 3, 3, hists, plottingfile);

    hists[0] =h_ptmc; hists[1]= h_pt; hists[2]= h_ptkari;
    makeSimpleMultiCanvas(1, 3, 3, hists, plottingfile);

    hists[0] = h_rdca; hists[1]= h_dcakari;
    makeSimpleMultiCanvas(1, 2, 2, hists, plottingfile);

    hists[0] = h_zdca; hists[1]= h_z0kari;
    makeSimpleMultiCanvas(1, 2, 2, hists, plottingfile);

    //kari stuff
    hists[0] = h_r3dkari; hists[1]= h_rinvkari;
    hists[2] = h_tchi2nkari; hists[3]= h_zchi2kari;
    makeSimpleMultiCanvas(1, 2, 2, hists, plottingfile);

    //graphs
    graphs[0] = g_pdev_z_dca; graphs[1]= g_ptkaridev_zdca;
    makeSimpleMultiCanvas(1, 2, 2, graphs, plottingfile);

    graphs[0] = g_pdev_dca; graphs[1]= g_ptkaridev_rdca;
    makeSimpleMultiCanvas(1, 2, 2, graphs, plottingfile);

    graphs[0] = g_pdev_phi; graphs[1]= g_ptkaridev_mcphi;
    makeSimpleMultiCanvas(1, 2, 2, graphs, plottingfile);

    graphs[0] = g_karizdca_mczdca; graphs[1]= g_ptdevkvsms_ptms;
    makeSimpleMultiCanvas(1, 2, 2, graphs, plottingfile);

    //mc info stuff
    hists[0] = h_rec_nhits; hists[1]= h_mctype;
    makeSimpleMultiCanvas(1, 2, 2, hists, plottingfile);

    //phi
    hists[0] =h_phimc; hists[1]= h_phi; hists[2]= h_phikari;
    makeSimpleMultiCanvas(1, 3, 3, hists, plottingfile);

    //theta
    hists[0] = h_theta; hists[1]= h_thetakari;
    makeSimpleMultiCanvas(1, 2, 2, hists, plottingfile + ")");

    ////###########################################################################

#endif

}


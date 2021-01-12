//
// Created by Konstantin Neureither on 08.07.20.
//

#define GET_DATA_FROM_MU3E false

#include <TFile.h>
#include <TROOT.h>
#include <algorithm>
#include <TTree.h>
#include <string>
#include <TH1.h>
#include <TGraph.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <filesystem>
#include <cmath>
#include <assert.h>
#include "utilityFunctions.h"
#include "plots.h"
#include "rootData.h"
#include "../karimaki/karimakiHelixfit.h"
#include "reconstruction_accuracy.h"
#include "basicDefines.h"
#include "SegsTreeRead.h"

using std::cout;
using std::endl;

void reconstructionAccuracyScript(int run, int FILTER) {

    const std::string pathtodata = "data/SimulationData/";
    const std::string pathtoplots = "output/0_RecAcc/";
    const bool RECONSTRUCTION_PRINTS = false;
    const bool HIT_PRINTS = DEBUG;
    const bool MAKE_PLOT = true;
    const int MAX_ENTRIES = 0;

    std::string filtertag;

    std::string runpadded = get_padded_string(run, 6, '0');
    std::string infile2 = pathtodata + "mu3e_run_" + runpadded + "_trirec_cosmic.root";
    std::string pathtorunplots = pathtoplots + "run_" + get_padded_string(run, 3, '0') + "/";

    check_create_directory(pathtodata);
    check_create_directory(pathtoplots);
    check_create_directory(pathtorunplots);

    TFile f2(infile2.c_str());
    TTree *t_segs;
    f2.GetObject("segs", t_segs);

#if GET_DATA_FROM_MU3E
    std::string infile1 = pathtodata + "mu3e_run_" + runpadded + ".root";
    TFile f1(infile1.c_str());
    TTree *t_mu3e;
    f1.GetObject("mu3e", t_mu3e);

    //data in mu3e tree (not used)
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

    t_mu3e->SetBranchAddress("hit_pixelid",&pixelid);
    t_mu3e->SetBranchAddress("traj_ID", &trajids);
    t_mu3e->SetBranchAddress("traj_px", &traj_px);
    t_mu3e->SetBranchAddress("traj_py", &traj_py);
    t_mu3e->SetBranchAddress("traj_pz", &traj_pz);

    cout << "Branches set for mu3e..." << endl;
#endif
    //class representation for segs tree and read functionality
    SegsTreeReadPlus Segs = SegsTreeReadPlus(t_segs);

    //For Helix Fit
    KariFit karires;


    //// further data for stats graphs and calculations
    int p_fail_count = 0;
    int rkari_swap_count = 0;
    int mckari_wrong_sign_count = 0;
    int mcrec_wrong_sign_count = 0;
    int reckari_wrong_sign_count = 0;
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

    //data is collected depending on entry
    std::vector<float> p_inv_rel_errors_hits[6];
    std::vector<float> p_inv_err_nhits[3];
    std::vector<float> pt_inv_err_nhits[3];
    std::vector<float> rec_rdca_nhits[3];
    std::vector<float> p_inv_err_r_dcas[4];
    std::vector<float> p_inv_err_z_dcas[4];

    //kari plots
    std::vector<float> pt_kari_inv_err_nhits[3];
    std::vector<float> pt_kari_inv_err_r_dcas[4];
    std::vector<float> pt_kari_inv_err_z_dcas[4];

    double RMS = 80 * 1e-3 / sqrt(12); //RMS in mm
    float BFIELD = 1.0;

    double tres[TRIPLET_HIT_ARRAY_LENGTH] = {RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS};
    double zres[TRIPLET_HIT_ARRAY_LENGTH] = {RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS, RMS};
    double rres[TRIPLET_HIT_ARRAY_LENGTH] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

#if GET_DATA_FROM_MU3E
    unsigned int mu3e_index = 1;
    t_mu3e->GetEntry(mu3e_index);
    mu3e_index++;
#endif


    for (unsigned int entryno = 0; entryno < (MAX_ENTRIES == 0 ? Segs.my_entries : MAX_ENTRIES); entryno++) {
        Segs.getEntry(entryno);

        if(HIT_PRINTS) printf("\nRECONSTRUCTION ACCURACY ENTRY\n------------------------------\n\n");

#if GET_DATA_FROM_MU3E
        //find corresponding entry in mu3e tree
        while (rec_event > header[0]) {
            t_mu3e->GetEntry(mu3e_index);
            mu3e_index++;
        }
        //if there is no corresponding entry found, skip
        if(rec_event != header[0]) continue;

        if(HIT_PRINTS) {
            //some status prints for debugging
            cout << "\tsegs_index=" << entryno << "\tmu3e_index=" << mu3e_index;
            cout << "\trec_event=" << rec_event << "\tmu3e_event= " << header[0];
            cout << "\trec_nhits=" << rec_nhit << "\tmu3e_nhit=" << mu3e_nhits << endl;
        }

        //get hit data from mu3e tree (not used, only for debugging)
        for (unsigned int hit = 0; hit < (unsigned int) mu3e_nhits; hit++) {
            PXID pixid = process_pixel_id((*pixelid)[hit]);
            unsigned int trajid = (*trajids)[hit];
            double trajpx = (*traj_px)[hit];

            if(HIT_PRINTS) {
                unsigned int layer = get_layer(pixid.sensor);
                cout << "\thit=" << hit << "\tlayer=" << layer << "\tsid=" << pixid.sensor << "\trow=" << pixid.row << "\tcol=" << pixid.column << endl;
            }
        }

        //get trajectory data from mu3e tree
        unsigned int trajid = (*trajids)[0];
        double trajpx = (*traj_px)[0];
        double trajpy = (*traj_py)[0];
        double trajpz = (*traj_pz)[0];
        double trajp = sqrt(pow(trajpx, 2) + pow(trajpy, 2) + pow(trajpz, 2));
#endif

        //print triplet data from segs tree
        if (HIT_PRINTS) {
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
        }

        //Theta, phi and traverse p of reconstruction
        Segs.calcAdditionalData();

        //Do Karimaki Helix fit
        std::vector<double> xp;
        std::vector<double> yp;
        std::vector<double> zp;
        std::vector<int> sids;
        std::vector<double> phi_hits;
        std::vector<double> thetas;

        int ncombinedhits = combineHits(xp, yp, zp, sids, phi_hits, thetas, Segs, false);

//        getResolutionArrs(ncombinedhits, RMS, tres, zres, rres, xp, yp, phi_hits, (double) rec_theta);
//        assert(tres.size() == ncombinedhits);

        if(HIT_PRINTS) {
            for(int i = 0; i < ncombinedhits; i++) {
                cout << "\tHits sorted for kari:  sid=" << sids[i] << " x=" << xp[i] << "  y=" << yp[i] << "  z=" << zp[i] <<" phi=" << phi_hits[i] << " theta=" << thetas[i] << endl;
            }
        }
        karimakiHelixfit(karires, ncombinedhits, &xp[0], &yp[0], &zp[0], &phi_hits[0], &thetas[0], &tres[0], &zres[0],
                         &rres[0]);
        correctKariDirection(karires);
        swapKariBField(karires);

        if (Segs.mc_p == 0 || Segs.mc_pt == 0 || Segs.rec_p == 0 || Segs.rec_pt == 0) {
            p_fail_count++;
        } else {

            float rec_p_corr = Segs.rec_p * sgn(Segs.rec_r);
            float rec_pt_corr = Segs.rec_pt;
            float kari_pt = 0.3 * karires.rad * BFIELD;


            float kari_inv_pt = 1. / kari_pt;
            float kari_inv_rad = 1. / karires.rad;
            float p_abs_error = (Segs.mc_p_corr - Segs.rec_p);
            float p_over_pmc = Segs.rec_p / Segs.mc_p_corr;
            float mc_z_dca = std::sin(Segs.mc_vpca_phi) * Segs.mc_vpca_offset;
            float pt_kari_inv_err = kari_inv_pt - Segs.mc_pt_inv_corr;
            float pt_kvsms_inv_err = kari_inv_pt - Segs.rec_pt_inv;


            //Filtering entries //FIXME filtertag is defined over and over again
            bool choice;
            switch(FILTER) {
                case 0: choice = Segs.mc_type % 10 == 3 || Segs.mc_type % 10 == 4;
                    filtertag = "muonsonly";
                    break;
                case 1: choice = pt_kari_inv_err / Segs.mc_pt_inv_corr < -1.8 && pt_kari_inv_err / Segs.mc_pt_inv_corr > -2.2;
                    filtertag = "kari-outliers";
                    break;
                case 2: choice = Segs.mc_type == 3;
                    filtertag = "mu+";
                    break;
                case 3: choice = Segs.mc_type == 4;
                    filtertag = "mu-";
                    break;
                case 4: choice = sgn(karires.rad) != sgn(Segs.mc_p_corr);
                    filtertag = "p-kari-mc-sgn";
                    break;
                case 5: choice = karires.phi > 0 && karires.phi < PI/2 && (pt_kari_inv_err / Segs.mc_pt_inv_corr < -1.8 && pt_kari_inv_err / Segs.mc_pt_inv_corr > -2.2);
                    filtertag = "special1";
                    break;
                case 6: choice = Segs.mc_type == 3 || Segs.mc_type == 4;
                    filtertag = "planemuonsonly";
                    break;
                case 7: choice = true;
                    filtertag = "all";
                    break;
                case 8: choice = (Segs.mc_type % 10 == 3 || Segs.mc_type % 10 == 4) && !(Segs.mc_type == 3 || Segs.mc_type == 4);
                    filtertag = "scatteredmuonsonly";
                    break;
                case 9: choice = !(Segs.mc_type == 3 || Segs.mc_type == 4);
                    filtertag = "exceptplanemuons";
                    break;
                default:choice = true;
                    filtertag = "";
            }

            if (choice) {


                cout << "event id " << Segs.rec_event << " index: " << entryno << " mc type: " << Segs.mc_type << endl;

                //calculated data
                rec_inv_ps.push_back(Segs.rec_p_inv);
                rec_p_corrs.push_back(rec_p_corr);
                rec_inv_pts.push_back(Segs.rec_p_inv);
                p_inv_rel_errors.push_back(Segs.p_inv_abs_error);
                p_rel_errors.push_back(p_abs_error);
                pt_inv_errors.push_back(Segs.pt_inv_abs_error);
                p_over_pmcs.push_back(p_over_pmc);
                rec_pts.push_back(Segs.rec_pt);
                rec_phis.push_back(Segs.rec_phi);
                rec_thetas.push_back(Segs.rec_theta);

                mc_p_corrs.push_back(Segs.mc_p_corr);
                mc_pt_corrs.push_back(Segs.mc_pt_corr);
                mc_inv_pts.push_back(Segs.mc_pt_inv_corr);
                mc_z_dcas.push_back(mc_z_dca);
                mc_inv_ps.push_back(Segs.mc_p_inv_corr);

                //reconstruction data
                rec_ps.push_back(Segs.rec_p);
                rec_rs.push_back(Segs.rec_r);
                rec_dca_rs.push_back(Segs.rec_zpca_r);
                rec_dca_xs.push_back(Segs.rec_zpca_x);
                rec_dca_ys.push_back(Segs.rec_zpca_y);
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
                kari_inv_rads.push_back(kari_inv_rad);
                pt_inv_kvsms_errors.push_back(pt_kvsms_inv_err);

                //monte carlo data
                mc_ps.push_back(Segs.mc_p);
                mc_pts.push_back(Segs.mc_pt);
                mc_phis.push_back(Segs.mc_phi);
                mc_dcas.push_back(Segs.mc_vpca_offset);
                mc_phi_dcas.push_back(Segs.mc_vpca_phi);
                mc_types.push_back(Segs.mc_type);

                //meta data
#if GET_DATA_FROM_MU3E
                sim_nhits.push_back(mu3e_nhits);
#endif
                rec_nhits.push_back(Segs.rec_nhit);
                rec_nhits_float.push_back((float)Segs.rec_nhit);

                if(Segs.rec_p / Segs.mc_p_corr < 0) mcrec_wrong_sign_count++;
                if(karires.rad / Segs.mc_p_corr < 0) mckari_wrong_sign_count++;
                if(karires.rad / Segs.rec_p < 0) reckari_wrong_sign_count++;
                processed_entries++;

                switch(Segs.rec_nhit) {
                    case 4:
                        p_inv_rel_errors_hits[0].push_back(Segs.p_inv_abs_error);
                        p_inv_err_nhits[0].push_back(Segs.p_inv_abs_error);
                        pt_inv_err_nhits[0].push_back(Segs.pt_inv_abs_error);
                        rec_rdca_nhits[0].push_back(Segs.rec_zpca_r);
                        pt_kari_inv_err_nhits[0].push_back(pt_kari_inv_err);
                        break;
                    case 5:
                        p_inv_rel_errors_hits[1].push_back(Segs.p_inv_abs_error);
                        p_inv_err_nhits[0].push_back(Segs.p_inv_abs_error);
                        pt_inv_err_nhits[0].push_back(Segs.pt_inv_abs_error);
                        rec_rdca_nhits[0].push_back(Segs.rec_zpca_r);
                        pt_kari_inv_err_nhits[0].push_back(pt_kari_inv_err);
                        break;
                    case 6:
                        p_inv_rel_errors_hits[2].push_back(Segs.p_inv_abs_error);
                        p_inv_err_nhits[1].push_back(Segs.p_inv_abs_error);
                        pt_inv_err_nhits[1].push_back(Segs.pt_inv_abs_error);
                        rec_rdca_nhits[1].push_back(Segs.rec_zpca_r);
                        pt_kari_inv_err_nhits[1].push_back(pt_kari_inv_err);
                        break;
                    case 7:
                        p_inv_rel_errors_hits[3].push_back(Segs.p_inv_abs_error);
                        p_inv_err_nhits[1].push_back(Segs.p_inv_abs_error);
                        pt_inv_err_nhits[1].push_back(Segs.pt_inv_abs_error);
                        rec_rdca_nhits[1].push_back(Segs.rec_zpca_r);
                        pt_kari_inv_err_nhits[1].push_back(pt_kari_inv_err);
                        break;
                    case 8:
                        p_inv_rel_errors_hits[4].push_back(Segs.p_inv_abs_error);
                        p_inv_err_nhits[2].push_back(Segs.p_inv_abs_error);
                        pt_inv_err_nhits[2].push_back(Segs.pt_inv_abs_error);
                        rec_rdca_nhits[2].push_back(Segs.rec_zpca_r);
                        pt_kari_inv_err_nhits[2].push_back(pt_kari_inv_err);
                        break;
                    default:
                        p_inv_rel_errors_hits[5].push_back(Segs.p_inv_abs_error);
                        p_inv_err_nhits[2].push_back(Segs.p_inv_abs_error);
                        pt_inv_err_nhits[2].push_back(Segs.pt_inv_abs_error);
                        rec_rdca_nhits[2].push_back(Segs.rec_zpca_r);
                        pt_kari_inv_err_nhits[2].push_back(pt_kari_inv_err);
                }

                //filling data for histograms of errors depending on dca r
                if(0 <= Segs.rec_zpca_r && Segs.rec_zpca_r < 40 ) {
                    p_inv_err_r_dcas[0].push_back(Segs.p_inv_abs_error);
                    pt_kari_inv_err_r_dcas[0].push_back(pt_kari_inv_err);
                } else if (40 <= Segs.rec_zpca_r && Segs.rec_zpca_r < 50) {
                    p_inv_err_r_dcas[1].push_back(Segs.p_inv_abs_error);
                    pt_kari_inv_err_r_dcas[1].push_back(pt_kari_inv_err);
                } else if (50 <= Segs.rec_zpca_r && Segs.rec_zpca_r < 60) {
                    p_inv_err_r_dcas[2].push_back(Segs.p_inv_abs_error);
                    pt_kari_inv_err_r_dcas[2].push_back(pt_kari_inv_err);
                } else {
                    p_inv_err_r_dcas[3].push_back(Segs.p_inv_abs_error);
                    pt_kari_inv_err_r_dcas[3].push_back(pt_kari_inv_err);
                }

                //filling data for histograms of errors depending on dca z
                int abs_dca_z = abs(Segs.rec_zpca_z);
                if(abs_dca_z < 50) {
                    p_inv_err_z_dcas[0].push_back(Segs.p_inv_abs_error);
                    pt_kari_inv_err_z_dcas[0].push_back(pt_kari_inv_err);
                } else if(50 <= abs_dca_z && abs_dca_z < 100) {
                    p_inv_err_z_dcas[1].push_back(Segs.p_inv_abs_error);
                    pt_kari_inv_err_z_dcas[1].push_back(pt_kari_inv_err);
                } else if(100 <= abs_dca_z && abs_dca_z < 200) {
                    p_inv_err_z_dcas[2].push_back(Segs.p_inv_abs_error);
                    pt_kari_inv_err_z_dcas[2].push_back(pt_kari_inv_err);
                } else {
                    p_inv_err_z_dcas[3].push_back(Segs.p_inv_abs_error);
                    pt_kari_inv_err_z_dcas[3].push_back(pt_kari_inv_err);
                }

                ////PRINT SECTION PER ENTRY IN TREE
                if(RECONSTRUCTION_PRINTS) {
                    //data from reconstruction
                    cout << "rec_event: " << Segs.rec_event;
#if GET_DATA_FROM_MU3E
                    //data from simulation
                    cout << " mu3e_event: " << header[0] << endl;
                    cout << "\t" << "traj id; (px, py, pz) " <<  trajid << "; (" << trajpx << ", " << trajpy << ", " << trajpz << ") ";
                    cout << "p_calc: " << trajp << "\t";
                    cout << "\tp_calc / mc_p: " << double(trajp / mc_p) << endl;
#endif
                    cout << " mc_p: " << Segs.mc_p << " rec_p *sgn(rec_r): " << Segs.rec_p*sgn(Segs.rec_r) << "rec_nhit: "<< Segs.rec_nhit <<"\t\t";
                    cout << "rec_r= " << Segs.rec_r << " mc_type=" << Segs.mc_type;

                    //mc and rec deviations
                    cout << "pt_mc= " << Segs.mc_pt << " pt_rec= " << Segs.rec_pt << " phi_mc= " << Segs.mc_phi << " rec_phi= " << Segs.rec_phi;
                    cout << " mc_theta= " << Segs.mc_theta << " rec_theta= " << Segs.rec_theta << "z-dca=" << mc_z_dca <<"\t\t";

                    cout << "(1/mc_p-1/rec_p)/(1/mc_p): " << Segs.p_inv_abs_error*100 << " %" << endl;
                }
            }
        }

    }

    if (MAKE_PLOT) {

        const float LEFT_BOUNDARY_P = -3e4;
        const float RIGHT_BOUNDARY_P = 3e4 ;
        const int BIN_COUNT = 50;

        const bool MAKE_468HIT_HISTOGRAMS = true;
        const bool MAKE_DCA_AREA_HISTOGRAMS = false;
        const bool MAKE_P_REL_PLOTS = false;
        const bool MAKE_STRANGE_BAND_CORR_PLOTS = true;

        std::string filename;
        std::string filename_template = pathtorunplots + "reconstruction-accuracy_run" +
                                        get_padded_string(run, 3, '0') + "_FILTER" + filtertag;
        std::string plottingfile = filename_template + "_plots.pdf";
        std::string plottingfilekari = filename_template + "_karifit_plots.pdf";
        std::string plottingfilerec = filename_template + "_msfit_plots.pdf";
        std::string plottingfilecomparison = filename_template + "_comparison_plots.pdf";


        ///FILLING THE HISTOGRAMS

        TH1F *h_nhits = new TH1F("h_nhits", "hits per frame", 20, 0., 20.);
        labelAxis(h_nhits, "number of hits", "count");
        fillHistWithVector(h_nhits, sim_nhits);

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

        TH1F *h_mctype = new TH1F("h_mctype", "Particle types of simulation", 100, -1000, 1000);
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

        //Calculated deviations and relations
        TH1F *h_p_relerror = new TH1F("h_p_relerror", "p_{rec} p_{mc} deviation", BIN_COUNT, -2.0, 1.0);
        labelAxis(h_p_relerror, "p - p_{mc}) [MeV]", "count");
        fillHistWithVector(h_p_relerror, p_rel_errors);


        TH1F *h_pinv_relerror = new TH1F("h_pinv_relerror", "p_{rec}^{-1} #minus p_{mc}^{-1} deviation",
                                         BIN_COUNT, -0.0003, 0.0003);
        labelAxis(h_pinv_relerror, "p_{rec}^{-1} #minus p_{mc}^{-1} [MeV^{-1}]", "count");
        fillHistWithVector(h_pinv_relerror, p_inv_rel_errors);

        TH1F *h_poverpmc = new TH1F("p / p_{mc} relation", "p / p_{mc}", 30, -0, 2);
        labelAxis(h_poverpmc, "#frac{p}{p_{mc}}", "count");
        fillHistWithVector(h_poverpmc, p_over_pmcs);

        //DCAs Monte Carlo
        TH1F *h_dcamc = new TH1F("h_dcamc", "dca monte carlo (var: mc_vpca_offset)", 30, -100, 100);
        labelAxis(h_dcamc, "dca [mm]", "count");
        fillHistWithVector(h_dcamc, mc_dcas);

        TH1F *h_zdcamc = new TH1F("h_zdcamc", "z-dca mc (var: sin(mc_vpca_phi) * mc_vpca_offset)[mm]", 40, -30, 5);
        labelAxis(h_zdcamc, "z dca", "count");
        fillHistWithVector(h_zdcamc, mc_z_dcas);

        TH1F *h_dcamc_phi = new TH1F("h_dcamc_phi", "dca phi monte carlo (var: mc_vpca_phi) [mm]", 30, -3.2, 1);
        labelAxis(h_dcamc_phi, "phi_{mc} of closest approach", "count");
        fillHistWithVector(h_dcamc_phi, mc_phi_dcas);

        //DCAs Reconstruction
        TH1F *h_rdca = new TH1F("h_rdca", "dca_{reconstruction} r", 30, -100, 100);
        labelAxis(h_rdca, "dca_{r} [mm]", "count");
        fillHistWithVector(h_rdca, rec_dca_rs);

        TH1F *h_xdca = new TH1F("h_xdca", "dca_{reconstruction} along x-axis", 30, -100, 100);
        labelAxis(h_xdca, "dca_{x} [mm]", "count");
        fillHistWithVector(h_xdca, rec_dca_xs);

        TH1F *h_ydca = new TH1F("h_ydca", "dca_{reconstruction} along y-axis", 30, -100, 100);
        labelAxis(h_ydca, "dca_{y} [mm]", "count");
        fillHistWithVector(h_ydca, rec_dca_ys);

        TH1F *h_zdca = new TH1F("h_zdca", "dca_{reconstruction} along z-axis", 30, -600, 600);
        labelAxis(h_zdca, "dca_{z} [mm]", "count");
        fillHistWithVector(h_zdca, rec_dca_zs);

        //pt nhit dependent error histograms
        std::vector<TH1F*> h_ptinv_errhits;
        for(int i=0; i<3; i++) {
            std::string histtitle = "pt_{rec}^{-1} #minus pt_{mc}^{-1} for " + get_string(2*i+4) + " and " +
                                    (i==2 ? " more" : get_string(2*i+5)) + " hits";
            std::string histhandle = "h_ptinv_" + get_string(2*i+4) + "hits";
            TH1F *h = new TH1F(histhandle.c_str(), histtitle.c_str(), BIN_COUNT, -0.0005, 0.0005);
            labelAxis(h, "pt_{rec}^{-1} #minus pt_{mc}^{-1} [Mev^{-1}]", "count");
            fillHistWithVector(h, pt_inv_err_nhits[i]);
            h_ptinv_errhits.push_back(h);
        }

        //r dca nhit dependent histograms
        std::vector<TH1F*> h_rdca_errhits;
        for(int i=0; i<3; i++) {
            std::string histtitle = "dca_{reconstruction} radius for " + get_string(2*i+4) + " and " +
                                    (i==2 ? " more" : get_string(2*i+5)) + " hits";
            std::string histhandle = "h_rdca_" + get_string(2*i+4) + "hits";
            TH1F *h = new TH1F(histhandle.c_str(), histtitle.c_str(), BIN_COUNT, -100, 100);
            labelAxis(h, "dca_{r} [mm]", "count");
            fillHistWithVector(h, rec_rdca_nhits[i]);
            h_rdca_errhits.push_back(h);
        }

        //p nhit dependent error histograms
        std::vector<TH1F*> h_pinv_errhits;
        for(int i=0; i<3; i++) {
            std::string histtitle = "p_{rec}^{-1} #minus p_{mc}^{-1} for " + get_string(2*i+4) + " and " +
                                    (i==2 ? " more" : get_string(2*i+5)) + " hits";
            std::string histhandle = "h_pinv_" + get_string(2*i+4) + "hits";
            TH1F *h = new TH1F(histhandle.c_str(), histtitle.c_str(), BIN_COUNT, -0.0005, 0.0005);
            labelAxis(h, "p_{rec}^{-1} #minus p_{mc}^{-1} [Mev^{-1}]", "count");
            fillHistWithVector(h, p_inv_err_nhits[i]);
            h_pinv_errhits.push_back(h);
        }

        //p r dca dependent error histograms
        std::vector<TH1F*> h_pinv_err_rdca;
        std::string radialintervals[4] = {"0-4", "4-5", "5-6", "> 6"};
        for(int i=0; i<4; i++) {
            std::string histtitle = "p_{rec}^{-1} #minus p_{mc}^{-1} for r-dca_{rec} " + radialintervals[i] + " cm";
            std::string histhandle = "h_pinv_r_" + radialintervals[i] + "cm";
            TH1F *h = new TH1F(histhandle.c_str(), histtitle.c_str(), BIN_COUNT, -0.0005, 0.0005);
            labelAxis(h, "p_{rec}^{-1} #minus p_{mc}^{-1} [Mev^{-1}]", "count");
            fillHistWithVector(h, p_inv_err_r_dcas[i]);
            h_pinv_err_rdca.push_back(h);
        }

        //p z dca dependent error histograms
        std::vector<TH1F*> h_pinv_err_zdca;
        std::string zdcaintervals[4] = {"0 - 5", "5 -10", "10-20", "> 20 "};
        for(int i=0; i<4; i++) {
            std::string histtitle = "p_{rec}^{-1} #minus p_{mc}^{-1} for abs(z-dca_{rec}) " + zdcaintervals[i] + " cm";
            std::string histhandle = "h_pinv_z_" + zdcaintervals[i] + "cm";
            TH1F *h = new TH1F(histhandle.c_str(), histtitle.c_str(), BIN_COUNT, -0.0005, 0.0005);
            labelAxis(h, "p_{rec}^{-1} #minus p_{mc}^{-1} [Mev^{-1}]", "count");
            fillHistWithVector(h, p_inv_err_z_dcas[i]);
            h_pinv_err_zdca.push_back(h);
        }

        //p nhit dependent error histograms (hits individually) (for overlay hist)
        std::vector<TH1F*> h_pinv_relerrorhits;
        for(int i = 0; i < 6; i++) {
            std::string histhandle = "h_pinv_relerror_hits_" + get_string(i);
            std::string histtitle = "p_{rec}^{-1} #minus p_{mc}^{-1} deviations, nhit=" + get_string(i+4);
            TH1F *h = new TH1F(histhandle.c_str(), histtitle.c_str(), BIN_COUNT, -0.0003, 0.0003);
            labelAxis(h, "p_{rec}^{-1} #minus p_{mc}^{-1} [MeV^{-1}]", "count");
            fillHistWithVector(h, p_inv_rel_errors_hits[i]);
            h_pinv_relerrorhits.push_back(h);
        }

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

        TH1F *h_pterr_kari = new TH1F("h_pterr_kari", "pt^{-1}_{kari} #minus pt^{-1}_{mc} transverse momentum error", 30, -0.0003, 0.0003);
        labelAxis(h_pterr_kari, "pt^{-1}_{kari} #minus pt^{-1}_{mc} [MeV^{-1}]", "count");
        fillHistWithVector(h_pterr_kari, p_inv_kari_errors);

        //KARIMAKI DCA, NHIT, Z0 HISTOGRAMS

        //kari pt nhit dependent error histograms
        std::vector<TH1F*> h_kari_ptinv_errhits;
        for(int i=0; i<3; i++) {
            std::string histtitle = "pt_{kari}^{-1} #minus pt_{mc}^{-1} for " + get_string(2*i+4) + " and " +
                                    (i==2 ? " more" : get_string(2*i+5)) + " hits";
            std::string histhandle = "h_k_ptinv_" + get_string(2*i+4) + "hits";
            TH1F *h = new TH1F(histhandle.c_str(), histtitle.c_str(), BIN_COUNT, -0.0005, 0.0005);
            labelAxis(h, "pt_{kari}^{-1} #minus pt_{mc}^{-1} [Mev^{-1}]", "count");
            fillHistWithVector(h, pt_kari_inv_err_nhits[i]);
            h_kari_ptinv_errhits.push_back(h);
        }

        //kari pt r dca dependent error histograms
        std::vector<TH1F*> h_kari_ptinv_err_rdca;
        for(int i=0; i<4; i++) {
            std::string histtitle = "p_{kari}^{-1} #minus p_{mc}^{-1} for r-dca_{kari} " + radialintervals[i] + " cm";
            std::string histhandle = "h_k_pinv_r_" + radialintervals[i] + "cm";
            TH1F *h = new TH1F(histhandle.c_str(), histtitle.c_str(), BIN_COUNT, -0.0005, 0.0005);
            labelAxis(h, "p_{rec}^{-1} #minus p_{mc}^{-1} [Mev^{-1}]", "count");
            fillHistWithVector(h, pt_kari_inv_err_r_dcas[i]);
            h_kari_ptinv_err_rdca.push_back(h);
        }

        //kari pt z dca dependent error histograms
        std::vector<TH1F*> h_kari_ptinv_err_zdca;
        for(int i=0; i<4; i++) {
            std::string histtitle = "p_{kari}^{-1} #minus p_{mc}^{-1} for abs(z-dca_{kari}) " + zdcaintervals[i] + " cm";
            std::string histhandle = "h_k_pinv_z_" + zdcaintervals[i] + "cm";
            TH1F *h = new TH1F(histhandle.c_str(), histtitle.c_str(), BIN_COUNT, -0.0005, 0.0005);
            labelAxis(h, "p_{kari}^{-1} #minus p_{mc}^{-1} [Mev^{-1}]", "count");
            fillHistWithVector(h, pt_kari_inv_err_z_dcas[i]);
            h_kari_ptinv_err_zdca.push_back(h);
        }

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

        TGraph *g_pdev_p = new TGraph(p_rel_errors.size(), &rec_ps[0],&p_inv_rel_errors[0]);
        g_pdev_p->SetTitle("p_{rec}^{-1} #minus p_{mc}^{-1} over p_{rec} correlation");
        labelAxis(g_pdev_p, "p_{rec} [MeV^{-1}]","p_{rec}^{-1} #minus p_{mc}^{-1} [MeV^{-1}]");
        setGraphRange(g_pdev_p, -1.5e5, 1.5e5, -4e-4, 4e-4);

        TGraph *g_prec_pmc = new TGraph(rec_ps.size(),&mc_p_corrs[0],&rec_ps[0]);
        g_prec_pmc->SetTitle("p_{rec} over p_{mc} correlation (p_{mc} charge corr)");
        labelAxis(g_prec_pmc, "p_{mc} * charge_{mc} [MeV]", "p_{rec} [MeV]");
        setGraphRange(g_prec_pmc, -1.5e5, 150e3, -1.5e5, 1.5e5);

        TGraph *g_prec_rrec = new TGraph(rec_ps.size(),&rec_ps[0],&rec_rs[0]);
        g_prec_rrec->SetTitle("r_{rec} over p_{rec} correlation");
        labelAxis(g_prec_rrec, "p_{rec} [MeV]", "r_{rec} [MeV]");
        setGraphRange(g_prec_rrec, -4e6,4e6,-10e6, 10e6);

        TGraph *g_invpt_invptmc = new TGraph(rec_inv_pts.size(),&mc_inv_pts[0],&rec_inv_pts[0]);
        g_invpt_invptmc->SetTitle("pt_{rec}^{-1} over pt_{mc}^{-1} correlation");
        labelAxis(g_invpt_invptmc, "pt_{mc}^{-1} [MeV^{-1}]", "pt_{rec}^{-1} [MeV^{-1}]");
        setGraphRange(g_invpt_invptmc, -5e-4,5e-4,-1e-3, 1e-3);

        TGraph *g_ptdev_ptmc = new TGraph(pt_inv_errors.size(),&mc_inv_pts[0],&pt_inv_errors[0]);
        g_ptdev_ptmc->SetTitle("pt_{rec}^{-1} #minus pt_{mc}^{-1} over pt_{mc}^{-1} correlation");
        labelAxis(g_ptdev_ptmc, "pt_{mc}^{-1} [MeV^{-1}]", "pt_{rec}^{-1} #minus pt_{mc}^{-1} [MeV^{-1}]");
        setGraphRange(g_ptdev_ptmc,-5e-4,5e-4, -1e-3,1e-3);

        //phi ms vs phi mc
        TGraph *g_phi_msvsmc = new TGraph(mc_phis.size(),&mc_phis[0],&rec_phis[0]);
        g_phi_msvsmc->SetTitle("#Phi_{msfit} vs #Phi_{mc} correlation");
        labelAxis(g_phi_msvsmc, "#Phi_{mc}", "#Phi_{msfit}");
        setGraphRange(g_phi_msvsmc,-3.2,0, -3.2,0);

        //KARI Scatter plots

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

        //theta kari vs rec theta
        TGraph *g_theta_kvsms = new TGraph(rec_thetas.size(),&rec_thetas[0],&kari_thetas[0]);
        g_theta_kvsms->SetTitle("#Theta_{kari} vs #Theta_{msfit} correlation");
        labelAxis(g_theta_kvsms, "#Theta_{msfit}", "#Theta_{kari}");
        setGraphRange(g_theta_kvsms,0,3.2, 0,3.2);

        //phi kari vs phi mc
        TGraph *g_phi_kvsmc = new TGraph(mc_phis.size(),&mc_phis[0],&kari_phis[0]);
        g_phi_kvsmc->SetTitle("#Phi_{kari} vs #Phi_{mc} correlation");
        labelAxis(g_phi_kvsmc, "#Phi_{mc}", "#Phi_{kari}");
        setGraphRange(g_phi_kvsmc,-3.2,0, -3.2,0);

        //phi kari vs phi ms
        TGraph *g_phi_kvsms = new TGraph(rec_phis.size(),&rec_phis[0],&kari_phis[0]);
        g_phi_kvsms->SetTitle("#Phi_{kari} vs #Phi_{msfit} correlation");
        labelAxis(g_phi_kvsms, "#Phi_{msfit}", "#Phi_{kari}");
        setGraphRange(g_phi_kvsms,-5,0, -3.2,0);

        //phi kari vs phi kari
        TGraph *g_phi_kvsk = new TGraph(rec_phis.size(),&kari_phis[0],&kari_phis[0]);
        g_phi_kvsk->SetTitle("#Phi_{kari} vs #Phi_{kari} correlation");
        labelAxis(g_phi_kvsk, "#Phi_{kari}", "#Phi_{kari}");
        setGraphRange(g_phi_kvsk,-3.2,0, -3.2,0);

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

        //pterr kari - rec over pt rec
        TGraph *g_ptdevkvsms_ptms = new TGraph(pt_inv_kvsms_errors.size(),&rec_inv_pts[0],&pt_inv_kvsms_errors[0]);
        g_ptdevkvsms_ptms->SetTitle("pt_{kari}^{-1} #minus pt_{msfit}^{-1} over #pt_{msfit}^{-1} correlation");
        labelAxis(g_ptdevkvsms_ptms, "#pt_{msfit}^{-1} [MeV^{-1}]", "pt_{kari}^{-1} #minus pt_{msfit}^{-1} [MeV^{-1}]");
        setGraphRange(g_ptdevkvsms_ptms,-5e-4,5e-4, -1e-3,1e-3);


        ////DCA, PHI, 1/PT CORR PLOTS

        // MC pt vs phi
        TGraph *g_mcpt_phi = new TGraph(mc_inv_pts.size(),&mc_phis[0],&mc_inv_pts[0]);
        g_mcpt_phi->SetTitle("pt_{mc}^{-1} over #Phi_{mc} correlation");
        labelAxis(g_mcpt_phi, " #Phi_{mc}", "pt_{mc}^{-1} [MeV^{-1}]");
        setGraphRange(g_mcpt_phi,-3.2, 3.2, -1e-3,1e-3);
        // REC pt vs phi
        TGraph *g_recpt_phi = new TGraph(rec_inv_pts.size(),&rec_phis[0],&rec_inv_pts[0]);
        g_recpt_phi->SetTitle("pt_{rec}^{-1} over #Phi_{rec} correlation");
        labelAxis(g_recpt_phi, " #Phi_{rec}", "pt_{rec}^{-1} [MeV^{-1}]");
        setGraphRange(g_recpt_phi,-3.2, 3.2, -1e-3,1e-3);
        // KARI pt vs phi
        TGraph *g_kaript_phi = new TGraph(kari_inv_pts.size(),&kari_phis[0],&kari_inv_pts[0]);
        g_kaript_phi->SetTitle("pt_{kari}^{-1} over #Phi_{kari} correlation");
        labelAxis(g_kaript_phi, " #Phi_{kari}", "pt_{kari}^{-1} [MeV^{-1}]");
        setGraphRange(g_kaript_phi,-3.2, 3.2, -1e-3,1e-3);

        // MC pt vs zdca
        TGraph *g_mcpt_zdca = new TGraph(mc_z_dcas.size(),&mc_z_dcas[0],&mc_inv_pts[0]);
        g_mcpt_zdca->SetTitle("pt_{mc}^{-1} over z-dca_{mc} correlation");
        labelAxis(g_mcpt_zdca, "z-dca_{mc}", "pt_{mc}^{-1} [MeV^{-1}]");
        setGraphRange(g_mcpt_zdca,-500, 500, -1e-3,1e-3);
        // REC pt vs zdca
        TGraph *g_recpt_zdca = new TGraph(rec_dca_zs.size(),&rec_dca_zs[0],&rec_inv_pts[0]);
        g_recpt_zdca->SetTitle("pt_{rec}^{-1} over z-dca_{rec} correlation");
        labelAxis(g_recpt_zdca, "z-dca_{rec}", "pt_{rec}^{-1} [MeV^{-1}]");
        setGraphRange(g_recpt_zdca,-500, 500, -1e-3,1e-3);
        // KARI pt vs zdca
        TGraph *g_kaript_zdca = new TGraph(kari_z0s.size(),&kari_z0s[0],&kari_inv_pts[0]);
        g_kaript_zdca->SetTitle("pt_{kari}^{-1} over z0_{kari} correlation");
        labelAxis(g_kaript_zdca, "z0_{kari}", "pt_{kari}^{-1} [MeV^{-1}]");
        setGraphRange(g_kaript_zdca,-500, 500, -1e-3,1e-3);

        // MC pt vs rdca
        TGraph *g_mcpt_rdca = new TGraph(mc_dcas.size(),&mc_dcas[0],&mc_inv_pts[0]);
        g_mcpt_rdca->SetTitle("pt_{mc}^{-1} over dca_{mc} correlation");
        labelAxis(g_mcpt_rdca, "dca_{mc}", "pt_{mc}^{-1} [MeV^{-1}]");
        setGraphRange(g_mcpt_rdca,-80, 80, -1e-3,1e-3);
        // REC pt vs rdca
        TGraph *g_recpt_rdca = new TGraph(rec_dca_rs.size(),&rec_dca_rs[0],&rec_inv_pts[0]);
        g_recpt_rdca->SetTitle("pt_{rec}^{-1} over r-dca_{rec} correlation");
        labelAxis(g_recpt_rdca, "r-dca_{rec}", "pt_{rec}^{-1} [MeV^{-1}]");
        setGraphRange(g_recpt_rdca,-80, 80, -1e-3,1e-3);
        // KARI pt vs rdca
        TGraph *g_kaript_rdca = new TGraph(kari_dcas.size(),&kari_dcas[0],&kari_inv_pts[0]);
        g_kaript_rdca->SetTitle("pt_{kari}^{-1} over dca_{kari} correlation");
        labelAxis(g_kaript_rdca, "dca_{kari}", "pt_{kari}^{-1} [MeV^{-1}]");
        setGraphRange(g_kaript_rdca,-500, 500, -1e-3,1e-3);

        // MC phi vs zdca
        TGraph *g_mcphi_zdca = new TGraph(mc_z_dcas.size(),&mc_z_dcas[0],&mc_phis[0]);
        g_mcphi_zdca->SetTitle("#phi_{mc} over z-dca_{mc} correlation");
        labelAxis(g_mcphi_zdca, "z-dca_{mc}", "#phi_{mc}");
        setGraphRange(g_mcphi_zdca,-500, 500, -3.2,3.2);
        // REC phi vs zdca
        TGraph *g_recphi_zdca = new TGraph(rec_dca_zs.size(),&rec_dca_zs[0],&rec_phis[0]);
        g_recphi_zdca->SetTitle("#phi_{rec} over z-dca_{rec} correlation");
        labelAxis(g_recphi_zdca, "z-dca_{rec}", "#phi_{rec}");
        setGraphRange(g_recphi_zdca,-500, 500, -3.2,3.2);
        // KARI phi vs zdca
        TGraph *g_kariphi_zdca = new TGraph(kari_z0s.size(),&kari_z0s[0],&kari_phis[0]);
        g_kariphi_zdca->SetTitle("#phi_{kari} over z0_{kari} correlation");
        labelAxis(g_kariphi_zdca, "z0_{kari}", "#phi_{kari}");
        setGraphRange(g_kariphi_zdca,-500, 500, -3.2,3.2);

        // MC phi vs rdca
        TGraph *g_mcphi_rdca = new TGraph(mc_dcas.size(),&mc_dcas[0],&mc_phis[0]);
        g_mcphi_rdca->SetTitle("#phi_{mc} over dca_{mc} correlation");
        labelAxis(g_mcphi_rdca, "dca_{mc}", "#phi_{mc}");
        setGraphRange(g_mcphi_rdca,-80, 80,  -3.2,3.2);
        // REC phi vs rdca
        TGraph *g_recphi_rdca = new TGraph(rec_dca_rs.size(),&rec_dca_rs[0],&rec_phis[0]);
        g_recphi_rdca->SetTitle("#phi_{rec} over r-dca_{rec} correlation");
        labelAxis(g_recphi_rdca, "r-dca_{rec}", "#phi_{rec}");
        setGraphRange(g_recphi_rdca,-80, 80,  -3.2,3.2);
        // KARI phi vs rdca
        TGraph *g_kariphi_rdca = new TGraph(kari_dcas.size(),&kari_dcas[0],&kari_phis[0]);
        g_kariphi_rdca->SetTitle("#phi_{kari} over dca_{kari} correlation");
        labelAxis(g_kariphi_rdca, "dca_{kari}", "#phi_{kari}");
        setGraphRange(g_kariphi_rdca,-500, 500, -3.2,3.2);


        //DCAz kari vs rec
        TGraph *g_karizdca_mczdca = new TGraph(kari_z0s.size(),&rec_dca_zs[0],&kari_z0s[0]);
        g_karizdca_mczdca->SetTitle("z-dca_{kari} over z-dca_{rec}  correlation");
        labelAxis(g_karizdca_mczdca, "z-dca_{rec}", "z-dca_{kari}");
        setGraphRange(g_karizdca_mczdca,-600, 600, -600, 600);
        //DCAr kari vs rec
        TGraph *g_karirdca_mcrdca = new TGraph(kari_dcas.size(),&rec_dca_rs[0],&kari_dcas[0]);
        g_karirdca_mcrdca->SetTitle("r-dca_{kari} over r-dca_{rec}  correlation");
        labelAxis(g_karirdca_mcrdca, "r-dca_{rec}", "r-dca_{kari}");
        setGraphRange(g_karirdca_mcrdca,-80, 80, -80, 80);
        //phi kari vs rec
        TGraph *g_kariphi_mcphi = new TGraph(kari_phis.size(),&rec_phis[0],&kari_phis[0]);
        g_kariphi_mcphi->SetTitle("#phi_{kari} over #phi_{rec}  correlation");
        labelAxis(g_kariphi_mcphi, "#phi_{rec}", "#phi_{kari}");
        setGraphRange(g_kariphi_mcphi,-3.2, 3.2,  -3.2,3.2);

        //rdca rec vs phi kari
        TGraph *g_recrdca_kariphi = new TGraph(kari_phis.size(),&kari_phis[0],&rec_dca_rs[0]);
        g_recrdca_kariphi->SetTitle("r-dca_{rec} over #phi_{kari}  correlation");
        labelAxis(g_recrdca_kariphi, "#phi_{kari}", "r-dca_{rec} [mm]");
        setGraphRange(g_recrdca_kariphi,-3.2, 3.2,  -80,80);

        //rdca kari vs phi kari
        TGraph *g_karirdca_kariphi = new TGraph(kari_phis.size(),&kari_phis[0],&kari_dcas[0]);
        g_karirdca_kariphi->SetTitle("dca_{kari} over #phi_{kari}  correlation");
        labelAxis(g_karirdca_kariphi, "#phi_{kari}", "dca_{kari} [mm]");
        setGraphRange(g_karirdca_kariphi,-3.2, 3.2,  -80,80);

        ////SOME MORE

        // x-axis p_rec / p_mc plots
        //vlt nicht nötig
        TGraph *g_ppmc_phi = new TGraph(p_over_pmcs.size(),&p_over_pmcs[0],&rec_phis[0]);
        g_ppmc_phi->SetTitle("#phi over p_{rec}/p_{mc} correlation");
        labelAxis(g_ppmc_phi, "#frac{p_{rec}}{p_{mc}}", "#phi_{rec}");
        setGraphRange(g_ppmc_phi, 0,2,-3, 0);

        TGraph *g_ppmc_dca = new TGraph(p_over_pmcs.size(),&p_over_pmcs[0],&rec_dca_rs[0]);
        g_ppmc_dca->SetTitle("r-dca_{rec} over p_{rec}/p_{mc} correlation");
        labelAxis(g_ppmc_dca, "#frac{p_{rec}}{p_{mc}}", "r-dca_{rec} [mm]");
        setGraphRange(g_ppmc_dca, 0,2,0, 60);

        TGraph *g_ppmc_p = new TGraph(p_over_pmcs.size(),&p_over_pmcs[0],&rec_ps[0]);
        g_ppmc_p->SetTitle("p_{rec} over p_{rec}/p_{mc} correlation");
        labelAxis(g_ppmc_p, "#frac{p_{rec}}{p_{mc}}", "p_{rec} [mm]");
        setGraphRange(g_ppmc_p, 0,2,-3e4, 3e4);

        //Nhits over rec r-dca
        TGraph *g_nhits_dca = new TGraph(rec_nhits_float.size(),&rec_nhits_float[0],&rec_dca_rs[0]);
        g_nhits_dca->SetTitle("nhits over r-dca_{rec} correlation");
        labelAxis(g_nhits_dca, "nhits", "r-dca_{rec} [mm]");
        setGraphRange(g_nhits_dca, 0,14,0, 60);


        ///######### PLOTTING PLOTS.PDF FILE ################################################

        //overview canvas containing 12 histograms
        TH1F * ov_hists[12];
        ov_hists[0] = h_p; ov_hists[1] = h_pmc; ov_hists[2] = h_pt; ov_hists[3] = h_ptmc;
        ov_hists[4] = h_dcamc; ov_hists[5] = h_phi; ov_hists[6] = h_nhits; ov_hists[7] = h_rec_nhits;
        ov_hists[8] = h_theta; ov_hists[9] = h_pinv_relerror; ov_hists[10] = h_poverpmc; ov_hists[11] = h_invptrec;
        makeSimpleMultiCanvas(3,4,12,ov_hists, true, false, plottingfile + "(");

        //p_rec and p_mc rel deviations 1x2 canvas
        auto *c_multi2 = new TCanvas("cmulti3", "cmulti3", 1200, 600);
        c_multi2->SetWindowPosition(0, 400);

        c_multi2->Divide(2,1);
        {
            c_multi2->cd(1);
            gPad->SetLeftMargin(0.15);
            g_prec_pmc->Draw("AP");

            c_multi2->cd(2);
            gPad->SetLeftMargin(0.15);
            g_pdev_p->Draw("AP");

            gPad->Update();
            TLine *l1=new TLine(gPad->GetUxmin(),0,gPad->GetUxmax(),0);
            l1->SetLineColor(kBlue);
            l1->Draw();

            c_multi2->Print(plottingfile.c_str(), "pdf");

        }

        //pt scatter plots (1x2 canvas)
        auto *c_multi5 = new TCanvas("cmulti5", "cmulti5", 900, 600);
        c_multi5->SetWindowPosition(0, 400);

        c_multi5->Divide(2,1);
        {
            c_multi5->cd(1);
            gPad->SetLeftMargin(0.15);
            g_invpt_invptmc->Draw("AP");

            gPad->Update();
            TLine *l1=new TLine(gPad->GetUxmin(),0,gPad->GetUxmax(),0);
            l1->SetLineColor(kBlue);
            l1->Draw();

            c_multi5->cd(2);
            gPad->SetLeftMargin(0.15);
            g_ptdev_ptmc->Draw("AP");

            c_multi5->Print(plottingfile.c_str(), "pdf");
        }

        TGraph* graph7[3] = {g_pdev_z_dca, g_pdev_dca, g_pdev_phi};
        makeSimpleMultiCanvas(1, 3, 3, graph7, plottingfile);

        //pt hist plots (1x3 canvas)
        TH1F * graph5[3] = {h_pt_inv_err, h_invptrec, h_invptmc};
        makeSimpleMultiCanvas(1,3,3,graph5, plottingfile);

        if(MAKE_468HIT_HISTOGRAMS) {
            //p error histograms for nhit=4,6,8 (1x3 canvas)
            makeSimpleMultiCanvas(1,3,3, &h_pinv_errhits[0], plottingfile);

            //pt error histograms for nhit=4,6,8 (1x3 canvas)
            makeSimpleMultiCanvas(1,3,3, &h_ptinv_errhits[0], plottingfile);

            //rdca error histograms for nhit=4,6,8 (1x3 canvas)
            makeSimpleMultiCanvas(1,3,3, &h_rdca_errhits[0], plottingfile);
        }

        if(MAKE_DCA_AREA_HISTOGRAMS) {
            //p error histograms for rdca intervals (2x2 canvas)
            makeSimpleMultiCanvas(2,2,4, &h_pinv_err_rdca[0], plottingfile);

            //p error histograms for zdca intervals (2x2 canvas)
            makeSimpleMultiCanvas(2,2,4, &h_pinv_err_zdca[0], plottingfile);
        }

        //rec dca plots (2x2 canvas)
        TH1F * graph6[4] = {h_xdca, h_ydca, h_zdca, h_rdca};
        makeSimpleMultiCanvas(2,2,4,graph6, plottingfile);

        if(MAKE_P_REL_PLOTS) {
            //p/p_mc plots with ref line (1x3 canvas)
            auto *c_multi3 = new TCanvas("cmulti2", "cmulti2", 1200, 600);
            c_multi3->SetWindowPosition(0, 400);

            c_multi3->Divide(3,1);
            {
                c_multi3->cd(1);
                gPad->SetLeftMargin(0.15);
                g_ppmc_phi->Draw("AP");

                gPad->Update();
                TLine *l=new TLine(1,gPad->GetUymin(),1,gPad->GetUymax());
                l->SetLineColor(kBlue);
                l->Draw();

                c_multi3->cd(2);
                gPad->SetLeftMargin(0.15);
                g_ppmc_dca->Draw("AP");

                gPad->Update();
                TLine *l1=new TLine(1,gPad->GetUymin(),1,gPad->GetUymax());
                l1->SetLineColor(kBlue);
                l1->Draw();

                c_multi3->cd(3);
                gPad->SetLeftMargin(0.15);
                g_ppmc_p->Draw("AP");

                gPad->Update();
                TLine *l2=new TLine(1,gPad->GetUymin(),1,gPad->GetUymax());
                l2->SetLineColor(kBlue);
                l2->Draw();

                c_multi3->Print(plottingfile.c_str(), "pdf");
            }
        }

        //karimaki histograms (2x2 canvas)
        TH1F* graph1[4] = {h_r3dkari, h_rinvkari, h_phikari, h_thetakari};
        makeSimpleMultiCanvas( 2, 2, 4, graph1, false, false, plottingfile);

        //karimaki histograms 2 (2x2 canvas)
        TH1F* graph2[4] = {h_dcakari, h_tchi2nkari, h_z0kari, h_zchi2kari};
        makeSimpleMultiCanvas( 2, 2, 4, graph2, false, false, plottingfile);

        // karimaki theta vs rec theta correlation plot
        auto *c_single4 = new TCanvas("cmulti4", "cmulti4", 900, 450);
        {
            c_single4->Divide(3,1);
            c_single4->cd(1);
            gPad->SetLeftMargin(0.15);
            h_thetakari->Draw();

            c_single4->cd(2);
            gPad->SetLeftMargin(0.15);
            h_theta->Draw();

            c_single4->cd(3);
            gPad->SetLeftMargin(0.15);
            g_theta_kvsms->Draw("AP");

            gPad->Update();
            TLine *l1=new TLine(0,0,PI,PI);
            TLine *l2=new TLine(0,PI,PI,0);
            l1->SetLineColor(kBlue);
            l2->SetLineColor(kBlue);
            l1->Draw();
            l2->Draw();

            c_single4->Print(plottingfile.c_str(), "pdf");
        }

        // karimaki perr vs kappa correlation plot
        auto *c_single5 = new TCanvas("cmulti5", "cmulti5", 900, 900);
        {
            c_single5->SetLeftMargin(0.15);
            g_ptkaridev_ptmc->Draw("AP");

            c_single5->Update();
            TLine *l1=new TLine(-0.25e-3,0.00045,0.25e-3,-0.00055);
            TLine *l2=new TLine(-0.25e-3,0.00055,0.25e-3,-0.00045);
            l1->SetLineColor(kBlue);
            l2->SetLineColor(kBlue);
            l1->Draw();
            l2->Draw();

            c_single5->Print(plottingfile.c_str(), "pdf");
        }

        // phi kari vs mc and ms vs mc
        TGraph * graph10[3] = {g_phi_kvsmc, g_phi_msvsmc};
        makeSimpleMultiCanvas( 1, 2, 2, graph10, plottingfile);

        //perr kari vs phi mc
        TGraph * graph11[2] = {g_phi_kvsms, g_ptkaridev_mcphi};
        makeSimpleMultiCanvas(1,2,2,graph11, plottingfile);

        //karimaki pt histograms 3 (3x1 canvas)
        TH1F* graph3[3] = {h_ptkari, h_ptkari_inv, h_pterr_kari};
        makeSimpleMultiCanvas(1, 3, 3, graph3, false, false, plottingfile);

        //karimaki scatter 1 (2x2 canvas)
        TGraph* graph4[4] = {g_ptkari_ptmc, g_ptkaridev_ptmc, g_invpt_invptmc, g_ptdev_ptmc};
        makeSimpleMultiCanvas(2, 2, 4, graph4, false, false, plottingfile);

        //karimaki nhit histograms
        if(MAKE_468HIT_HISTOGRAMS) makeSimpleMultiCanvas(1,3,3, &h_kari_ptinv_errhits[0], false, false, plottingfile);

        //karimaki r dca histograms
        if (MAKE_DCA_AREA_HISTOGRAMS) makeSimpleMultiCanvas(2, 2, 4, &h_kari_ptinv_err_rdca[0], false, false, plottingfile);

        //karimaki z dca histograms
        if( MAKE_DCA_AREA_HISTOGRAMS) makeSimpleMultiCanvas(2, 2, 4, &h_kari_ptinv_err_zdca[0], false, false, plottingfile);

        //karimaki pterr over rdca and zdca scatter plots (1x2 canvas)
        TGraph* graph9[2] = {g_ptkaridev_zdca, g_ptkaridev_rdca};
        makeSimpleMultiCanvas(1, 2, 2, graph9, plottingfile);

        ////DCA, PHI, 1/PT CORR PLOTS
        if(MAKE_STRANGE_BAND_CORR_PLOTS) {
            TGraph * graphs1[6];
            TH1F * hists1[6];

            graphs1[0] = g_karirdca_mcrdca; graphs1[1] = g_karizdca_mczdca;
            makeSimpleMultiCanvas(1, 2, 2, graphs1, plottingfile);

            graphs1[0] = g_recrdca_kariphi; graphs1[1] = g_karirdca_kariphi;
            makeSimpleMultiCanvas(1, 2, 2, graphs1, plottingfile);

            graphs1[0] = g_theta_kvsms; graphs1[1] = g_phi_kvsms;
            makeSimpleMultiCanvas(1, 2, 2, graphs1, plottingfile);

            hists1[0] = h_dcamc; hists1[1] = h_rdca; hists1[2] = h_dcakari;
            makeSimpleMultiCanvas(1, 3, 3, hists1, plottingfile);

            hists1[0] = h_zdcamc; hists1[1] = h_zdca; hists1[2] = h_z0kari;
            makeSimpleMultiCanvas(1, 3, 3, hists1, plottingfile);

            // pt vs phi rec and mc
            graphs1[0] = g_mcpt_phi; graphs1[1]= g_recpt_phi; graphs1[2]= g_kaript_phi;
            makeSimpleMultiCanvas(1, 3, 3, graphs1, plottingfile);

            // pt vs zdca and rdca rec and mc
            graphs1[0] = g_mcpt_rdca; graphs1[1]= g_recpt_rdca; graphs1[2]= g_kaript_rdca;
            graphs1[3] = g_mcpt_zdca; graphs1[4]= g_recpt_zdca; graphs1[5]= g_kaript_zdca;
            makeSimpleMultiCanvas(2, 3, 6, graphs1, plottingfile);

            // phi vs zdca and rdca rec and mc
            graphs1[0] = g_mcphi_rdca; graphs1[1]= g_recphi_rdca; graphs1[2]= g_kariphi_rdca;
            graphs1[3] = g_mcphi_zdca; graphs1[4]= g_recphi_zdca; graphs1[5]= g_kariphi_zdca;
            makeSimpleMultiCanvas(2, 3, 6, graphs1, plottingfile);

            makeSimpleSingleCanvas(g_ptdevkvsms_ptms, plottingfile);
        }

        makeSimpleSingleCanvas(h_mctype, plottingfile);

        g_nhits_dca->SetMarkerColor(4);
        makeSimpleSingleCanvas(g_nhits_dca, plottingfile);

        makeSimpleSingleCanvas(h_phimc, plottingfile);

        //Single histogram for different nhits
        auto  *c_single10 = new TCanvas();
        {
            c_single10->SetWindowPosition(0, 500 );
            c_single10->SetLogy(1);
            auto legend1 = new TLegend(0.2,0.2);

            h_pinv_relerror->SetFillStyle(3001);
            h_pinv_relerror->SetFillColor(16);
            h_pinv_relerror->SetLineColor(4);
            h_pinv_relerror->SetName("p_{rec}^{-1} reconstruction error");

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
            legend1->Draw();

            c_single10->Print(plottingfile.c_str(), "pdf");

//            filename = filename_template + "_hist_perr_nhits.pdf";
//            c_single10->SaveAs(filename.c_str());
        }

        //empty plot closes _plots.pdf file
        auto *c_final = new TCanvas("c_final", "c_final");
        filename = plottingfile + ")";
        c_final->Print(filename.c_str(), "pdf");


        ////###### COMPARISON PLOTS ##################################################

        TGraph * graphs[4];
        TH1F * hists[4];

        graphs[0] = g_ptdev_ptmc; graphs[1]= g_ptkaridev_ptmc;
        makeSimpleMultiCanvas(1, 2, 2, graphs, plottingfilecomparison + "(");

        graphs[0] = g_invpt_invptmc; graphs[1]= g_ptkari_ptmc;
        makeSimpleMultiCanvas(1, 2, 2, graphs, plottingfilecomparison);

        hists[0] = h_pt_inv_err; hists[1]= h_pterr_kari;
        makeSimpleMultiCanvas(1, 2, 2, hists, plottingfilecomparison);

        hists[0] =h_invptmc; hists[1]= h_invptrec; hists[2]= h_ptkari_inv;
        makeSimpleMultiCanvas(1, 3, 3, hists, plottingfilecomparison);

        hists[0] = h_rdca; hists[1]= h_dcakari;
        makeSimpleMultiCanvas(1, 2, 2, hists, plottingfilecomparison);

        graphs[0] = g_pdev_z_dca; graphs[1]= g_ptkaridev_zdca;
        makeSimpleMultiCanvas(1, 2, 2, graphs, plottingfilecomparison);

        graphs[0] = g_pdev_dca; graphs[1]= g_ptkaridev_rdca;
        makeSimpleMultiCanvas(1, 2, 2, graphs, plottingfilecomparison);

        hists[0] = h_theta; hists[1]= h_thetakari;
        makeSimpleMultiCanvas(1, 2, 2, hists, plottingfilecomparison + ")");


        ////###########################################################################

    }

    ////PRINT THE STATS
//    for(int i= 0; i<mc_types.size(); i++) {
//        cout << mc_types[i] << " ";
//    }

    cout << endl << endl << "---General Stats---\n" << endl;
    cout << "Filter set to : \""<< filtertag << "\", " << processed_entries << " entries of " << Segs.my_entries << " processed." << endl << endl;
    cout << "p fail analysis:" << endl;
    cout << "(p_mc == 0 || pt_mc == 0 || p_rec == 0 || pt_rec == 0) --- count: " << p_fail_count ;
    cout << " fails of total: " << processed_entries << ", fail rate: ";
    cout << (p_fail_count / (float) Segs.my_entries) * 100 << " %" << endl;
    printf("MC vs Kari sign fails: %d of %d\n", mckari_wrong_sign_count, processed_entries);
    printf("MC vs REC sign fails: %d of %d\n", mcrec_wrong_sign_count, processed_entries);
    printf("KARI vs REC sign fails: %d of %d\n", reckari_wrong_sign_count, processed_entries);
    //cout << "1 / p reconstruction error mean: "  << p_inv_rec_error_mean * 100 << "% " << endl;

    for(int i=0; i < 6; i++) {
        if(i < 5) {
            cout << i+4 << " hits count: " << p_inv_rel_errors_hits[i].size() << endl;
        } else {
            cout << "other: " << p_inv_rel_errors_hits[i].size() << endl;
        }
    }

}

//
// Created by Konstantin Neureither on 27.04.20.
//
#define TRIPLET_HIT_ARRAY_LENGTH 1024

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
#include <fstream>
#include <filesystem>
#include <cmath>
#include "../util/utility_functions.h"
#include "../util/plots.h"
#include "reconstruction_accuracy.h"

using std::cout;
using std::endl;


void reconstruction_accuracy() {

    const std::string pathtodata = "../data/";
    const std::string pathtoplots = "../plots/";
    const int run = 10;
    const bool RECONSTRUCTION_PRINTS = false;
    const bool HIT_PRINTS = true;
    const bool MAKE_PLOT = true;
    const bool ADDITIONAL_PLOTS = false;
    const int MAX_ENTRIES = 0;

    //TODO Get run number as cmd argument

    check_create_directory(pathtodata);
    check_create_directory(pathtoplots);
    std::string pathtorunplots = pathtoplots + "run_" + get_padded_string(run, 3, '0') + "/";
    check_create_directory(pathtorunplots);

    std::string runpadded = get_padded_string(run, 6, '0');
    std::string infile1 = pathtodata + "mu3e_run_" + runpadded + ".root";
    std::string infile2 = pathtodata + "mu3e_run_" + runpadded + "_trirec_cosmic.root";

    TFile f1(infile1.c_str());
    TFile f2(infile2.c_str());

//    gStyle->SetLegendBorderSize(1);
//    gStyle->SetLegendFillColor(0);
//    gStyle->SetLegendFont(42);
//    gStyle->SetLegendTextSize(0.03);

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

    t_mu3e->SetBranchAddress("hit_pixelid",&pixelid);
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
    int mc_type;
    int mc_pid;

    float rec_p;
    float rec_r;
    float rec_tan01[TRIPLET_HIT_ARRAY_LENGTH];
    float rec_lam01[TRIPLET_HIT_ARRAY_LENGTH];
    float rec_zpca_x;
    float rec_zpca_y;
    float rec_zpca_z;
    float rec_zpca_r;

    float x00[1024];
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

    //for calculation
    float rec_phi;
    float rec_theta;
    float rec_pt;

    t_segs->SetBranchAddress("eventId", &rec_event);
    t_segs->SetBranchAddress("nhit", &rec_nhit);

    t_segs->SetBranchAddress("p", &rec_p);
    t_segs->SetBranchAddress("r", &rec_r);
    t_segs->SetBranchAddress("tan01", &rec_tan01);
    t_segs->SetBranchAddress("lam01", &rec_lam01);
    t_segs->SetBranchAddress("zpca_z", &rec_zpca_z);
    t_segs->SetBranchAddress("zpca_x", &rec_zpca_x);
    t_segs->SetBranchAddress("zpca_y", &rec_zpca_y);
    t_segs->SetBranchAddress("zpca_r", &rec_zpca_r);

    t_segs->SetBranchAddress("mc_tid", &rec_trajid);
    t_segs->SetBranchAddress("mc_p", &mc_p);
    t_segs->SetBranchAddress("mc_pt", &mc_pt);
    t_segs->SetBranchAddress("mc_theta", &mc_theta);
    t_segs->SetBranchAddress("mc_phi", &mc_phi);
    t_segs->SetBranchAddress("mc_lam", &mc_lam);
    t_segs->SetBranchAddress("mc_vpca_offset", &mc_vpca_offset);
    t_segs->SetBranchAddress("mc_vpca_phi", &mc_vpca_phi);
    t_segs->SetBranchAddress("mc_type", &mc_type);
    t_segs->SetBranchAddress("mc_pid", &mc_pid);

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

    //data for frames tree
    unsigned int frames_entries = t_frames->GetEntries();
    //TODO Maybe switch from segs tree to frames tree
    cout << "Branches set for frames..." << endl;

    //stats data definitions
    int p_fail_count = 0;
    std::vector<float> p_rel_errors;
    std::vector<float> p_inv_rel_errors;
    std::vector<float> pt_inv_errors;
    std::vector<float> p_over_pmcs;
    std::vector<float> rec_rs;
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
    std::vector<float> p_inv_rel_errors_hits[6];
    std::vector<float> p_inv_err_nhits[3];
    std::vector<float> pt_inv_err_nhits[3];
    std::vector<float> rec_rdca_nhits[3];
    int rec_hits_count[6] = {0};

    unsigned int mu3e_index = 1;
    t_mu3e->GetEntry(mu3e_index);
    mu3e_index++;

//    for (unsigned int i = 0; i < segs_entries; i++) {
    for (unsigned int i = 1; i < (MAX_ENTRIES == 0 ? segs_entries : MAX_ENTRIES); i++) {
        t_segs->GetEntry(i);

        //find corresponding entry in mu3e tree
        while (rec_event > header[0]) {
            t_mu3e->GetEntry(mu3e_index);
            mu3e_index++;
        }

        if(HIT_PRINTS) {
            cout << "segs_index=" << i << "\tmu3e_index=" << mu3e_index;
            cout << "\trec_event=" << rec_event << "\tmu3e_event= " << header[0];
            cout << "\trec_nhits=" << rec_nhit << "\tmu3e_nhit=" << mu3e_nhits << endl;
        }

        for (unsigned int hit = 0; hit < (unsigned int) mu3e_nhits; hit++) {
            PXID pixid = process_pixel_id((*pixelid)[hit]);
            unsigned int trajid = (*trajids)[hit];
            double trajpx = (*traj_px)[hit];

            if(HIT_PRINTS) {
                unsigned int layer = get_layer(pixid.sensor);
                cout << "\thit=" << hit << "\tlayer=" << layer << "\tsid=" << pixid.sensor << "\trow=" << pixid.row << "\tcol=" << pixid.column << endl;
            }
        }

        if (HIT_PRINTS) {
            cout << "\tsegs data" << endl;
            cout << "\ttriplet 00: " << "\tsid=" << sid00[0] << "\tx=" << x00[0] << "\tx=" << x00[1] << "\tx=" << x00[4] << endl;
            cout << "\ttriplet 10: " << "\tsid=" << sid10[0] << "\tx=" << x10[0] << endl;
            cout << "\ttriplet 20: " << "\tsid=" << sid20[0] << "\tx=" << x20[0] << endl;
            cout << "\ttriplet 01: " << "\tsid=" << sid01[0] << "\tx=" << x01[0] << endl;
            cout << "\ttriplet 11: " << "\tsid=" << sid11[0] << "\tx=" << x11[0] << endl;
            cout << "\ttriplet 21: " << "\tsid=" << sid21[0] << "\tx=" << x21[0] << endl;

            for(int i=0; i<rec_nhit-2; i++) {
                cout << "\tx=" << x00[i];
            }
            cout << endl;
        }

        //store hit data
        std::vector<float> xp;
        std::vector<float> yp;
        std::vector<float> zp;
        combinehits(xp, yp, zp, rec_nhit, x00, x11, y00, y11, z00, z11);

        if(HIT_PRINTS) {
            for(int i = 0; i < rec_nhit; i++) {
                cout << "Gathered hits in array: (" << xp[i] << ", " << yp[i] << ", " << zp[i] << ")" << endl;
            }
        }

        //get trajectory data from mu3e tree
        unsigned int trajid = (*trajids)[0];
        double trajpx = (*traj_px)[0];
        double trajpy = (*traj_py)[0];
        double trajpz = (*traj_pz)[0];
        double trajp = sqrt(pow(trajpx, 2) + pow(trajpy, 2) + pow(trajpz, 2));




        //Theta, phi and traverse p of reconstruction
        float rec_pt = rec_p * std::cos((rec_lam01)[0]);
        float rec_phi = (rec_tan01)[0];
        float rec_theta = 3.14159*0.5 - (rec_lam01)[0];

        if (mc_p == 0 || mc_pt == 0 || rec_p == 0 || rec_pt == 0) {
            p_fail_count++;
        } else {
            if (true) {

                //TODO Use direction provided by phi to correct sign of momentum
                //The charge of the particle is given by the sign of the traj radius (rec_r)
                //This is used to correct the p_mc_corr = sgn(r) * p_mc, because the monte carlo
                //momenta are only given as absolutes.

                float mc_p_corr = mc_p * sgn(mc_pid);
                float mc_pt_corr = mc_pt * sgn(mc_pid);
                float rec_p_corr = rec_p * sgn(rec_r); //TODO rec_p anders korrigieren (phi)
                float rec_pt_corr = rec_pt;

                float mc_inv_p = 1. / mc_p_corr;
                float mc_inv_pt = 1. / mc_pt_corr;
                float rec_inv_p = 1. / rec_p;
                float rec_inv_pt = 1. / rec_pt; //TODO switch to rec_pt_corr
                float p_inv_rel_error = (rec_inv_p - mc_inv_p); //this will mostly be used as estimator for deviation
                float p_rel_error = (mc_p_corr - rec_p);
                float p_over_pmc = rec_p / mc_p_corr;
                float mc_z_dca = std::sin(mc_vpca_phi) * mc_vpca_offset;
                float pt_inv_error = (1./rec_pt_corr) - (1./mc_pt_corr);

                //calculated data
                rec_inv_ps.push_back(rec_inv_p);
                rec_p_corrs.push_back(rec_p_corr);
                rec_inv_pts.push_back(rec_inv_p);
                p_inv_rel_errors.push_back(p_inv_rel_error);
                p_rel_errors.push_back(p_rel_error);
                pt_inv_errors.push_back(pt_inv_error);
                p_over_pmcs.push_back(p_over_pmc);
                rec_pts.push_back(rec_pt);
                rec_phis.push_back(rec_phi);
                rec_thetas.push_back(rec_theta);

                mc_p_corrs.push_back(mc_p_corr);
                mc_pt_corrs.push_back(mc_pt_corr);
                mc_inv_pts.push_back(mc_inv_pt);
                mc_z_dcas.push_back(mc_z_dca);
                mc_inv_ps.push_back(mc_inv_p);

                //reconstruction data
                rec_ps.push_back(rec_p);
                rec_rs.push_back(rec_r);
                rec_dca_rs.push_back(rec_zpca_r);
                rec_dca_xs.push_back(rec_zpca_x);
                rec_dca_ys.push_back(rec_zpca_y);
                rec_dca_zs.push_back(rec_zpca_z);

                //monte carlo data
                mc_ps.push_back(mc_p);
                mc_pts.push_back(mc_pt);
                mc_phis.push_back(mc_phi);
                mc_dcas.push_back(mc_vpca_offset);
                mc_phi_dcas.push_back(mc_vpca_phi);
                mc_types.push_back(mc_type);

                //meta data
                sim_nhits.push_back(mu3e_nhits);
                rec_nhits.push_back(rec_nhit);
                rec_nhits_float.push_back((float)rec_nhit);


                switch(rec_nhit) {
                    case 4:
                        p_inv_rel_errors_hits[0].push_back(p_inv_rel_error);
                        p_inv_err_nhits[0].push_back(p_inv_rel_error);
                        pt_inv_err_nhits[0].push_back(pt_inv_error);
                        rec_rdca_nhits[0].push_back(rec_zpca_r);
                        break;
                    case 5:
                        p_inv_rel_errors_hits[1].push_back(p_inv_rel_error);
                        p_inv_err_nhits[0].push_back(p_inv_rel_error);
                        pt_inv_err_nhits[0].push_back(pt_inv_error);
                        rec_rdca_nhits[0].push_back(rec_zpca_r);
                        break;
                    case 6:
                        p_inv_rel_errors_hits[2].push_back(p_inv_rel_error);
                        p_inv_err_nhits[1].push_back(p_inv_rel_error);
                        pt_inv_err_nhits[1].push_back(pt_inv_error);
                        rec_rdca_nhits[1].push_back(rec_zpca_r);
                        break;
                    case 7:
                        p_inv_rel_errors_hits[3].push_back(p_inv_rel_error);
                        p_inv_err_nhits[1].push_back(p_inv_rel_error);
                        pt_inv_err_nhits[1].push_back(pt_inv_error);
                        rec_rdca_nhits[1].push_back(rec_zpca_r);
                        break;
                    case 8:
                        p_inv_rel_errors_hits[4].push_back(p_inv_rel_error);
                        p_inv_err_nhits[2].push_back(p_inv_rel_error);
                        pt_inv_err_nhits[2].push_back(pt_inv_error);
                        rec_rdca_nhits[2].push_back(rec_zpca_r);
                        break;
                    default:
                        p_inv_rel_errors_hits[5].push_back(p_inv_rel_error);
                        p_inv_err_nhits[2].push_back(p_inv_rel_error);
                        pt_inv_err_nhits[2].push_back(pt_inv_error);
                        rec_rdca_nhits[2].push_back(rec_zpca_r);
                }

                ////PRINT SECTION PER ENTRY IN TREE
                if(RECONSTRUCTION_PRINTS) {
                    //data from reconstruction
                    cout << "rec_event: " << rec_event << " mu3e_event: " << header[0];
                    cout << " mc_p: " << mc_p << " rec_p: " << rec_p*sgn(rec_r) << "rec_nhit: "<< rec_nhit <<"\t\t";

                    cout << "rec_r= " << rec_r << " mc_type=" << mc_type;
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
        //momenta
        const float LEFT_BOUNDARY = -3e4;
        const float RIGHT_BOUNDARY = 3e4 ;
        const int BIN_COUNT = 50;

        std::string filename;
        std::string filename_template = pathtorunplots + "reconstruction-accuracy_run" +
                get_padded_string(run, 3, '0');
        std::string plottingfile = filename_template + "_plots.pdf";



        ///FILLING THE HISTOGRAMS

        TH1F *h_nhits = new TH1F("h_nhits", "hits per frame", 20, 0., 20.);
        labelAxis(h_nhits, "number of hits", "count");
        fillHistWithVector(h_nhits, sim_nhits);

        TH1F *h_rec_nhits = new TH1F("h_rec_nhits", "hits per frame of trirec output", 20, 0., 20.);
        labelAxis(h_rec_nhits, "number of hits", "count");
        fillHistWithVector(h_rec_nhits, rec_nhits);

        //Trajectory data Monte Carlo
        TH1F *h_pmc = new TH1F("h_pmc", "muon momenta monte carlo (charge corrrected)", 30, LEFT_BOUNDARY, RIGHT_BOUNDARY);
        labelAxis(h_pmc, "p_{mc} * charge_{mc} [MeV]", "count");
        fillHistWithVector(h_pmc, mc_p_corrs);

        TH1F *h_ptmc = new TH1F("h_ptmc", "muon traverse momenta monte carlo (charge corrrected)", 30, LEFT_BOUNDARY, RIGHT_BOUNDARY);
        labelAxis(h_ptmc, "p_{t-mc} * charge_{mc} [MeV]", "count");
        fillHistWithVector(h_ptmc, mc_pt_corrs);

        TH1F *h_phimc = new TH1F("h_phimc", "#Phi monte carlo (charge corrrected)", 30, -3.2, 1);
        labelAxis(h_phimc, "#Phi", "count");
        fillHistWithVector(h_phimc, mc_phis);

        //Trajectory data reconstructed
        TH1F *h_p = new TH1F("h_p", "muon momenta", 30, LEFT_BOUNDARY, RIGHT_BOUNDARY);
        labelAxis(h_p, "p_{rec} [MeV]", "count");
        fillHistWithVector(h_p, rec_ps);

        TH1F *h_p_corr = new TH1F("h_p_corr", "muon momenta (charge corrected)", 30, 0, RIGHT_BOUNDARY);
        labelAxis(h_p_corr, "p_{rec-corr} [MeV^{-1}]", "count");
        fillHistWithVector(h_p_corr, rec_p_corrs);

        TH1F *h_pt = new TH1F("h_pt", "muon transverse momentum", 30, LEFT_BOUNDARY, RIGHT_BOUNDARY);
        labelAxis(h_pt, "p_{t} [MeV]", "count");
        fillHistWithVector(h_pt, rec_pts);

        TH1F *h_phi = new TH1F("h_phi", "reconstruction #phi (var: tan01)", 30, -3.2, 3.2);
        labelAxis(h_phi, "#phi", "count");
        fillHistWithVector(h_phi, rec_phis);

        TH1F *h_theta = new TH1F("h_theta", "reconstruction #Theta", 30, -3.2, 3.2 );
        labelAxis(h_theta, "#Theta", "count");
        fillHistWithVector(h_theta, rec_thetas);

        //rec pt
        TH1F *h_pt_inv_err = new TH1F("h_pt_inv_err", "pt_{rec}^{-1} #minus pt_{mc}^{-1} histogram", 30, -5e-4,5e-4);
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
        TH1F *h_dcamc = new TH1F("h_dcamc", "dca monte carlo (var: mc_vpca_offset)", 30, 0, 100);
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

        TH1F *h_zdca = new TH1F("h_zdca", "dca_{reconstruction} along z-axis", 30, -150, 150);
        labelAxis(h_zdca, "dca_{z} [mm]", "count");
        fillHistWithVector(h_zdca, rec_dca_zs);

        //pt nhit dependend error histograms
        std::vector<TH1F*> h_ptinv_errhits;
        for(int i=0; i<3; i++) {
            std::string histtitle = "pt_{rec}^{-1} #minus pt_{mc}^{-1} for " + get_string(2*i+4) + " and " +
                    (i==2 ? " more" : get_string(2*i+5)) + " hits";
            std::string histhandle = "h_ptinv_" + get_string(2*i+4) + "hits";
            TH1F *h = new TH1F(histhandle.c_str(), histtitle.c_str(), BIN_COUNT, -0.0003, 0.0003);
            labelAxis(h, "pt_{rec}^{-1} #minus pt_{mc}^{-1} [Mev^{-1}]", "count");
            fillHistWithVector(h, pt_inv_err_nhits[i]);
            h_ptinv_errhits.push_back(h);
        }

        //r dca nhit dependend histograms
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

        //pt nhit dependend error histograms
        std::vector<TH1F*> h_pinv_errhits;
        for(int i=0; i<3; i++) {
            std::string histtitle = "p_{rec}^{-1} #minus p_{mc}^{-1} for " + get_string(2*i+4) + " and " +
                                    (i==2 ? " more" : get_string(2*i+5)) + " hits";
            std::string histhandle = "h_pinv_" + get_string(2*i+4) + "hits";
            TH1F *h = new TH1F(histhandle.c_str(), histtitle.c_str(), BIN_COUNT, -0.0003, 0.0003);
            labelAxis(h, "p_{rec}^{-1} #minus p_{mc}^{-1} [Mev^{-1}]", "count");
            fillHistWithVector(h, p_inv_err_nhits[i]);
            h_pinv_errhits.push_back(h);
        }

        //p nhit dependend error histograms (hits individually) (for overlay hist)
        std::vector<TH1F*> h_pinv_relerrorhits;
        for(int i = 0; i < 6; i++) {
            std::string histhandle = "h_pinv_relerror_hits_" + get_string(i);
            std::string histtitle = "p_{rec}^{-1} #minus p_{mc}^{-1} deviations, nhit=" + get_string(i+4);
            TH1F *h = new TH1F(histhandle.c_str(), histtitle.c_str(), BIN_COUNT, -0.0003, 0.0003);
            labelAxis(h, "p_{rec}^{-1} #minus p_{mc}^{-1} [MeV^{-1}]", "count");
            fillHistWithVector(h, p_inv_rel_errors_hits[i]);
            h_pinv_relerrorhits.push_back(h);
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
        g_prec_rrec->SetTitle("r_{rec} over p_{mc} correlation");
        labelAxis(g_prec_rrec, "p_{rec} [MeV]", "r_{rec} [MeV]");
        setGraphRange(g_prec_rrec, -4e6,4e6,-10e6, 10e6);

        TGraph *g_invpt_invptmc = new TGraph(rec_inv_pts.size(),&mc_inv_pts[0],&rec_inv_pts[0]);
        g_invpt_invptmc->SetTitle("pt_{rec}^{-1} over pt_{mc}^{-1} correlation");
        labelAxis(g_invpt_invptmc, "pt_{mc}^{-1} [MeV^{-1}]", "pt_{rec}^{-1} [MeV^{-1}]");
        setGraphRange(g_invpt_invptmc, -5e-4,5e-4,-1e-3, 1e-3);

        TGraph *g_ptdev_ptmc = new TGraph(pt_inv_errors.size(),&mc_inv_pts[0],&pt_inv_errors[0]);
        g_ptdev_ptmc->SetTitle("pt_{rec}^{-1} #minus pt_{mc}^{-1} over pt_{mc}^{-1} correlation");
        labelAxis(g_ptdev_ptmc, "pt_{mc}^{-1} [MeV^{-1}]", "pt_{rec}^{-1} #minus pt_{mc}^{-1} [MeV^{-1}]");
        setGraphRange(g_ptdev_ptmc,-5e-4,5e-4, -5e-4,5e-4);


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

        ///PLOTTING THE HISTOGRAMS

        //overview canvas containing 16 plots
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

            filename = plottingfile + "(";
            c_multi->Print(filename.c_str(), "pdf");
        }

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

        //pt hist plots (1x3 canvas)
        auto *c_multi6 = new TCanvas("cmulti6", "cmulti6", 1200, 600);
        c_multi6->SetWindowPosition(0, 400);

        c_multi6->Divide(3,1);
        {
            c_multi6->cd(1);
            gPad->SetLeftMargin(0.15);
            h_pt_inv_err->Draw();

            c_multi6->cd(2);
            gPad->SetLeftMargin(0.15);
            h_invptrec->Draw();

            c_multi6->cd(3);
            gPad->SetLeftMargin(0.15);
            h_invptmc->Draw();

            c_multi6->Print(plottingfile.c_str(), "pdf");
        }

        //p error histograms for nhit=4,6,8 (1x3 canvas)
        auto *c_multi9 = new TCanvas("cmulti9", "cmulti9", 1200, 600);
        c_multi9->SetWindowPosition(0, 400);

        c_multi9->Divide(3,1);
        {
            for(int i=0; i<3; i++) {
                c_multi9->cd(i+1);
                gPad->SetLeftMargin(0.15);
                h_pinv_errhits[i]->Draw();
            }
            c_multi9->Print(plottingfile.c_str(), "pdf");
        }

        //pt error histograms for nhit=4,6,8 (1x3 canvas)
        auto *c_multi7 = new TCanvas("cmulti7", "cmulti7", 1200, 600);
        c_multi7->SetWindowPosition(0, 400);

        c_multi7->Divide(3,1);
        {
            for(int i=0; i<3; i++) {
                c_multi7->cd(i+1);
                gPad->SetLeftMargin(0.15);
                h_ptinv_errhits[i]->Draw();
            }
            c_multi7->Print(plottingfile.c_str(), "pdf");
        }

        //rdca error histograms for nhit=4,6,8 (1x3 canvas)
        auto *c_multi8 = new TCanvas("cmulti8", "cmulti8", 1200, 600);
        c_multi8->SetWindowPosition(0, 400);

        c_multi8->Divide(3,1);
        {
            for(int i=0; i<3; i++) {
                c_multi8->cd(i+1);
                gPad->SetLeftMargin(0.15);
                h_rdca_errhits[i]->Draw();
            }
            c_multi8->Print(plottingfile.c_str(), "pdf");
        }

        //rec dca plots (2x2 canvas)
        auto *c_multi4 = new TCanvas("cmulti4", "cmulti4", 900, 600);
        c_multi4->SetWindowPosition(0, 400);

        c_multi4->Divide(2,2);
        {
            c_multi4->cd(1);
            gPad->SetLeftMargin(0.15);
            h_xdca->Draw();

            c_multi4->cd(2);
            gPad->SetLeftMargin(0.15);
            h_ydca->Draw();

            c_multi4->cd(3);
            gPad->SetLeftMargin(0.15);
            h_zdca->Draw();

            c_multi4->cd(4);
            gPad->SetLeftMargin(0.15);
            h_rdca->Draw();

            c_multi4->Print(plottingfile.c_str(), "pdf");

        }

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



        //     ####### SINGLE PLOTS #######

        auto *c_single1 = new TCanvas("csinlge1", "csinlge1");
        c_single1->SetLeftMargin(0.15);
        g_pdev_dca->Draw("AP");
        c_single1->Print(plottingfile.c_str(), "pdf");

        auto *c_single2 = new TCanvas("csinlge2", "csinlge2");
        c_single2->SetLeftMargin(0.15);
        g_pdev_z_dca->Draw("AP");
        c_single2->Print(plottingfile.c_str(), "pdf");

        auto *c_single3 = new TCanvas("csinlge3", "csinlge3");
        c_single3->SetLeftMargin(0.15);
        g_nhits_dca->SetMarkerColor(4);
        g_nhits_dca->Draw("AP");
        c_single3->Print(plottingfile.c_str(), "pdf");


        ////Monte carlo dca
//        auto *c_single4 = new TCanvas("csinlge4", "csinlge4");
//        c_single4->SetLeftMargin(0.15);
//        h_dcamc_phi->Draw();
//        c_single4->Print(plottingfile.c_str(), "pdf");
//
//        auto *c_single5 = new TCanvas("csinlge5", "csinlge5");
//        c_single5->SetLeftMargin(0.15);
//        h_zdcamc->Draw();
//        c_single5->Print(plottingfile.c_str(), "pdf");

        ////Einzel rec dca
//        auto *c_single6 = new TCanvas("csinlge6", "csingle6");
//        c_single6->SetLeftMargin(0.15);
//        h_xdca->Draw();
//        c_single6->Print(plottingfile.c_str(), "pdf");
//
//        auto *c_single7 = new TCanvas("csinlge7", "csingle7");
//        c_single7->SetLeftMargin(0.15);
//        h_ydca->Draw();
//        c_single7->Print(plottingfile.c_str(), "pdf");
//
//        auto *c_single8 = new TCanvas("csinlge8", "csingle8");
//        c_single8->SetLeftMargin(0.15);
//        h_zdca->Draw();
//        c_single8->Print(plottingfile.c_str(), "pdf");
//
//        auto *c_single9 = new TCanvas("csinlge9", "csingle9");
//        c_single9->SetLeftMargin(0.15);
//        h_rdca->Draw();
//        c_single9->Print(plottingfile.c_str(), "pdf");

        auto *c_single11 = new TCanvas();
        c_single11->SetLeftMargin(0.15);
        h_phimc->Draw();
        c_single11->Print(plottingfile.c_str(), "pdf");

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
        }
        filename = filename_template + "_hist_perr_nhits.pdf";
        c_single10->SaveAs(filename.c_str());

        ////########## FINAL PLOT CLOSES FILE

        auto *c_final = new TCanvas("c_final", "c_final");
        filename = plottingfile + ")";
        c_final->Print(filename.c_str(), "pdf");

        ////###########################################################################

        if(ADDITIONAL_PLOTS) {

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

            filename = filename_template + "_hist_3plots_err_nhits.pdf";
            c1->SaveAs(filename.c_str());

            //Single histogram for nhits
            auto  *c3 = new TCanvas();
            c3->SetWindowPosition(0, 500 );
            c3->SetLogy(1);
            auto legend2 = new TLegend();

            h_nhits->SetLineColor(4);
            h_nhits->SetFillStyle(16);
            h_nhits->SetName("nhits for reconstructed tracks");
            h_nhits->DrawClone("");

            filename = filename_template + "_hist_nhits.pdf";
            c3->SaveAs(filename.c_str());

            auto *c4 = new TCanvas();
            c4->SetWindowPosition(0, 500 );
            auto legend4 = new TLegend();
            g_pdev_phi->SetTitle("#phi - p_{err} correlation");
            //"#frac{1/p_{err} #minus 1/p_{mc}}{1/p_{mc}}"
            g_pdev_phi->GetXaxis()->SetTitle("p_{inv_err}");
            g_pdev_phi->GetXaxis()->SetRangeUser(-10.,10.);
            g_pdev_phi->GetYaxis()->SetTitle("#phi_{mc}");
            g_pdev_phi->Draw("ap");
            filename = filename_template + "_perr-phi.pdf";
            c4->SaveAs(filename.c_str());

            auto *c5 = new TCanvas();
            c5->SetWindowPosition(0, 500 );
            auto legend5 = new TLegend();
            g_pdev_dca->SetTitle("Dca - p_{err} correlation");
            g_pdev_dca->GetXaxis()->SetTitle("p_{inv_err}");
            g_pdev_dca->GetYaxis()->SetTitle("#phi_{mc}");
            g_pdev_dca->GetXaxis()->SetRangeUser(-7.,7.);
            g_pdev_dca->GetYaxis()->SetRangeUser(0.,60.);
            g_pdev_dca->Draw("ap");
            filename = filename_template + "_perr-dca.pdf";
            c5->SaveAs(filename.c_str());

            auto *c6 = new TCanvas();
            c6->SetWindowPosition(0, 500 );
            auto legend6 = new TLegend();
            g_pdev_p->SetTitle("p_{rec} - p_{err} correlation");
            g_pdev_p->GetYaxis()->SetTitle("p_{inv_err}");
            g_pdev_p->GetXaxis()->SetTitle("#p_{rec}");
            g_pdev_p->GetXaxis()->SetRangeUser(-1.e5,1.e5);
            g_pdev_p->GetYaxis()->SetRangeUser(-2.,2.);
            g_pdev_p->Draw("ap");
            filename = filename_template + "_p-perr.pdf";
            c6->SaveAs(filename.c_str());

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
            filename = filename_template + "_p-pmc.pdf";
            c7->SaveAs(filename.c_str());
        }
    }

    ////PRINT THE STATS
    //float p_inv_rec_error_mean = vector_mean(p_inv_rel_errors);

    cout << endl << endl << "---General Stats---\n" << endl;
    cout << "p fail analysis:" << endl;
    cout << "(p_mc == 0 || pt_mc == 0 || p_rec == 0 || pt_rec == 0) --- count: " << p_fail_count ;
    cout << " fails of total: " << segs_entries << ", fail rate: ";
    cout << (p_fail_count / (float) segs_entries) * 100 << " %" << endl;
    //cout << "1 / p reconstruction error mean: "  << p_inv_rec_error_mean * 100 << "% " << endl;

    for(int i=0; i < 6; i++) {
        if(i < 5) {
            cout << i+4 << " hits count: " << p_inv_rel_errors_hits[i].size() << endl;
        } else {
            cout << "other: " << p_inv_rel_errors_hits[i].size() << endl;
        }
    }

}

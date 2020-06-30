//
// Created by Konstantin Neureither on 30.06.20.
//
#ifndef PI
#define PI 3.1415926535
#endif //PI

#define TRIPLET_HIT_ARRAY_LENGTH 1024

#include "SegsRepresentation.h"
#include <iostream>
#include "utilityFunctions.h"

SegsRepresentation::SegsRepresentation(TTree *t_segs) {
    this->tr_segs = t_segs;
    this->segs_entries = this->tr_segs->GetEntries();
    setBranches();
}

void SegsRepresentation::setBranches() {
    tr_segs->SetBranchAddress("eventId", &rec_event);
    tr_segs->SetBranchAddress("nhit", &rec_nhit);
    tr_segs->SetBranchAddress("n", &rec_ntriplet);

    tr_segs->SetBranchAddress("p", &rec_p);
    tr_segs->SetBranchAddress("r", &rec_r);
    tr_segs->SetBranchAddress("rt", &rec_rt);
    tr_segs->SetBranchAddress("tan01", &rec_tan01);
    tr_segs->SetBranchAddress("tan12", &rec_tan12);
    tr_segs->SetBranchAddress("lam01", &rec_lam01);
    tr_segs->SetBranchAddress("lam12", &rec_lam12);
    tr_segs->SetBranchAddress("zpca_z", &rec_zpca_z);
    tr_segs->SetBranchAddress("zpca_r", &rec_zpca_r);

    tr_segs->SetBranchAddress("mc_tid", &mc_tid);
    tr_segs->SetBranchAddress("mc_p", &mc_p);
    tr_segs->SetBranchAddress("mc_pt", &mc_pt);
    tr_segs->SetBranchAddress("mc_theta", &mc_theta);
    tr_segs->SetBranchAddress("mc_phi", &mc_phi);
    tr_segs->SetBranchAddress("mc_lam", &mc_lam);
    tr_segs->SetBranchAddress("mc_type", &mc_type);
    tr_segs->SetBranchAddress("mc_pid", &mc_pid);

    tr_segs->SetBranchAddress("x00", &x00);
    tr_segs->SetBranchAddress("x10", &x10);
    tr_segs->SetBranchAddress("x20", &x20);
    tr_segs->SetBranchAddress("x01", &x01);
    tr_segs->SetBranchAddress("x11", &x11);
    tr_segs->SetBranchAddress("x21", &x21);

    tr_segs->SetBranchAddress("y00", &y00);
    tr_segs->SetBranchAddress("y10", &y10);
    tr_segs->SetBranchAddress("y20", &y20);
    tr_segs->SetBranchAddress("y01", &y01);
    tr_segs->SetBranchAddress("y11", &y11);
    tr_segs->SetBranchAddress("y21", &y21);

    tr_segs->SetBranchAddress("z00", &z00);
    tr_segs->SetBranchAddress("z10", &z10);
    tr_segs->SetBranchAddress("z20", &z20);
    tr_segs->SetBranchAddress("z01", &z01);
    tr_segs->SetBranchAddress("z11", &z11);
    tr_segs->SetBranchAddress("z21", &z21);

    tr_segs->SetBranchAddress("sid00", &sid00);
    tr_segs->SetBranchAddress("sid10", &sid10);
    tr_segs->SetBranchAddress("sid20", &sid20);
    tr_segs->SetBranchAddress("sid01", &sid01);
    tr_segs->SetBranchAddress("sid11", &sid11);
    tr_segs->SetBranchAddress("sid21", &sid21);

    std::cout << "Branches set for segs..." << std::endl;
}

void SegsRepresentationAndCalc::calcAdditionalData() {
    //basic data: Theta, phi and traverse p of reconstruction
    this->rec_pt = this->rec_p * std::cos((this->rec_lam01)[0]);
    this->rec_phi = (this->rec_tan01)[0];
    this->rec_theta = PI*0.5 - (this->rec_lam01)[0];

    mc_p_corr = this->mc_p * (this->mc_type % 10 == 3 ? -1 : (this->mc_type % 10 == 4 ? 1 : sgn(this->mc_pid)));
    mc_pt_corr = this->mc_pt * (this->mc_type % 10 == 3 ? -1 : (this->mc_type % 10 == 4 ? 1 : sgn(this->mc_pid)));

    this->mc_p_inv_corr = 1. / this->mc_p_corr;
    this->mc_pt_inv_corr = 1. / this->mc_pt_corr;

    this->rec_p_inv = 1. / rec_p;
    this->rec_pt_inv = 1. / rec_pt;

    this->p_inv_abs_error = this->rec_p_inv - this->mc_p_inv_corr;
    this->pt_inv_abs_error = this->rec_pt_inv - this-> mc_pt_inv_corr;
}

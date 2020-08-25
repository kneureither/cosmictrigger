//
// Created by Konstantin Neureither on 25.08.20.
//

#include "mu3eTree.h"

mu3eTree::mu3eTree(TTree *mu3e) {
    this->mu3e = mu3e;
    this->setBranches();
    this->my_entries = this->mu3e->GetEntries();
}

void mu3eTree::setBranches() {
    mu3e->SetBranchAddress("Header", &header);
    mu3e->SetBranchAddress("Nhit", &Nhit);
    mu3e->SetBranchAddress("hit_mc_i", &hit_mc_i);
    mu3e->SetBranchAddress("hit_mc_n", &hit_mc_n);
    mu3e->SetBranchAddress("traj_ID", &traj_ID);
    mu3e->SetBranchAddress("traj_PID", &traj_PID);
}

void mu3eTree::getEntry(const int &index) {
    mu3e->GetEntry(index);
    event = header[0];
    run = header[1];
}


////###########################################################################



mu3eMChitsTree::mu3eMChitsTree(TTree *mu3e_mchits) {
    this->mu3e_mchits = mu3e_mchits;
    this->setBranches();
    this->my_entries = this->mu3e_mchits->GetEntries();
}

void mu3eMChitsTree::setBranches() {
    mu3e_mchits->SetBranchAddress("tid", &tid);
    mu3e_mchits->SetBranchAddress("pdg", &pdg);
    mu3e_mchits->SetBranchAddress("hid", &hid);
    mu3e_mchits->SetBranchAddress("hid_g", &hid_g);
    mu3e_mchits->SetBranchAddress("pos_g_x", &pos_g_x);
    mu3e_mchits->SetBranchAddress("pos_g_y", &pos_g_y);
    mu3e_mchits->SetBranchAddress("pos_g_z", &pos_g_z);
}

void mu3eMChitsTree::getEntry(const int &index) {

}

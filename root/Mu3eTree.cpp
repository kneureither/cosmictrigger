//
// Created by Konstantin Neureither on 25.08.20.
//

#include "Mu3eTree.h"
#include <iostream>

Mu3eTree::Mu3eTree(TTree *mu3e) {
    this->mu3e = mu3e;
    this->setBranches();
    this->my_entries = this->mu3e->GetEntries();
}

void Mu3eTree::setBranches() {
    mu3e->SetBranchAddress("Header", &header);
    mu3e->SetBranchAddress("Nhit", &Nhit);
    mu3e->SetBranchAddress("hit_mc_i", &hit_mc_i);
    mu3e->SetBranchAddress("hit_mc_n", &hit_mc_n);
    mu3e->SetBranchAddress("traj_ID", &traj_ID);
    mu3e->SetBranchAddress("traj_PID", &traj_PID);
}

void Mu3eTree::getEntry(const int &index) {
    mu3e->GetEntry(index);
    event = header[0];
    run = header[1];
}

void Mu3eTree::Print(int entry) {
    std::cout << "Mu3e entry=" << entry << " event=" << event << " run=" << run << "  NHIT=" << Nhit << std::endl;
    std::cout << "\thit_mc_i\t";
    for(int i=0; i<hit_mc_i->size(); i++) std::cout << (*hit_mc_i)[i] << "  ";
    std::cout << std::endl << "\thit_mc_n\t";
    for(int i=0; i<hit_mc_n->size(); i++) std::cout << (*hit_mc_n)[i] << "  ";
    std::cout << std::endl << "\ttraj_ID\t";
    for(int i=0; i<traj_ID->size(); i++) std::cout << (*traj_ID)[i] << "  ";
    std::cout << std::endl << "\ttraj_PID\t";
    for(int i=0; i<traj_PID->size(); i++) std::cout << (*traj_PID)[i] << "  ";
    std::cout << std::endl << std::endl;
}


////###########################################################################



Mu3eMChitsTree::Mu3eMChitsTree(TTree *mu3e_mchits) {
    this->mu3e_mchits = mu3e_mchits;
    this->setBranches();
    this->my_entries = this->mu3e_mchits->GetEntries();
}

void Mu3eMChitsTree::setBranches() {
    mu3e_mchits->SetBranchAddress("tid", &tid);
    mu3e_mchits->SetBranchAddress("pdg", &pdg);
    mu3e_mchits->SetBranchAddress("hid", &hid);
    mu3e_mchits->SetBranchAddress("hid_g", &hid_g);
    mu3e_mchits->SetBranchAddress("pos_g_x", &pos_g_x);
    mu3e_mchits->SetBranchAddress("pos_g_y", &pos_g_y);
    mu3e_mchits->SetBranchAddress("pos_g_z", &pos_g_z);
}

void Mu3eMChitsTree::getEntry(const int &index) {
    this->mu3e_mchits->GetEntry(index);
}

void Mu3eMChitsTree::Print(const int entry) {
    std::cout << "  --> Mchits entry=" << entry << "  tid=" << tid << " pdg=" << pdg;
    std::cout << " hid=" << hid << " hid_g=" << hid_g << " pos_g=(" << pos_g_x << " / ";
    std::cout << pos_g_y << " / " << pos_g_z << " )" << std::endl;
}

//
// Created by Konstantin Neureither on 25.08.20.
//

#ifndef COSMICTRIGGER_MU3ETREE_H
#define COSMICTRIGGER_MU3ETREE_H

#include <vector>
#include "TTree.h"


class mu3eTree {
private:
public:
    TTree *mu3e;
    unsigned int my_entries;
    int header[4]; //event, run, type (empty), setup (emtpy)
    int event;
    int run;
    int Nhit;

    std::vector<int> *hit_mc_i = nullptr;
    std::vector<int> *hit_mc_n = nullptr;
    std::vector<unsigned int> *traj_ID = nullptr;
    std::vector<unsigned int> *traj_PID = nullptr;

    mu3eTree(TTree *mu3e);
    void setBranches();
    void getEntry(const int &index);
};

class mu3eMChitsTree {
private:
public:
    TTree *mu3e_mchits;
    unsigned int my_entries;

    int tid; //Unique trajectory ID
    int pdg; //Particle type encoded using the PDG scheme
    int hid; //Hit number (of corresponding detector) along trajectory, negative if trajectory is radially inwards moving
    int hid_g; //Global hit id (pixel hid)

    float pos_g_x;
    float pos_g_y;
    float pos_g_z;

    mu3eMChitsTree(TTree *mu3e_mchits);
    void setBranches();
    void getEntry(const int &index);
};


#endif //COSMICTRIGGER_MU3ETREE_H

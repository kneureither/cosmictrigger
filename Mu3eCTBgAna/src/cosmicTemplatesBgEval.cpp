//
// Created by Konstantin Neureither on 25.08.20.
//

#include "TFile.h"

#include "../inc/cosmicTemplatesBgEval.h"
#include "utilityFunctions.h"
#include "Mu3eTree.h"
#include <cassert>

void cosmicTemplatesBgEval(const int run) {
    /*
     * Read and analyse the mu3e mc hits
     * get the hits in xyz
     * assign sps to these hits
     * initialize template bank
     * check the frequency
     */

    const int MAX_ENTRIES = 10;

    const std::string pathtoBGdata = "data/SimulationData/";
    const std::string pathtoTemplateData = "data/TemplateData/";
    const std::string pathtoplots = "plots/Mu3eCosPatBgEval/";
    const std::string infile = pathtoBGdata + "mu3e_run_" + get_padded_string(run, 6, '0') + ".root";
//    std::string pathtorunplots = pathtoplots + "dataset_" + get_padded_string(dataset, 3, '0') + "/";
//    std::string pathtorundata = pathtodata + "dataset_" + get_padded_string(dataset, 3, '0') + "/";

    check_create_directory(pathtoBGdata);
    check_create_directory(pathtoTemplateData);
    check_create_directory(pathtoplots);

    // FILE FOR READING
    TFile tinF(infile.c_str());
    if (!tinF.IsOpen()) {
        std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
        exit(0);
    }

    TTree *t_mu3e;
    tinF.GetObject("mu3e", t_mu3e);
    TTree *t_mu3e_mchits;
    tinF.GetObject("mu3e_mchits", t_mu3e_mchits);

    Mu3eTree Mu3e = Mu3eTree(t_mu3e);
    Mu3eMChitsTree Mu3eMChits = Mu3eMChitsTree(t_mu3e_mchits);

    for(int frame=0; frame <= (MAX_ENTRIES == 0 ? Mu3e.my_entries : MAX_ENTRIES); frame++) {
        Mu3e.getEntry(frame);
        Mu3e.Print(frame);

        std::vector<bghit> bgframehits;

        assert(Mu3e.Nhit == Mu3e.hit_mc_i->size());
        for(int hit=0; hit< Mu3e.Nhit; hit++) {
            Mu3eMChits.getEntry((*Mu3e.hit_mc_i)[hit]);
            Mu3eMChits.Print((*Mu3e.hit_mc_i)[hit]);
            BGHIT.fill(Mu3eMChits.pos_g_x, Mu3eMChits.pos_g_y, Mu3eMChits.pos_g_z, 0);
            bgframehits.push_back(BGHIT);
        }
    }

}

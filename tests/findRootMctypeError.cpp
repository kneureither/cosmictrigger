//
// Created by Konstantin Neureither on 08.07.20.
//

//basics
#include <string>
#include <iostream>

//root
#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>

//custom
#include "utilityFunctions.h"
#include "rootData.h"
#include "SegsTreeRead.h"

void testSimpleRun() {
    const int run = 10;
    const std::string pathtodata = "data/SimulationData/";
    const std::string infile = pathtodata + "mu3e_run_" + get_padded_string(run, 6, '0') + "_trirec_cosmic.root";
    const int MAX_ENTRIES = 0;

    TFile f(infile.c_str());
    TTree *segs;

    f.GetObject("segs", segs);

    //data for trirec result in "segs" tree
    unsigned int segs_entries = segs->GetEntries();
    int rec_event;
    int rec_nhit;
    int mc_type;

    segs->SetBranchAddress("eventId", &rec_event);
    segs->SetBranchAddress("nhit", &rec_nhit);
    segs->SetBranchAddress("mc_type", &mc_type);

    cout << "Branches set for segs..." << endl;

    for (unsigned int entry = 0; entry < (MAX_ENTRIES == 0 ? segs_entries : MAX_ENTRIES); entry++) {
        segs->GetEntry(entry);

        if(!(mc_type == 3 || mc_type == 4)) {
            cout << "event id " << rec_event << " index: " << entry << " mc type: " << mc_type << endl;
        }
    }
}

void testRunWithSegsClass() {
    const int run = 10;
    const std::string pathtodata = "data/SimulationData/";
    const std::string infile = pathtodata + "mu3e_run_" + get_padded_string(run, 6, '0') + "_trirec_cosmic.root";
    const int MAX_ENTRIES = 0;

    // FILE FOR READING
    TFile tinF(infile.c_str());
    if (!tinF.IsOpen()) {
        std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
        exit(0);
    }

    TTree *t_segs = (TTree*)tinF.Get("segs");
//  tinF.GetObject("segs", t_segs);

    //class representation for segs tree and read functionality
    SegsTreeRead *Segs = new SegsTreeRead(t_segs);



    for (unsigned int entry = 0; entry < (MAX_ENTRIES == 0 ? Segs->my_entries : MAX_ENTRIES); entry++) {
        Segs->getEntry(entry);

        if((Segs->mc_type == 3 || Segs->mc_type == 4)) {
            cout << "event id " << Segs->rec_event << " index: " << entry << " mc type: " << Segs->mc_type << endl;
        }
    }


    delete Segs;
}



int main (int argc, char *argv[]) {
    testSimpleRun();

    testRunWithSegsClass();
}
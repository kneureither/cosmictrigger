//
// Created by Konstantin Neureither on 26.06.20.
//
#include "slimSegsData.h"
#include <iostream>
#include "utilityFunctions.h"

int main(int argc, char *argv[]) {
    slimSegsData(std::string("mu3e_slimmed_segs_000000.root"), 10, false);
}



//int run, outnum;
//if(argc < 2) {
//std::cout << "ERROR: Error in argument! Usage: Mu3eSlimSegs <run number> <slimmed file number> " << std::endl;
//exit(0);
//} else {
//run = atoi(argv[1]);
//outnum = atoi(argv[2]);
//}
//
////    run = 13;
////    outnum = 1;
//
//std::string outfile = "mu3e_slimmed_segs_" + get_padded_string(outnum, 6, '0') + ".root";
//
//std::cout << outfile << endl;
//
//std::cout << "Running Mu3ePreSlimSegs for run " << run << " with outfile " << outnum << "..." << std::endl;
//slimSegsData(outfile, run, false);
//
// Created by Konstantin Neureither on 26.06.20.
//
#include "slimSegsDataScript.h"
#include <iostream>
#include <vector>
#include "utilityFunctions.h"

int main(int argc, char *argv[]) {
    std::vector<int> runs = {14,16,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37};
    const int dataset = 6;

    for(auto it = std::begin(runs); it != runs.end(); ++it) {
        slimSegsDataScript(dataset, *it, true);
    }

    //data set 4 contains 14,16,19,20,21,22,23,24,25
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
//slimSegsDataScript(outfile, run, false);
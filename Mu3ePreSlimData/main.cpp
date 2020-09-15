//
// Created by Konstantin Neureither on 26.06.20.
//
#include "slimSegsDataScript.h"
#include <iostream>
#include <vector>
#include "utilityFunctions.h"

int main(int argc, char *argv[]) {
//    std::vector<int> runs = {14,16,19,20,21,22,23,24,25,27,28,29,30,31,32,33,34,35,36,37};
    std::vector<int> runs = {1031, 1032, 1033, 1034, 1035};
//    std::vector<int> runs;

    //add the bunch data
//    for(int i=1000; i<=1025; i++) {
//        runs.push_back(i);
//    }

    const int dataset = 10;

    bool append = true;

    for(auto it = std::begin(runs); it != runs.end(); ++it) {
        std::cout << "STATUS : adding run " << *it << " to dataset " << dataset << std::endl;
        slimSegsDataScript(dataset, *it, append);
        append = true;
    }

    //data set 4 contains 14,16,19,20,21,22,23,24,25
    //data set 6 contains 14,16,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37
    //run 26 is broken
    //data set 7 contains 14,16,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37 and 1000-1025
    //data set 8 contains 14,16,19,20,21,22,23,24,25,27,28,29,30,31,32,33,34,35,36,37 and 1000-1025
    //data set 9 contains 14,16,19,20,21,22,23,24,25,27,28,29,30,31,32,33,34,35,36,37 and 1000-1025 and 1026, 1036
    //data set 10 contains 14,16,19,20,21,22,23,24,25,27,28,29,30,31,32,33,34,35,36,37 and 1000-1025 and 1026, 1036, 1027-1030
        //to be added 1028-1035

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
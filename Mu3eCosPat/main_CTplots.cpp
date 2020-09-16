//
// Created by Konstantin Neureither on 11.08.20.
//
#include "include/makePlots.h"
#include <iostream>
#include <vector>

int main(int argc, char *argv[]) {
    int dataset;
    int combination_id;

    /*dataset 09 id 02*/ std::vector<int> cycle_plotting_order = {1,2,3,4};
//    /*dataset 09 id 01*/std::vector<int> cycle_plotting_order = {20, 16, 15, 14, 13, 19, 12, 11, 10, 9, 18, 8, 7, 6, 5, 17, 4, 3, 2 ,1}; //order dataset=9, id=001
//    std::vector<int> cycle_plotting_order {11, 9, 7, 6, 5, 4, 3, 2, 1, 8, 12}; //dataset 6, id=1
//    std::vector<int> cycle_plotting_order = {20, 19, 18, 17, 16, 12, 8, 4, 11, 7, 3, 10, 6, 2, 9, 5, 1}; //order dataset=9, id=001, eff sorted
//    std::vector<int> cycle_plotting_order = {16,15,14,13}; //order dataset=9, id=001

    if(argc < 2) {
        std::cout << "ERROR: Error in argument! Usage: "
                     "Mu3eCosPatPlots <dataset number> <sp ratio combination id>" << std::endl;
        exit(0);
    } else {
        dataset = atoi(argv[1]);
        combination_id = atoi(argv[2]);


    }

    std::cout << "Producing combined plots for dataset " << dataset << " with id " << combination_id << "..." << std::endl;
    makeCosPatPlots(dataset, combination_id, cycle_plotting_order);
    return 0;
}
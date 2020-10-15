//
// Created by Konstantin Neureither on 11.08.20.
//
#include "inc/cosmicTemplatesTrainingPlots.h"
#include <iostream>
#include <vector>
#include "Configuration.h"

int main(int argc, char *argv[]) {
    Configuration CONFIG;
    CONFIG.BUILDTB_TRAIN_PLOT();
    int dataset = CONFIG.dataset;
    int combination_id = CONFIG.TrainPlots.combination_id;

    std::vector<int> cycle_plotting_order = CONFIG.TrainPlots.cycle_plotting_order;
//    /*dataset 09 id 02*/ std::vector<int> cycle_plotting_order = {57,56,55,54,53,52};
//    /*dataset 09 id 01*/std::vector<int> cycle_plotting_order = {20, 16, 15, 14, 13, 19, 12, 11, 10, 9, 18, 8, 7, 6, 5, 17, 4, 3, 2 ,1}; //order dataset=9, id=001
//    std::vector<int> cycle_plotting_order {11, 9, 7, 6, 5, 4, 3, 2, 1, 8, 12}; //dataset 6, id=1
//    std::vector<int> cycle_plotting_order = {20, 19, 18, 17, 16, 12, 8, 4, 11, 7, 3, 10, 6, 2, 9, 5, 1}; //order dataset=9, id=001, eff sorted
//    std::vector<int> cycle_plotting_order = {16,15,14,13}; //order dataset=9, id=001

//    for(int i=1; i<=54; i++) {
//        cycle_plotting_order.push_back(i);
//    }

    std::cout << "Producing combined plots for dataset " << dataset << " with id " << combination_id << "..." << std::endl;
    makeCosPatPlots(dataset, combination_id, cycle_plotting_order, std::string());
    return 0;
}
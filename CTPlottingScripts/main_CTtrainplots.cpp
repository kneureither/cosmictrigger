//
// Created by Konstantin Neureither on 11.08.20.
//
#include "inc/ctTrainingPlots.h"
#include <iostream>
#include <vector>
#include "../CTCoreModules/Configuration.h"

int main(int argc, char *argv[]) {
    Configuration CONFIG;
    CONFIG.PLOT_BUILDTB();
    int dataset = CONFIG.dataset;
    int combination_id = CONFIG.TrainPlots.trainplotoutput_id;

    std::vector<int> cycle_plotting_order = CONFIG.TrainPlots.cycle_plotting_order;

    std::cout << "Producing combined plots for dataset " << dataset << " with id " << combination_id << "..." << std::endl;
    makeCosPatPlots(dataset, combination_id, cycle_plotting_order, std::string());
    return 0;
}
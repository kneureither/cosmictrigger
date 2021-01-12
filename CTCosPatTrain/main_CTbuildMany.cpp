//
// Created by Konstantin Neureither on 16.06.20.
//

#include "TemplateBank.h"
#include "cosmicTemplatesBuild.h"
#include "../CTPlottingScripts/inc/ctTrainingPlots.h"
#include <vector>
#include "../CTCoreModules/Configuration.h"

int main(int argc, char *argv[]) {

    Configuration CONFIG;
    CONFIG.BUILDTB();

    std::vector<float> sp_ratios = CONFIG.sp_ratios;
    std::vector<int> sp_count = CONFIG.sp_cnt;
    std::vector<float> stopping_effs = CONFIG.stopping_effs;
    int combination_id = CONFIG.TrainPlots.trainplotoutput_id; // produces a separate file for each trainplotoutput_id
    int dataset = CONFIG.dataset;
    bool append_to_outfile = true;

    std::vector<int> cycle_plotting_order;


    for(int i=0; i < sp_ratios.size(); i++) {
        for(int j=0; j < sp_count.size(); j++) {
            for(int n=0; n < stopping_effs.size(); n++) {
                std::cout << "Building Template Database for dataset " << dataset << " SPratio=" << sp_ratios[i]
                          << " SPcout=" << sp_count[j] << ".." << std::endl;
                cycle_plotting_order.push_back(cycle_plotting_order.size()+1);
                cosmicTemplatesBuild(dataset, sp_count[j], sp_ratios[i], combination_id, stopping_effs[n], append_to_outfile);
                append_to_outfile = true;
            }
        }
    }

    std::cout << "Producing combined plots for dataset " << dataset << "..." << std::endl;
    makeCosPatPlots(dataset, combination_id, cycle_plotting_order, std::string());

    return 0;
}
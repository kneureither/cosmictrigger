//
// Created by Konstantin Neureither on 27.10.20.
//

#include "cosmicTemplatesBuild.h"
#include "cosmicTemplatesEfficiency.h"
#include "../CTPlottingScripts/inc/cosmicTemplatesTrainingPlots.h"
#include <vector>
#include "Configuration.h"

int main(int argc, char *argv[]) {

    Configuration CONFIG;
    CONFIG.BUILDTB_DATAPOINTS();

    int combination_id = CONFIG.TrainPlots.combination_id; //will produce a separate file
    int dataset = CONFIG.dataset;
    bool append_to_outfile = true;

    std::vector<int> cycle_plotting_order;

    for (auto &curves : CONFIG.DBconfigDatapoints) {
        for (auto &config : curves) {
            std::cout << "(STATUS) : Building Template Database for dataset " << dataset
                      << " CONFIG: WxZ bins " << config.wbins << " x " << config.zbins
                      << " | SPR " << config.spr << " | SPC " << config.spc << " | EFF " << config.stopp_eff
                      << std::endl;

            cycle_plotting_order.push_back(cycle_plotting_order.size() + 1);
            cosmicTemplatesBuild(dataset, config.spc, config.spr, combination_id, config.stopp_eff, append_to_outfile);
            append_to_outfile = true;
        }
    }

    std::cout << "Producing combined plots for dataset " << dataset << "..." << std::endl;
    makeCosPatPlots(dataset, combination_id, cycle_plotting_order, CONFIG.set_description);

}
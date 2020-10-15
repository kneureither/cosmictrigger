//
// Created by Konstantin Neureither on 15.10.20.
//

#ifndef COSMICTRIGGER_CONFIGURATION_H
#define COSMICTRIGGER_CONFIGURATION_H

#include <vector>
#include "TemplateBank.h"

struct TltBnkTrnSet {
    std::vector<int> cycle_plotting_order = {1};
    combination_id = 2;
    std::string filelabel = "";

    void setFullCyclePlottingOrder(int plotcount) {
        for(in i = 1; i<=plotcount; i++) {
            cycle_plotting_order.push_back(i);
        }
    }
};

struct TltBnkFltr {
    std::vector<TIDLoadingFilter> filters = {ALL, NO_CENTER, CUT_ON_FREQ};

    void setDefault() {
        filters.clear();
        filters.push_back(ALL);
        filters.push_back(NO_CENTER);
        filters.push_back(CUT_ON_FREQ);
    }

    void setAllFilters() {
        filters.clear();
        filters.push_back(ALL);
        filters.push_back(NO_CENTER);
        filters.push_back(RECURL_ONLY);
        filters.push_back(CENTER_ONLY);
        filters.push_back(MIXED_ONLY);
        filters.push_back(CUT_ON_FREQ);
    }
};


class Configuration {
public:
    int dataset = 12;
    std::vector<float> sp_ratios;
    std::vector<int> sp_res;
    std::vector<float> stopping_effs;

    int cosmic_testing_dataset = 30;
    int max_cosmic_events = 0;
    int background_run = 107;
    int max_bg_events = 0;

    void michelBackground() {background_run = 107;};
    void fullBackground() {background_run = 109;}

    TltBnkTrnSet TrainPlots;
    TltBnkFltr TmplBankFilter;

    void set1() {
        this->resetMembers();

        sp_ratios.push_back(1);
        sp_res.push_back(400);
        stopping_effs.push_back(0.6);
        this->TrainPlots.setFullCyclePlottingOrder(numOfSettings());
    }
    void set2();
    void set3();

private:
    void resetMembers() {
        sp_ratios.clear();
        sp_res.clear();
        stopping_effs.clear();
        TrainPlots.cycle_plotting_order.clear();
    }

    int numOfSettings() {
        return sp_ratios.size() * sp_res.size() * stopping_effs.size();
    }
};


#endif //COSMICTRIGGER_CONFIGURATION_H

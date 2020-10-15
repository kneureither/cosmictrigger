//
// Created by Konstantin Neureither on 15.10.20.
//

#ifndef COSMICTRIGGER_CONFIGURATION_H
#define COSMICTRIGGER_CONFIGURATION_H

#include <vector>
#include "TemplateBank.h"

struct TltBnkTrnSet {
    std::vector<int> cycle_plotting_order = {1};
    int combination_id = 0;
    std::string filelabel = "";

    void setFullCyclePlottingOrder(int plotcount) {
        for(int i = 1; i<=plotcount; i++) {
            cycle_plotting_order.push_back(i);
        }
    }
};

struct TltBnkFltr {
    std::vector<TIDLoadingFilter> filters;

    TltBnkFltr() {
        setDefault();
    }

    void setAll() {
        filters.clear();
        filters.push_back(ALL);
    }

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
    int mode = 0;
    bool prints = false;
    std::vector<float> sp_ratios;
    std::vector<int> sp_res;
    std::vector<float> stopping_effs;

    int cosmic_testing_dataset = 30;
    int max_cosmic_events = 0;
    int background_run = 107;
    int max_bg_frames = 0;
    int max_bg_frame_nhits = 0;

    std::string set_description;

    Configuration &michelBackground() {background_run = 107; return *this;}
    Configuration &fullBackground() {background_run = 109; return *this;}

    TltBnkTrnSet TrainPlots;
    TltBnkFltr TmplBankFilter;

    Configuration &set1() {
        this->resetMembers();

        sp_ratios.push_back(1);
        sp_res.push_back(400);
        stopping_effs.push_back(0.6);
        this->TrainPlots.setFullCyclePlottingOrder(numOfSettings());

        return *this;
    }

    Configuration &set2_spcount() {
        this->resetMembers();
        dataset = 12;
        background_run = 107;
        max_bg_frames = 99900;
        max_bg_frame_nhits = 100;

        sp_ratios.push_back(32);
        sp_ratios.push_back(128);
        sp_res.push_back(1152);
        sp_res.push_back(2048);
        sp_res.push_back(3200);
        stopping_effs.push_back(0.5);

        TmplBankFilter.filters.push_back(ALL);

        set_description = "stripes_128_32";

        return *this;
    }

    Configuration &set2_roc() {
        this->resetMembers();
        dataset = 12;
        background_run = 107;
        max_bg_frames = 99900;
        max_bg_frame_nhits = 100;

        sp_ratios.push_back(128);
        sp_res.push_back(3200);
        stopping_effs.push_back(0.5);
        stopping_effs.push_back(0.6);

        TmplBankFilter.filters.push_back(ALL);

        set_description = "stripes_128";

        return *this;
    }

    Configuration &set3_build() {
        this->resetMembers();
        dataset = 12;
        background_run = 107;
        max_bg_frames = 99900;
        max_bg_frame_nhits = 100;

        TrainPlots.combination_id = 0;
        TrainPlots.setFullCyclePlottingOrder(numOfSettings());

//        sp_ratios.push_back(1);
//        sp_ratios.push_back(4);
//        sp_ratios.push_back(16);
        sp_ratios.push_back(64);
//        sp_ratios.push_back(256);


//        sp_res.push_back(576);
        sp_res.push_back(1024);
//        sp_res.push_back(1600);
//        sp_res.push_back(4096);
        stopping_effs.push_back(0.6);
        stopping_effs.push_back(0.8);

        TmplBankFilter.filters.push_back(ALL);
        TmplBankFilter.filters.push_back(NO_CENTER);

        set_description = "fltr_ready";

        return *this;
    };

    Configuration &set3_spc() {
        this->resetMembers();
        dataset = 12;
        background_run = 107;
        max_bg_frames = 99900;
        max_bg_frame_nhits = 100;

        TrainPlots.combination_id = 0;
        TrainPlots.setFullCyclePlottingOrder(numOfSettings());

//        sp_ratios.push_back(1);
//        sp_ratios.push_back(4);
//        sp_ratios.push_back(16);
        sp_ratios.push_back(64);
//        sp_ratios.push_back(256);


        sp_res.push_back(576);
        sp_res.push_back(1024);
        sp_res.push_back(1600);
//        sp_res.push_back(4096);
        stopping_effs.push_back(0.6);

        set_description = "fltr_ready";

        return *this;
    };

    Configuration &set3_roc() {
        this->resetMembers();
        dataset = 12;
        background_run = 107;
        max_bg_frames = 99900;
        max_bg_frame_nhits = 100;

        sp_ratios.push_back(64);

//        sp_res.push_back(576);
        sp_res.push_back(1024);
//        sp_res.push_back(1600);
//        sp_res.push_back(4096);
        stopping_effs.push_back(0.6);
        stopping_effs.push_back(0.8);

        set_description = "fltr_ready";

        return *this;
    };

    std::string RatiosToString() {
        std::string ratios = "";
        for(int i=0; i< sp_ratios.size(); i++) {
            ratios = ratios + get_string(sp_ratios[i]) + (i < sp_ratios.size() - 1 ? "_" : "");
        }
        return ratios;
    }

    std::string StoppingEffsToString() {
        std::string effs = "";
        for(int i=0; i< stopping_effs.size(); i++) {
            effs = effs + get_string(stopping_effs[i]) + (i < stopping_effs.size() - 1 ? "_" : "");
        }
        return effs;
    }

    std::string ResolutionsToString() {
        std::string res = "";
        for(int i=0; i< sp_res.size(); i++) {
            res = res + get_string(sp_res[i]) + (i < sp_res.size() - 1 ? "_" : "");
        }
        return res;
    }


    //// SETTING CALLS OF DIFFERENT SCRIPTS

    void BUILDTB() {
        set3_build();
    }

    void BUILDTB_TRAIN_PLOT(){
        set3_build();
        TrainPlots.setFullCyclePlottingOrder(2);

    }
    void BUILDTB_COS_EFF() {
        set3_build();
    }

    void BGANA() {
        set3_build();

    }

    void BGANA_PLOT_ROC() {
//        set2_roc();
        set3_roc();
    }
    void BGANA_PLOT_SPC() {
        set3_spc();
//        set2_spcount();
    }


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

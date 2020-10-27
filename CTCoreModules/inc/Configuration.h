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
    }

    TltBnkFltr &CLR() {
        filters.clear();
        return *this;
    }

    TltBnkFltr &addALL() {
        filters.push_back(ALL);
        return *this;
    }

    TltBnkFltr &addNO_CENTER() {
        filters.push_back(NO_CENTER);
        return *this;
    }

    TltBnkFltr &addCUT_ON_FREQ() {
        filters.push_back(CUT_ON_FREQ);
        return *this;
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

struct BkgRateAna {
    // Final Plot
//    std::vector<int> bg_runs = {120, 121, 122, 126, 123};
//    std::vector<float> beam_rates = {2e7, 6e7, 1e8, 2e8, 1e9};
//    std::vector<float> nhit_mean = {5.3, 14.5, 24.3, 0, 242};
//    std::vector<float> nhit_sigma = {10.55, 17.65, 22.9, 0, 71.9};
//    std::vector<int> max_bkg_frames = {0, 0, 0, 0, 0};  //to reduce claculation time


    // Testing
    std::vector<int> bg_runs = {107, 122};
    std::vector<float> beam_rates = {6e7, 1e8};
    std::vector<float> nhit_mean = {14.5, 24.3};
    std::vector<float> nhit_sigma = {17.65, 22.9};
    std::vector<int> max_bkg_frames = {99900, 99900};  //to reduce claculation time


    std::vector<int> max_nhits;

    std::vector<int> make_max_nhit_cut = {1, 0}; // yes no
    BkgRateAna(){
        calc_max_nhits(1);
    };

    void calc_max_nhits(float sigmas) {
      for(int i=0; i<nhit_mean.size(); i++) {
          max_nhits.push_back((int) (nhit_mean[i] + (float) sigmas * nhit_sigma[i]));
      }
    }

};


class Configuration {
public:
    int dataset = 12;
    int mode = 0;
    bool prints = false;
    std::vector<float> sp_ratios;
    std::vector<int> sp_cnt;
    std::vector<float> stopping_effs;

    int cosmic_testing_dataset = 30;
    int max_cosmic_events = 0;
    int background_run = 107;
    int max_bg_frames = 0;
    int max_bg_frame_nhits = 0;

    std::string pathtosimfiles = "../../../../../../../../../Volumes/Extreme SSD KN/MAC EXTENSION/Mu3e_data/ServerData/data/";

    std::string set_description;

    Configuration &michelBackground() {background_run = 107; return *this;}
    Configuration &fullBackground() {background_run = 109; return *this;}

    TltBnkTrnSet TrainPlots;
    TltBnkFltr TmplBankFilter;
    BkgRateAna BkgFiles;

    Configuration &set1() {
        this->resetMembers();

        sp_ratios.push_back(1);
        sp_cnt.push_back(400);
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
        sp_cnt.push_back(1152);
        sp_cnt.push_back(2048);
        sp_cnt.push_back(3200);
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
        sp_cnt.push_back(3200);
        stopping_effs.push_back(0.5);
        stopping_effs.push_back(0.6);

        TmplBankFilter.filters.push_back(ALL);

        set_description = "stripes_128";

        return *this;
    }

    Configuration &set3_base() {
        this->resetMembers();
        dataset = 12;
        background_run = 107;
        max_bg_frames = 99900;
        max_bg_frame_nhits = 100;
        TrainPlots.combination_id = 0;
        TmplBankFilter.CLR().addALL().addNO_CENTER().addCUT_ON_FREQ();
        set_description = "fltr_ready";
        return *this;
    }

    Configuration &set3_trainplot() {
        //only SPratio=64, tb_eff=0.6, 0.8
        std::vector<int> cycles = {3,6,4,1,5,2};
        TrainPlots.cycle_plotting_order = cycles;
    }

    Configuration &set3_build() {

//        sp_ratios.push_back(1);
//        sp_ratios.push_back(4);
//        sp_ratios.push_back(16);
        sp_ratios.push_back(64);
//        sp_ratios.push_back(256);

//        sp_cnt.push_back(576);
        sp_cnt.push_back(1024);
//        sp_cnt.push_back(1600);
//        sp_cnt.push_back(4096);
        stopping_effs.push_back(0.6);
        stopping_effs.push_back(0.8);

        return *this;
    };

    Configuration &set3_spc(float ratio, float train_eff) {
        sp_ratios.clear();
        sp_ratios.push_back(ratio);

        sp_cnt.push_back(576);
        sp_cnt.push_back(1024);
        sp_cnt.push_back(1600);
//        sp_cnt.push_back(4096);
        stopping_effs.push_back(train_eff);
        return *this;
    };

    Configuration &set3_roc() {

        sp_ratios.push_back(64);

//        sp_cnt.push_back(576);
        sp_cnt.push_back(1024);
//        sp_cnt.push_back(1600);
//        sp_cnt.push_back(4096);
        stopping_effs.push_back(0.6);
        stopping_effs.push_back(0.8);
        return *this;
    };

    Configuration &bkg_rates() {
        sp_cnt.push_back(1024);
        sp_ratios.push_back(64);
        stopping_effs.push_back(0.6);


        int datapoint_bg_run = 0;
        max_bg_frames = 99900;
//        max_bg_frame_nhits = 0;
        max_bg_frame_nhits = BkgFiles.max_nhits[datapoint_bg_run];
        background_run = BkgFiles.bg_runs[datapoint_bg_run];

        TmplBankFilter.filters.clear();
        TmplBankFilter.addALL().addNO_CENTER();

        return *this;
    }


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
        for(int i=0; i < sp_cnt.size(); i++) {
            res = res + get_string(sp_cnt[i]) + (i < sp_cnt.size() - 1 ? "_" : "");
        }
        return res;
    }

    std::string FiltersToString() {
        std::string filter = "";
        for(int i=0; i< TmplBankFilter.filters.size(); i++) {
            filter = filter + enum_to_string(TmplBankFilter.filters[i]) + (i < TmplBankFilter.filters.size() - 1 ? "_" : "");
        }
        return filter;
    }


    //// SETTING CALLS OF DIFFERENT SCRIPTS

    void BUILDTB() {
//        set3_base().set3_build();
        dataset=12;
        stopping_effs.push_back(0.9);
        sp_cnt.push_back(400);
        sp_ratios.push_back(25);
        sp_ratios.push_back(0.04);
        sp_ratios.push_back(1);
        TrainPlots.combination_id = 3;
    }

    void BUILDTB_TRAIN_PLOT(){
        set3_base().set3_trainplot();
    }

    void BUILDTB_COS_EFF() {
        set3_base().set3_build();
    }

    void BGANA() {
//        bkg_rates();
        set3_base().set3_build();
    }

    void BGANA_PLOT_ROC() {
        set3_base().set3_roc();
    }
    void BGANA_PLOT_SPC() {
        set3_base().set3_spc(64, 0.6);
    }


private:
    void resetMembers() {
        sp_ratios.clear();
        sp_cnt.clear();
        stopping_effs.clear();
        TrainPlots.cycle_plotting_order.clear();
    }

    int numOfSettings() {
        return sp_ratios.size() * sp_cnt.size() * stopping_effs.size();
    }
};


#endif //COSMICTRIGGER_CONFIGURATION_H

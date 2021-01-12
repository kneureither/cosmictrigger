//
// Created by Konstantin Neureither on 15.10.20.
//

#ifndef COSMICTRIGGER_CONFIGURATION_H
#define COSMICTRIGGER_CONFIGURATION_H

#include <vector>
#include "TemplateBank.h"
#include <cmath>

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

struct InfoForBkgRejVsBeamRatePlot {
    // Final Plot
//    std::vector<int> bg_runs =        {120,   121,    122,    126,    128,    123};
//    std::vector<float> beam_rates =   {2e7,   6e7,    1e8,    2e8,    4.932e8,1e9};
//    std::vector<float> nhit_mean =    {5.3,   14.5,   24.3,   48.3,   119.5,  242};
//    std::vector<float> nhit_sigma =   {10.55, 17.65,  22.9,   32.1,   50.6,   71.9};
//    std::vector<int> max_bkg_frames = {0,     0,      0,      0,      0,      0};  //to reduce claculation time

    // michel: 107; 6e7; 14.5; 17.65
    // plotting

    std::vector<int> bg_runs =        {125,   122,    126,   127};
    std::vector<float> beam_rates =   {5e7,   1e8,    2e8,   2.466e8};
    std::vector<float> nhit_mean =    {14.5,  24.3,   48.3,  59.86};
    std::vector<float> nhit_sigma =   {17.65, 22.9,   32.1,  36.07};
    std::vector<int> max_bkg_frames = {249900, 150000,  99900, 80000 /* 70000 for 50% */};

//    std::vector<int> bg_runs =        {127};
//    std::vector<float> beam_rates =   {2.466e8};
//    std::vector<float> nhit_mean =    {59.86};
//    std::vector<float> nhit_sigma =   {36.07};
//    std::vector<int> max_bkg_frames = {70000};

    // Testing
//    std::vector<int> bg_runs = {127};
//    std::vector<float> beam_rates = {2.466e8};
//    std::vector<float> nhit_mean = {59.86};
//    std::vector<float> nhit_sigma = {36.07};
//    std::vector<int> max_bkg_frames = {99900};  //to reduce calculation time

    std::vector<int> max_nhits; //calculated automatically

    std::vector<int> make_max_nhit_cut = {0,1}; // yes no

    int ndatapoints;

    InfoForBkgRejVsBeamRatePlot(){
        calc_max_nhits(1);
        ndatapoints = bg_runs.size();
    };

    void calc_max_nhits(float sigmas) {
      for(int i=0; i<nhit_mean.size(); i++) {
          max_nhits.push_back((int) (nhit_mean[i] + (float) sigmas * nhit_sigma[i]));
      }
    }

};

struct DatabaseConfigBuild {
    int spc;
    float spr;
    int wbins;
    int zbins;
    float stopp_eff;

    DatabaseConfigBuild(int wbins, int zbins, float stopp_eff) {
        this->wbins = wbins;
        this->zbins = zbins;
        this->stopp_eff = stopp_eff;
        this->spc = wbins*zbins;
        this->spr = wbins / (float) zbins;
    }

    DatabaseConfigBuild(int spc, float spr, float stopp_eff) {
        this->wbins = (int) sqrt((float) spr * (float) spc);
        this->zbins = (int) sqrt((float) spc / (float) spr);
        this->spr = spr;
        this->spc = spc;
        this->stopp_eff = stopp_eff;
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

    Configuration &michelBackground() {background_run = 107; max_bg_frames=99900; return *this;}
    Configuration &fullBackground() {background_run = 109; return *this;}

    TltBnkTrnSet TrainPlots;
    TltBnkFltr TmplBankFilter;
    InfoForBkgRejVsBeamRatePlot BkgFiles;
    std::vector<std::vector<DatabaseConfigBuild>> DBconfigDatapoints;

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
        return *this;
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

    Configuration &set13_build_spc_final() {
        this->resetMembers();
        dataset = 13;
        TrainPlots.combination_id = 0;
        sp_ratios.push_back(4);
//        sp_cnt.push_back(512);
        sp_cnt.push_back(784);
//        sp_cnt.push_back(1024);
        stopping_effs.push_back(0.8);

        return *this;
    };

    Configuration &set12_plot_traineff_spc_final() {
        this->resetMembers();
        dataset = 12;
        TrainPlots.combination_id = 0;
        TrainPlots.cycle_plotting_order.clear();
        TrainPlots.cycle_plotting_order.push_back(3);
        TrainPlots.cycle_plotting_order.push_back(6);
        TrainPlots.cycle_plotting_order.push_back(4);
        return *this;
    };

    Configuration &set12_plot_traineff_spc_appendix() {
        this->resetMembers();
        dataset = 12;
        TrainPlots.combination_id = 5;
        TrainPlots.cycle_plotting_order.clear();
        TrainPlots.cycle_plotting_order.push_back(3);
        TrainPlots.cycle_plotting_order.push_back(16);
        TrainPlots.cycle_plotting_order.push_back(6);
        TrainPlots.cycle_plotting_order.push_back(4);
        TrainPlots.cycle_plotting_order.push_back(15);

        return *this;
    };

    Configuration &set13_plot_traineff_software() {
        this->resetMembers();
        dataset = 13;
        TrainPlots.combination_id = 9;
        TrainPlots.cycle_plotting_order.clear();
        TrainPlots.cycle_plotting_order.push_back(20);
        TrainPlots.cycle_plotting_order.push_back(48);

        return *this;
    };

    Configuration &set12_plot_traineff_spr_final() {
        this->resetMembers();
        dataset=12;
        TrainPlots.cycle_plotting_order.clear();
        TrainPlots.cycle_plotting_order.push_back(4);
        TrainPlots.cycle_plotting_order.push_back(1);
        TrainPlots.cycle_plotting_order.push_back(3);
        TrainPlots.cycle_plotting_order.push_back(2);
        TrainPlots.cycle_plotting_order.push_back(5);
        TrainPlots.combination_id=3;
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

//    Configuration &set13_base() {
//        this->resetMembers();
//        dataset = 13;
//        background_run = 107;
//        max_bg_frames = 99900;
//        max_bg_frame_nhits = 0;
//        TrainPlots.combination_id = 0;
//        TmplBankFilter.CLR().addALL().addNO_CENTER().addCUT_ON_FREQ();
//        set_description = "final";
//        return *this;
//    }

    Configuration &set13_spc_tmpl_plot() {
        /**
         * This creates the datapoints for the #templ / spc plot
         */
        this->resetMembers();
        dataset = 13;
        TrainPlots.combination_id=0;
        this->michelBackground();
        TmplBankFilter.CLR().addALL().addNO_CENTER();

        float stopping_eff = 0.6;
        std::vector<DatabaseConfigBuild> plot_datapoints;
        plot_datapoints.push_back(DatabaseConfigBuild(1024, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(768, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(512, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(384, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(128, 4, stopping_eff));

        DBconfigDatapoints.push_back(plot_datapoints); plot_datapoints.clear();
        stopping_eff = 0.8;
        plot_datapoints.push_back(DatabaseConfigBuild(512, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(384, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(128, 4, stopping_eff));
        DBconfigDatapoints.push_back(plot_datapoints);

        return *this;
    };

    Configuration &set13_high_spr_check_build_data() {
        /**
         * This creates the datapoints for an appendix plot
         */
        this->resetMembers();
        dataset = 12;
        float stopping_eff = 0.8;
        TrainPlots.combination_id=0;
        std::vector<DatabaseConfigBuild> plot_datapoints;
        plot_datapoints.push_back(DatabaseConfigBuild(400, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(144, 4, stopping_eff));

        DBconfigDatapoints.push_back(plot_datapoints);

        return *this;
    };



    Configuration &set13_plot_traineff_final() {
        this->resetMembers();
        dataset = 13;
        TrainPlots.combination_id = 0;
        TrainPlots.cycle_plotting_order.clear();
        TrainPlots.cycle_plotting_order.push_back(10);
        TrainPlots.cycle_plotting_order.push_back(11);
        TrainPlots.cycle_plotting_order.push_back(12);
        TrainPlots.cycle_plotting_order.push_back(13);
        return *this;
    };

    Configuration &set13_plot_traineff_1024_final() {
        this->resetMembers();
        dataset = 13;
        TrainPlots.combination_id = 6;
        TrainPlots.cycle_plotting_order.clear();
        TrainPlots.cycle_plotting_order.push_back(17);
        return *this;
    };

    Configuration &set13_bkg_spr_plot() {
        /**
         * This creates the datapoints for the #templ / spc plot
         */
        this->resetMembers();
        this->michelBackground();
        dataset = 13;
        TrainPlots.combination_id=0;
        TmplBankFilter.CLR().addALL();

        float stopping_eff = 0.6;
        std::vector<DatabaseConfigBuild> plot_datapoints;
        plot_datapoints.push_back(DatabaseConfigBuild(256, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(128, 8, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(64, 16, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(32, 32, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(16, 64, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(8, 128, stopping_eff));

        DBconfigDatapoints.push_back(plot_datapoints);
        plot_datapoints.clear();

        stopping_eff = 0.8;
        plot_datapoints.push_back(DatabaseConfigBuild(256, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(128, 8, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(64, 16, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(32, 32, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(16, 64, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(8, 128, stopping_eff));
        DBconfigDatapoints.push_back(plot_datapoints);


        return *this;
    };


    Configuration &set13_bkg_spc_plot() {
        /**
         * This creates the datapoints for the bkg SPC eval
         */
        this->resetMembers();
        this->michelBackground();
        dataset = 13;
        TrainPlots.combination_id=0;
        TmplBankFilter.CLR().addALL();
        max_bg_frame_nhits=0;

        float stopping_eff = 0.6;
        std::vector<DatabaseConfigBuild> plot_datapoints;
        plot_datapoints.push_back(DatabaseConfigBuild(1024, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(768, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(512, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(384, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(128, 4, stopping_eff));

        DBconfigDatapoints.push_back(plot_datapoints); plot_datapoints.clear();
        stopping_eff = 0.8;
        plot_datapoints.push_back(DatabaseConfigBuild(512, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(384, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(128, 4, stopping_eff));

        DBconfigDatapoints.push_back(plot_datapoints); plot_datapoints.clear();
        stopping_eff = 0.6;
        plot_datapoints.push_back(DatabaseConfigBuild(64, 8, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(128, 8, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 8, stopping_eff));

        DBconfigDatapoints.push_back(plot_datapoints); plot_datapoints.clear();
        stopping_eff = 0.8;
        plot_datapoints.push_back(DatabaseConfigBuild(64, 8, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(128, 8, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 8, stopping_eff));
        DBconfigDatapoints.push_back(plot_datapoints);

        return *this;
    };

    Configuration &set13_bkg_spc_zbins_plot() {
        /**
         * This creates the datapoints for the bkg SPC eval
         */

        this->resetMembers();
        this->michelBackground();
        dataset = 13;
        TrainPlots.combination_id=0;
        TmplBankFilter.CLR().addALL();
        max_bg_frame_nhits=0;

        float stopping_eff = 0.6;
        std::vector<DatabaseConfigBuild> plot_datapoints;
        plot_datapoints.push_back(DatabaseConfigBuild(256, 2, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 8, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 12, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 16, stopping_eff));

        DBconfigDatapoints.push_back(plot_datapoints); plot_datapoints.clear();
        stopping_eff = 0.8;
        plot_datapoints.push_back(DatabaseConfigBuild(256, 2, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 8, stopping_eff));
//        plot_datapoints.push_back(DatabaseConfigBuild(256, 12, stopping_eff));
//        plot_datapoints.push_back(DatabaseConfigBuild(256, 16, stopping_eff));

        DBconfigDatapoints.push_back(plot_datapoints); plot_datapoints.clear();

        return *this;
    };

    Configuration &set13_bkg_roc_plot() {
        /**
         * This creates the datapoints for the bkg SPC eval
         */
        this->resetMembers();
        this->michelBackground();
        dataset = 13;
        TrainPlots.combination_id=0;
        TmplBankFilter.CLR().addALL().addNO_CENTER().addCUT_ON_FREQ();
        max_bg_frame_nhits=0;

        std::vector<DatabaseConfigBuild> plot_datapoints;
//        DBconfigDatapoints.push_back(plot_datapoints);
        DBconfigDatapoints.push_back(plot_datapoints);
//        DBconfigDatapoints.push_back(plot_datapoints);

        for(int i=1; i<10; i++ ) {
            std::cout << "adding datapoints " << i*0.1 << std::endl;
//            DBconfigDatapoints[0].push_back(DatabaseConfigBuild(128, 4, 0.1*i));
            DBconfigDatapoints[0].push_back(DatabaseConfigBuild(256, 4, 0.1*i));
//            if(i<9) DBconfigDatapoints[2].push_back(DatabaseConfigBuild(256, 8, 0.1*i));
        }

        return *this;
    };


    Configuration &bkg_rates() {
        dataset=13;

        sp_cnt.push_back(2048);
        sp_ratios.push_back(32);
        stopping_effs.push_back(0.8);

        BkgFiles.make_max_nhit_cut.clear();
        BkgFiles.make_max_nhit_cut.push_back(0);
        BkgFiles.make_max_nhit_cut.push_back(1);

        TmplBankFilter.filters.clear();
        TmplBankFilter.addALL().addNO_CENTER();
//        TmplBankFilter.addALL();

        return *this;
    }


    //// SETTING CALLS OF DIFFERENT SCRIPTS

    void BUILDTB() {
//        set3_base().set3_build();
//        set13_build_spc_final();

        dataset = 13;
        stopping_effs.push_back(1);
        sp_ratios.push_back(64);
        sp_cnt.push_back(1024);
        TrainPlots.combination_id = 0;


//        dataset=12;
//        stopping_effs.push_back(0.9);
//        sp_cnt.push_back(400);
//        sp_ratios.push_back(400);
//        sp_ratios.push_back(0.0025);
//        TrainPlots.combination_id = 3;
    }

    void BUILDTB_DATAPOINTS() {
//        set13_spc_tmpl_plot();
//        set13_high_spr_check_build_data(); //final plot fitted in thesis
//        set13_bkg_spr_plot();
//        set13_bkg_spc_zbins_plot();
        dataset=13;
        TrainPlots.combination_id = 0;

        std::vector<DatabaseConfigBuild> plot_datapoints;
        for(int i=6; i<=6; i++ ) {
//            plot_datapoints.push_back(DatabaseConfigBuild(8, 32, 0.1*i));
            plot_datapoints.push_back(DatabaseConfigBuild(24, 24, 0.1*i));
//            plot_datapoints.push_back(DatabaseConfigBuild(32, 8, 0.1*i));
        }
        DBconfigDatapoints.push_back(plot_datapoints);
    }

    void PLOT_BUILDTB(){
//        set3_base().set3_trainplot();
//        set12_plot_traineff_spr_final(); //check
//        set12_plot_traineff_spc_final(); //
        set13_plot_traineff_final(); //final plot in thesis // check
//        set13_plot_train eff_1024_final(); //final plot in thesis //check
//        set12_plot_traineff_spc_appendix(); // final plot / check
//        set13_plot_traineff_software(); // check
    }

    void BUILDTB_COS_EFF() {
        set3_base().set3_build();

    }

    void BGANA_DATAPOINTS() {
//        set13_bkg_spr_plot();
//        set13_bkg_roc_plot();
        dataset=13;
        TrainPlots.combination_id = 0;
        this->michelBackground();
        TmplBankFilter.CLR().addCUT_ON_FREQ();

        std::vector<DatabaseConfigBuild> plot_datapoints;
        for(int i=6; i<10; i++ ) {
//            plot_datapoints.push_back(DatabaseConfigBuild(128, 4, 0.1*i));
            plot_datapoints.push_back(DatabaseConfigBuild(256, 4, 0.1*i));
//            plot_datapoints.push_back(DatabaseConfigBuild(256, 8, 0.1*i));
        }
        DBconfigDatapoints.push_back(plot_datapoints);
    }

    void BGANA() {
        set3_base().set3_build();
    }

    void BGANA_MULTI() {
        bkg_rates();
    }

    void BGANA_PLOT_ROC() {
        set3_base().set3_roc();
    }
    void BGANA_PLOT_SPC() {
        set3_base().set3_spc(64, 0.6);
    }
    void PLOT_BGANA_SPC_DATAPOINTS() {
//        set13_bkg_spc_plot();
        set13_bkg_spc_zbins_plot();
    }

    void PLOT_BGANA_SPR_DATAPOINTS() {
        set13_bkg_spr_plot();
    }

    void PLOT_BGANA_ROC_DATAPOINTS() {
        set13_bkg_roc_plot();
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

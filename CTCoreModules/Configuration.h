//
// Created by Konstantin Neureither on 15.10.20.
//

#ifndef COSMICTRIGGER_CONFIGURATION_H
#define COSMICTRIGGER_CONFIGURATION_H

#include <vector>
#include "TemplateBank.h"
#include <cmath>

/**
 * This file and Classes are the first approach to a centralised parameter handling by the use of config files.
 * Almost every parameter that might need to be changed for the analysis was moved to this class.
 * (For the plotting scripts, it might still be required to perform some code apations to produce the desired plots)
 *
 * The idea is, that the Configuration class has all parameters as members, and different scripts call their parameter
 * set by calling the classes member functions.
 *
 * Especially for the training phase and background evaluation phase, it was necessary to simulate very different
 * template databases with many different configs. Also, plots needed to be produced, that include single numbers from
 * each of these simulated banks and their performances. To cope with this task, a central 2d vector was introduced:
 *
 *     'std::vector<std::vector<DatabaseConfigBuild>> DBconfigCurveDatapoints;'
 *
 * Its inner vector includes elements of the struct DatabaseConfigBuild which handles the very basic template bank
 * paramters which are mostly examined, namely SPC, SPR and STOPP_EFF (cosmic efficiency) and WBINS and ZBINS respectively
 * The inner vector therefore includes a list of different template configs, which can produce a curve in a plot such as
 * a ROC Curve (traineff vs. bkg discrimination) or a simple SPC vs. template cnt plot.
 *
 * The outer vector combines several of such curves into one datastructure. For most plotting scripts, this was the
 * easiest way to easily set these parameters.
 *
 * The Analysis Scripts such as the actual Template Bank Building or the Background Analysis simply omit the second
 * dimension of this plot and just go trough every curve vector to process each configuration step by step. In that way,
 * the same datastructure can be used to produce nice plots, but also simulate the data needed for exactly those plots.
 *
 *
 *
 * Unfortunately I was not able to move the whole configuration and parameter framework to config files during the
 * development of the code. In case the framework is going to be used again, I highly encourage you, the developer, to
 * do so if applicable. It's going to save a lot of building and compilation time, as well as helps really keeping
 * track of all of these parameters ;-)
 *
 */


struct TltBnkTrnPlotSet {
    /**
     * Stores a vector of cycles that are to be included in the training overview plots
     */
    std::vector<int> cycle_plotting_order = {1};
    int trainplotoutput_id = 0;
    std::string filelabel = "";

    void setFullCyclePlottingOrder(int plotcount) {
        for(int i = 1; i<=plotcount; i++)  cycle_plotting_order.push_back(i);
    }
};

struct TltBnkFltr {
    /**
     * Class that handles the template filters, that can be easily set by calling
     * TltBnkFlt.addALL().addNO_CENTER().[ETC]
     */
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

    TltBnkFltr &addCENTER_ONLY() {
        filters.push_back(CENTER_ONLY);
        return *this;
    }

    TltBnkFltr &addRECULR_ONLY() {
        filters.push_back(RECURL_ONLY);
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
    /**
     * This function contains the info from several background simulation files. These paramters were obtained manually
     * from the files via the root browser.
     *
     * It is also possible to add a cut, that excludes the frames with the highest nhit numbers in the bkg analysis.
     */

    // Final Plot -- values for beam rate, nhit_mean and nhit_sigma was obtained manually
    std::vector<int> bg_runs =        {125,   122,    126,   127};
    std::vector<float> beam_rates =   {5e7,   1e8,    2e8,   2.466e8};
    std::vector<float> nhit_mean =    {14.5,  24.3,   48.3,  59.86};
    std::vector<float> nhit_sigma =   {17.65, 22.9,   32.1,  36.07};
    std::vector<int> max_bkg_frames = {249900, 150000,  99900, 80000 /* 70000 for 50% */};

    std::vector<int> make_max_nhit_cut = {0,1}; // yes no

    std::vector<int> max_nhits; //calculated automatically
    int ndatapoints;

    InfoForBkgRejVsBeamRatePlot(){
        calc_max_nhits(1); // sigmas = where to cut of the upper nhit frames
        ndatapoints = bg_runs.size();
    };

    void calc_max_nhits(float sigmas) {
      for(int i=0; i<nhit_mean.size(); i++) {
          max_nhits.push_back((int) (nhit_mean[i] + (float) sigmas * nhit_sigma[i]));
      }
    }

};

struct DatabaseConfigBuild {
    /**
     * This class can be initialized by either defining SPR and SPC or by defining Wbins and Zbins.
     * The additional parameters will be calculated from it.
     */

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

    // the basic meta parameters are set to a default value
    int dataset = 0;
    int mode = 0;
    bool prints = false;

    // these vectors are used in main_CTbuildMany.cpp
    std::vector<float> sp_ratios;
    std::vector<int> sp_cnt;
    std::vector<float> stopping_effs;

    // some more basic paramters with default value
    int cosmic_testing_dataset = 30; // for testing cosmic eff with external data
    int max_cosmic_events = 0; // number of events taken from cosmic_testing_dataset (0 = complete file)
    int background_run = 107; // run for bkgAna
    int max_bg_frames = 0; // frames to use from background_run (0 = complete file)
    int max_bg_frame_nhits = 0; // frames from background_run file with more than these nhits, are excluded. (0 = use complete)

    // change to your own filesystem path, where simulation files are stored
    std::string pathtosimfiles = "../../../../../../../../../Volumes/Extreme SSD KN/MAC EXTENSION/Mu3e_data/ServerData/data/";

    // give a comment on the trainplotoutput_id (optional)
    std::string set_description;

    // two methods that manage background file settings. Update values to your own!!
    Configuration &michelBackground() {background_run = 107; max_bg_frames=99900; return *this;}
    Configuration &fullBackground() {background_run = 109; return *this;}

    // Some more variables that store more complex data.
    TltBnkTrnPlotSet TrainPlots; // This function holds a vector containing the cycles of training plots that are then comobined in the training plots
    TltBnkFltr TmplBankFilter; // Template bank filters, that are used
    InfoForBkgRejVsBeamRatePlot BkgFiles; // Data and files that are used for the eff vs. beam rate plot

    // very important structure, vector (curves) of vectors (datapoints per curve) used to train and eval specific
    // configs but also to produce the plots.
    std::vector<std::vector<DatabaseConfigBuild>> DBconfigCurveDatapoints;


    //// SETTING CALLS OF DIFFERENT SCRIPTS
    /* These calls and parameters settings should be moved to a config file, which can then be modified, so that
     * the code must not be compiled every time. The name of the member function could be used as a tag, that would then
     * also run the corresponding script like:
     *
     * [BUILDDB_DATAPOINTS]
     * dataset=13
     * id=0
     *
     * CURVE
     * wbins=64 zbins=4 stopp_eff=0.5
     * wbins=64 zbins=4 stopp_eff=0.6
     *
     * CURVE
     * wbins=64 zbins=8 stopp_eff=0.5
     * wbins=64 zbins=8 stopp_eff=0.6
     *
     * END
     *
     *
     */


    void BUILDTB() {
        /**
         * This implementation and the corresponding script main_CTbuildMany.cpp (BuildDBParamCombinations) was used
         * in the beginning. It is highly recommended to used the approach of main_CTbuildDatapoints.cpp
         * (BuildDBMultiConfigs).
         *
         * This implementation is still useful, when a lot of different configs need to be simulated.
         */

        //Cosmic Dataset
        dataset = 13;
        TrainPlots.trainplotoutput_id = 0;

        //SPM parameters, that are to be used to train a template database for each of their combinations
        stopping_effs.push_back(1); // simulate with every training data present in dataset (possible: train_eff < 1)
        sp_ratios.push_back(64);
        sp_cnt.push_back(1024);
    }

    void BUILDTB_DATAPOINTS() {
        /**
         * This function is called by main_CTbuildDatapoints.cpp (BuildDBMultiConfigs)
         * As it uses the DBconfigCurveDatapoints vector and performs a training for every config setting in it, one
         * can also just call BGANA_DATAPOINTS() in here.
         *
         * By using DatabaseConfigBuild() it is possible to either specify wbins and zbins or spc and spr.
         */

        //Cosmic Dataset
        dataset=13;
        TrainPlots.trainplotoutput_id = 0;

        std::vector<DatabaseConfigBuild> plot_datapoints;

        //using wbin and zbin infos by using (int, int, float) params
        plot_datapoints.push_back(DatabaseConfigBuild(64,4,0.5));
        plot_datapoints.push_back(DatabaseConfigBuild(64,4,0.6));

        //same config, but using spr and spc by using (int, float, float) params
        plot_datapoints.push_back(DatabaseConfigBuild(256, (float) 16.0, 0.7));

        //add more configs that are to be build here ....

        //push back the added datapoints as a curve
        //note: for training the curve / datapoint hierarchy is not relevant. Everything inside DBconfigCurveDatapoints is going to be build.
        DBconfigCurveDatapoints.push_back(plot_datapoints);


        // It is also possible instead of the above code, just to call the settings for a specific plotting script,
        // so that every data necessary for that script is just simulated, for example just call one of these:

//        this->PLOT_BGANA_ROC_DATAPOINTS();
//        this->PLOT_BGANA_SPC_DATAPOINTS();
//        this->PLOT_BGANA_SPR_DATAPOINTS();
    }

    void BGANA_DATAPOINTS() {
        /**
         * Called by 'main_CTBkgAnaDatapoints.cpp' (EvalBkgDBMulti)
         * It uses the same parameter format as this->BUILDTB_DATAPOINTS()
         */

        //basics
        dataset=13;
        TrainPlots.trainplotoutput_id = 0;

        //sets background run number
        this->michelBackground();

        //activates three filters, that are simulated for each of the below combination
        TmplBankFilter.addALL().addNO_CENTER().addCUT_ON_FREQ();

        this->BUILDTB_DATAPOINTS();
    }

    void BGANA_MULTI() {
        /**
         * Called by 'main_CTBkgAnaMany.cpp' (EvalBkgDBParamCombinations)
         * This uses the same parameter format as this->BUILDTB()
         * Therefore, this function can be called to set them.
         */

        dataset=13;
        TrainPlots.trainplotoutput_id = 0;
        this->michelBackground();
        TmplBankFilter.addALL().addNO_CENTER().addCUT_ON_FREQ();

        this->BUILDTB();
    }


    //// TRAINING PLOTS #################################

    void PLOT_BUILDTB(){
        /**
         * Called by main_CTtrainplots.cpp, which produces a couple of train overview plots, such as
         * traineff vs. trainevents / traineff vs. templcnt / templcnt vs. trainevents
         *
         * Variables to be specified:
         * - dataset
         * - trainplotoutput_id
         * - root file cycles to include in the plot -> check file first in TBrowser to produce this plot
         */

        dataset=13;
        TrainPlots.trainplotoutput_id=10;

        TrainPlots.cycle_plotting_order.push_back(1);
        // etc...
    }

    void PLOT_TRAIN_TEMPL_SPC() {
        /**
         * Called by 'main_CTtrainTemplSPC.cpp' (PlotsBuildDBTemplSPC), which produces templ cnt vs. spc plot and
         * fits it with a quadratic function.
         */
        dataset = 13;

        float stopping_eff = 0.6;
        std::vector<DatabaseConfigBuild> plot_datapoints;
        plot_datapoints.push_back(DatabaseConfigBuild(256, 2, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 8, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 12, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 16, stopping_eff));
        DBconfigCurveDatapoints.push_back(plot_datapoints); plot_datapoints.clear();

        stopping_eff = 0.8;
        plot_datapoints.push_back(DatabaseConfigBuild(256, 2, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 8, stopping_eff));
        DBconfigCurveDatapoints.push_back(plot_datapoints); plot_datapoints.clear();
    }

    //// BACKGROUND EVAL PLOTS ########################

    void PLOT_BGANA_SPC_DATAPOINTS() {
        /**
         * Called by 'bgAnaPlots_EffSPC.cpp'
         * This creates the datapoints for the bkg rejection vs. SPC and #templ vs SPR
         */

        this->resetMembers();
        this->michelBackground();
        dataset = 13;
        TrainPlots.trainplotoutput_id=0;
        TmplBankFilter.CLR().addALL();
        max_bg_frame_nhits=0;

        float stopping_eff = 0.6;
        std::vector<DatabaseConfigBuild> plot_datapoints;
        plot_datapoints.push_back(DatabaseConfigBuild(256, 2, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 8, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 12, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 16, stopping_eff));
        DBconfigCurveDatapoints.push_back(plot_datapoints); plot_datapoints.clear();

        stopping_eff = 0.8;
        plot_datapoints.push_back(DatabaseConfigBuild(256, 2, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(256, 8, stopping_eff));
//        plot_datapoints.push_back(DatabaseConfigBuild(256, 12, stopping_eff));
//        plot_datapoints.push_back(DatabaseConfigBuild(256, 16, stopping_eff));
        DBconfigCurveDatapoints.push_back(plot_datapoints); plot_datapoints.clear();

    }

    void PLOT_BGANA_SPR_DATAPOINTS() {
        /**
         * Called by 'bgAnaPlots_EffSPR.cpp'
         * This creates the datapoints for the bkg rejection vs. spr and #templ / spr plot
         */

        this->resetMembers();
        this->michelBackground();
        dataset = 13;
        TrainPlots.trainplotoutput_id=0;
        TmplBankFilter.CLR().addALL();

        // set up the data for two curves, one at traineff=0.6 and one at traineff=0.8
        float stopping_eff = 0.6;
        std::vector<DatabaseConfigBuild> plot_datapoints;
        plot_datapoints.push_back(DatabaseConfigBuild(256, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(128, 8, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(64, 16, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(32, 32, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(16, 64, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(8, 128, stopping_eff));
        DBconfigCurveDatapoints.push_back(plot_datapoints);
        plot_datapoints.clear();

        stopping_eff = 0.8;
        plot_datapoints.push_back(DatabaseConfigBuild(256, 4, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(128, 8, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(64, 16, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(32, 32, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(16, 64, stopping_eff));
        plot_datapoints.push_back(DatabaseConfigBuild(8, 128, stopping_eff));
        DBconfigCurveDatapoints.push_back(plot_datapoints);
    }

    void PLOT_BGANA_ROC_DATAPOINTS() {
        /**
         * Called by 'bgAnaPlots_EffROC.cpp'
         * This creates the datapoints for the bkg rejection vs. spr and #templ / spr plot
         */

        this->resetMembers();
        this->michelBackground();
        dataset = 13;
        TrainPlots.trainplotoutput_id=0;
        TmplBankFilter.CLR().addALL();
        max_bg_frame_nhits=0;

        // set up three curves, wbins x zbins @ 128x4 256x4 256x8
        // traineff reaching from 10% to 90%
        std::vector<DatabaseConfigBuild> plot_datapoints;
        DBconfigCurveDatapoints.push_back(plot_datapoints);
        DBconfigCurveDatapoints.push_back(plot_datapoints);
        DBconfigCurveDatapoints.push_back(plot_datapoints);

        for(int i=1; i<10; i++ ) {
            std::cout << "adding datapoints " << i*0.1 << std::endl;
            DBconfigCurveDatapoints[0].push_back(DatabaseConfigBuild(128, 4, 0.1*i));
            DBconfigCurveDatapoints[0].push_back(DatabaseConfigBuild(256, 4, 0.1*i));
            DBconfigCurveDatapoints[2].push_back(DatabaseConfigBuild(256, 8, 0.1*i));
        }
    }

    void PLOT_BGANA_BEAMRATES() {
        /**
         * Called by 'bgAnaPlots_EffBeamRate.cpp'
         * Uses the data manually added to BkgFiles struct.
         * Produces bkg rejection with higher beamrates plot.
         */

        dataset=13;

        sp_cnt.push_back(2048);
        sp_ratios.push_back(32);
        stopping_effs.push_back(0.8);

        BkgFiles.make_max_nhit_cut.clear();
        BkgFiles.make_max_nhit_cut.push_back(0);
        BkgFiles.make_max_nhit_cut.push_back(1);

        TmplBankFilter.filters.clear();
        TmplBankFilter.addALL().addNO_CENTER();
    }


    void BGANA_PLOT_ROC() {
        /**
         * WARNING: This code uses the old (non-datapoint oriented) parameter format
         * Called by 'bgAnaPlots_legacyROC.cpp'.
         * Using the data format in spr, spc and stopp_eff vectors
         */

        this->resetMembers();
        dataset = 13;
        background_run = 107;
        max_bg_frames = 99900;
        max_bg_frame_nhits = 100;
        TrainPlots.trainplotoutput_id = 0;
        TmplBankFilter.CLR().addALL().addNO_CENTER().addCUT_ON_FREQ();
        set_description = "fltr_ready";

        sp_ratios.push_back(64);
        sp_cnt.push_back(1024);
        stopping_effs.push_back(0.6);
        stopping_effs.push_back(0.8);
    }


    void BGANA_PLOT_SPC() {
        /**
         * WARNING: This code uses the old (non-datapoint oriented) parameter format
         * Called by 'bgAnaPlots_legacySPC.cpp'.
         * Using the data format in spr, spc and stopp_eff vectors
         */

        this->resetMembers();
        dataset = 12;
        background_run = 107;
        max_bg_frames = 99900;
        max_bg_frame_nhits = 100;
        TmplBankFilter.addALL();
        set_description = "stripes_128_32";

        //two curves for spr=32, 128
        sp_ratios.push_back(32);
        sp_ratios.push_back(128);

        //each curve has got three datapoints at stopping_eff=50%
        sp_cnt.push_back(1152);
        sp_cnt.push_back(2048);
        sp_cnt.push_back(3200);
        stopping_effs.push_back(0.5);

    }

    //// These function are used in many scripts to compute consistent file names for PDF plots.
    // Warning: These do not take into accout the datapoint vectors, when they were used. So its partially legacy code.

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

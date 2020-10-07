//
// Created by Konstantin Neureither on 23.06.20.
//

#ifndef COSMICTRIGGER_TEMPLATEBANK_H
#define COSMICTRIGGER_TEMPLATEBANK_H

#include <map>
#include <queue>
#include "TemplateData.h"
#include "SPCalculations.h"

//typedef unsigned long long temid;
typedef TemplateID temid;
typedef std::map<temid, TemplateData> AssociativeMemory;
typedef std::pair<temid, TemplateData> AssociativePair;
typedef std::map<temid, int> CheckedMemory;

//to distinquish between TID types (area dependency)
enum TrackType {RLRL, RLCE, CECE, RRCE, RRRR};
static std::string enum_to_string(TrackType tracktype) {
    switch(tracktype) {
        case RLRL: return "RLRL";
        case RLCE: return "RLCE";
        case CECE: return "CECE";
        case RRCE: return "RRCE";
        case RRRR: return "RRRR";
        default: return "Invalid TrackType!";
    }
}

//filter when loading tids from file
enum TIDLoadingFilter {CENTER_ONLY, RECURL_ONLY, MIXED_ONLY, NO_CENTER, CUT_ON_FREQ, ALL};
static std::string enum_to_string(TIDLoadingFilter filter) {
    switch(filter) {
        case CENTER_ONLY: return "CENTER_ONLY";
        case RECURL_ONLY: return "RECURL_ONLY";
        case MIXED_ONLY: return "MIXED_ONLY";
        case NO_CENTER: return "NO_CENTER";
        case CUT_ON_FREQ: return "CUT_ON_FREQ";
        case ALL: return "ALL";
        default: return "Invalid TIDLoadingFilter!";
    }
}

struct tidQueueNode {
    temid TID; //the template
    unsigned int frequency; //population
    bool operator<(const tidQueueNode& rhs) const {
        return frequency < rhs.frequency;
    }
};

class TemplateBank {
private:
    //stores cosmic templates
    AssociativeMemory AMem;
    //stores found matches with data that is checked
    CheckedMemory CMem;
    //some data structure that stores population of templates: THF1 ?
    SPCalculations SPC;

    //for creating templates
    unsigned int newtemplatecount = 0;
    unsigned int eventcount = 0;
    unsigned int matchedtemplatecount = 0;

    // vectors contain data for training stat graphs.
    // WARNING: when loading templates from file and using a filter, these will still give the "original" data
    //          (wihtout taking into account the filter)
    std::vector<float> Nevents;
    std::vector<float>Ntemplates;
    std::vector<float> efficiency;
    std::vector<int> hitorder;

    //keep track of template types when loading from file with filter
    unsigned int count_loaded_template_types[5] = {0,0,0,0,0};
    TIDLoadingFilter loading_filter = ALL;

    //for checking templates (background eval)
    unsigned int rejectedcount = 0;
    unsigned int acceptedcount = 0;

    std::string plottingpath;

    void fillTemplateFromDB(temid TID, int frequency);
    void initializeMembers(int dataset, int mode, int wBins, int zBins);

    float stopping_efficiency;
    bool PRINTS = false;

public:

    //FIXME Make private again (usage in CosmicTemplatesBgEval)
    //for some meta information
    int mywbins;
    int myzbins;
    int mydataset;
    int mymode;

    TemplateBank(std::string plottingpath, int dataset, int mode, int wBins, int zBins);
    TemplateBank(std::string plottingpath, float stopping_efficiency, int dataset, int mode, int wBins,
                 int zBins);
    ~TemplateBank();
    void SetPrints(bool opt);

    void PlotTemplatePopulationHistogram(); //how many with one, two, three, ..., n roads stored
    void PlotTemplateMatchedFreqHistogram(std::string filetag);
    void PlotTemplatePopHistSortedbyFreq();
    void PlotEfficiency();
    void PlotFreqTimesTemplatecount(); //not finished
    void PlotTemplateTypeDistribution();

    bool fillTemplate(unsigned int * SPIDs, int hitcount, float p, float dca, float phi, float theta);
    bool checkTemplate(temid &TID);
    int getRejectedCount();

    int getAcceptedCount();

    void rmSinglePopulatedTemplates();
    std::vector<temid>  getMostPopulatedTemplates(int howmany);
    std::vector<temid>  getMostMatchedTemplates(int howmany);
    template <typename spidtype>
    temid getTemplateID(spidtype *SPIDs, int hitcount);
    unsigned int getSPIDfromTemplateID(temid TID, int index);
    TrackType GetTypeOfTID(temid TID);
    float getEfficiency();

    int getTemplateCount();
    void writeAMtoFile(std::string path, const int *zBins, const int *wBins, char areaDescript[3][8],
                       std::string mode_description);

    bool readAMfromFile(std::string path, float stopping_efficiency, TIDLoadingFilter filter);

    std::string getfileidtag(int format);

    void resetStats();
    //tests
    void testTemplateID();
    void testFill();
    void testCheck();

    void testGetMostPopTemplates();

    int getTrainingEventCount();

    void fillTemplateFromDB(temid TID, int frequency, TIDLoadingFilter filter);
};


#endif //COSMICTRIGGER_TEMPLATEBANK_H

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
    std::vector<float> Nevents;
    std::vector<float>Ntemplates;
    std::vector<float> efficiency;
    std::vector<int> hitorder;

    //for checking templates
    unsigned int rejectedcount = 0;
    unsigned int acceptedcount = 0;

    std::string plottingpath;

    void fillTemplateFromDB(temid TID, int frequency);
    void initializeMembers();

    float stopping_efficiency;

public:

    //FIXME MAke private again
    //for some meta information
    int mywbins;
    int myzbins;
    int mydataset;
    int mymode;

    TemplateBank(std::string plottingpath);
    TemplateBank(std::string plottingpath, float stopping_efficiency);
    ~TemplateBank();
    bool PRINTS = false;

    void displayTemplatePopulationHistogram(std::string filetag); //how many with one, two, three, ..., n roads stored
    void displayTemplateMatchedFreqHistogram(std::string filetag);
    void displayTemplatePopHistSortedbyFreq(std::string filetag);
    void displayEfficiency(std::string filetag);
    void plotFreqTimesTemplatecount(std::string filetag);

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
    float getEfficiency();

    int getTemplateCount();
    void writeAMtoFile(std::string path, const int *zBins, const int *wBins, char areaDescript[3][8],
                       const int &dataset, const int &mode, std::string mode_description);

    bool readAMfromFile(std::string path, int wbins, int zbins, int mode, int dataset, float stopping_efficiency);

    std::string getfileidtag(int format);

    void resetStats();
    //tests
    void testTemplateID();
    void testFill();
    void testCheck();

    void testGetMostPopTemplates();

    int getTrainingEventCount();
};


#endif //COSMICTRIGGER_TEMPLATEBANK_H

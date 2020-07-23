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
    unsigned int templatecount = 0;
    unsigned int matchedtemplatecount = 0;
    std::vector<float> Nevents;
    std::vector<float>Ntemplates;
    std::vector<float> efficiency;
    std::vector<int> hitorder;

    //for checking templates
    unsigned int rejectedcount = 0;
    unsigned int accepetedcount = 0;

    std::string plottingpath;

public:
    TemplateBank(std::string plottingpath);
    ~TemplateBank();
    bool PRINTS = false;

    void displayTemplatePopulationHistogram(std::string filetag); //how many with one, two, three, ..., n roads stored
    void displayEfficiency(std::string filetag);

    void fillTemplate(unsigned int * SPIDs, int count, float p, float dca, float phi, float theta);
    bool checkTemplate(temid &TID);

    void rmSinglePopulatedTemplates();
    std::vector<temid>  getMostPopulatedTemplates(int howmany);

    temid getTemplateID(unsigned int *SPIDs, int count);
    unsigned int getSPIDfromTemplateID(temid TID, int index);

    void writeAMtoFile(std::string path, int *zBins, int *wBins, char areaDescript[3][8],
                       const int &dataset, const int &mode, std::string mode_description);
    void readAMfromFile(std::string path);

    void resetStats();

    //tests
    void testTemplateID();
    void testFill();
    void testCheck();
    void testGetMostPopTemplates();
};


#endif //COSMICTRIGGER_TEMPLATEBANK_H

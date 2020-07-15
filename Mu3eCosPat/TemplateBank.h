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
typedef std::map<temid, std::vector<TemplateData>> AssociativeMemory;

struct tidQueueNode {
    temid TID; //the template
    unsigned int frequency; //population
    bool operator<(const tidQueueNode& rhs) const {
        return frequency < rhs.frequency;
    }
};

class TemplateBank {
private:
    AssociativeMemory AMem;
    //some data structure that stores population of templates: THF1 ?
    SPCalculations SPC;

    unsigned int newtemplatecount = 0;
    unsigned int templatecount = 0;
    unsigned int matchedtemplatecount = 0;
    std::vector<unsigned int> Nevents;
    std::vector<unsigned int>Ntemplates;
    std::vector<float> efficiency;
    std::vector<int> hitorder;

    bool PRINTS = false;

public:
    TemplateBank();
    ~TemplateBank();

    void displayTemplatePopulationHistogram(); //how many with one, two, three, ..., n roads stored
    void fillTemplate(unsigned int * SPIDs, int count, float p, float dca, float phi, float theta);

    void rmSinglePopulatedTemplates();
    std::vector<temid>  getMostPopulatedTemplates(int howmany);

    temid getTemplateID(unsigned int *SPIDs, int count);
    unsigned int getSPIDfromTemplateID(temid TID, int index);

    //tests
    void testTemplateID();
    void testFill();
    void testGetMostPopTemplates();
};


#endif //COSMICTRIGGER_TEMPLATEBANK_H

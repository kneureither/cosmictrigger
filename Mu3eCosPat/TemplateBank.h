//
// Created by Konstantin Neureither on 23.06.20.
//

#ifndef COSMICTRIGGER_TEMPLATEBANK_H
#define COSMICTRIGGER_TEMPLATEBANK_H

#include <map>
#include <queue>
#include "TemplateData.h"

typedef std::map<unsigned long long, std::vector<TemplateData>> AssociativeMemory;

struct tidQueueNode {
    unsigned long long TID; //the template
    unsigned int frequency; //population
    bool operator<(const tidQueueNode& rhs) const {
        return frequency < rhs.frequency;
    }
};

class TemplateBank {
private:
    AssociativeMemory AMem;
    //some data structure that stores population of templates: THF1 ?

    unsigned int newtemplatecount = 0;
    unsigned int templatecount = 0;
    unsigned int matchedtemplatecount = 0;
    std::vector<unsigned int> Nevents;
    std::vector<unsigned int>Ntemplates;
    std::vector<float> efficiency;

public:
    TemplateBank();
    ~TemplateBank();

    void displayTemplatePopulationHistogram(); //how many with one, two, three, ..., n roads stored
    void fillTemplate(unsigned int * SPIDs, int count, float p, float dca, float phi, float theta);
    void fillTemplate(unsigned int * SPIDs, float * xps, float* yps, float* zps, int count, float p, float dca, float phi, float theta);

    void rmSinglePopulatedTemplates();
    std::vector<unsigned long long>  getMostPopulatedTemplates(int howmany);

    unsigned long long getTemplateID(unsigned int *SPIDs, int count);
    unsigned int getSPIDfromTemplateID(unsigned long long TemplateID, int index);

    //tests
    void testTemplateID();
    void testFill();
    void testGetMostPopTemplates();
};


#endif //COSMICTRIGGER_TEMPLATEBANK_H

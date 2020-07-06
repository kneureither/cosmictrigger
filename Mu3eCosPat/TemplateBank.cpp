//
// Created by Konstantin Neureither on 23.06.20.
//
#include <assert.h>

#include "TemplateBank.h"
#include "TemplateData.h"
#include <queue>
#include <TROOT.h>
#include <TH1F.h>
#include "plots.h"


void TemplateBank::fillTemplate(unsigned int *SPIDs, const int count, const float p, const float dca, const float phi, const float theta) {

    //assume that SPIDs are already cleaned (only 4 entries)
    unsigned long long TID = getTemplateID(SPIDs, count);
    TemplateData TD(SPIDs, count, p, dca, phi, theta);

    if(AMem.count(TID) > 0) {
        AMem[TID].push_back(TD);
        printf(" -- appended to entry TID=%#018llX, TID count=%lu\n", TID, AMem[TID].size());
        matchedtemplatecount++;
    } else {
        std::vector<TemplateData> TDlist;
        TDlist.push_back(TD);
        AMem[TID] = TDlist;
        printf(" -- added new entry TID=%#018llX\n", TID);
        newtemplatecount++;
    }

    if(templatecount / 100 == 0) {
        float efficiency = 1.0 - (newtemplatecount - Nevents[Nevents.size() - 1]) / (float) (templatecount - Ntemplates[Ntemplates.size()-1]);

        this->Nevents.push_back(matchedtemplatecount + newtemplatecount);
        this->Ntemplates.push_back(newtemplatecount);
        this->efficiency.push_back(efficiency);
    }
}

void TemplateBank::fillTemplate(unsigned int *SPIDs, float *xps, float *yps, float *zps, int count, float p, float dca,
                                float phi, float theta) {

}

unsigned long long TemplateBank::getTemplateID(unsigned int *SPIDs, const int count) {
    assert(count == 4);
    unsigned long long templateID; //uint64_t
    templateID = 0x00;

    for(int i = 0; i<count; i++) {
        templateID |= (unsigned long long) (SPIDs[i] & 0xFFFF) << (16 * ((count - 1) - i));
    }
    return templateID;
}

unsigned int TemplateBank::getSPIDfromTemplateID(unsigned long long TemplateID, int index) {
    assert(0 <= index && index < 4);
    return (unsigned int) (TemplateID>>(16*(3-index))) &0xFFFF;
}

void TemplateBank::testTemplateID() {
//    printf("1 & 0xFFFF = %#018llX\n", (unsigned long long) (1 & 0xFFFF) << 48);
    unsigned int SPIDs[4] = {10465, 20345, 20845, 10245};
    for(int i=0; i<4; i++) printf("SPIDs=%d", SPIDs[i]);
    unsigned long long TID = getTemplateID(SPIDs, 4);
    printf("Template ID dec=%llu hex=%#018llX\n", TID, TID);

    for(int i = 0; i<4; i++) {
        unsigned int SPID = getSPIDfromTemplateID(TID, i);
        printf("SID[%d]=%d", i, SPID);
        assert(SPID == SPIDs[i]);
    }
    printf("\n");
}

TemplateBank::TemplateBank() {
    this->Nevents.push_back(0);
    this->Ntemplates.push_back(0);
    this->efficiency.push_back(1.0);
}

TemplateBank::~TemplateBank() {
    //TODO call destructors of class members
}

void TemplateBank::testFill() {
    unsigned int SPIDs1[4] = {1, 2, 3, 4};
    unsigned int SPIDs2[4] = {10465, 20345, 20845, 10245};
    unsigned int SPIDs3[4] = {234, 1231, 464, 3455};
    unsigned int SPIDs4[4] = {3, 4, 5, 6};

    fillTemplate(SPIDs1, 4, 1.0, 1.0, 2.0, 2.0);
    fillTemplate(SPIDs1, 4, 1.0, 1.0, 2.0, 2.0);
    fillTemplate(SPIDs1, 4, 1.0, 1.0, 2.0, 2.0);
    fillTemplate(SPIDs1, 4, 1.0, 1.0, 2.0, 2.0);
    fillTemplate(SPIDs2, 4, 1.0, 1.0, 2.0, 2.0);
    fillTemplate(SPIDs2, 4, 1.0, 1.0, 2.0, 2.0);
    fillTemplate(SPIDs2, 4, 1.0, 1.0, 2.0, 2.0);
    fillTemplate(SPIDs3, 4, 1.0, 1.0, 2.0, 2.0);
    fillTemplate(SPIDs4, 4, 1.0, 1.0, 2.0, 2.0);
    fillTemplate(SPIDs4, 4, 1.0, 1.0, 2.0, 2.0);
    printf("- AMem filled with test data!\n");
}

std::vector<unsigned long long> TemplateBank::getMostPopulatedTemplates(int howmany) {
    assert(howmany <= this->AMem.size());
    AssociativeMemory::iterator it;
    std::priority_queue<tidQueueNode> templQueue;

    unsigned long long TID;
    unsigned int frequency;

    for(it = AMem.begin(); it != AMem.end(); it++){
        TID = it->first;
        frequency = (it->second).size();
        templQueue.push(tidQueueNode{TID, frequency});
    }

    std::vector<unsigned long long> priorityTemplates;
    for(int i=0; i<howmany; i++) {
        priorityTemplates.push_back(templQueue.top().TID);
        templQueue.pop();
    }
    return priorityTemplates;
}

void TemplateBank::displayTemplatePopulationHistogram() {
    auto *canvas = new TCanvas("template frequency", "template frequency", 1200, 900);
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.15);
    canvas->SetGrid(1,1);
    canvas->SetTicks(1, 1);

    TH1F * h_templfreq = new TH1F("h", "Template frequency distribution", 30, 0, 10);
    h_templfreq->SetStats(true);
    labelAxis(h_templfreq, "frequency of templates", "count");
    AssociativeMemory::iterator it;


    unsigned int frequency;

    for(it = AMem.begin(); it != AMem.end(); it++){
        frequency = (it->second).size();
        h_templfreq->Fill(frequency);
    }
    h_templfreq->Draw();
    saveCanvas(canvas, "templateFrequency", "plots/Mu3eCosPat");
}

void TemplateBank::testGetMostPopTemplates() {
    testFill();
    std::vector<unsigned long long> priorityTemplates;

    int howmany=3;
    priorityTemplates = getMostPopulatedTemplates(howmany);
    for(int i=0; i<priorityTemplates.size(); i++) {
        printf("First run (howmany=%d) priorityTemplates[%d]=%#018llX\n", howmany, i, priorityTemplates[i]);
    }

    howmany=1;
    priorityTemplates = getMostPopulatedTemplates(howmany);
    for(int i=0; i<priorityTemplates.size(); i++) {
        printf("second run (howmany=%d) priorityTemplates[%d]=%#018llX\n", howmany, i, priorityTemplates[i]);
    }

//    howmany=10; //should throw an error
//    priorityTemplates = getMostPopulatedTemplates(howmany);
}



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

    //assume that SPIDs are already cleaned (only two hits per entry)
    temid TID = getTemplateID(SPIDs, count);
    TemplateData TD(count, p, dca, phi, theta);

    if(AMem.count(TID) > 0) {
        AMem[TID].push_back(TD);
        std::cout << " -- appended to entry TID=" + TID.toString() + " count=" << AMem[TID].size() << std::endl;
        matchedtemplatecount++;
    } else {
        std::vector<TemplateData> TDlist;
        TDlist.push_back(TD);
        AMem[TID] = TDlist;
        std::cout << " -- added new entry TID=" + TID.toString() << std::endl;
        newtemplatecount++;
    }

    if(templatecount / 100 == 0) {
        float efficiency = 1.0 - (newtemplatecount - Nevents[Nevents.size() - 1]) / (float) (templatecount - Ntemplates[Ntemplates.size()-1]);

        this->Nevents.push_back(matchedtemplatecount + newtemplatecount);
        this->Ntemplates.push_back(newtemplatecount);
        this->efficiency.push_back(efficiency);
    }
}

temid TemplateBank::getTemplateID(unsigned int *SPIDs, const int count) {
    assert(count <= TID_LEN);
    temid TID;
    unsigned short SPID;
    int sidindex=0;

    for(int i = 0; i<TID_LEN; i++) {
        SPID = SPIDs[sidindex];
        if(PRINTS) printf("SPID to be added=%d ", SPID);
        if(SPC.getLayerFromSPID(SPID) == hitorder[i]) {
            if(PRINTS) printf(" -- added sidindex=%d as TID index=%d\n", sidindex, i);
            TID.HIDS[i] = SPID;
            sidindex++;
        } else {
            if(PRINTS) printf(" -- set to 0\n");
            TID.HIDS[i] = 0;
        }
    }
    if(PRINTS) printf("\n count=%d, sidindex=%d \n", count, sidindex);
    assert(count == sidindex);
    return TID;
}

unsigned int TemplateBank::getSPIDfromTemplateID(temid TID, int index) {
    assert(0 <= index && index < TID_LEN);
    return (unsigned int) TID.HIDS[index];
}

void TemplateBank::testTemplateID() {
//    printf("1 & 0xFFFF = %#018llX\n", (unsigned long long) (1 & 0xFFFF) << 48);
    unsigned int SPIDs[4] = {10465, 20345, 20845, 10245};
    for(int i=0; i<4; i++) printf("SPIDs=%d", SPIDs[i]);
    temid TID = getTemplateID(SPIDs, 4);
    std::cout << "Template ID hex=" + TID.toString() << std::endl;

    for(int i = 0; i<4; i++) {
        unsigned int SPID = getSPIDfromTemplateID(TID, i);
        printf("SID[%d]=%d  ", i, SPID);
        assert(SPID == SPIDs[i]);
    }
    printf("\n");
}

TemplateBank::TemplateBank() {
    this->Nevents.push_back(0);
    this->Ntemplates.push_back(0);
    this->efficiency.push_back(1.0);

    assert(TID_LEN % 2 == 0);
    hitorder.push_back(3);
    hitorder.push_back(2);
    if(TID_LEN > 4) {
        hitorder.push_back(1);
        if(TID_LEN > 6) {
            hitorder.push_back(0);
            hitorder.push_back(0);
        }
        hitorder.push_back(1);
    }
    hitorder.push_back(2);
    hitorder.push_back(3);

    printf("asserting hitorder.size() == TID_LEN : %d == %d\n", hitorder.size(), TID_LEN);
    assert(hitorder.size() == TID_LEN);
}

TemplateBank::~TemplateBank() {
    //TODO call destructors of class members
}

void TemplateBank::testFill() {
    unsigned int SPIDs1[4] = {3, 2, 2, 3};
    unsigned int SPIDs2[4] = {0x1243, 0x2342, 0x3462, 0x2343};
    unsigned int SPIDs3[8] = {0x0017, 0x0026, 0x0035, 0x0044, 0x0054, 0x0065, 0x0076, 0x0087};
    unsigned int SPIDs4[6] = {7, 6,5, 5, 6, 7};
    unsigned int SPIDs5[4] = {0x0017, 0x0026, 0x0076, 0x0087};
    unsigned int SPIDs6[4] = {7, 6, 6, 7};

    fillTemplate(SPIDs1, 4, 1.0, 1.0, 2.0, 2.0);
    fillTemplate(SPIDs1, 4, 1.0, 1.0, 2.0, 2.0);
    fillTemplate(SPIDs1, 4, 1.0, 1.0, 2.0, 2.0);
    fillTemplate(SPIDs1, 4, 1.0, 1.0, 2.0, 2.0);
    fillTemplate(SPIDs2, 4, 1.0, 1.0, 2.0, 2.0);
    fillTemplate(SPIDs2, 4, 1.0, 1.0, 2.0, 2.0);
    fillTemplate(SPIDs2, 4, 1.0, 1.0, 2.0, 2.0);
    if(TID_LEN == 8) {
        fillTemplate(SPIDs3, 8, 1.0, 1.0, 2.0, 2.0);
        fillTemplate(SPIDs4, 6, 1.0, 1.0, 2.0, 2.0);
        fillTemplate(SPIDs4, 6, 1.0, 1.0, 2.0, 2.0);
    } else {
        fillTemplate(SPIDs5, 4, 1.0, 1.0, 2.0, 2.0);
        fillTemplate(SPIDs6, 4, 1.0, 1.0, 2.0, 2.0);
        fillTemplate(SPIDs6, 4, 1.0, 1.0, 2.0, 2.0);
    }
    printf("- AMem filled with test data!\n");
}

std::vector<temid> TemplateBank::getMostPopulatedTemplates(int howmany) {
    assert(howmany <= this->AMem.size());
    AssociativeMemory::iterator it;
    std::priority_queue<tidQueueNode> templQueue;

    temid TID;
    unsigned int frequency;

    for(it = AMem.begin(); it != AMem.end(); it++){
        TID = it->first;
        frequency = (it->second).size();
        templQueue.push(tidQueueNode{TID, frequency});
    }

    std::vector<temid> priorityTemplates;
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
    std::vector<temid> priorityTemplates;

    int howmany=3;
    priorityTemplates = getMostPopulatedTemplates(howmany);
    for(int i=0; i<priorityTemplates.size(); i++) {
        printf("First run (howmany=%d) priorityTemplates[%d]=%s\n", howmany, i, priorityTemplates[i].toString().c_str());
    }

    howmany=1;
    priorityTemplates = getMostPopulatedTemplates(howmany);
    for(int i=0; i<priorityTemplates.size(); i++) {
        printf("second run (howmany=%d) priorityTemplates[%d]=%s\n", howmany, i, priorityTemplates[i].toString().c_str());
    }

//    howmany=10; //should throw an error
//    priorityTemplates = getMostPopulatedTemplates(howmany);
}



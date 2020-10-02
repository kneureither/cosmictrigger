//
// Created by Konstantin Neureither on 23.06.20.
//
#include <assert.h>
#include <numeric>

#include "TemplateBank.h"
#include "TemplateData.h"
#include "utilityFunctions.h"
#include <queue>
#include <TROOT.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TTree.h>
#include "plots.h"
#include "TemplateDatabaseFile.h"


bool TemplateBank::fillTemplate(unsigned int *SPIDs, const int hitcount, const float p, const float dca, const float phi, const float theta) {
    //returns 0 if stopping point reached, else 1

    //assume that SPIDs are already cleaned (only two hits per entry)
    temid TID = getTemplateID(SPIDs, hitcount);

    if(AMem.count(TID) > 0) {
        AMem[TID].add_track(hitcount, p, dca, phi, theta);
//        std::cout << " -- appended to entry TID=" + TID.toString() + " hitcount=" << AMem[TID].size() << std::endl;
        matchedtemplatecount++;
    } else {
        TemplateData TD = TemplateData(hitcount, p, dca, phi, theta);
        AMem[TID] = TD;
//        std::cout << " -- added new entry TID=" + TID.toString() << std::endl;
        newtemplatecount++;
    }

    eventcount = matchedtemplatecount + newtemplatecount;
    if(((eventcount % (int) pow((float) 10, (float) std::floor(log10(eventcount))) == 0) && (eventcount >= 1000))
            || (eventcount >= 1000000 && eventcount % 100000 == 0)) {

        float efficiency = 1.0 - (newtemplatecount - Ntemplates[Ntemplates.size() - 1]) / (float) (eventcount - Nevents[Nevents.size() - 1]);

        this->Nevents.push_back((float) eventcount);
        this->Ntemplates.push_back((float) newtemplatecount);
        this->efficiency.push_back(efficiency);

        if(efficiency >= this->stopping_efficiency) return true; //stopping efficiency reached
    }

    return false; //efficiency not reached or not checked
}

void TemplateBank::fillTemplateFromDB(temid TID, int frequency) {
    TemplateData TD = TemplateData(frequency);
    AMem[TID] = TD;
    newtemplatecount++;
    matchedtemplatecount += (frequency - 1);
    eventcount += frequency;
}

TrackType TemplateBank::GetTypeOfTID(temid TID) {
    //return type of TID (recurlL, recurlL), (recurlL, center), (center, center) (center, rcurlR) (recurlR, recurlR)
    int tidareas[TID_LEN];
    int type = 0;
    int area = 0;
    for(int i=0; i<TID_LEN; i++) {
        area = SPC.getAreaFromSPID(TID.HIDS[i]);
        type += pow(10, area);
    }

    if(type == TID_LEN) {
        return CECE;
    } else if(type == 10 * TID_LEN) {
        return RRRR;
    } else if(type == 100*TID_LEN) {
        return RLRL;
    } else if(10 < type && type < 100) {
        return RRCE;
    } else if(100 < type && type < 1000) {
        return RLCE;
    } else {
        std::cout << "(WARING): Something went wrong it TID type classification" << std::endl;
        exit(0);
    }
}

void TemplateBank::fillTemplateFromDB(temid TID, int frequency, TIDLoadingFilter filter) {
    TrackType tidAreaType = this->GetTypeOfTID(TID); //get area type of track

    switch(filter) {
        case ALL:
            break;
        case CENTER_ONLY:
            if(!(tidAreaType == CECE)) return;
            break;
        case RECURL_ONLY:
            if(!(tidAreaType == RRRR || tidAreaType == RLRL)) return;
            break;
        case RECURL_CENTER:
            if(tidAreaType == CECE) return;
            break;
        case CUT_ON_FREQ:
            if(frequency == 1) return;
            break;
        default:
            break;
    }

//    std::cout << "  -> TID: " << TID.toString() << "tid type " << tidAreaType << std::endl;

    TemplateData TD = TemplateData(frequency);
    AMem[TID] = TD;
    newtemplatecount++;
    matchedtemplatecount += (frequency - 1);
    eventcount += frequency;
    this->count_loaded_template_types[tidAreaType]++;
}

bool TemplateBank::checkTemplate(temid &TID) {
    if(AMem.count(TID) == 0) {
        rejectedcount++;
        return false;
    } else {
        acceptedcount++;
        if(CMem.count(TID) > 0) {
            CMem[TID]++;
            if(PRINTS) std::cout << " -- found multi occurence TID=" + TID.toString() + " count=" << CMem[TID] << std::endl;
        } else {
            CMem[TID] = 1;
            if(PRINTS) std::cout << " -- first occurence of TID=" + TID.toString() << std::endl;
        }
        return true;
    }
}

template <typename spidtype>
temid TemplateBank::getTemplateID(spidtype *SPIDs, const int hitcount) {
    assert(hitcount <= TID_LEN);
    temid TID;
    unsigned short SPID;
    int sidindex=0;

    for(int i = 0; i<TID_LEN; i++) {
        SPID = (unsigned short) SPIDs[sidindex];
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
    if(PRINTS) printf("\nasserting... hitcount=%d, sidindex=%d \n", hitcount, sidindex);
    assert(hitcount == sidindex);
    return TID;
}

unsigned int TemplateBank::getSPIDfromTemplateID(temid TID, int index) {
    assert(0 <= index && index < TID_LEN);
    return (unsigned int) TID.HIDS[index];
}

float TemplateBank::getEfficiency() {
    return this->efficiency[this->efficiency.size() - 1];
}

int TemplateBank::getTrainingEventCount() {
    return this->eventcount;
}

int TemplateBank::getTemplateCount() {
    return this->newtemplatecount;
}

void TemplateBank::initializeMembers(int dataset, int mode, int wBins, int zBins) {
    this->Nevents.push_back(0);
    this->Ntemplates.push_back(0);
    this->efficiency.push_back(0.0);

    mywbins = wBins;
    myzbins = zBins;
    mymode = mode;
    mydataset = dataset;

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

//    printf("asserting hitorder.size() == TID_LEN : %lu == %d\n", hitorder.size(), TID_LEN);
    assert(hitorder.size() == TID_LEN);
}

TemplateBank::TemplateBank(std::string plottingpath, int dataset, int mode, int wBins, int zBins) {
    initializeMembers(dataset, mode, wBins, zBins);
    this->plottingpath = plottingpath;
}

TemplateBank::TemplateBank(std::string plottingpath, float stopping_efficiency, int dataset, int mode, int wBins,
                           int zBins) {
    initializeMembers(dataset, mode, wBins, zBins);
    this->plottingpath = plottingpath;
    this->stopping_efficiency = stopping_efficiency;
}

TemplateBank::~TemplateBank() {
    //TODO call destructors of class members
}

void TemplateBank::testTemplateID() {
//    printf("1 & 0xFFFF = %#018llX\n", (unsigned long long) (1 & 0xFFFF) << 48);
    unsigned short SPIDs[4] = {0x1243, 0x2342, 0x3462, 0x2343};
    for(int i=0; i<4; i++) printf("SPIDs=%d ", SPIDs[i]);
    temid TID = getTemplateID(SPIDs, 4);
    std::cout << "Template ID hex=" + TID.toString() << std::endl;

    for(int i = 0; i<4; i++) {
        unsigned int SPID = getSPIDfromTemplateID(TID, i);
        printf("SID[%d]=%d  ", i, SPID);
        assert(SPID == SPIDs[i]);
    }
    printf("\n");
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

void TemplateBank::testCheck() {
    assert(AMem.size() > 0);
    unsigned int SPIDs5[4] = {0x0017, 0x0026, 0x0076, 0x0087};
    unsigned int SPIDs1[4] = {3, 2, 2, 3};
    unsigned int SPIDs2[4] = {0xB, 0x2, 0x2, 0x3};

    temid TID1 = getTemplateID(SPIDs1, 4);
    temid TID2 = getTemplateID(SPIDs2, 4);
    temid TID3 = getTemplateID(SPIDs5, 4);


    assert(checkTemplate(TID1) == true);
    assert(checkTemplate(TID1) == true);
    assert(checkTemplate(TID2) == false);
    assert(checkTemplate(TID3) == true);
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

std::vector<temid> TemplateBank::getMostPopulatedTemplates(int howmany) {
    assert(howmany <= this->AMem.size());
    AssociativeMemory::iterator it;
    std::priority_queue<tidQueueNode> templQueue;

    temid TID;
    unsigned int frequency;

    for(it = AMem.begin(); it != AMem.end(); it++){
        TID = it->first;
        frequency = (it->second).frequency;
        templQueue.push(tidQueueNode{TID, frequency});
    }

    std::cout << " -- Getting the " << howmany << " most populated templates: " << std::endl;

    std::vector<temid> priorityTemplates;
    for(int i=0; i<howmany; i++) {
        priorityTemplates.push_back(templQueue.top().TID);
        frequency = templQueue.top().frequency;
        std::cout << "  *rank["  << i+1 << "] " << priorityTemplates[i].toString() << "  frequency: " << frequency << std::endl;
        templQueue.pop();
    }
    return priorityTemplates;
}

//TODO: Would be nicer to call a generic function e.g. getTopElementsFromMap(int howmany, T* Mem) to avoid code duplication
std::vector<temid>  TemplateBank::getMostMatchedTemplates(int howmany) {
    if(howmany > this->CMem.size()) {
        howmany = this->CMem.size();
        std::cout << " WARNING : howmany > CMem.size(). Setting howmany to CMem.size()=" << howmany << std::endl;
    }

    CheckedMemory::iterator it;
    std::priority_queue<tidQueueNode> templQueue;

    temid TID;
    unsigned int frequency;

    for(it = CMem.begin(); it != CMem.end(); it++){
        TID = it->first;
        frequency = (it->second);
        templQueue.push(tidQueueNode{TID, frequency});
    }

    std::cout << " -- Getting the " << howmany << " most matched templates: " << std::endl;

    std::vector<temid> priorityTemplates;
    for(int i=0; i<howmany; i++) {
        priorityTemplates.push_back(templQueue.top().TID);
        frequency = templQueue.top().frequency;
        std::cout << "  *rank["  << i+1 << "] " << priorityTemplates[i].toString() << "  frequency: " << frequency << std::endl;
        templQueue.pop();
    }
    return priorityTemplates;
}

void TemplateBank::displayTemplateMatchedFreqHistogram(std::string filetag) {
    auto *canvas = new TCanvas("template matching frequency", "template matching frequency", 1200, 900);
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.15);
    canvas->SetGrid(1,1);
    canvas->SetTicks(1, 1);
    canvas->SetLogy(1);

    TH1F * h_templfreq = new TH1F("h", "Template matching freq distribution", 101, -0.5, 100.5);
    h_templfreq->SetName("h_templmatchfreq");
    h_templfreq->SetStats(true);
    labelAxis(h_templfreq, "number of matches per template", "# templates");
    CheckedMemory::iterator it;


    unsigned int frequency;

    for(it = CMem.begin(); it != CMem.end(); it++){
        frequency = (it->second);
        h_templfreq->Fill(frequency);
    }
    h_templfreq->Draw();
    h_templfreq->Write();
    saveCanvas(canvas, ("templateMatchingFrequency_" + filetag).c_str(), plottingpath);
}



void TemplateBank::displayTemplatePopulationHistogram() {
    auto *canvas = new TCanvas("template frequency", "template frequency", 1200, 900);
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.15);
    canvas->SetGrid(1,1);
    canvas->SetTicks(1, 1);
    canvas->SetLogy(1);

    TH1F * h_templfreq = new TH1F("h", "Template frequency distribution", 100, 0, 100);
    h_templfreq->SetName("h_templfreq");
    h_templfreq->SetStats(true);
    labelAxis(h_templfreq, "frequency of templates", "# templates");
    AssociativeMemory::iterator it;


    unsigned int frequency;

    for(it = AMem.begin(); it != AMem.end(); it++){
        frequency = (it->second).frequency;
        h_templfreq->Fill(frequency);
    }
    h_templfreq->Draw();
    h_templfreq->Write();
    saveCanvas(canvas, ("templateFrequency_" + getfileidtag(0)).c_str(), plottingpath);
}

void TemplateBank::plotFreqTimesTemplatecount() {
    auto *canvas = new TCanvas("template weights", "template weights", 1200, 900);
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.15);
    canvas->SetGrid(1,1);
    canvas->SetTicks(1, 1);
    canvas->SetLogy(0);

    TH1F * h_templweights = new TH1F("h", "Template generation weights", 200, 0, 20000);
    h_templweights->SetName("h_templweight");
    h_templweights->SetStats(true);
    labelAxis(h_templweights, "freq * # templates", "# templates");

    TH1F * h_templfreq = new TH1F("h", "Template frequency distribution", 100, 0, 100);
    h_templfreq->SetName("h_templfreq");
    h_templfreq->SetStats(true);
    labelAxis(h_templfreq, "frequency of templates", "# templates");


    AssociativeMemory::iterator it;

    it = std::max_element(AMem.begin(),AMem.end(),[] (const AssociativePair& a, const AssociativePair& b)->bool{ return (a.second).frequency < (b.second).frequency; } );

    if(PRINTS){
        std::cout << "INFO   : most frequent checked TID=";
        std::cout << (it->first).HIDS[0]  << (it->first).HIDS[1] << (it->first).HIDS[2] << (it->first).HIDS[3];
        std::cout << " , freq=" << (it->second).frequency << std::endl;
    }

    unsigned int max_frequency = (it->second).frequency;
    unsigned int frequency;

    AssociativeMemory::iterator it2;
    for(it2 = AMem.begin(); it2 != AMem.end(); it2++){
        frequency = (it2->second).frequency;
        h_templfreq->Fill(frequency);
    }

    unsigned int freqtemplcount;

    for(int i=0; i <= max_frequency; i++) {
        freqtemplcount = (*h_templfreq)[h_templfreq->GetBin(i)];
        if(this->PRINTS) std::cout << "freq: " << i << " bin: " << (*h_templfreq)[h_templfreq->GetBin(i)] << " bin index: " << h_templfreq->GetBin(i) << std::endl;
        h_templweights->Fill(freqtemplcount*i);
    }

//    if(this->PRINTS) std::cout << h_templweights->GetMaximum() << std::endl;

    h_templweights->Draw();
    h_templweights->Write();
    saveCanvas(canvas, ("templateFrequencyTimesCount_" + getfileidtag(0)).c_str(), plottingpath);
}

void TemplateBank::displayEfficiency() {
    auto *canvas = new TCanvas("template bank stats", "cosmic template bank stats", 1200, 900);

    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.15);
    canvas->SetGrid(1,1);
    canvas->SetTicks(1, 1);
    canvas->SetLogx(1);

    auto *pad1 = new TPad("template efficiency", "template efficiency", 0, 0.5, 1, 1);
    pad1->Draw();
    auto *pad2 = new TPad("template count", "template count", 0, 0, 1, 0.5);
    pad2->Draw();

    //Efficiency
    TGraph *g_efficiency = new TGraph( Nevents.size(),&Nevents[0],&efficiency[0]);
    g_efficiency->SetName("g_efficiency");
    g_efficiency->SetTitle("template efficiency");
    labelAxis(g_efficiency, "N events", "efficiency");
    setGraphRange(g_efficiency,100, Nevents[Nevents.size()-1], 0, 1);

    g_efficiency->SetLineColor(kBlue);
    g_efficiency->SetMarkerStyle(23);
    g_efficiency->SetMarkerSize(1);
    pad1->cd();
    g_efficiency->Draw("ALP");
    g_efficiency->Write();

    //template count
    TGraph *g_tnumber = new TGraph( Nevents.size(),&Nevents[0],&Ntemplates[0]);
    g_tnumber->SetName("g_tnumber");
    g_tnumber->SetTitle("template count");
    labelAxis(g_tnumber, "N events", "number of template");
    setGraphRange(g_tnumber,100, Nevents[Nevents.size()-1], 0, Ntemplates[Ntemplates.size()-1]);

    g_tnumber->SetLineColor(kRed);
    g_tnumber->SetMarkerStyle(23);
    g_tnumber->SetMarkerSize(1);
    pad2->cd();
    g_tnumber->Draw("ALP");
    g_tnumber->Write();

    canvas->Update();
    saveCanvas(canvas, ("templateBankStats" + getfileidtag(0)).c_str(), plottingpath);
}

void TemplateBank::writeAMtoFile(std::string path, const int *zBins, const int *wBins, char areaDescript[3][8],
                                 std::string mode_description) {
    //iterate over the Associative Memory map this->AM and write data to root file

    assert(myzbins == zBins[0]);
    assert(mywbins == wBins[0]);

    std::cout << "(INFO)   : Writing AM Template Database to file..." << std::endl;

    std::string customnametag = getfileidtag(0);
    TFile tF((path + "CosmicPatternDatabase_" + customnametag + ".root").c_str(), "recreate");
    if (!tF.IsOpen()) {
        std::cout << "[ERROR] File " << tF.GetName() << " is not open!" << std::endl;
    }

    TTree tT_spconfig("ConfigTree","Tree with Superpixel configuration information");
    TTree tT_tids("TIDTree","Tree with Template IDentification (TID) number");

    TemplateDatabaseWrite TDB = TemplateDatabaseWrite(&tT_spconfig, &tT_tids, mydataset, zBins, wBins, areaDescript,
                                                      mymode, this->efficiency[this->efficiency.size() - 1], eventcount,
                                                      mode_description, this->newtemplatecount, this->stopping_efficiency);

    int tid_len = TID_LEN;
    AssociativeMemory::iterator it;

    int i = 0;

    for(it = AMem.begin(); it != AMem.end(); it++){
        if(i % 1000 == 0) print_status_bar(i, AMem.size(), "writing templates", "");
        i++;
        temid TID(it->first);
        TDB.fillTIDData(TID.HIDS, tid_len, TID.toString(), (it->second).frequency);
    }

    std::string filename = tF.GetName();

    tF.Write();
    tF.Close();

    std::cout << "(INFO)   : CHECK -> Wrote AM Template Database to file " << filename << std::endl;
    std::cout << "(CONFIG) : wBins " << mywbins << " | zBins " << myzbins << " | T count " << newtemplatecount;
    std::cout << " | training eff " << this->getEfficiency() << " | training events " << eventcount << std::endl;
}

bool
TemplateBank::readAMfromFile(std::string path, float stopping_efficiency, TIDLoadingFilter filter) {
    this->stopping_efficiency = stopping_efficiency;

    std::string customnametag = getfileidtag(0);
    std::string filename = path + "CosmicPatternDatabase_" + customnametag + ".root";
    std::cout << "(INFO)   : START -> Getting AM Template Database from file " << filename << std::endl;

    // FILE FOR READING
    TFile tF(filename.c_str());
    if (!tF.IsOpen()) {
        std::cout << "[ERROR] File " << tF.GetName() << " is not open!" << std::endl;
        return 0;
    }

    TTree* tT_tids;
    TTree* tT_spconfig;
    tF.GetObject("TIDTree", tT_tids);
    tF.GetObject("ConfigTree", tT_spconfig);

    TemplateDatabaseRead TDB = TemplateDatabaseRead(tT_spconfig, tT_tids);

    //check if the meta data is okay
    assert(TDB.tid_len == TID_LEN);
    assert(TDB.wBins[0] == mywbins);
    assert(TDB.zBins[0] == myzbins);
    assert(TDB.mode == mymode);
    assert(TDB.dataset == mydataset);

    unsigned int entries = TDB.tT_tid->GetEntries();
    this->efficiency.push_back(TDB.training_efficiency); //not exactly correct, as the efficiency is calculated upt to the last 10e6 value.
    this->Nevents.push_back(TDB.training_events); //not exactly correct, as the Nevents is calculated upt to the last 10e6 value.
    this->Ntemplates.push_back(entries); //correct!

    for(int i=0; i<entries; i++) {
        if(i % 1000 == 0) print_status_bar(i, entries, "reading templates", "");
        TDB.getEntry(i);
        temid TID = getTemplateID(TDB.tid, TID_LEN);
        this->fillTemplateFromDB(TID, TDB.frequency, filter);
    }

    if(filter == ALL) assert(newtemplatecount == entries);

    tF.Close();

    std::cout << "(INFO)   : CHECK -> Got AM Template Database from file " << filename << std::endl;
    std::cout << "(INFO)   : loaded Templ types: RLRL " << count_loaded_template_types[RLRL] << string_perc(count_loaded_template_types[RLRL], entries);
    std::cout << " | RLCE " << count_loaded_template_types[RLCE] << string_perc(count_loaded_template_types[RLCE], entries);
    std::cout << " | CECE " << count_loaded_template_types[CECE] << string_perc(count_loaded_template_types[CECE], entries);
    std::cout << " | RRCE " << count_loaded_template_types[RRCE] << string_perc(count_loaded_template_types[RRCE], entries);
    std::cout << " | RRRR " << count_loaded_template_types[RRRR] << string_perc(count_loaded_template_types[RRRR], entries);
    std::cout << " | used " << string_perc(newtemplatecount, entries) << "of db templates" << std::endl;
    std::cout << "(CONFIG) : wBins " << mywbins << " | zBins " << myzbins << " | T count " << newtemplatecount;
    std::cout << " | training eff " << this->getEfficiency()  << " | train events " << eventcount << std::endl;
}

int TemplateBank::getRejectedCount() {
    return rejectedcount;
}

int TemplateBank::getAcceptedCount() {
    return acceptedcount;
}

std::string TemplateBank::getfileidtag(int format) {
    if(format == 1) {
        //no max eff
        return ::getfileidtag(mydataset, mymode, mywbins, myzbins);
    } else {
        //default
        return ::getfileidtag(mydataset, mymode, mywbins, myzbins, stopping_efficiency);
    }
}


void TemplateBank::displayTemplatePopHistSortedbyFreq() {
    std::vector<temid> priorityTemplates = getMostPopulatedTemplates(10000);
    std::vector<int> templateIndex(priorityTemplates.size());
    std::iota (std::begin(templateIndex), std::end(templateIndex), 0);

    auto *canvas = new TCanvas("templates sorted by population", "templates sorted by population", 1200, 900);

    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.15);
    canvas->SetGrid(1,1);
    canvas->SetTicks(1, 1);

    TGraph *g_efficiency = new TGraph( priorityTemplates.size(),&Nevents[0],&efficiency[0]);
//    TGraph *g_templpop = new TGraph(priorityTemplates.size(), &templateIndex[0], &priorityTemplates[0]);
//    g_templpop->SetName("g_efficiency");
//    g_templpop->SetTitle("template efficiency");
//    labelAxis(g_templpop, "N events", "efficiency");
//    setGraphRange(g_templpop,0, templateIndex[templateIndex.size()-1], 0, 4000);
//
//    g_templpop->Draw("ALP");
//    saveCanvas(canvas, ("TemplatesByPopulation_" + filetag).c_str(), plottingpath);
}


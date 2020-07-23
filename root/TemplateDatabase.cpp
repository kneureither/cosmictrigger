//
// Created by Konstantin Neureither on 15.07.20.
//

#include "TemplateDatabase.h"
#include "../Mu3eCosPat/TemplateData.h"

void TemplateDatabase::reinitializeData() {
    tid_len = -1;
    tid_repr = "default";

    frequency = -1.1;

    for(int i=0; i<this->tid_len; i++) this->tid[i] = 0;

    nhit.clear();
    pt.clear();
    phi.clear();
    theta.clear();
    dca.clear();
}

TemplateDatabaseWrite::TemplateDatabaseWrite(TTree *tT_meta, TTree *tT_tid, const int dataset, const int *zBins, const int *wBins, char areaDescript[3][8],
                                             const int mode, const float efficiency, std::string mode_description) {
    this->tT_meta = tT_meta;
    this->tT_tid = tT_tid;
    this->dataset = dataset;
    for(int i=0; i<this->tid_len; i++) this->tid[i] = tid[i];
    for(int i=0; i<3; i++) {
        this->zBins[i] = zBins[i];
        this->wBins[i] = wBins[i];
        for(int j=0; j<8; j++) this->areaDescript[i][j] = areaDescript[i][j];
    }

    this->mode = mode;
    this->efficiency = efficiency;
    this->mode_description = mode_description;

    this->tid_len = TID_LEN;

    this->tT_meta->Branch("dataset", &this->dataset, "dataset/I");
    this->tT_meta->Branch("area0Description", &this->areaDescript[0], "area0Description/C");
    this->tT_meta->Branch("area1Description", &this->areaDescript[1], "area1Description/C");
    this->tT_meta->Branch("area2Description", &this->areaDescript[2], "area2Description/C");
    this->tT_meta->Branch("wBins0", &this->wBins[0], "wBins0/I");
    this->tT_meta->Branch("wBins1", &this->wBins[1], "wBins1/I");
    this->tT_meta->Branch("wBins2", &this->wBins[2], "wBins2/I");
    this->tT_meta->Branch("zBins0", &this->zBins[0], "zBins0/I");
    this->tT_meta->Branch("zBins1", &this->zBins[1], "zBins1/I");
    this->tT_meta->Branch("zBins2", &this->zBins[2], "zBins2/I");
    this->tT_meta->Branch("mode", &this->mode, "mode/I");
    this->tT_meta->Branch("mode_description", &this->mode_description, "mode_description/C");
    this->tT_meta->Branch("efficiency", &this->efficiency, "efficiency/F");
    this->tT_meta->Fill();


    this->tT_tid->Branch("tid_len", &this->tid_len, "tid_len/I");
    this->tT_tid->Branch("tid", tid, "tid[tid_len]/s");
    this->tT_tid->Branch("tid_repr", &this->tid_repr, "tid_repr/C");
    this->tT_tid->Branch("freq", &this->frequency, "frequency/I");

//really needed?
//    this->tT_tid->Branch("nhit", &this->nhit);
//    this->tT_tid->Branch("p", &this->pt);
//    this->tT_tid->Branch("phi", &this->phi);
//    this->tT_tid->Branch("theta", &this->theta);
//    this->tT_tid->Branch("dca", &this->dca);
}


void TemplateDatabaseWrite::fillTIDData(unsigned short *tid, const int tid_len, std::string tid_repr, const int &freq, std::vector<int> &nhit,
                                     std::vector<float> &p, std::vector<float> &phi, std::vector<float> &theta,
                                     std::vector<float> dca, std::vector<unsigned int> &uEventIDs) {
    this->reinitializeData();
    this->tid_len = tid_len;
    for(int i=0; i<this->tid_len; i++) this->tid[i] = tid[i];
    this->frequency = freq;

    this->pt = p;
    this->phi = phi;
    this->theta = theta;
    this->dca = dca;

    this->tT_tid->Fill();
}

void TemplateDatabaseWrite::fillTIDData(unsigned short *tid, const int &tid_len,std::string tid_repr, const int &freq) {
    this->reinitializeData();
    this->tid_len = tid_len;
    for(int i=0; i<this->tid_len; i++) this->tid[i] = tid[i];
    this->frequency = freq;
    this->tid_repr = tid_repr;

    this->tT_tid->Fill();
}

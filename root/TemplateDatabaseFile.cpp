//
// Created by Konstantin Neureither on 15.07.20.
//

#include "TemplateDatabaseFile.h"
#include "../CTCoreModules/inc/TemplateData.h"

void TemplateDatabaseFile::reinitializeData() {
    tid_len = -1;
    tid_repr = "0000000000000000";

    frequency = -1;

    for(int i=0; i<this->tid_len; i++) this->tid[i] = 0;

    nhit.clear();
    pt.clear();
    phi.clear();
    theta.clear();
    dca.clear();
}

TemplateDatabaseWrite::TemplateDatabaseWrite(TTree *tT_meta, TTree *tT_tid, const int dataset, const int *zBins,
                                             const int *wBins,
                                             char areaDescript[3][8], const int mode, const float efficiency,
                                             const int eventcount,
                                             std::string mode_description, unsigned int template_count,
                                             const float stopping_efficiency) {

    this->tT_meta = tT_meta;
    this->tT_tid = tT_tid;
    this->dataset = dataset;
    this->tid_len = TID_LEN;
    for(int i=0; i<3; i++) {
        this->zBins[i] = zBins[i];
        this->wBins[i] = wBins[i];
        for(int j=0; j<8; j++) this->areaDescript[i][j] = areaDescript[i][j];
    }

    this->mode = mode;
    this->training_efficiency = efficiency;
    this->training_events = eventcount;
    this->stopping_efficiency = stopping_efficiency;
    this->template_count = template_count;
    this->mode_description = mode_description;

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
    this->tT_meta->Branch("mode_description", &this->mode_description);
    this->tT_meta->Branch("training_events", &this->training_events, "training_events/I");
    this->tT_meta->Branch("template_count", &this->template_count, "template_count/i");
    this->tT_meta->Branch("training_efficiency", &this->training_efficiency, "training_efficiency/F");
    this->tT_meta->Branch("stopping_efficiency", &this->stopping_efficiency, "stopping_efficiency/F");
    this->tT_meta->Branch("tid_len", &this->tid_len, "tid_len/I");
    this->tT_meta->Fill();


    this->tT_tid->Branch("tid_len", &this->tid_len, "tid_len/I");
    this->tT_tid->Branch("tid", tid, "tid[tid_len]/s");
    this->tT_tid->Branch("tid_repr", &this->tid_repr, "tid_repr/C");
    this->tT_tid->Branch("freq", &this->frequency, "frequency/I");
}


void TemplateDatabaseWrite::fillTIDData(unsigned short *tid, const int tid_len, std::string tid_repr, const int &freq, std::vector<int> &nhit,
                                        std::vector<float> &p, std::vector<float> &phi, std::vector<float> &theta,
                                        std::vector<float> dca, std::vector<unsigned int> &uEventIDs) {
    //this method also stores the traj parameters to the database file

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

TemplateDatabaseRead::TemplateDatabaseRead(TTree *tT_meta, TTree *tT_tid) {
    this->tT_tid = tT_tid;
    this->tT_meta = tT_meta;
    this->setBranches();
}

void TemplateDatabaseRead::setBranches() {

    //meta data
    tT_meta->SetBranchAddress("dataset", &dataset);
    tT_meta->SetBranchAddress("area0Description", &this->areaDescript[0]);
    tT_meta->SetBranchAddress("area1Description", &this->areaDescript[1]);
    tT_meta->SetBranchAddress("area2Description", &this->areaDescript[2]);
    tT_meta->SetBranchAddress("wBins0", &this->wBins[0]);
    tT_meta->SetBranchAddress("wBins1", &this->wBins[1]);
    tT_meta->SetBranchAddress("wBins2", &this->wBins[2]);
    tT_meta->SetBranchAddress("zBins0", &this->zBins[0]);
    tT_meta->SetBranchAddress("zBins1", &this->zBins[1]);
    tT_meta->SetBranchAddress("zBins2", &this->zBins[2]);
    tT_meta->SetBranchAddress("mode", &this->mode);
    tT_meta->SetBranchAddress("mode_description", &this->mode_description_ptr);
    tT_meta->SetBranchAddress("training_events", &this->training_events);
    tT_meta->SetBranchAddress("template_count", &this->template_count);
    tT_meta->SetBranchAddress("training_efficiency", &this->training_efficiency);
    tT_meta->SetBranchAddress("stopping_efficiency", &this->stopping_efficiency);
    tT_meta->SetBranchAddress("tid_len", &this->tid_len);
    tT_meta->GetEntry(0);

    //template data
    this->tT_tid->SetBranchAddress("tid_len", &this->tid_len);
    this->tT_tid->SetBranchAddress("tid", &this->tid);
    this->tT_tid->SetBranchAddress("tid_repr", &this->tid_repr_char);
    this->tT_tid->SetBranchAddress("freq", &this->frequency);
}

void TemplateDatabaseRead::getEntry(const int &index) {
    this->tT_tid->GetEntry(index);
    this->mode_description = *mode_description_ptr;
    this->tid_repr = std::string(tid_repr_char);
}

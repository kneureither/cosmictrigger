//
// Created by Konstantin Neureither on 11.09.20.
//

#include "BackgroundDataFile.h"

void BackgroundDataFile::reinitializeData() {
    bgtids.clear();
    tidtypes.clear();
    max_cosmic_hits = 0;
    nhit = 0;
    for(int i=0; i<TID_LEN; i++) cosmic_track[i] = 0;
}

BackgroundDataWrite::BackgroundDataWrite(TTree *tT_meta, TTree *tT_frames, const int *zBins,
                                         const int *wBins, char (*areaDescript)[8], const int mode, const int eventcount, std::string mode_description) {
    this->tT_meta = tT_meta;
    this->tT_frames = tT_frames;
    this->bgtid_len = TID_LEN;
    for(int i=0; i<3; i++) {
        this->zBins[i] = zBins[i];
        this->wBins[i] = wBins[i];
        for(int j=0; j<8; j++) this->areaDescript[i][j] = areaDescript[i][j];
    }

    this->mode = mode;
    this->training_efficiency = training_efficiency;
    this->training_events = eventcount;
    this->mode_description = mode_description;

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
    this->tT_meta->Branch("training_efficiency", &this->training_efficiency, "training_efficiency/F");
    this->tT_meta->Branch("training_events", &this->training_events, "training_events/I");
    this->tT_meta->Branch("bgtid_len", &this->bgtid_len, "bgtid_len/I");
    this->tT_meta->Branch("max_cosmic_hits", &this->max_cosmic_hits, "max_cosmic_hits/I");
    this->tT_meta->Fill();

    this->tT_frames->Branch("bgtids", &bgtids);
    this->tT_frames->Branch("bgtypes", &tidtypes);
    this->tT_frames->Branch("nhit", &nhit, "nhit/I");
    this->tT_frames->Branch("cosmic_track", &cosmic_track, "cosmic_track[bgtid_len]/s");
}

void
BackgroundDataWrite::fillBGTIDData(std::vector<temidarr> &bgtids, std::vector<int> types, temidarr cosmic_track, int nhits, int max_cosmics_hits) {
    this->reinitializeData();
    this->bgtids=bgtids;
    this->tidtypes=types;
    for(int i=0; i<this->bgtid_len; i++) this->cosmic_track[i] = cosmic_track[i];
    this->nhit = nhit;
    this->tT_frames->Fill();
}

BackgroundDataRead::BackgroundDataRead(TTree *tT_meta, TTree *tT_frames) {
    this->tT_frames = tT_frames;
    this->tT_meta = tT_meta;
    this->setBranches();
}

void BackgroundDataRead::setBranches() {
    //meta data
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
    tT_meta->SetBranchAddress("training_efficiency", &this->training_efficiency);
    tT_meta->SetBranchAddress("training_events", &this->training_events);
    tT_meta->SetBranchAddress("bgtid_len", &this->bgtid_len);
    tT_meta->SetBranchAddress("max_cosmic_hits", &this->max_cosmic_hits);
    tT_meta->GetEntry(0);

    //template data
    tT_frames->SetBranchAddress("bgtids", &this->bgtids_ptr);
    tT_frames->SetBranchAddress("bgtypes", &this->tidtype_ptr);
    tT_frames->SetBranchAddress("nhit", &this->nhit);
    tT_frames->SetBranchAddress("cosmic_track", &this->cosmic_track);
}

void BackgroundDataRead::getEntry(const int &index) {
    tT_frames->GetEntry(index);
    this->bgtids=*bgtids_ptr;
    this->tidtypes=*tidtype_ptr;
}

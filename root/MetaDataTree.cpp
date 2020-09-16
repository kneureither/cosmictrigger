//
// Created by Konstantin Neureither on 14.09.20.
//

#include "MetaDataTree.h"
#include "TFile.h"
#include "TTree.h"

void MetaDataTreeFile::reinitializeData() {

}

MetaDataTreeWrite::MetaDataTreeWrite(TTree *tT_meta, const int dataset, const int *zBins, const int *wBins,
                                     char areaDescript[3][8], const int mode, const float training_efficiency,
                                     const int bg_events, std::string mode_description, int bg_run,
                                     int max_muon_hits, int max_frame_nhits, int processed_frames,
                                     float tb_stopping_eff, unsigned int sp_count, float sp_target_ratio,
                                     int tb_training_eventcount, unsigned int template_count) {

    this->tT_meta = tT_meta;
    this->dataset = dataset;
    for(int i=0; i<3; i++) {
        this->zBins[i] = zBins[i];
        this->wBins[i] = wBins[i];
        for(int j=0; j<8; j++) this->areaDescript[i][j] = areaDescript[i][j];
    }

    this->mode = mode;
    this->efficiency = training_efficiency;
    this->training_eventcount = bg_events;
    this->tb_training_eventcount = tb_training_eventcount;
    this->mode_description = mode_description;

    this->bg_run = bg_run;
    this->max_frame_nhits = max_frame_nhits;
    this->max_muon_hits = max_muon_hits;
    this->bg_events = processed_frames;
    this->tb_stopping_eff = tb_stopping_eff;
    this->sp_count = sp_count;
    this->sp_target_ratio = sp_target_ratio;
    this->template_count = template_count;

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

    this->tT_meta->Branch("training_efficiency", &this->efficiency, "training_efficiency/F");
    this->tT_meta->Branch("stopping_efficiency", &this->tb_stopping_eff, "stopping_efficiency/F");
    this->tT_meta->Branch("training_events", &this->training_eventcount, "training_events/I");
    this->tT_meta->Branch("bg_events", &this->bg_events, "bg_events/I");
    this->tT_meta->Branch("bg_run", &this->bg_run, "bg_run/I");
//    this->tT_meta->Branch("processed_frames", &this->processed_frames, "processed_frames/I");
//    this->tT_meta->Branch("tb_stopping_eff", &this->tb_stopping_eff, "tb_stopping_eff/F");
    this->tT_meta->Branch("sp_count", &this->sp_count, "spcount/i");
    this->tT_meta->Branch("template_count", &this->template_count, "template_count/i");
    this->tT_meta->Branch("sp_target_ratio", &this->sp_target_ratio, "sp_target_ratio/F");
    this->tT_meta->Branch("max_muon_hits", &this->max_muon_hits, "max_muon_hits/I");
    this->tT_meta->Branch("max_frame_nhits", &this->max_frame_nhits, "max_frame_nhits/I");

    this->tT_meta->Fill();
}

MetaDataTreeRead::MetaDataTreeRead(TTree *tT_meta) {
    this->tT_meta = tT_meta;
    this->setBranches();
}

void MetaDataTreeRead::setBranches() {

    //meta data
//    tT_meta->SetBranchAddress("dataset", &dataset);
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
//    tT_meta->SetBranchAddress("mode_description", &this->mode_description_ptr);
//    tT_meta->SetBranchAddress("training_efficiency", &this->training_efficiency);

    tT_meta->SetBranchAddress("training_efficiency", &this->efficiency);
    tT_meta->SetBranchAddress("stopping_efficiency", &this->tb_stopping_eff);
    tT_meta->SetBranchAddress("training_events", &this->training_eventcount);
    tT_meta->SetBranchAddress("bg_events", &bg_events);
    tT_meta->SetBranchAddress("bg_run", &bg_run);
    tT_meta->SetBranchAddress("tb_training_eventcount", &tb_training_eventcount);
    tT_meta->SetBranchAddress("tb_stopping_eff", &tb_stopping_eff);
    tT_meta->SetBranchAddress("sp_count", &sp_count);
    tT_meta->SetBranchAddress("sp_target_ratio", &sp_target_ratio);
    tT_meta->SetBranchAddress("max_muon_hits", &max_muon_hits);
    tT_meta->SetBranchAddress("max_frame_nhits", &max_frame_nhits);
    tT_meta->GetEntry(0);

}

void MetaDataTreeRead::getEntry(const int &index) {
    this->tT_meta->GetEntry(index);
    this->mode_description = *mode_description_ptr;
}


BGAnaResTreeRead::BGAnaResTreeRead(TTree *tTMeta, TTree *t_bgeff) : MetaDataTreeRead(tTMeta) {
    this->tT_bgeff=t_bgeff;
    this->setBranches();
}

void BGAnaResTreeRead::setBranches() {

    tT_bgeff->SetBranchAddress("background_efficiency", &bg_discr_eff);
    tT_bgeff->SetBranchAddress("tb_training_eff", &tb_training_eff);
    tT_bgeff->SetBranchAddress("frame_eff", &frame_effs_ptr);
    tT_bgeff->SetBranchAddress("frame_bghits", &frame_hits_ptr);
    tT_bgeff->GetEntry(0);

    frame_effs = *frame_effs_ptr;
    frame_hits = *frame_hits_ptr;
}

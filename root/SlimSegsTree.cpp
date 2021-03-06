//
// Created by Konstantin Neureither on 30.06.20.
//

#include "SlimSegsTree.h"

SlimSegsTreeWrite::SlimSegsTreeWrite(TTree *slimSegs) {
    //get the tree pointer
    this->t_slimSegs = slimSegs;

    //meta data
    this->t_slimSegs->Branch("eventID", &this->eventID, "eventID/i");
    this->t_slimSegs->Branch("runID", &this->runID, "runID/i");
    this->t_slimSegs->Branch("uEventID", &this->uEventID, "uEventID/i");
    this->t_slimSegs->Branch("segsIndex", &this->segsIndex, "segsIndex/i");

    //monte carlo data
    this->t_slimSegs->Branch("mc_type", &this->mc_type, "mc_type/I");
    this->t_slimSegs->Branch("mc_tid", &this->mc_tid, "mc_tid/I");
    this->t_slimSegs->Branch("mc_pid", &this->mc_pid, "mc_pid/I");
    this->t_slimSegs->Branch("mc_p", &this->mc_p, "mc_p/F");
    this->t_slimSegs->Branch("mc_p_corr", &this->mc_p_corr, "mc_p_corr/F");
    this->t_slimSegs->Branch("mc_pt", &this->mc_pt, "mc_pt/F");
    //maybe also add mc_r

    this->t_slimSegs->Branch("nhit", &this->rec_nhit, "nhit/I");
    this->t_slimSegs->Branch("nsegs", &this->rec_ntriplet, "nsegs/I");

#if GET_ORIGINAL_HITS
    //actual hit data for each segment
    this->t_slimSegs->Branch("x00", &this->x00);
    this->t_slimSegs->Branch("x01", &this->x01);
    this->t_slimSegs->Branch("x10", &this->x10);
    this->t_slimSegs->Branch("x11", &this->x11);
    this->t_slimSegs->Branch("x20", &this->x20);
    this->t_slimSegs->Branch("x21", &this->x21);

    this->t_slimSegs->Branch("y00", &this->y00);
    this->t_slimSegs->Branch("y01", &this->y01);
    this->t_slimSegs->Branch("y10", &this->y10);
    this->t_slimSegs->Branch("y11", &this->y11);
    this->t_slimSegs->Branch("y20", &this->y20);
    this->t_slimSegs->Branch("y21", &this->y21);

    this->t_slimSegs->Branch("z00", &this->z00);
    this->t_slimSegs->Branch("z01", &this->z01);
    this->t_slimSegs->Branch("z10", &this->z10);
    this->t_slimSegs->Branch("z11", &this->z11);
    this->t_slimSegs->Branch("z20", &this->z20);
    this->t_slimSegs->Branch("z21", &this->z21);

    this->t_slimSegs->Branch("sid00", &this->sid00);
    this->t_slimSegs->Branch("sid01", &this->sid01);
    this->t_slimSegs->Branch("sid10", &this->sid10);
    this->t_slimSegs->Branch("sid11", &this->sid11);
    this->t_slimSegs->Branch("sid20", &this->sid20);
    this->t_slimSegs->Branch("sid21", &this->sid21);

    //angles from triplets
    this->t_slimSegs->Branch("tan01", &this->rec_tan01);
    this->t_slimSegs->Branch("tan12", &this->rec_tan12);
    this->t_slimSegs->Branch("lam01", &this->rec_lam01);
    this->t_slimSegs->Branch("lam12", &this->rec_lam12);
#endif

    //the combined hits from x00 x10 and x20 (y,z respectively) hits.
    this->t_slimSegs->Branch("ncombinedhits", &this->ncombinedhits, "ncombinedhits/i");
    this->t_slimSegs->Branch("xp", &this->xp);
    this->t_slimSegs->Branch("yp", &this->yp);
    this->t_slimSegs->Branch("zp", &this->zp);
    this->t_slimSegs->Branch("layerp", &this->layerp);

    //reconstruction (multi-scattering) fit data
    this->t_slimSegs->Branch("rec_p", &this->rec_p, "rec_p/F");
    this->t_slimSegs->Branch("rec_pt", &this->rec_pt, "rec_pt/F");
    this->t_slimSegs->Branch("rec_r", &this->rec_r, "rec_r/F");
    this->t_slimSegs->Branch("rec_phi", &this->rec_phi, "rec_phi/F");
    this->t_slimSegs->Branch("rec_theta", &this->rec_theta, "rec_theta/F");
    this->t_slimSegs->Branch("rec_dca_r", &this->rec_dca_r, "rec_dca_r/F");
    this->t_slimSegs->Branch("rec_dca_z", &this->rec_dca_z, "rec_dca_z/F");

    //karimaki helix fit data
    this->t_slimSegs->Branch("kari_r3d", &this->kari_r3d, "kari_r3d/F");
    this->t_slimSegs->Branch("kari_rad", &this->kari_rad, "kari_rad/F");
    this->t_slimSegs->Branch("kari_dca", &this->kari_dca, "kari_dca/F");
    this->t_slimSegs->Branch("kari_p", &this->kari_p, "kari_p/F");
    this->t_slimSegs->Branch("kari_pt", &this->kari_pt, "kari_pt/F");
    this->t_slimSegs->Branch("kari_phi", &this->kari_phi, "kari_phi/F");
    this->t_slimSegs->Branch("kari_theta", &this->kari_theta, "kari_theta/F");
    this->t_slimSegs->Branch("kari_z0", &this->kari_z0, "kari_z0/F");
    this->t_slimSegs->Branch("kari_tchi2", &this->kari_tchi2, "kari_tchi2/F");
    this->t_slimSegs->Branch("kari_lchi2", &this->kari_lchi2, "kari_lchi2/F");
}

/* combines all necessary data and fills the members of the slim segs tree representation class. */
/* After acquiring the data, it is written to the TTree */
void SlimSegsTreeWrite::fillData(const SegsTreeReadPlus &Segs,
                                 const SlimSegsMeta &Meta,
                                 const KariFitCalc &Karires,
                                 const unsigned int &ncombinedhits,
                                 const std::vector<double> &xps,
                                 const std::vector<double> &zps,
                                 const std::vector<double> &yps,
                                 const std::vector<int> &layerps) {

    //set everything to 0 and clear vectors
    this->reInitializeData();

    //meta data
    this->eventID = (unsigned int) Segs.rec_event;
    this->runID = (unsigned int) Meta.runID;
    this->uEventID = Meta.uEventID;
    this->segsIndex = (unsigned int) Meta.segsIndex;

    //traj meta data
    this->rec_nhit = Segs.rec_nhit;
    this->rec_ntriplet = Segs.rec_ntriplet;

    //monte carlo data
    this->mc_tid = Segs.mc_tid;
    this->mc_p = Segs.mc_p;
    this->mc_p_corr = Segs.mc_p_corr;
    this->mc_pt = Segs.mc_p_corr;
    this->mc_type = Segs.mc_type;
    this->mc_pid = Segs.mc_pid;

#if GET_ORIGINAL_HITS
    //convert hit arrays to vectors (standardize data types)
    for(int i = 0; i<Segs.rec_ntriplet; i++) {
        this->x00.push_back(Segs.x00[i]);
        this->x01.push_back(Segs.x01[i]);
        this->x10.push_back(Segs.x10[i]);
        this->x11.push_back(Segs.x11[i]);
        this->x20.push_back(Segs.x20[i]);
        this->x21.push_back(Segs.x21[i]);

        this->y00.push_back(Segs.y00[i]);
        this->y01.push_back(Segs.y01[i]);
        this->y10.push_back(Segs.y10[i]);
        this->y11.push_back(Segs.y11[i]);
        this->y20.push_back(Segs.y20[i]);
        this->y21.push_back(Segs.y21[i]);

        this->z00.push_back(Segs.z00[i]);
        this->z01.push_back(Segs.z01[i]);
        this->z10.push_back(Segs.z10[i]);
        this->z11.push_back(Segs.z11[i]);
        this->z20.push_back(Segs.z20[i]);
        this->z21.push_back(Segs.z21[i]);

        this->sid00.push_back(Segs.sid00[i]);
        this->sid01.push_back(Segs.sid01[i]);
        this->sid10.push_back(Segs.sid10[i]);
        this->sid11.push_back(Segs.sid11[i]);
        this->sid20.push_back(Segs.sid20[i]);
        this->sid21.push_back(Segs.sid21[i]);

        this->rec_tan01.push_back(Segs.rec_tan01[i]);
        this->rec_tan12.push_back(Segs.rec_tan12[i]);
        this->rec_lam01.push_back(Segs.rec_lam01[i]);
        this->rec_lam12.push_back(Segs.rec_lam12[i]);
    }
#endif

    //reconstruction fit data from segs
    this->rec_p = Segs.rec_p;
    this->rec_pt = Segs.rec_pt;
    //perr
    this->rec_r = Segs.rec_r;
    //rerr
    this->rec_phi = Segs.rec_phi;
    this->rec_theta = Segs.rec_theta;
    this->rec_pt = Segs.rec_pt;
    this->rec_dca_r = Segs.rec_zpca_r;
    this->rec_dca_z = Segs.rec_zpca_z;

    this->kari_rad = Karires.rad;
    this->kari_r3d = Karires.r3d;
    this->kari_dca = Karires.dca;
    this->kari_z0 = Karires.z0;
    this->kari_theta = Karires.theta;
    this->kari_phi = Karires.phi;
    this->kari_p = Karires.p;
    this->kari_pt = Karires.pt;
    this->kari_tchi2 = Karires.tchi2n;
    this->kari_lchi2 = Karires.zchi2n;

    this->xp = xps;
    this->yp = yps;
    this->zp = zps;
    this->layerp = layerps;
    this->ncombinedhits = ncombinedhits;

    this->fillTree();
}

void SlimSegsTree::reInitializeData() {
/* clears the vectors and sets everything else to 0 in order to receive another entry from the tree*/

    //meta data
    this->eventID = 0;
    this->runID = 0;
    this->uEventID = 0;
    this->segsIndex = 0;

    this->rec_nhit = 0;
    this->rec_ntriplet = 0;

    //monte carlo data
    this->mc_tid = 0;
    this->mc_p = 0;
    this->mc_p_corr = 0;
    this->mc_pt = 0;
    this->mc_type = 0;
    this->mc_pid = 0;

#if GET_ORIGINAL_HITS
    this->x00.clear();
    this->x01.clear();
    this->x10.clear();
    this->x11.clear();
    this->x20.clear();
    this->x21.clear();

    this->y00.clear();
    this->y01.clear();
    this->y10.clear();
    this->y11.clear();
    this->y20.clear();
    this->y21.clear();

    this->z00.clear();
    this->z01.clear();
    this->z10.clear();
    this->z11.clear();
    this->z20.clear();
    this->z21.clear();

    this->sid00.clear();
    this->sid01.clear();
    this->sid10.clear();
    this->sid11.clear();
    this->sid20.clear();
    this->sid21.clear();

    this->rec_tan01.clear();
    this->rec_tan12.clear();
    this->rec_lam01.clear();
    this->rec_lam12.clear();
#endif

    //ms fit data
    this->rec_p = 0;
    this->rec_perr = 0; //vlt
    this->rec_r = 0;
    this->rec_rerr = 0; //vlt
    this->rec_phi = 0;
    this->rec_theta = 0;
    this->rec_pt = 0;
    this->rec_dca_z = 0;
    this->rec_dca_r = 0;

    //kari fit data
    this->kari_rad = 0;
    this->kari_r3d = 0;
    this->kari_dca = 0;
    this->kari_z0 = 0;
    this->kari_theta = 0;
    this->kari_phi = 0;
    this->kari_p = 0;
    this->kari_pt = 0;
    this->kari_tchi2 = 0;
    this->kari_lchi2 = 0;

    //combined data
    this->xp.clear();
    this->yp.clear();
    this->zp.clear();
    this->layerp.clear();
    this->ncombinedhits = 0;
}

void SlimSegsTree::fillTree() {
    this->t_slimSegs->Fill();
}

SlimSegsTreeRead::SlimSegsTreeRead(TTree *slimSegs) {
    this->t_slimSegs = slimSegs;
    this->entries = this->t_slimSegs->GetEntries();
    this->setBranches();
}

void SlimSegsTreeRead::setBranches() {
    //meta data
    t_slimSegs->SetBranchAddress("eventID", &eventID);
    t_slimSegs->SetBranchAddress("runID", &runID);
    t_slimSegs->SetBranchAddress("uEventID", &uEventID);
    t_slimSegs->SetBranchAddress("segsIndex", &segsIndex);

    //monte carlo data
    t_slimSegs->SetBranchAddress("mc_type", &mc_type);
    t_slimSegs->SetBranchAddress("mc_tid", &mc_tid);
    t_slimSegs->SetBranchAddress("mc_pid", &mc_pid);
    t_slimSegs->SetBranchAddress("mc_p", &mc_p);
    t_slimSegs->SetBranchAddress("mc_p_corr", &mc_p_corr);
    t_slimSegs->SetBranchAddress("mc_pt", &mc_pt);

    t_slimSegs->SetBranchAddress("nhit", &rec_nhit);
    t_slimSegs->SetBranchAddress("nsegs", &rec_ntriplet);

    //the combined hits from x00, x10 and x20
    t_slimSegs->SetBranchAddress("ncombinedhits", &ncombinedhits);
    t_slimSegs->SetBranchAddress("xp", &xpp);
    t_slimSegs->SetBranchAddress("yp", &ypp);
    t_slimSegs->SetBranchAddress("zp", &zpp);
    t_slimSegs->SetBranchAddress("layerp", &layerpp);

    //reconstruction (multi-scattering) fit data
    t_slimSegs->SetBranchAddress("rec_p", &rec_p);
    t_slimSegs->SetBranchAddress("rec_pt", &rec_pt);
    t_slimSegs->SetBranchAddress("rec_r", &rec_r);
    t_slimSegs->SetBranchAddress("rec_phi", &rec_phi);
    t_slimSegs->SetBranchAddress("rec_theta", &rec_theta);
    t_slimSegs->SetBranchAddress("rec_dca_r", &rec_dca_r);
    t_slimSegs->SetBranchAddress("rec_dca_z", &rec_dca_z);

    //karimaki helix fit data
    t_slimSegs->SetBranchAddress("kari_r3d", &kari_r3d);
    t_slimSegs->SetBranchAddress("kari_rad", &kari_rad);
    t_slimSegs->SetBranchAddress("kari_dca", &kari_dca);
    t_slimSegs->SetBranchAddress("kari_p", &kari_p);
    t_slimSegs->SetBranchAddress("kari_pt", &kari_pt);
    t_slimSegs->SetBranchAddress("kari_phi", &kari_phi);
    t_slimSegs->SetBranchAddress("kari_theta", &kari_theta);
    t_slimSegs->SetBranchAddress("kari_z0", &kari_z0);
    t_slimSegs->SetBranchAddress("kari_tchi2", &kari_tchi2);
    t_slimSegs->SetBranchAddress("kari_lchi2", &kari_lchi2);

}

void SlimSegsTreeRead::getEntry(const int &index) {
//    this->t_slimSegs->Print();
    this->t_slimSegs->GetEntry(index);

    //copy these vectors
    this->xp = (*this->xpp);
    this->yp = (*this->ypp);
    this->zp = (*this->zpp);
    this->layerp = (*this->layerpp);

    assert(xp[0] == (*xpp)[0]);
    assert(yp[1] == (*ypp)[1]);
    assert(zp[2] == (*zpp)[2]);
    assert(layerp[3] == (*layerpp)[3]);
}

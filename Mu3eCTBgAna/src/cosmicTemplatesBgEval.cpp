//
// Created by Konstantin Neureither on 25.08.20.
//

#include "TFile.h"
#include "TH1F.h"

#include <cassert>
#include <map>
#include <stdlib.h>

#include "../inc/cosmicTemplatesBgEval.h"
#include "utilityFunctions.h"
#include "Mu3eTree.h"
#include "../Mu3eCosPat/include/PatternEngine.h"
#include "../Mu3eCosPat/include/TemplateBank.h"
#include "../Mu3eCosPat/include/TemplateData.h"
#include "MetaDataTree.h"


void cosmicTemplatesBgEval(const int run, int dataset, unsigned int centralTPcount, float spWZratio,
                           const float tb_stopping_efficiency, const bool append_outfile) {
    /*
     * Read and analyse the mu3e mc hits
     * get the hits in xyz
     * assign sps to these hits
     * initialize template bank
     * check the frequency
     */

    int MAX_ENTRIES = 0;
    int MAX_MUON_HITS = 0;
    int MAX_NHITS = 100;
    float TB_STOPPING_EFF = tb_stopping_efficiency;
    const bool RECREATE_FILE = !append_outfile;
    const int MUONTYPE = 1;
    const int PRINTS = false;
    const int mode = 0;

    const std::string pathtoBGdata = "data/SimulationData/";
    const std::string pathtoTemplateData = "data/TemplateData/";
    const std::string pathtoplots = "output/Mu3eCosPatBgEval/dataset_" + get_padded_string(dataset, 3, '0') + "/";
    const std::string pathtooutfile = pathtoplots + "bgrun_" + get_padded_string(run, 3, '0') + "/"; //this is where the root file is stored
    const std::string pathtorunplots = pathtooutfile +"/PDF/"; //this is where the pdf files are stored
    const std::string pathtodatasettemplatedata = pathtoTemplateData + "dataset_" + get_padded_string(dataset, 3, '0') + "/";

    check_create_directory(pathtoBGdata);
    check_create_directory(pathtoTemplateData);
    check_create_directory(pathtoplots);
    check_create_directory(pathtooutfile);
    check_create_directory(pathtorunplots);

    // FILE FOR READING BACKGROUND DATA
    const std::string infile = pathtoBGdata + "mu3e_run_" + get_padded_string(run, 6, '0') + ".root";
    TFile tinF(infile.c_str());
    if (!tinF.IsOpen()) {
        std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
        exit(0);
    }

    TTree *t_mu3e;
    tinF.GetObject("mu3e", t_mu3e);
    TTree *t_mu3e_mchits;
    tinF.GetObject("mu3e_mchits", t_mu3e_mchits);

    //Get the Pattern Engine and Template Manager
    const int spWbins = (int) sqrt((float) spWZratio * (float) centralTPcount);
    const int spZbins = (int) sqrt((float) centralTPcount / (float) spWZratio);

    Mu3eTree Mu3e = Mu3eTree(t_mu3e);
    Mu3eMChitsTree Mu3eMChits = Mu3eMChitsTree(t_mu3e_mchits);

    PatternEngine PE(spWbins, spZbins, pathtorunplots);
    TemplateBank TB(pathtorunplots, dataset, mode, spWbins, spZbins);
    TB.PRINTS = PRINTS;
    TB.readAMfromFile(pathtodatasettemplatedata, TB_STOPPING_EFF, ALL);

    //make some analysis plots
    TH1F h_bgeff("h_bgeff", "background match efficiency", 100, 0, 0.01);
    TH1F h_discreff(("h_discreff_cosmax" + get_string(MAX_MUON_HITS)).c_str(), ("h_discreff_cosmax" + get_string(MAX_MUON_HITS)).c_str(), 100, 0, 0.1);
    TH1F h_oheff(("h_oheff_cosmax" + get_string(MAX_MUON_HITS)).c_str(), "h_oheff", 100, 0, 1);
    TH1F h_heff(("h_heff_cosmax" + get_string(MAX_MUON_HITS)).c_str(), "h_heff", 100, 0, 1);


    //get some cosmic hits to play with
    srand((unsigned) time(0));
    int cosmicpool = 50;
//    std::vector<temid> COSMICTIDs = TB.getMostPopulatedTemplates(cosmicpool);


    int bg_events = (MAX_ENTRIES == 0 ? Mu3e.my_entries : MAX_ENTRIES);
    int processed_frames = 0;
    int rejected_frames = 0;

    //for some plots with bg eff / hits in bg frame
    std::vector<float> frame_eff;
    std::vector<int> frame_bghits;

    std::cout << "(CONFIG) : TDB: mywbins " << TB.mywbins << " | myzbins " << TB.myzbins << " | PE: wbins " << PE.WBins[0] << " zbins "<< PE.ZBins[0] << std::endl;

    for(int frame=0; frame <= bg_events; frame++) {
        Mu3e.getEntry(frame);
        if(PRINTS) Mu3e.Print(frame);

        if(Mu3e.Nhit > MAX_NHITS) {
            continue;
        }
        processed_frames++;




        std::vector<BGhit> bgframehits;
        BGSortedHits hits;
        std::vector<TemplateID> TIDS;
        std::vector<TemplateID> TIDSnocosmics;
        std::map<unsigned short, int> SIDMem;
        BGhit BGHIT;
        unsigned short SID;
        int innerhits=0;
        int outerhits=0;
        int doublesid=0;

        //get all hits from one frame
        assert(Mu3e.Nhit == Mu3e.hit_mc_i->size());
        for(int hitno=0; hitno < Mu3e.Nhit; hitno++) {
            Mu3eMChits.getEntry((*Mu3e.hit_mc_i)[hitno]);
            if(PRINTS) Mu3eMChits.Print((*Mu3e.hit_mc_i)[hitno]);
            BGHIT.fill(Mu3eMChits.pos_g_x, Mu3eMChits.pos_g_y, Mu3eMChits.pos_g_z, 0);
            bgframehits.push_back(BGHIT);

            SID = (unsigned short) PE.getSuperPixel(BGHIT.x, BGHIT.y, BGHIT.z);
//            if(SIDMem.count(SID) == 0) {
//                SIDMem[SID] = 1;
//            } else {
//                doublesid++;
//                continue;
//            }

            int layer=PE.getLayerFromSPID(SID);

            if(layer == 3 && BGHIT.y >= 0) {
                hits.h0.push_back(SIDtype(SID, BGHIT.type));
            } else if (layer == 2 && BGHIT.y >= 0) {
                hits.h1.push_back(SIDtype(SID, BGHIT.type));
            } else if (layer == 2 && BGHIT.y < 0) {
                hits.h2.push_back(SIDtype(SID, BGHIT.type));
            } else if (layer == 3 && BGHIT.y < 0) {
                hits.h3.push_back(SIDtype(SID, BGHIT.type));
            } else {
                innerhits++;
            }
        }

        int cosmiccount[4] = {0,0,0,0};
        int toomanycosmicscount = 0;
        outerhits = bgframehits.size() - innerhits;

//        int randindex = (rand() % cosmicpool);
//        temid COSMICTID = COSMICTIDs[randindex];

//        if(PRINTS) std::cout << "   -> got cosmic at index " << randindex << " TID=" << COSMICTID.toString() << std::endl;

//        if(SIDMem.count(COSMICTID.HIDS[0]) == 0) hits.h0.push_back(SIDtype(COSMICTID.HIDS[0], MUONTYPE));
//        if(SIDMem.count(COSMICTID.HIDS[1]) == 0) hits.h1.push_back(SIDtype(COSMICTID.HIDS[1], MUONTYPE));
//        if(SIDMem.count(COSMICTID.HIDS[2]) == 0) hits.h2.push_back(SIDtype(COSMICTID.HIDS[2], MUONTYPE));
//        if(SIDMem.count(COSMICTID.HIDS[3]) == 0) hits.h3.push_back(SIDtype(COSMICTID.HIDS[3], MUONTYPE));

//        hits.h0.push_back(SIDtype(COSMICTID.HIDS[0], MUONTYPE));
//        hits.h1.push_back(SIDtype(COSMICTID.HIDS[1], MUONTYPE));
//        hits.h2.push_back(SIDtype(COSMICTID.HIDS[2], MUONTYPE));
//        hits.h3.push_back(SIDtype(COSMICTID.HIDS[3], MUONTYPE));

        //Go through all possible combinations of hits and create corresponding Template IDs
        //for all hits in upper layer 3
        for(const auto &h0 : hits.h0) {
            cosmiccount[0] = (h0.type == MUONTYPE ? 1 : 0);

            //for all hits in upper layer 2
            for (const auto &h1 : hits.h1) {
                cosmiccount[1] = (h1.type == MUONTYPE ? 1 : 0);

                //for all hits in lower layer 2
                for (const auto &h2 : hits.h2) {
                    cosmiccount[2] = (h2.type == MUONTYPE ? 1 : 0);

                    //for all hits in lower layer 3
                    for (const auto &h3 : hits.h3) {
                        cosmiccount[3] = (h3.type == MUONTYPE ? 1 : 0);

                        //when cosmic hits are also added to data -> check if the traj does not consist of too many cosmic hits
                        if(cosmiccount[0] + cosmiccount[1] + cosmiccount[2] + cosmiccount[3] == 0) {
                            TIDSnocosmics.push_back(TemplateID(h0.SID, h1.SID, h2.SID, h3.SID));
                        }

                        if(cosmiccount[0] + cosmiccount[1] + cosmiccount[2] + cosmiccount[3] <= MAX_MUON_HITS) {
                            TIDS.push_back(TemplateID(h0.SID, h1.SID, h2.SID, h3.SID));
                        } else {
                            toomanycosmicscount++;
                        }

                    }
                }
            }
        }

//        assert(bgframehits.size() == (hits.h0.size() + hits.h1.size() + hits.h2.size() + hits.h3.size() + innerhits) - 4);

        //calculate the discrimination efficiency in a frame ---->
        int acceptedcount=0;
        int acceptedcountnocos=0;

        //check if the templates match to templates in the database
        for(auto &tid : TIDS) {
            if (TB.checkTemplate(tid)) acceptedcount++;
        }
        for(auto &tid : TIDSnocosmics) {
            if (TB.checkTemplate(tid)) acceptedcountnocos++;
        }

        float eff_templatecount = (float) acceptedcount / (float) TIDS.size();
        float eff_templatecountnc = (float) acceptedcountnocos / (float) TIDSnocosmics.size();
        float eff_outerhits = (float) acceptedcount / (float) outerhits;
        float eff_hits = (float) acceptedcount / (float) bgframehits.size();
        // <---- end of calculating the efficiencies in a frame


        if(PRINTS) {

            //show some stats
            std::cout << std::endl;
            std::cout << "----frame number " << frame << std::endl;
            std::cout << "  > total hits            --- " << bgframehits.size() << std::endl;
            std::cout << "  > inner hits            --- " << innerhits << std::endl;
            std::cout << "  > hits in upper layer 3 --- " << hits.h0.size() << std::endl;
            std::cout << "  > hits in upper layer 2 --- " << hits.h1.size() << std::endl;
            std::cout << "  > hits in lower layer 2 --- " << hits.h2.size() << std::endl;
            std::cout << "  > hits in lower layer 3 --- " << hits.h3.size() << std::endl;
            std::cout << "  > double sids           --- " << doublesid << std::endl;

            std::cout << "  > templates created: " << TIDS.size() << "    (templates excluded: " << toomanycosmicscount
                      << ")" << std::endl;
            std::cout << "  > from which were accepted: " << acceptedcount << std::endl;
            std::cout << "  > template efficiency (accepted tmpl/tested tmpl) : " << eff_templatecount << std::endl;
            std::cout << "  > template efficiency no cosmics                  : " << eff_templatecountnc << std::endl;
            std::cout << "  > outerhit efficiency (accepted tmpl/outerhits)   : " << eff_outerhits << std::endl;
            std::cout << "  > hit efficiency      (accepted tmpl/hits)        : " << eff_hits << std::endl;
            std::cout << "  > TB accepted " << TB.getAcceptedCount() << "    TB rejected " << TB.getRejectedCount();
            std::cout << "    TB ratio " << TB.getAcceptedCount() / (float) TB.getRejectedCount() << std::endl
                      << std::endl;
        }


        h_bgeff.Fill(eff_templatecount);
        h_discreff.Fill(eff_templatecount);
        h_oheff.Fill(eff_outerhits);
        h_heff.Fill(eff_hits);

        frame_eff.push_back(eff_templatecount);
        frame_bghits.push_back(bgframehits.size());

        if(acceptedcount == 0) {
            rejected_frames++;
        }

        float bg_discr_eff = (float) rejected_frames / (float) processed_frames;
//        if(frame % 1000 == 0) std::cout << "STATUS : processed background frame " << frame << " of " << bg_events << " | bg discr eff: " << bg_discr_eff*100 << " %" << std::endl;

        if(frame % 1000 == 0){
            int MAX_LEN = 50;
            float prog_perc = frame /  (float) bg_events;
            std::string prog_bar_fill((int) (MAX_LEN * prog_perc), '=');
            std::string prog_bar_empty((int) (MAX_LEN * (1-prog_perc)), ' ');
            std::cout << "\r(STATUS) : " << "[" << prog_bar_fill << ">" << prog_bar_empty << "] ";
            std::cout << frame /  (float) bg_events * 100 << "% | BG discr eff " << bg_discr_eff * 100 << "%" << std::flush;
        }

        //calculate rations:
            /*
             * - accepted templates / tested templates
             * - accepted templates / outer hits
             *
             * Also do this for a cosmic in the frame (with adding 1,2,3 and 4 hits of it to the combinatorics
             */
    }
    std::cout << std::endl;

    float background_efficiency = (float) rejected_frames / (float) processed_frames;

    std::cout << "INFO    : rejected frames: " << rejected_frames << " processed frames: " << processed_frames << std::endl;
    std::cout << "INFO    : --  BG efficiency: " << background_efficiency << std::endl;

    //open new TFile for plots
    TFile * tF = new TFile((pathtooutfile + "CosmicBackgroundEval_bgevents_" + get_padded_string(bg_events, 6, '0') +
                            "_run_" + get_padded_string(run, 6, '0') + "_" +
            TB.getfileidtag(1) + "_plots.root").c_str(), (RECREATE_FILE ? "recreate" : "update"));
    if (!tF->IsOpen()) {
        std::cout << "[ERROR] File " << tF->GetName() << " is not open!" << std::endl;
    }

    std::string mydirectory = "trainingEff" + get_string(TB_STOPPING_EFF);
    tF->mkdir(mydirectory.c_str());
    tF->cd(mydirectory.c_str());

    //add some meta data for the output file
    int bg_run = run;
    float tb_max_efficiency = TB.getEfficiency();

    std::cout << "tb event count " << TB.getTrainingEventCount() << " templ count " << TB.getTemplateCount() << std::endl;

    TTree tT_met("MetadataTree","Metadata associated with these plots (SID config and dataset)");
    MetaDataTreeWrite Meta = MetaDataTreeWrite(&tT_met, dataset, PE.ZBins, PE.WBins, PE.areaDescript,
                                               mode, tb_max_efficiency, bg_events, "default", bg_run, MAX_MUON_HITS,
                                               MAX_NHITS, processed_frames, tb_stopping_efficiency,
                                               (unsigned int) centralTPcount, spWZratio, TB.getTrainingEventCount(), TB.getTemplateCount());
    tT_met.Write();


//    tT_met.Branch("bg_run", &bg_run, "bg_run/I");
//    tT_met.Branch("bg_events", &bg_events, "bg_events/I");
//    tT_met.Branch("max_muon_hits", &MAX_MUON_HITS, "max_muon_hits/I");
//    tT_met.Branch("max_frame_nhits", &MAX_NHITS, "max_frame_nhits/I");
//    tT_met.Branch("processed_frames", &processed_frames, "processed_frames/I");
//    tT_met.Branch("tb_stopping_eff", &TB_STOPPING_EFF, "tb_stopping_eff/F");
//    tT_met.Branch("sp_count", &centralTPcount, "spcount/i");
//    tT_met.Branch("sp_target_ratio", &spWZratio, "sp_target_ratio/i");
//
//
//    tT_met.Branch("area0Description", &PE.areaDescript[0], "area0Description/C");
//    tT_met.Branch("area1Description", &PE.areaDescript[1], "area1Description/C");
//    tT_met.Branch("area2Description", &PE.areaDescript[2], "area2Description/C");
//    tT_met.Branch("wBins0", &PE.WBins[0], "wBins0/I");
//    tT_met.Branch("wBins1", &PE.WBins[1], "wBins1/I");
//    tT_met.Branch("wBins2", &PE.WBins[2], "wBins2/I");
//    tT_met.Branch("zBins0", &PE.ZBins[0], "zBins0/I");
//    tT_met.Branch("zBins1", &PE.ZBins[1], "zBins1/I");
//    tT_met.Branch("zBins2", &PE.ZBins[2], "zBins2/I");
//    tT_met.Branch("mode", &PE.mode, "mode/I");
//    tT_met.Fill();
//    tT_met.Write();


    TTree tT_efficiencies("BackgroundEfficiency", "Contains the Bg efficiencies per Frame");
    tT_efficiencies.Branch("background_efficiency", &background_efficiency, "background_efficiency/F");
    tT_efficiencies.Branch("tb_training_eff", &tb_max_efficiency, "tb_training_eff/F");
    tT_efficiencies.Branch("frame_eff", &frame_eff);
    tT_efficiencies.Branch("frame_bghits", &frame_bghits);
    tT_efficiencies.Fill();
    tT_efficiencies.Write();

    PE.displayBinBoundaries();
    PE.displayBinWeightDistribution();
    PE.closePlot();

//    TB.getMostPopulatedTemplates(50);
//    TB.getMostMatchedTemplates(50);
    TB.PlotTemplatePopulationHistogram();
    TB.PlotTemplateMatchedFreqHistogram(TB.getfileidtag(0));
//    TB.PlotTemplatePopHistSortedbyFreq(TB.getfileidtag());

    // Norm
    h_bgeff.Scale(1.0 / h_bgeff.Integral());
    h_discreff.Scale(1.0 / h_discreff.Integral());
    h_oheff.Scale(1.0 / h_oheff.Integral());
    h_heff.Scale(1.0 / h_heff.Integral());

    //Write
    h_bgeff.Write();
    h_discreff.Write();
    h_oheff.Write();
    h_heff.Write();
    tF->Close();
}

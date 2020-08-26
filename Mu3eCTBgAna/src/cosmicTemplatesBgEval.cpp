//
// Created by Konstantin Neureither on 25.08.20.
//

#include "TFile.h"
#include "TH1F.h"

#include <cassert>
#include <stdlib.h>

#include "../inc/cosmicTemplatesBgEval.h"
#include "utilityFunctions.h"
#include "Mu3eTree.h"
#include "../Mu3eCosPat/include/PatternEngine.h"
#include "../Mu3eCosPat/include/TemplateBank.h"
#include "../Mu3eCosPat/include/TemplateData.h"


void cosmicTemplatesBgEval(const int run, unsigned int centralTPcount, float spWZratio) {
    /*
     * Read and analyse the mu3e mc hits
     * get the hits in xyz
     * assign sps to these hits
     * initialize template bank
     * check the frequency
     */

    const int MAX_ENTRIES = 1000;
    const int PRINTS = false;
    const int dataset = 6; //determines which pretrained database will be used
    const int mode = 0;
    const int MUONTYPE = 1;
    const int MAX_MUON_HITS = 0;

    const std::string pathtoBGdata = "data/SimulationData/";
    const std::string pathtoTemplateData = "data/TemplateData/";
    const std::string pathtoplots = "plots/Mu3eCosPatBgEval/";
    const std::string infile = pathtoBGdata + "mu3e_run_" + get_padded_string(run, 6, '0') + ".root";
    std::string pathtorunplots = pathtoplots + "run_" + get_padded_string(run, 3, '0') + "/";
    std::string pathtodatasettemplatedata = pathtoTemplateData + "dataset_" + get_padded_string(dataset, 3, '0') + "/";

    check_create_directory(pathtoBGdata);
    check_create_directory(pathtoTemplateData);
    check_create_directory(pathtoplots);
    check_create_directory(pathtorunplots);

    // FILE FOR READING
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
    std::cout << "\n -- PE config data:" << std::endl << "  wbins=" << spWbins << std::endl << "  zbins=" << spZbins << std::endl << std::endl;

    Mu3eTree Mu3e = Mu3eTree(t_mu3e);
    Mu3eMChitsTree Mu3eMChits = Mu3eMChitsTree(t_mu3e_mchits);

    PatternEngine PE(spWbins, spZbins, pathtorunplots);
    TemplateBank TB(pathtoplots);
    TB.PRINTS = PRINTS;
    TB.readAMfromFile(pathtodatasettemplatedata, spWbins, spZbins, mode, dataset);



    //make some analysis plots
    TH1F h_discreff(("h_discreff_cosmax" + get_string(MAX_MUON_HITS)).c_str(), ("h_discreff_cosmax" + get_string(MAX_MUON_HITS)).c_str(), 100, 0, 0.01);
    TH1F h_oheff(("h_oheff_cosmax" + get_string(MAX_MUON_HITS)).c_str(), "h_oheff", 100, 0, 1);
    TH1F h_heff(("h_heff_cosmax" + get_string(MAX_MUON_HITS)).c_str(), "h_heff", 100, 0, 1);


    //get some cosmic hits to play with
    srand((unsigned) time(0));
    int cosmicpool = 50;
    std::vector<temid> COSMICTIDs = TB.getMostPopulatedTemplates(cosmicpool);



    for(int frame=0; frame <= (MAX_ENTRIES == 0 ? Mu3e.my_entries : MAX_ENTRIES); frame++) {
        Mu3e.getEntry(frame);
        if(PRINTS) Mu3e.Print(frame);

        if(Mu3e.Nhit > 50) {
            continue;
        }

        std::vector<BGhit> bgframehits;
        BGSortedHits hits;
        std::vector<TemplateID> TIDS;
        std::vector<TemplateID> TIDSnocosmics;
        BGhit BGHIT;
        unsigned short SID;
        int innerhits=0;
        int outerhits=0;

        //get all hits from one frame
        assert(Mu3e.Nhit == Mu3e.hit_mc_i->size());
        for(int hitno=0; hitno < Mu3e.Nhit; hitno++) {
            Mu3eMChits.getEntry((*Mu3e.hit_mc_i)[hitno]);
            if(PRINTS) Mu3eMChits.Print((*Mu3e.hit_mc_i)[hitno]);
            BGHIT.fill(Mu3eMChits.pos_g_x, Mu3eMChits.pos_g_y, Mu3eMChits.pos_g_z, 0);
            bgframehits.push_back(BGHIT);

            SID = (unsigned short) PE.getSuperPixel(BGHIT.x, BGHIT.y, BGHIT.z);
            int layer=PE.getLayerFromSPID(SID);

            //FIXME: Now double SIDs are not ignored. Should be checked before pushing back into hits vector
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

        int randindex = (rand() % cosmicpool);
        temid COSMICTID = COSMICTIDs[randindex];

        std::cout << "   -> got cosmic at index " << randindex << " TID=" << COSMICTID.toString() << std::endl;

        hits.h0.push_back(SIDtype(COSMICTID.HIDS[0], MUONTYPE));
        hits.h1.push_back(SIDtype(COSMICTID.HIDS[1], MUONTYPE));
        hits.h2.push_back(SIDtype(COSMICTID.HIDS[2], MUONTYPE));
        hits.h3.push_back(SIDtype(COSMICTID.HIDS[3], MUONTYPE));


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

        assert(bgframehits.size() == (hits.h0.size() + hits.h1.size() + hits.h2.size() + hits.h3.size() + innerhits) - 4);

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

        //show some stats
        std::cout << std::endl;
        std::cout << "----frame number " << frame << std::endl;
        std::cout << "  > total hits            --- " << bgframehits.size() << std::endl;
        std::cout << "  > inner hits            --- " << innerhits << std::endl;
        std::cout << "  > hits in upper layer 3 --- " << hits.h0.size() << std::endl;
        std::cout << "  > hits in upper layer 2 --- " << hits.h1.size() << std::endl;
        std::cout << "  > hits in lower layer 2 --- " << hits.h2.size() << std::endl;
        std::cout << "  > hits in lower layer 3 --- " << hits.h3.size() << std::endl;

        std::cout << "  > templates created: " << TIDS.size() << "    (templates excluded: " << toomanycosmicscount << ")" << std::endl;
        std::cout << "  > from which were accepted: " << acceptedcount << std::endl;
        std::cout << "  > template efficiency (accepted tmpl/tested tmpl) : " << eff_templatecount << std::endl;
        std::cout << "  > template efficiency no cosmics                  : " << eff_templatecountnc << std::endl;
        std::cout << "  > outerhit efficiency (accepted tmpl/outerhits)   : " << eff_outerhits << std::endl;
        std::cout << "  > hit efficiency      (accepted tmpl/hits)        : " << eff_hits << std::endl;
        std::cout << "  > TB accepted " << TB.getAcceptedCount() << "    TB rejected " << TB.getRejectedCount();
        std::cout << "    TB ratio " << TB.getAcceptedCount() / (float) TB.getRejectedCount() << std::endl << std::endl;

        h_discreff.Fill(eff_templatecount);
        h_oheff.Fill(eff_outerhits);
        h_heff.Fill(eff_hits);


        //calculate rations:
            /*
             * - accepted templates / tested templates
             * - accepted templates / outer hits
             *
             * Also do this for a cosmic in the frame (with adding 1,2,3 and 4 hits of it to the combinatorics
             */
    }

    //open new TFile for plots
    TFile * tF = new TFile((pathtorunplots +"CosmicBackgroundEval_dataset_" + get_padded_string(dataset, 3, '0') +
                            "_run_" + get_padded_string(run, 6, '0') + "_" +
                            TB.getcustomnamestring() + "_plots.root").c_str(), "update");
    if (!tF->IsOpen()) {
        std::cout << "[ERROR] File " << tF->GetName() << " is not open!" << std::endl;
    }

    PE.displayBinBoundaries();
    PE.displayBinWeightDistribution();
    PE.closePlot();

    // Norm
    h_discreff.Scale(1.0 / h_discreff.Integral());
    h_oheff.Scale(1.0 / h_oheff.Integral());
    h_heff.Scale(1.0 / h_heff.Integral());

    //Write
    h_discreff.Write();
    h_oheff.Write();
    h_heff.Write();
    tF->Close();
}

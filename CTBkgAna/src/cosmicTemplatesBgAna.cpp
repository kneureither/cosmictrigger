//
// Created by Konstantin Neureither on 12.10.20.
//

#include "TFile.h"
#include "TH1F.h"

#include <cassert>
#include <map>
#include <stdlib.h>

#include "cosmicTemplatesBgAna.h"
#include "utilityFunctions.h"
#include "Mu3eTree.h"
#include "../../CTCoreModules/inc/PatternEngine.h"
#include "../../CTCoreModules/inc/TemplateBank.h"
#include "../../CTCoreModules/inc/TemplateData.h"
#include "../../CTCosPatTrain/inc/getCosmicSIDtracks.h"
#include "MetaDataTree.h"
#include "../../CTCoreModules/Configuration.h"

/**
 * This file contains the background analysis.
 * To save computation time, it get the background hits from the file once per sp config
 * and stores it inside a struct BGcombinatorics.
 * For different filter settings, this struct is reused, which is a lot faster, as the data
 * is in the RAM.
 *
 * However, depending on the beam rate of the simulation, this data grows very fast, so that
 * the RAM size can get a critical parameter.
 *
 * For higher beam rates, this must be changed, so that the data ist always read from the root file and
 * never completely stored inside the RAM.
 */

struct BGcombinatorics {
    std::vector<std::vector<TemplateID>> frames_TIDS;
    std::vector<int> frames_outerhits;
    std::vector<int> frames_innerhits;
    std::vector<int> frames_bgframehits;
    std::vector<BGSortedSIDs> frames_sortedhits;
    int total_nhits = 0;
    int processed_frames = 0;
    int nhit_excluded_frames = 0;
    int bg_events = 0;
};

BGcombinatorics produceBGcmbncsTID(const int bg_run, PatternEngine &PEbg, int max_frame_hits, const int MAX_BG_ENTRIES, std::string pathtoBGdata, std::string pathtoplots);


void cosmicTemplatesBgAna(const int run, int dataset, unsigned int centralTPcount, float spWZratio,
                          const float tb_stopping_efficiency, const bool append_outfile,
                          std::vector<TIDLoadingFilter> filters, int max_bg_entries, int max_cosmic_entries,
                          int max_bg_frame_nhits, int cosmic_testing_dst) {
    /**
     * Read and analyse the mu3e mc hits of a background run
     * get the hits in xyz
     * assign sids to these hits
     * make full combinatorics of these hits into templates
     * initialize a pretrained template bank
     * check the frequency of false-positive matches
     */

    int MAX_BG_ENTRIES = max_bg_entries;
    int MAX_COSMIC_ENTRIES = max_cosmic_entries;
    int MAX_MUON_HITS = 0;
    int MAX_NHITS = max_bg_frame_nhits;
    float TB_STOPPING_EFF = tb_stopping_efficiency;
    const bool RECREATE_FILE = !append_outfile;
    const int PRINTS = false;
    const int mode = 0;
    const int cosmic_testing_dataset = cosmic_testing_dst;

    Configuration CONFIG;

    const std::string pathtoBGdata = CONFIG.pathtosimfiles;
    const std::string pathtoTemplateData = "data/TemplateData/";
    const std::string pathtobkgeval = "output/3_BKGEvaluation/";
    const std::string pathtoplots = pathtobkgeval + "dataset_" + get_padded_string(dataset, 3, '0') + "/";
    const std::string pathtooutfile =
            pathtoplots + "bgrun_" + get_padded_string(run, 3, '0') + "/"; //this is where the root file is stored
    const std::string pathtorunplots = pathtooutfile + "/PDF/"; //this is where the pdf files are stored
    const std::string pathtodatasettemplatedata =
            pathtoTemplateData + "dataset_" + get_padded_string(dataset, 3, '0') + "/";

    check_create_directory(pathtoBGdata);
    check_create_directory(pathtoTemplateData);
    check_create_directory(pathtobkgeval);
    check_create_directory(pathtoplots);
    check_create_directory(pathtooutfile);
    check_create_directory(pathtorunplots);

    //make some analysis plots
    TH1F h_bgeff("h_bgeff", "background match efficiency", 100, 0, 0.01);
    TH1F h_oheff(("h_oheff_cosmax" + get_string(MAX_MUON_HITS)).c_str(), "h_oheff", 100, 0, 1);
    TH1F h_heff(("h_heff_cosmax" + get_string(MAX_MUON_HITS)).c_str(), "h_heff", 100, 0, 1);

    //Get the Pattern Engine for background data eval.
    const int spWbins = (int) sqrt((float) spWZratio * (float) centralTPcount);
    const int spZbins = (int) sqrt((float) centralTPcount / (float) spWZratio);
    //for acuiring the background data
    PatternEngine PEbg(spWbins, spZbins, pathtorunplots);


    //store all background TIDs as a vector per frame in these vectors, so that the combinatoric
    //must not be repeated when calling another filter setting.
    BGcombinatorics BGcmbncsResult = produceBGcmbncsTID(run, PEbg, MAX_NHITS, MAX_BG_ENTRIES, pathtoBGdata, pathtoplots);

    for (auto &filter : filters) {
        std::cout << "\n(STATUS) : running analysis for filter " << enum_to_string(filter) << std::endl;

        //// Calculate the cosmic efficiency

        //check the templates collected from file above
        std::vector<std::vector<unsigned int>> cosmic_spid_tracks;
        cosmic_spid_tracks = getCosmicSIDtracks(cosmic_testing_dataset, centralTPcount, spWZratio, mode);
        TemplateBank TBC = TemplateBank(pathtoplots, dataset, mode, spWbins, spZbins);
        TBC.readAMfromFile(pathtodatasettemplatedata,tb_stopping_efficiency, filter);
        for(auto &SPIDs : cosmic_spid_tracks) {
            TBC.checkCosmicTemplate(&SPIDs[0], SPIDs.size(), filter);
        }

        //store data
        float train_eff_total = TBC.GetTrainEffTotal();
        float train_eff_relative = TBC.GetTrainEffRelative();
        std::cout << "\n(INFO)   : [COS-EFF] total eff " << TBC.GetTrainEffTotal() << " | relative eff " << TBC.GetTrainEffRelative() << std::endl;


        //// Calculate the bg discrimination

        TemplateBank TB(pathtorunplots, dataset, mode, spWbins, spZbins);
        TB.SetPrints(false);
        TB.readAMfromFile(pathtodatasettemplatedata, TB_STOPPING_EFF, filter);
        TB.PlotTemplateTypeDistribution();
        std::cout << "(CONFIG) : [BKG-EFF] wBins " << TB.mywbins << " | zBins " << TB.myzbins << " | max nhits " << MAX_NHITS << std::endl;

        int rejected_frames = 0;
        int accepted_frames = 0;
        //for some plots with bg eff / hits in bg frame
        std::vector<float> frame_eff;
        std::vector<int> frame_bghits;

        for (int frame = 0; frame < BGcmbncsResult.frames_TIDS.size(); frame++) {

            //calculate the discrimination efficiency in a frame ---->
            int acceptedcount = 0;
            int tid_count = BGcmbncsResult.frames_TIDS[frame].size();
            int bg_hits = BGcmbncsResult.frames_bgframehits[frame];

            //check if the templates match to templates in the database
            for (auto &tid : BGcmbncsResult.frames_TIDS[frame]) {
                if (TB.checkTemplate(tid)) acceptedcount++;
            }

            float eff_templatecount = (float) acceptedcount / (float) tid_count;
            float eff_outerhits = (float) acceptedcount / (float) BGcmbncsResult.frames_outerhits[frame];
            float eff_hits = (float) acceptedcount / (float) bg_hits;
            // <---- end of calculating the efficiencies in a frame


            if (PRINTS) {

                //show some stats
                std::cout << std::endl;
                std::cout << "----frame number " << frame << std::endl;
                std::cout << "  > total hits            --- " << bg_hits << std::endl;
                std::cout << "  > inner hits            --- " << BGcmbncsResult.frames_innerhits[frame] << std::endl;
//                std::cout << "  > hits in upper layer 3 --- " << BGcmbncsResult.frames_sortedhits[frame].h0.size() << std::endl;
//                std::cout << "  > hits in upper layer 2 --- " << BGcmbncsResult.frames_sortedhits[frame].h1.size() << std::endl;
//                std::cout << "  > hits in lower layer 2 --- " << BGcmbncsResult.frames_sortedhits[frame].h2.size() << std::endl;
//                std::cout << "  > hits in lower layer 3 --- " << BGcmbncsResult.frames_sortedhits[frame].h3.size() << std::endl;
//                std::cout << "  > double sids           --- " << doublesid << std::endl;

                std::cout << "  > templates created: " << tid_count << std::endl;
                std::cout << "  > from which were accepted: " << accepted_frames << std::endl;
                std::cout << "  > template efficiency (accepted tmpl/tested tmpl) : " << eff_templatecount << std::endl;
                std::cout << "  > outerhit efficiency (accepted tmpl/outerhits)   : " << eff_outerhits << std::endl;
                std::cout << "  > hit efficiency      (accepted tmpl/hits)        : " << eff_hits << std::endl;
                std::cout << "  > TB accepted " << TB.getAcceptedCount() << "    TB rejected " << TB.getRejectedCount();
                std::cout << "    TB ratio " << TB.getAcceptedCount() / (float) TB.getRejectedCount() << std::endl
                          << std::endl;
            }

            h_bgeff.Fill(eff_templatecount);
            h_oheff.Fill(eff_outerhits);
            h_heff.Fill(eff_hits);

            frame_eff.push_back(eff_templatecount);
            frame_bghits.push_back(bg_hits);

            if (acceptedcount == 0) {
                rejected_frames++;
            } else {
                accepted_frames++;
            }

            float bg_discr_eff = (float) rejected_frames / (float) (accepted_frames + rejected_frames);
            print_status_bar(frame, BGcmbncsResult.frames_TIDS.size(), "checking frames ", "BG discr eff " + get_string(bg_discr_eff * 100) + "%");

        }


        std::cout << std::endl;

        float background_efficiency = (float) rejected_frames /  (float) (accepted_frames + rejected_frames);
        float mean_frame_nhits = BGcmbncsResult.total_nhits / (float) (accepted_frames + rejected_frames);

        std::cout << "(INFO)    : [BKG-EFF] rejected frames: " << rejected_frames
                  << " accepted frames: " << accepted_frames
                  << " processed frames: " << BGcmbncsResult.processed_frames
                  << " cut off frames (nhit): " << BGcmbncsResult.nhit_excluded_frames
                  << std::endl;
        std::cout << "(INFO)    : [BKG-EFF] BG efficiency: " << background_efficiency << "   (avg. nhits: " << mean_frame_nhits << ")"
                  << std::endl;

        //open new TFile for plots
        TFile *tF = new TFile((pathtooutfile + get_bgevalfile(BGcmbncsResult.bg_events, MAX_COSMIC_ENTRIES,
                MAX_NHITS, run, cosmic_testing_dataset, dataset, mode, spWbins, spZbins)).c_str(),
                        (RECREATE_FILE ? "recreate" : "update"));

//        TFile *tF = new TFile((pathtooutfile + "CosmicTBGAna" +
//                               "_bkgEv" + get_padded_string(BGcmbncsResult.bg_events, 6, '0') +
//                               "_cosEv" + get_padded_string(MAX_COSMIC_ENTRIES, 6, '0') +
//                               "_bkgrun_" + get_padded_string(run, 6, '0') +
//                               "_cosdst_" + get_padded_string(cosmic_testing_dataset, 6, '0') + "_" +
//                               TB.getfileidtag(1) + "_plots.root").c_str(), (RECREATE_FILE ? "recreate" : "update"));

        if (!tF->IsOpen()) {
            std::cout << "[ERROR] File " << tF->GetName() << " is not open!" << std::endl;
        }

        //make directories in file for different train effs
        std::string dirlevelone = enum_to_string(filter);
        tF->mkdir(dirlevelone.c_str());
        tF->cd(dirlevelone.c_str());
        std::string dirleveltwo = "trainingEff" + get_string(TB_STOPPING_EFF);
        tF->mkdir((dirlevelone + "/" + dirleveltwo).c_str());
        tF->cd((dirlevelone + "/" + dirleveltwo).c_str());

        //add some meta data for the output file
        int bg_run = run;
        float tb_max_efficiency = TB.getEfficiency();

        std::cout << "tb event count " << TB.getTrainingEventCount() << " templ count " << TB.getTemplateCount()
                  << std::endl;

        TTree tT_met("MetadataTree", "Metadata associated with these plots (SID config and dataset)");

        MetaDataTreeWrite Meta = MetaDataTreeWrite(&tT_met, dataset, PEbg.ZBins, PEbg.WBins, PEbg.areaDescript,
                                                   mode, tb_max_efficiency, BGcmbncsResult.bg_events, "default", bg_run, MAX_MUON_HITS,
                                                   MAX_NHITS, BGcmbncsResult.processed_frames, tb_stopping_efficiency,
                                                   (unsigned int) centralTPcount, spWZratio, TB.getTrainingEventCount(), TB.getTemplateCount());
        tT_met.Write();

        std::string filter_string = enum_to_string(filter);


        TTree tT_efficiencies("BackgroundEfficiency", "Contains the Bg efficiencies per Frame");
        tT_efficiencies.Branch("background_efficiency", &background_efficiency, "background_efficiency/F");
        tT_efficiencies.Branch("tb_training_eff", &tb_max_efficiency, "tb_training_eff/F");
        tT_efficiencies.Branch("tb_train_eff_total", &train_eff_total, "tb_train_eff_total/F");
        tT_efficiencies.Branch("tb_train_eff_relative", &train_eff_relative, "tb_train_eff_relative/F");
        tT_efficiencies.Branch("tb_filter_str", &filter_string);
        tT_efficiencies.Branch("tb_filter", &filter, "tb_filter/I");

        tT_efficiencies.Branch("frame_eff", &frame_eff);
        tT_efficiencies.Branch("frame_bghits", &frame_bghits);
        tT_efficiencies.Branch("mean_frame_nhits", &mean_frame_nhits, "mean_frame_nhits/F");
        tT_efficiencies.Branch("max_bg_frame_nhits", &max_bg_frame_nhits, "max_bg_frame_nhits/I");
        tT_efficiencies.Branch("nhit_cut_off_frames", &BGcmbncsResult.nhit_excluded_frames, "nhit_cut_off_frames/I");
        tT_efficiencies.Branch("processed_frames", &BGcmbncsResult.processed_frames, "processed_frames/I");
        tT_efficiencies.Fill();
        tT_efficiencies.Write();

        TB.PlotTemplatePopulationHistogram();
        TB.PlotTemplateMatchedFreqHistogram(TB.getfileidtag(0));
    //    TB.PlotTemplatePopHistSortedbyFreq(TB.getfileidtag());

        // Norm
        h_bgeff.Scale(1.0 / h_bgeff.Integral());
        h_oheff.Scale(1.0 / h_oheff.Integral());
        h_heff.Scale(1.0 / h_heff.Integral());

        //Write
        h_bgeff.Write();
        h_oheff.Write();
        h_heff.Write();
        tF->Close();
    }
}


BGcombinatorics produceBGcmbncsTID(const int bg_run, PatternEngine &PEbg, int max_frame_hits, const int MAX_BG_ENTRIES, std::string pathtoBGdata, std::string pathtoplots) {
//    std::string pathtorunplots = pathtoplots + "BGcmbncs/";
//    check_create_directory(pathtorunplots);
    const int MAX_NHITS = max_frame_hits;
    bool PRINTS = false;

    // FILE FOR READING BACKGROUND DATA
    const std::string infile = pathtoBGdata + "mu3e_run_" + get_padded_string(bg_run, 6, '0') + ".root";
    TFile tinF(infile.c_str());
    if (!tinF.IsOpen()) {
        std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
        exit(0);
    }

    TTree *t_mu3e;
    tinF.GetObject("mu3e", t_mu3e);
    TTree *t_mu3e_mchits;
    tinF.GetObject("mu3e_mchits", t_mu3e_mchits);

    Mu3eTree Mu3e = Mu3eTree(t_mu3e);
    Mu3eMChitsTree Mu3eMChits = Mu3eMChitsTree(t_mu3e_mchits);

    int bg_events = (MAX_BG_ENTRIES == 0 ? Mu3e.my_entries : MAX_BG_ENTRIES);

    //store all background TIDs as a vector per frame in these vectors, so that the combinatoric
    //must not be repeated when calling another filter setting.
    BGcombinatorics BGcmbncsResult;


    //bg analysis.
    for (int frame = 0; frame <= bg_events; frame++) {
        Mu3e.getEntry(frame);
        if (PRINTS) Mu3e.Print(frame);

        if (MAX_NHITS != 0 && Mu3e.Nhit > MAX_NHITS) {
            BGcmbncsResult.nhit_excluded_frames++;
            continue;
        } else if(Mu3e.Nhit == 0) {
            continue;
        }

        BGcmbncsResult.processed_frames++;
        print_status_bar(frame, bg_events, "bg event combinatorics", "");

        std::vector<BGhit> bgframehits;
        BGSortedSIDs hits;
        std::vector<TemplateID> TIDS;
        std::map<unsigned short, int> SIDMem;
        BGhit BGHIT;
        unsigned short SID;
        int innerhits = 0;
        int outerhits = 0;
//            int doublesid = 0;

        //get all hits from one frame
        assert(Mu3e.Nhit == Mu3e.hit_mc_i->size());
        for (int hitno = 0; hitno < Mu3e.Nhit; hitno++) {
            Mu3eMChits.getEntry((*Mu3e.hit_mc_i)[hitno]);
            if (PRINTS) Mu3eMChits.Print((*Mu3e.hit_mc_i)[hitno]);
            BGHIT.fill(Mu3eMChits.pos_g_x, Mu3eMChits.pos_g_y, Mu3eMChits.pos_g_z, 0);
            bgframehits.push_back(BGHIT);

            SID = (unsigned short) PEbg.getSuperPixel(BGHIT.x, BGHIT.y, BGHIT.z);
            //            if(SIDMem.count(SID) == 0) {
            //                SIDMem[SID] = 1;
            //            } else {
            //                doublesid++;
            //                continue;
            //            }

            int layer = PEbg.getLayerFromSPID(SID);

            //sort the hits into their layers
            if (layer == 3 && BGHIT.y >= 0) {
                hits.h0.push_back(SIDtype(SID, BGHIT.type));
                outerhits++;
            } else if (layer == 2 && BGHIT.y >= 0) {
                hits.h1.push_back(SIDtype(SID, BGHIT.type));
                outerhits++;
            } else if (layer == 2 && BGHIT.y < 0) {
                hits.h2.push_back(SIDtype(SID, BGHIT.type));
                outerhits++;
            } else if (layer == 3 && BGHIT.y < 0) {
                hits.h3.push_back(SIDtype(SID, BGHIT.type));
                outerhits++;
            } else {
                innerhits++;
            }
        }


        //Go through all possible combinations of hits and create corresponding Template IDs
        //for all hits in upper layer 3
        for (const auto &h0 : hits.h0) {
            //for all hits in upper layer 2
            for (const auto &h1 : hits.h1) {
                //for all hits in lower layer 2
                for (const auto &h2 : hits.h2) {
                    //for all hits in lower layer 3
                    for (const auto &h3 : hits.h3) {
                        TIDS.push_back(TemplateID(h0.SID, h1.SID, h2.SID, h3.SID));
                    }
                }
            }
        }

        assert(outerhits == bgframehits.size() - innerhits);

        //save combinatorics result per frame
        BGcmbncsResult.frames_TIDS.push_back(TIDS); //calculated tids
        BGcmbncsResult.frames_outerhits.push_back(outerhits); //int number of hits in outer layers per frame
        BGcmbncsResult.frames_innerhits.push_back(innerhits); //int number of hits in inner layers per frame
        BGcmbncsResult.frames_bgframehits.push_back(bgframehits.size()); // int number of hits per frame
//        BGcmbncsResult.frames_sortedhits.push_back(hits); //SIDs sorted in layers per frame
        BGcmbncsResult.total_nhits += bgframehits.size();
    }

    BGcmbncsResult.bg_events = bg_events;

    return BGcmbncsResult;
}


//void produceBackgroundSIDs(const int bg_run, unsigned int centralTPcount, float spWZratio, int mode, int max_frame_hits) {
//    /**
//     * This function creates TIDs of background frames
//     */
//
//    const bool PRINTS = false;
//
//    const std::string path_to_inputdata = "data/SimulationData/";
//    const std::string path_to_bkg_tids = "data/BackgroundSIDData/";
//    const std::string path_to_plots = "output/produceBackgroundSIDs/";
//    const std::string path_to_plots_dataset = path_to_plots + "bgrun_" + get_padded_string(bg_run, 3, '0') + "/";
//    const std::string path_to_bkg_tids_dataset = path_to_bkg_tids + "bgrun_" + get_padded_string(bg_run, 3, '0') + "/";
//
//    check_create_directory(path_to_inputdata);
//    check_create_directory(path_to_bkg_tids);
//    check_create_directory(path_to_plots);
//    check_create_directory(path_to_plots_dataset);
//    check_create_directory(path_to_bkg_tids_dataset);
//
//    //Get the Pattern Engine and Template Manager
//    const int spWbins = (int) sqrt((float) spWZratio * (float) centralTPcount);
//    const int spZbins = (int) sqrt((float) centralTPcount / (float) spWZratio);
//    PatternEngine PE(spWbins, spZbins, path_to_plots_dataset);
//
//    std::string infile = path_to_inputdata + "mu3e_run_" + get_padded_string(bg_run, 6, '0') + ".root";
//    std::string background_sid_file = path_to_bkg_tids_dataset + "BackgroundSIDtracks_" + getfileidtag(bg_run, mode, spWbins, spZbins) + ".root";
//    std::vector<std::vector<unsigned int>> cosmic_spid_tracks;
//
//    //First check if a file exists already
//    TFile tbkgTIDF(background_sid_file.c_str());
//    if (!tbkgTIDF.IsOpen()) {
//        std::cout << "(INFO)   : No TID File " << background_sid_file << " exists. Starting to create data!" << std::endl;
//
//        //open file again with writing access
//        tbkgTIDF.Close();
//        TFile tbkgTIDF(background_sid_file.c_str(), "recreate");
//
//        // FILE FOR READING BACKGROUND DATA
//        const std::string infile = pathtoBGdata + "mu3e_run_" + get_padded_string(run, 6, '0') + ".root";
//        TFile tinF(infile.c_str());
//        if (!tinF.IsOpen()) {
//            std::cout << "[ERROR] File " << tinF.GetName() << " is not open!" << std::endl;
//            exit(0);
//        }
//
//        //Get the Pattern Engine and Template Manager
//        const int spWbins = (int) sqrt((float) spWZratio * (float) centralTPcount);
//        const int spZbins = (int) sqrt((float) centralTPcount / (float) spWZratio);
//        //for acuiring the background data
//        PatternEngine PEbg(spWbins, spZbins, pathtorunplots);
//
//        TTree *t_mu3e;
//        tinF.GetObject("mu3e", t_mu3e);
//        TTree *t_mu3e_mchits;
//        tinF.GetObject("mu3e_mchits", t_mu3e_mchits);
//
//        Mu3eTree Mu3e = Mu3eTree(t_mu3e);
//        Mu3eMChitsTree Mu3eMChits = Mu3eMChitsTree(t_mu3e_mchits);
//
//        int bg_events = (MAX_BG_ENTRIES == 0 ? Mu3e.my_entries : MAX_BG_ENTRIES);
//        int processed_frames = 0;
//
//        //store all background TIDs as a vector per frame in these vectors, so that the combinatoric
//        //must not be repeated when calling another filter setting.
//        BGcombinatorics BGcmbncsResult;
//
//        //for some plots with bg eff / hits in bg frame
//        std::vector<float> frame_eff;
//        std::vector<int> frame_bghits;
//
//        //bg analysis.
//        for (int frame = 0; frame <= bg_events; frame++) {
//            Mu3e.getEntry(frame);
//            if (PRINTS) Mu3e.Print(frame);
//
//            if (Mu3e.Nhit > MAX_NHITS) {
//                continue;
//            }
//            processed_frames++;
//            print_status_bar(frame, bg_events, "bg event combinatorics", "");
//
//            std::vector<BGhit> bgframehits;
//            BGSortedSIDs hits;
//            std::vector<TemplateID> TIDS;
//            std::map<unsigned short, int> SIDMem;
//            BGhit BGHIT;
//            unsigned short SID;
//            int innerhits = 0;
//            int outerhits = 0;
////            int doublesid = 0;
//
//            //get all hits from one frame
//            assert(Mu3e.Nhit == Mu3e.hit_mc_i->size());
//            for (int hitno = 0; hitno < Mu3e.Nhit; hitno++) {
//                Mu3eMChits.getEntry((*Mu3e.hit_mc_i)[hitno]);
//                if (PRINTS) Mu3eMChits.Print((*Mu3e.hit_mc_i)[hitno]);
//                BGHIT.fill(Mu3eMChits.pos_g_x, Mu3eMChits.pos_g_y, Mu3eMChits.pos_g_z, 0);
//                bgframehits.push_back(BGHIT);
//
//                SID = (unsigned short) PEbg.getSuperPixel(BGHIT.x, BGHIT.y, BGHIT.z);
//                //            if(SIDMem.count(SID) == 0) {
//                //                SIDMem[SID] = 1;
//                //            } else {
//                //                doublesid++;
//                //                continue;
//                //            }
//
//                int layer = PEbg.getLayerFromSPID(SID);
//
//                //sort the hits into their layers
//                if (layer == 3 && BGHIT.y >= 0) {
//                    hits.h0.push_back(SIDtype(SID, BGHIT.type));
//                    outerhits++;
//                } else if (layer == 2 && BGHIT.y >= 0) {
//                    hits.h1.push_back(SIDtype(SID, BGHIT.type));
//                    outerhits++;
//                } else if (layer == 2 && BGHIT.y < 0) {
//                    hits.h2.push_back(SIDtype(SID, BGHIT.type));
//                    outerhits++;
//                } else if (layer == 3 && BGHIT.y < 0) {
//                    hits.h3.push_back(SIDtype(SID, BGHIT.type));
//                    outerhits++;
//                } else {
//                    innerhits++;
//                }
//            }
//
//
//            //Go through all possible combinations of hits and create corresponding Template IDs
//            //for all hits in upper layer 3
//            for (const auto &h0 : hits.h0) {
//                //for all hits in upper layer 2
//                for (const auto &h1 : hits.h1) {
//                    //for all hits in lower layer 2
//                    for (const auto &h2 : hits.h2) {
//                        //for all hits in lower layer 3
//                        for (const auto &h3 : hits.h3) {
//                            TIDS.push_back(TemplateID(h0.SID, h1.SID, h2.SID, h3.SID));
//                        }
//                    }
//                }
//            }
//
//            assert(outerhits == bgframehits.size() - innerhits);
//
//            //save combinatorics result per frame
//            BGcmbncsResult.frames_TIDS.push_back(TIDS); //calculated tids
//            BGcmbncsResult.frames_outerhits.push_back(outerhits); //int number of hits in outer layers per frame
//            BGcmbncsResult.frames_innerhits.push_back(innerhits); //int number of hits in inner layers per frame
//            BGcmbncsResult.frames_bgframehits.push_back(bgframehits.size()); // int number of hits per frame
//            BGcmbncsResult.frames_sortedhits.push_back(hits); //SIDs sorted in layers per frame
//            BGcmbncsResult.total_nhits += bgframehits.size();
//        }
//
//
//        PEbg.displayBinBoundaries();
//        PEbg.displayBinWeightDistribution();
//        PEbg.closePlot();
//
//        std::cout << "(STATUS) : Got background data. processed bg events: " << bg_events << std::endl;
//
//        tbkgTIDF.cd();
//
//        //add some meta data for the
//        int wbins = spWbins;
//        int zbins = spZbins;
//        TTree tT_met("CosmicSIDMeta","Metadata associated with these plots (PE config and dataset)");
//        tT_met.Branch("cosmic_dataset", &cosmic_testing_dataset, "cosmic_dataset/I");
//        tT_met.Branch("cosmic_testing_processed_entries", &processed_entries, "cosmic_processed_entries/I");
//        tT_met.Branch("area0Description", &PE.areaDescript[0], "area0Description/C");
//        tT_met.Branch("area1Description", &PE.areaDescript[1], "area1Description/C");
//        tT_met.Branch("area2Description", &PE.areaDescript[2], "area2Description/C");
//        tT_met.Branch("wBins0", &PE.WBins[0], "wBins0/I");
//        tT_met.Branch("wBins1", &PE.WBins[1], "wBins1/I");
//        tT_met.Branch("wBins2", &PE.WBins[2], "wBins2/I");
//        tT_met.Branch("zBins0", &PE.ZBins[0], "zBins0/I");
//        tT_met.Branch("zBins1", &PE.ZBins[1], "zBins1/I");
//        tT_met.Branch("zBins2", &PE.ZBins[2], "zBins2/I");
//        tT_met.Branch("mode", &PE.mode, "mode/I");
//        tT_met.Branch("sp_count", &centralTPcount, "spcount/i");
//        tT_met.Branch("sp_target_ratio", &spWZratio, "sp_target_ratio/F");
//        tT_met.Fill();
//        tT_met.Write();
//
//        TTree tT_sids("CosmicSIDtracks","SIDs of cosmic muon tracks (one cosmic track per entry)");
//        std::vector<unsigned int> SPIDs;
//        tT_sids.Branch("cosmic_track_sids", &SPIDs);
//
//        for(int track=0; track < cosmic_spid_tracks.size(); track++) {
//            SPIDs = cosmic_spid_tracks[track];
//            tT_sids.Fill();
//        }
//
//        tT_sids.Write();
//        tbkgTIDF.Close();
//
//        return cosmic_spid_tracks;
//
//    } else {
//
//        std::cout << "(INFO)   : Opened file " << background_sid_file << std::endl;
//
//        TTree *tT_met;
//        TTree *tT_sids;
//
//        tbkgTIDF.GetObject("CosmicSIDMeta", tT_met);
//        tbkgTIDF.GetObject("CosmicSIDtracks", tT_sids);
//
//        int tf_cosmic_dataset;
//        int tf_cosmic_testing_processed_entries;
//        int tf_mode;
//        int wBins;
//        int zBins;
//        unsigned int tf_sp_count;
//        float tf_sp_target_ratio;
//
//        std::vector<unsigned int> *SPIDs = nullptr;
//        cosmic_spid_tracks.clear();
//
//
//        tT_met->SetBranchAddress("cosmic_dataset", &tf_cosmic_dataset);
//        tT_met->SetBranchAddress("cosmic_testing_processed_entries", &tf_cosmic_testing_processed_entries);
//        tT_met->SetBranchAddress("wBins0", &wBins);
//        tT_met->SetBranchAddress("zBins0", &zBins);
//        tT_met->SetBranchAddress("sp_count", &tf_sp_count);
//        tT_met->SetBranchAddress("sp_target_ratio", &tf_sp_target_ratio);
//        tT_met->GetEntry(0);
//
//        tT_sids->SetBranchAddress("cosmic_track_sids",&SPIDs);
//        int entries = tT_sids->GetEntries();
//
//        std::cout << "(INFO)   : Got data from configuration: wbins " << wBins << " | zbins " << zBins << " | mode " << mode << std::endl;
//        std::cout << "(INFO)   : Processed entries " << tf_cosmic_testing_processed_entries << " entries in file " <<  entries << std::endl;
//
//        for(int i=0; i<entries; i++) {
//            tT_sids->GetEntry(i);
//            cosmic_spid_tracks.push_back(*SPIDs);
//        }
//
//
//        assert(cosmic_spid_tracks.size() == tf_cosmic_testing_processed_entries);
//
//        return cosmic_spid_tracks;
//    }
//
//}

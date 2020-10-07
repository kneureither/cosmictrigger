//
// Created by Konstantin Neureither on 02.10.20.
//

#ifndef COSMICTRIGGER_TESTTEMPLATEBANK_H
#define COSMICTRIGGER_TESTTEMPLATEBANK_H

#include "../Mu3eCosPat/include/TemplateBank.h"

void testTemplateBank() {
    int dataset = 11;
    int mode = 0;
    int spWbins = 224;
    int spZbins = 7;

    std::string pathtorunplots = "output/tests/testTemplateBank/";
    std::string pathtodatasettemplatedata = "data/TemplateData/dataset_011/";
    float TB_STOPPING_EFF = 0.6;

    TemplateBank TB(pathtorunplots, dataset, mode, spWbins, spZbins);
    TB.PRINTS = false;
//    TB.readAMfromFile(pathtodatasettemplatedata, TB_STOPPING_EFF);

    //test Area type function

    TemplateID TID;
    //center
    TID.HIDS[0] = 0x3FF1; TID.HIDS[1] = 0x2AA2; TID.HIDS[2] = 0x2BB3; TID.HIDS[3] = 0x3233;
    assert(TB.GetTypeOfTID(TID) == CECE);

    //recurl R
    TID.HIDS[0] = 0x7FF7; TID.HIDS[1] = 0x6AA6; TID.HIDS[2] = 0x6BB6; TID.HIDS[3] = 0x7237;
    assert(TB.GetTypeOfTID(TID) == RRRR);

    //recurl L
    TID.HIDS[0] = 0xBFFB; TID.HIDS[1] = 0xAAAA; TID.HIDS[2] = 0xABBA; TID.HIDS[3] = 0xB23B;
    assert(TB.GetTypeOfTID(TID) == RLRL);

    //recurl L center
    TID.HIDS[0] = 0xBFFB; TID.HIDS[1] = 0xAAAA; TID.HIDS[2] = 0x2BB3; TID.HIDS[3] = 0x3232;
    assert(TB.GetTypeOfTID(TID) == RLCE);

    //recurl R center
    TID.HIDS[0] = 0x7FF7; TID.HIDS[1] = 0x6AA2; TID.HIDS[2] = 0x2BB2; TID.HIDS[3] = 0x3233;
    assert(TB.GetTypeOfTID(TID) == RRCE);

    std::cout << "(STATUS) : Finished Template Bank Area Type test. Success." << std::endl;

    // second test (just check console out)
    // options: {ALL, CENTER_ONLY, RECURL_ONLY, MIXED_ONLY, , NO_CENTER, CUT_ON_FREQ};
    TB.readAMfromFile(pathtodatasettemplatedata, TB_STOPPING_EFF, NO_CENTER);

    // test histogram template type
    TB.PlotTemplateTypeDistribution();



    std::cout << "\n(STATUS) : Finished all tests for Template Bank!" << std::endl;

}

#endif //COSMICTRIGGER_TESTTEMPLATEBANK_H

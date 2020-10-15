//
// Created by Konstantin Neureither on 15.07.20.
//
// This code is not properly debugged. Especially the compiler directives were quickly implemented
//

#ifndef COSMICTRIGGER_TESTTEMPLATEID_H
#define COSMICTRIGGER_TESTTEMPLATEID_H

#include "../CTCoreModules/inc/TemplateData.h"
#include <cassert>

void testTemplateID() {
    TemplateID TID1;
    TemplateID TID2;
    TemplateID TID3;
    TemplateID TID4;
    TemplateID TID5;

    TID1.HIDS[0] = 0x111B;
    TID1.HIDS[1] = 0x139A;
    TID1.HIDS[2] = 0;
    TID1.HIDS[3] = 0;
#if TID_LEN == 8
    TID1.HIDS[4] = 0;
    TID1.HIDS[5] = 0;
    TID1.HIDS[6] = 0x2BEA;
    TID1.HIDS[7] = 0x2E6B;
#endif

    TID2.HIDS[0] = 0x111B;
    TID2.HIDS[1] = 0x139A;
    TID2.HIDS[2] = 0;
    TID2.HIDS[3] = 0;
#if TID_LEN == 8
    TID2.HIDS[4] = 0;
    TID2.HIDS[5] = 0;
    TID2.HIDS[6] = 0x2BEA;
    TID2.HIDS[7] = 0x2E6B;
#endif

    TID3.HIDS[0] = 0x111A;
    TID3.HIDS[1] = 0x139A;
    TID3.HIDS[2] = 0;
    TID3.HIDS[3] = 0;
#if TID_LEN == 8
    TID3.HIDS[4] = 0;
    TID3.HIDS[5] = 0;
    TID3.HIDS[6] = 0x2BEA;
    TID3.HIDS[7] = 0x2E6B;
#endif

    TID4.HIDS[0] = 0x111A;
    TID4.HIDS[1] = 0x139A;
    TID4.HIDS[2] = 0x139A;
    TID4.HIDS[3] = 0;
#if TID_LEN == 8
    TID4.HIDS[4] = 0;
    TID4.HIDS[5] = 0;
    TID4.HIDS[6] = 0x2BEA;
    TID4.HIDS[7] = 0x2E6B;
#endif

    TID5.HIDS[0] = 0x111A;
    TID5.HIDS[1] = 0x139A;
    TID5.HIDS[2] = 0x139A;
    TID5.HIDS[3] = 0;
#if TID_LEN == 8
    TID5.HIDS[4] = 0;
    TID5.HIDS[5] = 0;
    TID5.HIDS[6] = 0x2BEA;
    TID5.HIDS[7] = 0x2E6B;
#endif

    std::cout << TID1.toString() << std::endl;
    std::cout << TID2.toString() << std::endl;
    std::cout << TID3.toString() << std::endl;

    assert(TID1 == TID2);
    assert(!(TID1 == TID3));
    assert(TID3 < TID1);
    assert(TID3 < TID2);
    assert(TID3 < TID4);
    assert(!(TID4 < TID5));

}

#endif //COSMICTRIGGER_TESTTEMPLATEID_H

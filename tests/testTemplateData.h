//
// Created by Konstantin Neureither on 09.07.20.
//

#ifndef COSMICTRIGGER_TESTTEMPLATEDATA_H
#define COSMICTRIGGER_TESTTEMPLATEDATA_H

#include "../Mu3eCosPat/include/TemplateData.h"

void TESTtemplateData() {
    TemplateID TID;
    TID.HIDS[0] = 0xFFF1;
    TID.HIDS[1] = 0xAAA2;
    TID.HIDS[2] = 0xBBB3;
    TID.HIDS[3] = 0xCCC4;

    short s;
    int i;
    long l;
    unsigned long long ll;

    std::cout << TID.toString() << std::endl;

    std::cout <<"short " << sizeof(s) << std::endl;
    std::cout <<"int " << sizeof(i) << std::endl;
    std::cout <<"long " << sizeof(l) << std::endl;
    std::cout <<"long long " << sizeof(ll) << std::endl;
}

#endif //COSMICTRIGGER_TESTTEMPLATEDATA_H

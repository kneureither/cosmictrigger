//
// Created by Konstantin Neureither on 23.06.20.
//

#ifndef COSMICTRIGGER_TEMPLATEDATA_H
#define COSMICTRIGGER_TEMPLATEDATA_H

#include <vector>

struct TemplateData {
    unsigned int SPIDs[16];
    float xp[16];
    float yp[16];
    float zp[16];

    int count;

    float p;
    float dca;
    float phi;
    float theta;

    TemplateData(unsigned int * SPIDs, float * xps, float * yps, float * zps, int count, float p, float dca, float phi, float theta);
    TemplateData(unsigned int * SPIDs, int count, float p, float dca, float phi, float theta);
    TemplateData(const TemplateData &other);
    TemplateData& operator=(const TemplateData& other);
};


#endif //COSMICTRIGGER_TEMPLATEDATA_H

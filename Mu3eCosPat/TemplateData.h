//
// Created by Konstantin Neureither on 23.06.20.
//

#ifndef COSMICTRIGGER_TEMPLATEDATA_H
#define COSMICTRIGGER_TEMPLATEDATA_H

#include <vector>
#include "basicDefines.h"

struct TemplateData {
    unsigned int SPIDs[TRIPLET_HIT_ARRAY_LENGTH];
    float xp[TRIPLET_HIT_ARRAY_LENGTH];
    float yp[TRIPLET_HIT_ARRAY_LENGTH];
    float zp[TRIPLET_HIT_ARRAY_LENGTH];

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

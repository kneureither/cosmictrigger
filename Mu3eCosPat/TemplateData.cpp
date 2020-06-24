//
// Created by Konstantin Neureither on 23.06.20.
//

#include "TemplateData.h"

TemplateData::TemplateData(unsigned int *SPIDs, float *xps, float *yps, float *zps, int count, float p, float dca, float phi,
                           float theta) {
    this->p = p;
    this->dca = dca;
    this->phi = phi;
    this->theta = theta;
    this->count = count;

    for(int i = 0; i<count; i++) {
        this->SPIDs[i] = SPIDs[i];
        this->xp[i] = xps[i];
        this->yp[i] = yps[i];
        this->zp[i] = zps[i];
    }
}

TemplateData::TemplateData(unsigned int *SPIDs, int count, float p, float dca, float phi, float theta) {
    this->p = p;
    this->dca = dca;
    this->phi = phi;
    this->theta = theta;
    this->count = count;

    for(int i = 0; i<count; i++) {
        this->SPIDs[i] = SPIDs[i];
    }
}

TemplateData::TemplateData(const TemplateData &other) {
    //copy constructor
    this->p = other.p;
    this->dca = other.dca;
    this->phi = other.phi;
    this->theta = other.theta;
    this->count = other.count;

    for(int i = 0; i < other.count; i++) {
        this->SPIDs[i] = other.SPIDs[i];
        this->xp[i] = other.xp[i];
        this->yp[i] = other.yp[i];
        this->zp[i] = other.zp[i];
    }

}

TemplateData &TemplateData::operator=(const TemplateData &other) {
    //copy assignment constructor
    this->p = other.p;
    this->dca = other.dca;
    this->phi = other.phi;
    this->theta = other.theta;
    this->count = other.count;

    for(int i = 0; i < other.count; i++) {
        this->SPIDs[i] = other.SPIDs[i];
        this->xp[i] = other.xp[i];
        this->yp[i] = other.yp[i];
        this->zp[i] = other.zp[i];
    }

}

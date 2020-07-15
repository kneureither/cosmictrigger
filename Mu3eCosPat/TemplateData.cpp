//
// Created by Konstantin Neureither on 23.06.20.
//

#include "TemplateData.h"

TemplateData::TemplateData(const int &count, const float &p, const float &dca, const float &phi,
                           const float &theta) {
    this->p = p;
    this->dca = dca;
    this->phi = phi;
    this->theta = theta;
    this->count = count;
    //uEventId
}

TemplateData::TemplateData(const TemplateData &other) {
    //copy constructor
    this->p = other.p;
    this->dca = other.dca;
    this->phi = other.phi;
    this->theta = other.theta;
    this->count = other.count;
    this->uEventID = other.uEventID;
}

TemplateData &TemplateData::operator=(const TemplateData &other) {
    //copy assignment constructor
    this->p = other.p;
    this->dca = other.dca;
    this->phi = other.phi;
    this->theta = other.theta;
    this->count = other.count;
    this->uEventID = other.uEventID;

    return *this;
}



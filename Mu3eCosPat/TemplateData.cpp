//
// Created by Konstantin Neureither on 23.06.20.
//

#include "TemplateData.h"

TemplateData::TemplateData(const int &count, const float &p, const float &dca, const float &phi,
                           const float &theta) {
    this->count.push_back(count);
    this->p.push_back(p);
    this->dca.push_back(dca);
    this->phi.push_back(phi);
    this->theta.push_back(theta);
    this->frequency = 1;
    //uEventId
}

TemplateData::TemplateData(const TemplateData &other) {
    //copy constructor
    this->count = other.count;
    this->p = other.p;
    this->dca = other.dca;
    this->phi = other.phi;
    this->theta = other.theta;
    this->frequency = other.frequency;
    this->uEventID = other.uEventID;
}

TemplateData &TemplateData::operator=(const TemplateData &other) {
    //copy assignment constructor
    this->count = other.count;
    this->p = other.p;
    this->dca = other.dca;
    this->phi = other.phi;
    this->theta = other.theta;
    this->frequency = other.frequency;
    this->uEventID = other.uEventID;

    return *this;
}

void TemplateData::add_track(const int &count, const float &p, const float &dca, const float &phi, const float &theta) {
    this->count.push_back(count);
    this->p.push_back(p);
    this->dca.push_back(dca);
    this->phi.push_back(phi);
    this->theta.push_back(theta);
    this->frequency++;
}

TemplateData::TemplateData() {
    this->frequency =0;
}





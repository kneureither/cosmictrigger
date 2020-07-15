//
// Created by Konstantin Neureither on 23.06.20.
//

#ifndef COSMICTRIGGER_TEMPLATEDATA_H
#define COSMICTRIGGER_TEMPLATEDATA_H

#define TID_LEN 4

#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include "basicDefines.h"

struct TemplateID {
    //SIDs ids in order of hits each hit 16bit
    unsigned short HIDS[TID_LEN];

    std::string toString() {
        char buffer[TID_LEN * 4]; //FIXME not working for 24 len. (when TID_LEN = 6)
        for(int i = 0; i<TID_LEN; i++) {
            sprintf(&(buffer[i*4]), "%04X", HIDS[i]);
        }
        return buffer;
    }

    bool operator<(const TemplateID& other) const {
        for(int i=0; i < TID_LEN; i++) {
            if(this->HIDS[i] < other.HIDS[i]) {
//                std::cout << "ordered after inspecting the " << i << " element" << std::endl;
                return true;
            } else if (this->HIDS[i] > other.HIDS[i]){
                return false;
            } else {
                if(i==TID_LEN-1) return false;
            }
        }
    }

    bool operator==(const TemplateID& other) const {
        return std::equal(std::begin(this->HIDS), std::end(this->HIDS), std::begin(other.HIDS));
    }

    TemplateID &operator=(const TemplateID &other) {
        for(int i=0; i < TID_LEN; i++) {
            this->HIDS[i] = other.HIDS[i];
        }
        return *this;
    }
};

struct TemplateData {
    int count;
    float p;
    float dca;
    float phi;
    float theta;
    unsigned int uEventID;

    TemplateData(const int &count, const float &p, const float &dca, const float &phi, const float &theta);
    TemplateData(const TemplateData &other);
    TemplateData& operator=(const TemplateData& other);
};


#endif //COSMICTRIGGER_TEMPLATEDATA_H

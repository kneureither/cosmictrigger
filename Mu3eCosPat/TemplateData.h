//
// Created by Konstantin Neureither on 23.06.20.
//

#ifndef COSMICTRIGGER_TEMPLATEDATA_H
#define COSMICTRIGGER_TEMPLATEDATA_H

#define TID_LEN 8

#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include "basicDefines.h"

struct TemplateID {
    //SIDs ids in order of hits each hit 16bit
    unsigned short HIDS[TID_LEN];

    std::string toString() {
        char buffer[4*TID_LEN];
        for(int i = 0; i<TID_LEN; i++) {
            sprintf(&(buffer[i*4]), "%04X", HIDS[i]);
        }
        return buffer;
    }

    bool operator<(const TemplateID& other) const {
        for(int i=0; i < TID_LEN; i++) {
            if(this->HIDS[i] < other.HIDS[i]) {
                std::cout << "ordered after inspecting the " << i << " element" << std::endl;
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
};

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

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
#include <sstream>
#include "basicDefines.h"

struct TemplateID {
    //SIDs ids in order of hits each hit 16bit
    unsigned short HIDS[TID_LEN];

    std::string toString() {
        std::stringstream ss;
        char buffer[4]; //FIXME not working for 24 len. (when TID_LEN = 6)
        for(int i = 0; i<TID_LEN; i++) {
            sprintf(buffer, "%04X", HIDS[i]);
            ss << buffer;
        }
        return ss.str();
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

        return false;
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
    int frequency = 0;
    std::vector<int> count;
    std::vector<float> p;
    std::vector<float> dca;
    std::vector<float> phi;
    std::vector<float> theta;
    std::vector<unsigned int> uEventID;

    TemplateData(const int &count, const float &p, const float &dca, const float &phi, const float &theta);
    TemplateData();
    TemplateData(const TemplateData &other);
    TemplateData& operator=(const TemplateData& other);
    void add_track(const int &count, const float &p, const float &dca, const float &phi, const float &theta);
};


#endif //COSMICTRIGGER_TEMPLATEDATA_H

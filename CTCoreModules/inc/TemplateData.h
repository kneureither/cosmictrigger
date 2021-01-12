//
// Created by Konstantin Neureither on 23.06.20.
//

#ifndef COSMICTRIGGER_TEMPLATEDATA_H
#define COSMICTRIGGER_TEMPLATEDATA_H

#define TID_LEN 4
/**
 * Theoretically, it is possible to use 8 instead of 4 hits for the templates. This means that also the inner detector
 * layers are taken into account. Tracks that do not pass the inner layers are then described as
 * [SID1, SID2, 0, 0 ,0 ,0 SID3, SID5].
 *
 * However, even if the data strture supports TID_LEN > 4, there is some routines that need to be updated manually.
 * So be cautious if you plan on doing so.
 */

#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include "basicDefines.h"

struct TemplateID {
    /**
     * Contains a template consisting of SIDs and some convenience function, such as (copy) constructors.
     * Theoretically supports TID_LEN>4  with some adaptions
     */

    //SIDs ids in order of hits each hit 16bit
    unsigned short HIDS[TID_LEN];

    std::string toString() {
        std::stringstream ss;
        char buffer[4]; //FIXME not working for 24 len (when TID_LEN = 6). Should be working for 32 len
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

    TemplateID(const unsigned short &h0, const unsigned short &h1, const unsigned short &h2, const unsigned short &h3) {
        HIDS[0] = h0;
        HIDS[1] = h1;
        HIDS[2] = h2;
        HIDS[3] = h3;
    }

    TemplateID() {}
};

struct TemplateData {
    /**
     * This is the data that is stored for every template in the database file.
     *
     * In the beginning of the study it was thinkable, that storing track parameters for the templates could be helpful,
     * but was eventually never used. Maybe it is needed in the future, so it was left in the code.
     */

    int frequency = 0;
    //FIXME z0 is missing!
    std::vector<int> count;
    std::vector<float> p;
    std::vector<float> dca;
    std::vector<float> phi;
    std::vector<float> theta;
    std::vector<unsigned int> uEventID;

    TemplateData(const int &count, const float &p, const float &dca, const float &phi, const float &theta);
    TemplateData(const int frequency);
    TemplateData();
    TemplateData(const TemplateData &other);
    TemplateData& operator=(const TemplateData& other);
    void add_track(const int &count, const float &p, const float &dca, const float &phi, const float &theta);
};


#endif //COSMICTRIGGER_TEMPLATEDATA_H

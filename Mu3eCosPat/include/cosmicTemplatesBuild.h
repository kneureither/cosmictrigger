//
// Created by Konstantin Neureither on 09.07.20.
//

#ifndef COSMICTRIGGER_COSMICTEMPLATESBUILD_H
#define COSMICTRIGGER_COSMICTEMPLATESBUILD_H
#include <vector>
#include "SlimSegsTree.h"
#include "PatternEngine.h"

void cosmicTemplatesBuild(const int dataset, unsigned int centralTPcount, float spWZratio, int combination_id,
                          float max_efficiency, bool append_to_outfile);

static bool getSymmetricRefHits(std::vector<float> &xpr, std::vector<float> &ypr,
        std::vector<float> &zpr, std::vector<int> &layerpr, SlimSegsTree &SlimSegs, const int &minlayer) {

    std::vector<int>::iterator itfwd = SlimSegs.layerp.begin();
    std::vector<int>::iterator itbackw = SlimSegs.layerp.end() - 1;
    int indexfwd = std::distance(SlimSegs.layerp.begin(), itfwd);
    int indexbackw = std::distance(SlimSegs.layerp.begin(), itbackw);
    int layercounter[4] = {0,0,0,0};
    int countfrontinserted = 0;
    unsigned int SPID;

    while(itfwd != itbackw) {
        while(!(layercounter[*itfwd] < 2) && itfwd != itbackw) {
            itfwd++;
            indexfwd = std::distance(SlimSegs.layerp.begin(), itfwd);
        }

        if (itfwd != itbackw) {
            if(*itfwd >= minlayer) {
                layercounter[*itfwd]++;
                xpr.insert(xpr.begin() + countfrontinserted, SlimSegs.xp[indexfwd]);
                ypr.insert(ypr.begin() + countfrontinserted, SlimSegs.yp[indexfwd]);
                zpr.insert(zpr.begin() + countfrontinserted, SlimSegs.zp[indexfwd]);
                layerpr.insert(layerpr.begin() + countfrontinserted, SlimSegs.layerp[indexfwd]);
                countfrontinserted++;
            }
        } else {
            break;
        }


        while(!(layercounter[*itbackw] < 2) && itfwd != itbackw) {
            itbackw--;
            indexbackw = std::distance(SlimSegs.layerp.begin(), itbackw);
        }

        if(itfwd != itbackw) {
            if(*itbackw >= minlayer) {
                layercounter[*itbackw]++;
                xpr.insert(xpr.begin() + countfrontinserted, SlimSegs.xp[indexbackw]);
                ypr.insert(ypr.begin() + countfrontinserted, SlimSegs.yp[indexbackw]);
                zpr.insert(zpr.begin() + countfrontinserted, SlimSegs.zp[indexbackw]);
                layerpr.insert(layerpr.begin() + countfrontinserted, SlimSegs.layerp[indexbackw]);
            }
        } else {
            break;
        }

        itfwd++;
        indexfwd = std::distance(SlimSegs.layerp.begin(), itfwd);

        if(itfwd != itbackw) {
            itbackw--;
            indexbackw = std::distance(SlimSegs.layerp.begin(), itbackw);
        } else {
            break;
        }
    }

    assert(xpr.size() == ypr.size());
    assert(ypr.size() == zpr.size());
    assert(zpr.size() == layerpr.size());

    if(layercounter[3] < 2 || layercounter[2] < 2) {
        return false;
    } else {
        return true;
    }
}



#endif //COSMICTRIGGER_COSMICTEMPLATESBUILD_H

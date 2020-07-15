//
// Created by Konstantin Neureither on 09.07.20.
//

#ifndef COSMICTRIGGER_BUILDCOSMICTRAJTEMPLATESSCRIPT_H
#define COSMICTRIGGER_BUILDCOSMICTRAJTEMPLATESSCRIPT_H
#include <vector>
#include "SlimSegsTree.h"
#include "PatternEngine.h"

void buildCosmicTemplatesScript(const int);

static void getSymmetricRefHits(std::vector<unsigned int> &SPIDs, SlimSegsTree &SlimSegs, PatternEngine &PE) {
    std::vector<int>::iterator itfwd = SlimSegs.layerp.begin();
    std::vector<int>::iterator itbackw = SlimSegs.layerp.end() - 1;
    int indexfwd = std::distance(SlimSegs.layerp.begin(), itfwd);
    int indexbackw = std::distance(SlimSegs.layerp.begin(), itbackw);
    int layercounter[4] = {0,0,0,0};
    int countfrontinserted = 0;
    unsigned int SPID;

    int minlayer = 0;

    while(itfwd != itbackw) {
        while(!(layercounter[*itfwd] < 2) && itfwd != itbackw) {
            itfwd++;
            indexfwd = std::distance(SlimSegs.layerp.begin(), itfwd);
        }

        if (itfwd != itbackw) {
            if(*itfwd >= minlayer) {
                layercounter[*itfwd]++;
                SPID = PE.getSuperPixel(SlimSegs.xp[indexfwd], SlimSegs.yp[indexfwd], SlimSegs.zp[indexfwd]);
                SPIDs.insert(SPIDs.begin() + countfrontinserted, SPID);
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
                SPID = PE.getSuperPixel(SlimSegs.xp[indexbackw], SlimSegs.yp[indexbackw],
                                        SlimSegs.zp[indexbackw]);
                SPIDs.insert(SPIDs.begin() + countfrontinserted, SPID);
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
}



#endif //COSMICTRIGGER_BUILDCOSMICTRAJTEMPLATESSCRIPT_H

//
// Created by Konstantin Neureither on 06.07.20.
//

#include "SPCalculations.h"
#include <assert.h>

int SPCalculations::getLayerFromSPID(const unsigned int SPID) {
    unsigned int zone = (SPID & 0xF);
    return zone % 4;
}

int SPCalculations::getAreaFromSPID(const unsigned int SPID) {
    unsigned int zone = (SPID & 0xF);
    return zone / 4;
}

int SPCalculations::getIndexFromSPID(const unsigned int SPID) {
    assert((SPID & 0xFFFF0000) == 0x00000000); //no more entries, than 16 bits

    unsigned int index = (unsigned int) ((SPID & 0xFFF0) >> 4);
    return (int) index;
}

int SPCalculations::computeIndex(const int zSPIndex, const int phiSPIndex, const int wBinCount) {
    return zSPIndex * wBinCount + phiSPIndex;
}

unsigned int SPCalculations::computeSPID(const int layer, const int area, const unsigned int index) {
    unsigned int zone = area * 4 + layer;
    unsigned int SPID = (unsigned int) (index & 0xFFFF) << 4;
    SPID |= (zone & 0xF);
    return SPID;
}

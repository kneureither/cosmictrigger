//
// Created by Konstantin Neureither on 06.07.20.
//

#ifndef COSMICTRIGGER_SPCALCULATIONS_H
#define COSMICTRIGGER_SPCALCULATIONS_H


class SPCalculations {
public:
    static const int areaCount = 3;
    bool PRINTS = false;

    ////Geometry
    //area index convention: 0=central, 1=recurl right, 2=recurl left
    const float DetectorZmin[areaCount] = {-200.00, 200, -600};
    const float DetectorZmax[areaCount] = {200.00, 600, -200};
    char areaDescript[areaCount][8] = {"central", "recurlR", "recurlL"};
    const float layerBoundaries[5] = {0.0, 26.00, 51.0, 78.5, 100.00}; //radial decision boundaries for layers

    ////Utility functions
    int getLayerFromSPID(const unsigned int SPID);
    int getAreaFromSPID(const unsigned int SPID);
    int getIndexFromSPID(const unsigned int SPID);
    int computeIndex(const int zSPIndex, const int phiSPIndex, const int wBinCount);
    unsigned int computeSPID(const int layer, const int area, const unsigned int index);
};

#endif //COSMICTRIGGER_SPCALCULATIONS_H

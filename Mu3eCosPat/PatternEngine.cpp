//
// Created by Konstantin Neureither on 06.07.20.
//

#include "include/PatternEngine.h"

//init class members and Area pattern engines
PatternEngine::PatternEngine(const int spWBins, const int spZBins, std::string plottingpath)
    : PES{{spWBins, spZBins, 0, 0, plottingpath},
          {spWBins, spZBins, 0, 1, plottingpath},
          {spWBins, spZBins, 0, 2, plottingpath}}
 {
    for(int area=0; area<this->areaCount; area++) {
        this->WBins[area] = spWBins;
        this->ZBins[area] = spZBins;
    }
    this->plottingpath = plottingpath;
    this->mode = 0;

 }

PatternEngine::PatternEngine(const int spWBinsCenter, const int spZBinsCenter,
                             const int spWBinsRecurl, const int spZBinsRecurl, std::string plottingpath)
        : PES{{spWBinsCenter, spZBinsCenter, 0, 0, plottingpath},
              {spWBinsRecurl, spZBinsRecurl, 0, 1, plottingpath},
              {spWBinsRecurl, spZBinsRecurl, 0, 2, plottingpath}}
{
    this->WBins[0] = spWBinsCenter;
    this->ZBins[0] = spZBinsCenter;
    for(int area=1; area<this->areaCount; area++) {
        this->WBins[area] = spWBinsRecurl;
        this->ZBins[area] = spZBinsRecurl;
    }
    this->plottingpath = plottingpath;
    this->mode = 0;
}

int PatternEngine::getArea(const float z) {
    for(int area=0; area<this->areaCount; area++){
        if(DetectorZmin[area] <= z && z <= DetectorZmax[area]) {
            return area;
        }
    }
    return 0;
}

unsigned int PatternEngine::getSuperPixel(const float x, const float y, const float z) {
    int area = getArea(z);
    return PES[area].getSuperPixel(x, y, z);
}

void PatternEngine::displayBinBoundaries() {
    for(int area=0; area<this->areaCount; area++) PES[area].displayBinBoundaries();
}

void PatternEngine::displayBinWeightDistribution() {
    for(int area=0; area < this->areaCount; area++) this->PES[area].displayBinWeightDistribution();
}

void PatternEngine::closePlot() {
    for(int area=0; area < this->areaCount; area++) this->PES[area].closePlot();
}

std::string PatternEngine::getModeTag() {
    return PES[0].getRunSpecs();
}



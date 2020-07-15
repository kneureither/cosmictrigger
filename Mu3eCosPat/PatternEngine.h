//
// Created by Konstantin Neureither on 06.07.20.
//

#ifndef COSMICTRIGGER_PATTERNENGINE_H
#define COSMICTRIGGER_PATTERNENGINE_H

#include "PatternEngineSingle.h"
#include <map>


class PatternEngine : public SPCalculations {
private:
    PatternEngineSingle PES[areaCount];

    int WBins[areaCount];
    int ZBins[areaCount];

    std::string plottingpath;

public:
    PatternEngine(const int spWBins, const int spZBins, std::string plottingpath);
    PatternEngine(const int spWBinsCenter, const int spZBinsCenter,
                  const int spWBinsRecurl, const int spZBinsRecurl, std::string plottingpath);
    int getArea(const float z);

    unsigned int getSuperPixel(const float x, const float y, const float z);
    void displayBinBoundaries();
    void displayBinWeightDistribution();
    std::string getModeTag();
    void closePlot();
};


#endif //COSMICTRIGGER_PATTERNENGINE_H

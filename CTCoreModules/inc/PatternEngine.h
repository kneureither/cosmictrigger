//
// Created by Konstantin Neureither on 06.07.20.
//

#ifndef COSMICTRIGGER_PATTERNENGINE_H
#define COSMICTRIGGER_PATTERNENGINE_H

#include "PatternEngineSingle.h"
#include <map>


class PatternEngine : public SPCalculations {
    /**
     * Handles pixel hit (cartesian coord) to super pixel assignment and super pixel mapping in general. For more
     * details see PatternEngineSingle.
     *
     * This class handles the different PatternEngineSingle modules that all manage a single detector station and
     * provides an interface that covers the whole detector and internally forwards requests to the PES responsible.
     *
     * Remember to call closePlot() at the end of your script, if you have used any of the plotting scripts.
     */

private:

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

    PatternEngineSingle PES[areaCount];
    int WBins[areaCount];
    int ZBins[areaCount];
    int mode;
};


#endif //COSMICTRIGGER_PATTERNENGINE_H

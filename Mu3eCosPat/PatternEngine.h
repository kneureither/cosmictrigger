//
// Created by Konstantin Neureither on 16.06.20.
//

#ifndef COSMICTRIGGER_PATTERNENGINE_H
#define COSMICTRIGGER_PATTERNENGINE_H


class PatternEngine {

private:
    int getLayer(float x, float y);
    int getZSP(float z);
    int getXSP(float x);
    int getYSP(float y);

    float getRadius(float x, float y);
    float getHitPhi(float x, float y);

    float layerFactor[4] = {0.0}; // scale down x boundaries with reference to outer layer
    const float layerBoundaries[4] = {0.0, 40.0, 50.0, 60.0}; //radial decision boundaries for layers

    int mode = 0; //default kartesian x,y,z mode
    float initialSPXwidth;
    float initialSPZwidth;
    float initialSPPhi; // for alternative definition of SP

public:
    PatternEngine(float spWidth, float spLength, int mode);
    float getSuperPixel(float x, float y, float z);
};
#endif //COSMICTRIGGER_PATTERNENGINE_H

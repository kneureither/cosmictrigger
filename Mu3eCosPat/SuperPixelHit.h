//
// Created by Konstantin Neureither on 16.06.20.
//

#ifndef COSMICTRIGGER_SUPERPIXELHIT_H
#define COSMICTRIGGER_SUPERPIXELHIT_H


class SuperPixelHit {
private:
    int layer;
    int spPhi;
    int spZ;
    int area;
    int binsZ;
    int binsPhi;
    float initX;
    float initY;
    float initZ;


public:
    unsigned int SPID;
    void testIntegrity();
    SuperPixelHit();
    SuperPixelHit(float x, float y, float z, int layer, int spPhi, int spZ, int area, int binsPhi, int binsZ);
};


#endif //COSMICTRIGGER_SUPERPIXELHIT_H

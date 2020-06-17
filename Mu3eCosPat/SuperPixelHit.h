//
// Created by Konstantin Neureither on 16.06.20.
//

#ifndef COSMICTRIGGER_SUPERPIXEL_H
#define COSMICTRIGGER_SUPERPIXEL_H


class SuperPixel {
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
    SuperPixel();
    SuperPixel(float x, float y, float z, int layer, int spPhi, int spZ, int area, int binsPhi, int binsZ);
};


#endif //COSMICTRIGGER_SUPERPIXEL_H

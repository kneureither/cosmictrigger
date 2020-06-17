//
// Created by Konstantin Neureither on 16.06.20.
//

#include "SuperPixel.h"

SuperPixel::SuperPixel() {
    this->layer=-1;
    this->initX=0.0;
    this->initY=0.0;
    this->initZ=0.0;
    this->spPhi=0;
    this->spZ=0;
    this->area=0;
}

SuperPixel::SuperPixel(float x, float y, float z, int layer, int spPhi, int spZ, int area, int binsPhi, int binsZ) {
    this->layer=layer;
    this->initX=x;
    this->initY=y;
    this->initZ=z;
    this->spPhi=spPhi;
    this->spZ=spZ;
    this->area=area;
    this->binsPhi = binsPhi;
    this->binsPhi = binsZ;
    this->SPID = layer + 10 * area + 100 * (spZ * binsPhi + spPhi);
}

void SuperPixel::testIntegrity() {

}

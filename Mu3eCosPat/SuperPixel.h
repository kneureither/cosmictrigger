//
// Created by Konstantin Neureither on 16.06.20.
//

#ifndef COSMICTRIGGER_SUPERPIXEL_H
#define COSMICTRIGGER_SUPERPIXEL_H


class SuperPixel {
private:
    int layer;
    int spx;
    int spy;
    int spz;
    int initx;
    int inity;
    int initz;

public:
    int SPID;
    void testIntegrity();

};


#endif //COSMICTRIGGER_SUPERPIXEL_H

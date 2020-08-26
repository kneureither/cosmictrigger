//
// Created by Konstantin Neureither on 25.08.20.
//

#ifndef COSMICTRIGGER_COSMICTEMPLATESBGEVAL_H
#define COSMICTRIGGER_COSMICTEMPLATESBGEVAL_H

#include <vector>

void cosmicTemplatesBgEval(const int, unsigned int, const float);

struct BGhit{
    double x;
    double y;
    double z;
    int type;

    void fill(double x, double y, double z, int type) {
        this->x = x;
        this->y = y;
        this->z = z;
        this->type = type;
    }
};

struct SIDtype{
    unsigned short SID;
    int type;
    SIDtype(unsigned short SID, int type) {
        this->SID = SID;
        this->type = type;
    }
};

struct BGSortedHits{
    std::vector<SIDtype> h0; //first hit in TID order
    std::vector<SIDtype> h1; //second hit (layer 2)
    std::vector<SIDtype> h2; //third hit (layer 2, y<0)
    std::vector<SIDtype> h3; //fourth hit (layer 3, y<0)
};


#endif //COSMICTRIGGER_COSMICTEMPLATESBGEVAL_H

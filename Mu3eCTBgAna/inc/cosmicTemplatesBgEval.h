//
// Created by Konstantin Neureither on 25.08.20.
//

#ifndef COSMICTRIGGER_COSMICTEMPLATESBGEVAL_H
#define COSMICTRIGGER_COSMICTEMPLATESBGEVAL_H

void cosmicTemplatesBgEval(const int);

struct bghit{
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
}BGHIT;


#endif //COSMICTRIGGER_COSMICTEMPLATESBGEVAL_H

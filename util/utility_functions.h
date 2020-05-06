#ifndef COSMICTRIGGER_UTILITY_FUNCTIONS_H
#define COSMICTRIGGER_UTILITY_FUNCTIONS_H

#endif //COSMICTRIGGER_UTILITY_FUNCTIONS_H

#include <string>
#include <iostream>
#include "custom_types.h"

struct PXID {
    unsigned int sensor;
    unsigned int column;
    unsigned int row;
    unsigned int columnaddress;
};


// pad an int number
std::string get_padded_string(int number, int n, char c) {
    stringstream ss;
    ss << number;
    string str = ss.str();
    str.insert(str.begin(), n-str.length(), c);
    return str;
}

PXID process_pixel_id(unsigned int pixel_id) {
    PXID pixid;
    pixid.sensor        = (pixel_id)>>16;
    pixid.column        = ((pixel_id)>>8)&0xFF;
    pixid.row           = ((pixel_id))&0xFF;
    pixid.columnaddress = pixid.column + (pixid.sensor << 8);

    return pixid;
}

template <class T>
float vector_mean(std::vector<T> vec) {
    float sum = 0.0;
    for(int i = 0; i < vec.size(); i++) {
        sum += vec[i] / (float) vec.size();
    }
    return sum;
}


void set_plotting_style() {
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(111111);
    gStyle->SetOptFit(1111);
    gStyle->SetPalette(57);
    gStyle->SetOptTitle(0);

    // gStyle->fCanvasImp->SetWindowPosition(0, 400);

    TF1 f1("f1", "sin(x)/x", 0., 10.);
    f1.Draw();

    TCanvas c1("c1", "<This is a canvas", 0,0,400,300);
//    c1.Divide(2,2);
//    c1.cd(1);
//
//    TF1 f1("f1", "sin(x)/x", 0., 10.);
//    f1.Draw();
}


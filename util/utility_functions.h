#ifndef COSMICTRIGGER_UTILITY_FUNCTIONS_H
#define COSMICTRIGGER_UTILITY_FUNCTIONS_H

#endif //COSMICTRIGGER_UTILITY_FUNCTIONS_H

#include <string>
#include <iostream>
#include <experimental/filesystem>
#include "custom_types.h"

namespace fs = std::experimental::filesystem;

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

template <class T>
std::string get_string(T number) {
    stringstream ss;
    ss << number;
    return ss.str();
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

template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void labelAxis(TH1 * h, const char * xtitle, const char * ytitle){
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle(ytitle);

    h->GetXaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetTitleSize(0.04);

    h->GetXaxis()->SetTitleOffset(0.9);
    h->GetYaxis()->SetTitleOffset(1.2);

    h->GetXaxis()->SetLabelSize(0.03);
    h->GetYaxis()->SetLabelSize(0.03);
}

void labelAxis(TGraph * g, const char * xtitle, const char * ytitle) {
    g->GetXaxis()->SetTitle(xtitle);
    g->GetYaxis()->SetTitle(ytitle);

    g->GetXaxis()->SetTitleSize(0.04);
    g->GetYaxis()->SetTitleSize(0.04);

    g->GetXaxis()->SetTitleOffset(0.9);
    g->GetYaxis()->SetTitleOffset(1.8);

    g->GetXaxis()->SetLabelSize(0.03);
    g->GetYaxis()->SetLabelSize(0.03);
}

void setGraphRange(TGraph * g, float xrange_in, float xrange_out, float yrange_in, float yrange_out) {
    g->GetXaxis()->SetRangeUser(xrange_in, xrange_out);
    g->GetYaxis()->SetRangeUser(yrange_in, yrange_out);
}

template <typename T>
void fillHistWithVector(TH1F * h, std::vector<T> &data) {
    for (int i = 0; i < data.size(); i++) {
        h->Fill(data[i]);
    }
}


void SaveCanvas(TCanvas* c1, const std::string& name, const std::string& path) {
    c1->SaveAs((path + "/" + name + ".eps").c_str());
    c1->SaveAs((path + "/" + name + ".png").c_str());
    c1->SaveAs((path + "/" + name + ".pdf").c_str());
}

void check_create_directory(std::string path) {
    if (fs::exists(path.c_str())) {
        cout << "STATUS : path exists: " << path << endl;
    } else {
        fs::create_directory(path);
        cout << "STATUS : no such path - created: " << path << endl;
    }
}

int get_charge_from_type(int type) {
    if (type == 4) {
        //mu-
        return 1;
    } else if (type == 3) {
        //mu+
        return -1;
    } else {
        //not a muon
        return 1;
    }
}

unsigned int get_layer(int sid) {
    if(0 <= sid && sid < 1024) {
        return 0;
    } else if (1024 <= sid && sid < 2048) {
        return 1;
    } else if ((2000 <= sid && sid < 3000) || (10000 <= sid && sid < 11500) || (14000 <= sid && sid < 15200)) {
        return 2;
    } else if ((3000 <= sid && sid < 4000) || (11500 <= sid && sid < 12500) || (15200 <= sid && sid < 16500)) {
        return 3;
    } else {
        return 100;
    }
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


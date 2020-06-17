//
// Created by Konstantin Neureither on 16.06.20.
//
#ifndef PI
#define PI 3.1415926535
#endif //PI

#include <vector>
#include <TFile.h>
#include <TROOT.h>
#include <algorithm>
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <cmath>
#include <assert.h>
#include <iostream>
#include "PatternEngine.h"
#include "plots.h"
#include "SuperPixelHit.h"

int PatternEngine::getXSP(float x) {
    return 0;
}

int PatternEngine::getZSP(float z) {
    return binSearch(this->zBins, z, 0, zBins.size() - 1);
}

int PatternEngine::getPhiSP(float phi) {
    return binSearch(this->phiBins, phi, 0, this->phiBins.size() - 1);
}


int PatternEngine::getLayer(float r) {
    if (layerBoundaries[0] < r && r <= layerBoundaries[1]) {
        return 0;
    } else if (layerBoundaries[1] < r && r <= layerBoundaries[2]) {
        return 1;
    } else if (layerBoundaries[2] < r && r <= layerBoundaries[3]) {
        return 2;
    } else if (layerBoundaries[3] < r && r <= layerBoundaries[4]) {
        return 3;
    } else {
        return -1;
    }
}

void PatternEngine::initializeSPbins() {
    float SPZlength = centralDetectorLength / zBinCount;
    for(int i=0; i <= zBinCount; i++) {
        zBins.push_back(i * SPZlength + centralDetectorZmin);
    }

    //wBinCount only for half detector
    for(int i = 0; i <= wBinCount; i++) {
        phiBins.push_back((float) i * 2 * PI / (float) wBinCount);
    }

    //TODO implement xBins
}

PatternEngine::PatternEngine(float spXpartition, float spZpartition, int mode) {
    this->wBinCount = spXpartition;
    this->zBinCount = spZpartition;
    this->centralDetectorLength = centralDetectorZmax - centralDetectorZmin;

    this->initializeSPbins();
}

PatternEngine::PatternEngine(float spXpartition, float spZpartition, int mode, std::string plottingpath) {
    PatternEngine(spXpartition, spZpartition, mode);
    this->plottingpath = plottingpath;
}

unsigned int PatternEngine::getSuperPixel(float x, float y, float z) {
    float r = getRadius(x, y);
    int layer = getLayer(r);
    int zSPIndex = getZSP(z);
    int phiSPIndex = getPhiSP(getHitPhi(x, y));
    int area = 0;

    printf("Hit SuperPixel params\t x=%f, y=%f, z=%f, layer=%d zIndex=%d, phiIndex=%d\n", x,y,z,layer, zSPIndex, phiSPIndex);
    return layer + 10 * area + 100 * (zSPIndex * wBinCount + phiSPIndex);
}

void PatternEngine::displayBinBoundaries() {
    float widthMax=2*PI;
    float widthMin=0.0;

    auto *canvas = new TCanvas("bin_boundaries", "bin_boundaries", 1200, 900);
    canvas->SetLeftMargin(0.15);
    canvas->SetGrid(1,1);
    canvas->SetTicks(1, 1);

    TH2F * grid = new TH2F("h", "Superpixel bin grid boundaries", 10, -200, 200, 10, 0, 2*PI);
    grid->SetStats(false);
    labelAxis(grid, "z [mm]", "#phi");
    grid->Draw();

    std::vector<TLine*>boundaries;

    for(int i=0; i < phiBins.size(); i++) {
        TLine * l = new TLine(this->centralDetectorZmin, this->phiBins[i], this->centralDetectorZmax, this->phiBins[i]);
        l->SetLineColor(kBlue);
        boundaries.push_back(l);
    }

    for(int i=0; i < zBins.size(); i++) {
        TLine * l = new TLine(this->zBins[i], widthMin, this->zBins[i], widthMax);
        l->SetLineColor(kRed);
        boundaries.push_back(l);
    }

    for(int i = 0; i< boundaries.size(); i++) boundaries[i]->Draw();

    canvas->SaveAs((this->plottingpath + "binBoundaries.pdf").c_str());
}

float PatternEngine::getRadius(float x, float y) {
    return sqrt(pow(x,2) + pow(y,2));
}

float PatternEngine::getHitPhi(float x, float y) {
    return atan2(y, x);
}

template<typename T>
int PatternEngine::binSearch(std::vector<T> targetVector, T value, int start, int end) {
    int size = end-start;
    if (size <= 0) {
        return -1;
    }
    if(size == 1) {
        return start;
    }
    int center = (start + end) / 2;
    if (targetVector[center] <= value && value <= targetVector[end]) {
        return binSearch(targetVector, value, center, end);
    } else if (targetVector[start] <= value && value < targetVector[center]){
        return binSearch(targetVector, value, start, center);
    } else {
        return -1;
    }
}

void PatternEngine::testbinSearch() {
    std::vector<float> bins;
    printf("Testing binSearch().... \n Added data: \n");
    for(int i = 0; i < 20; i++) {
        bins.push_back(i*PI);
        printf("  %f", i*PI);
    }
    int result;
    float data;

    printf("\n\nIn-Between test data results \n");
    for(int i = 0; i<19; i++) {
        data = PI/2 + i*PI;
        result = binSearch(bins, data, 0, bins.size() - 1);
        printf("Testing data=%f, index result=%d\n", data, result);
        assert(i==result);
    }

    printf("\n\nSpot on test data results \n");
    for(int i = 0; i<19; i++) {
        data = i*PI;
        result = binSearch(bins, data, 0, bins.size() - 1);
        printf("Testing data=%f, index result=%d\n", data, result);

        assert(i==result);
    }

    printf("\n\nOut of bounds test \n");
    data = -PI;
    result = binSearch(bins, data, 0, bins.size() - 1);
    printf("Testing data=%f, index result=%d\n", data, result);
    assert(result == -1);

    data = 21 * PI;
    result = binSearch(bins, data, 0, bins.size() - 1);
    printf("Testing data=%f, index result=%d\n", data, result);
    assert(result == -1);
}



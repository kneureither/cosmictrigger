//
// Created by Konstantin Neureither on 16.06.20.
//

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

#include "utilityFunctions.h"
#include "plots.h"
#include "PatternEngineSingle.h"
#include "SuperPixelHit.h"
#include "basicDefines.h"

PatternEngineSingle::PatternEngineSingle() = default;

PatternEngineSingle::PatternEngineSingle(const int spXpartition, const int spZpartition, const int mode) {
    /**
     * Initalize as single PE for only one area.
     *
     * @param spWpartition number of super pixel bins in X or Phi (=W for width) direction
     * @param spZpartition number of super pixel bins in Z direction
     * @param mode type of sp mapping (0=uniform)
     */

    printf("[ ---- INIT SINGLE PATTERN ENGINE ---- ]\n");
    this->initializeMembers(spXpartition, spZpartition, 0, mode);
    this->initializeSPbins();
    printf("wBins=%d, zBins=%d\n", wBinCount, zBinCount);
    printf("----SUCCESS! PE INITIATED---\n\n");

}


PatternEngineSingle::PatternEngineSingle(const int spWpartition, const int spZpartition, const int mode, std::string plottingpath) {
    /**
     * initializes as single PE for only one area.
     *
     * @param spWpartition number of super pixel bins in X or Phi (=W for width) direction
     * @param spZpartition number of super pixel bins in Z direction
     * @param mode type of sp mapping (0=uniform)
     * @param path to plot pdfs to.
     */

    printf("[ ---- INIT SINGLE PATTERN ENGINE ---- ]\n");
    this->initializeMembers(spWpartition, spZpartition, 0, mode);
    this->initializeSPbins();
    this->plottingpath = plottingpath;
    this->runspecs = ::getfileidtag(this->mode, this->wBinCount, this->zBinCount);
    printf("wBins=%d, zBins=%d\n", wBinCount, zBinCount);
    printf("\n");
}


PatternEngineSingle::PatternEngineSingle(const int spWpartition, const int spZpartition, const int mode, const int area, std::string plottingpath) {
    /**
     * initializes as one detector segment of compound PE.
     *
     * @param spWpartition number of super pixel bins in X or Phi (=W for width) direction
     * @param spZpartition number of super pixel bins in Z direction
     * @param mode type of sp mapping (0=uniform)
     * @param area detector area that is handled by this single pattern engine.
     * @param path to plot pdfs to.
     */

    std::cout << "[ ---- INIT SINGLE PATTERN ENGINE AS AREA " << area << " (" << areaDescript[area] << ") ---- ]" << std::endl;
    std::cout << "[ CONFIG :  Z Coord (min, max)=(" << DetectorZmin[area] << " , " <<  DetectorZmax[area] << ")" << " mode=" << mode << std::endl;
    this->centralDetectorZmin = DetectorZmin[area];
    this->centralDetectorZmax = DetectorZmax[area];
    this->initializeMembers(spWpartition, spZpartition, area, mode);
    this->initializeSPbins();
    this->plottingpath = plottingpath;
    this->plottingfile = "PEplots_" + getAreaTag() + ".pdf";
    this->runspecs = ::getfileidtag(this->mode, this->wBinCount, this->zBinCount);
}

//init members (called by constructors)
void PatternEngineSingle::initializeMembers(const int spWpartition, const int spZpartition, const int area, const int mode) {

    this->wBinCount = spWpartition;
    this->zBinCount = spZpartition;
    this->totalBinCount = spWpartition * spZpartition;
    this->centralDetectorLength = this->centralDetectorZmax - this->centralDetectorZmin;
    this->widthMin = -PI;
    this->widthMax = PI;
    this->widthDimLength = this->widthMax - this->widthMin;
    this->SPZsize = centralDetectorLength / zBinCount;
    this->SPWsize = widthDimLength / wBinCount;
    this->mode = mode;
    this->area = area;

    std::cout << "[ CONFIG : ZBins=" << zBinCount << " (size=" << SPZsize << ") ";
    std::cout << "WBins=" << wBinCount << " (size=" << SPWsize << ")" << " TOTAL BINS=" << totalBinCount << std::endl;
    std::cout << std::endl;

    for(int i=0; i<4; i++) {
//        this->spWeights[i].reserve(this->totalBinCount);
        for(int j=0; j<(this->totalBinCount); j++) this->spWeights[i].push_back(0);
//        printf("initialized spWeights[%d]\n", i);
    }
}


void PatternEngineSingle::initializeSPbins() {
    /**
     * init the superpixel bin boundary arrays for uniform equidistant mapping.
     */

    float SPZlength = centralDetectorLength / zBinCount;
    for(int i=0; i <= zBinCount; i++) {
        zBins.push_back(i * SPZlength + centralDetectorZmin);
    }

    float SPWlength = widthDimLength / wBinCount;
    for(int i = 0; i <= wBinCount; i++) {
        phiBins.push_back(i * SPWlength + widthMin);
    }

    //TODO implement xBins
}

////SOME UTILITY FUNCTIONS ------------------->
int PatternEngineSingle::getZSP(const float z) {
    return binSearch(this->zBins, z, 0, zBins.size() - 1);
//    return seqSearch(this->zBins, z);

}

int PatternEngineSingle::getPhiSP(const float phi) {
    return binSearch(this->phiBins, phi, 0, this->phiBins.size() - 1);
//    return seqSearch(this->phiBins, phi);
}

int PatternEngineSingle::getXSP(const float x) {
    return binSearch(this->xBins, x, 0, this->xBins.size() - 1);
}


int PatternEngineSingle::getLayer(const float r) {
    if (layerBoundaries[0] <= r && r <= layerBoundaries[1]) {
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

float PatternEngineSingle::getRadius(const float x, const float y) {
    return sqrt(pow(x,2) + pow(y,2));
}

float PatternEngineSingle::getHitPhi(const float x, const float y) {
    return atan2(y, x);
}

int PatternEngineSingle::getSPWcoord(const int sp2DID) {
    /**
     * computes 2D w-coordinate from continuous super pixel index.
     *
     * @param sp2DID index in sp grid (12 bit used, without layer+area info)
     * @return W (x or phi) coordinate in SP grid
     */

    return sp2DID % this->wBinCount;
}

int PatternEngineSingle::getSPZcoord(const int sp2DID) {
    /**
     * computes 2D z-coordinate from continuous super pixel index.
     *
     * @param sp2DID index in sp grid (12 bit used, without layer+area info)
     * @return Z coordinate in SP grid
     */

    return sp2DID / this->wBinCount;
}

template<typename T>
int PatternEngineSingle::binSearch(std::vector<T> targetVector, T value, int start, int end) {
    /**
     * Recursive bin search on sorted vector. Used to get super pixel index from bin boundary array.
     * This is necessary, when using a non-uniform bin mapping.
     * Computation time on average is better than a sequential search, as it is O(log(n)) instead O(n).
     *
     * @param targetVector vector to be searched
     * @param value that is searched
     * @param start index of area to be searched
     * @param end index of area to be searched
     * @return index of bin where `value` lies in.
     */

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

template<typename T>
int PatternEngineSingle::seqSearch(std::vector<T> targetVector, T value) {
    for(int i = 0; i < targetVector.size() - 1; i++) {
        if(targetVector[i] <= value && value < targetVector[i+1]){
            return i;
        }
    }
    if(value == targetVector[targetVector.size()-1]) {
        return targetVector.size()-2;
    } else {
        return -1;
    }
}

//// <-------------------- SOME UTILITY FUNCTIONS



unsigned int PatternEngineSingle::getSuperPixel(const float x, const float y, const float z) {
    /**
     * Computes the Super Pixel ID for a given hit on the detector area. In case the hit can't be assigned to a
     * Super Pixel, it returns 0.
     *
     * @param x,y,z coordinates of hit
     * @return 16bit SuperPixel ID. bit 1-4 encode area+layer, bit 5-16 encode super pixel index
     */

    float r = getRadius(x, y);
    int layer = getLayer(r);
    int phiSPIndex = getPhiSP(getHitPhi(x, y));
    int zSPIndex = getZSP(z);

    // data valid?
    assert(0 <= area && area < 3);
    assert(0 <= layer && layer < 4);
    assert(totalBinCount < 4096); // fatal because SID needs format of 0xFFF

    int sp2DID = computeIndex(zSPIndex, phiSPIndex, this->wBinCount);
    if(sp2DID < 0 || totalBinCount < sp2DID) sp2DID = 0;
    unsigned int SPID = computeSPID(layer, area, sp2DID);

    if(PRINTS) printf("Hit SuperPixel params\t x=%f, y=%f, z=%f, layer=%d zIndex=%d, phiIndex=%d, Sp2D=%d, SPID=%#X\n", x,y,z,layer, zSPIndex, phiSPIndex, sp2DID, SPID);

    //    printf("SP2D = %d, weights = %d \n", sp2DID, this->spWeights[0][sp2DID]);
    if(SPID != 0) this->spWeights[layer].at(sp2DID) = this->spWeights[layer].at(sp2DID) + 1;

    return SPID;
}

//// TESTS ( ONLY WORK ON CENTRAL DETECTOR SO FAR) FIXME ------------------->
void PatternEngineSingle::testbinSearch() {
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

void PatternEngineSingle::testCoordImpl() {
    int z = -100;
    int y = 10;
    int x = 10;

    float phi = getHitPhi(x, y);
    int phiSP = getPhiSP(phi);
    int zSP = getZSP(z);

    int SP2D = zSP * wBinCount + phiSP;
    assert(getLayer(getRadius(x,y)) == 0);
    assert(phiSP == getSPWcoord(SP2D));
    assert(zSP == getSPZcoord(SP2D));
}

void PatternEngineSingle::testSPID() {
    unsigned int sp1 = getSuperPixel(1, 0, -100);
    printf("layer=%d, area=%d, index=%d\n", getLayerFromSPID(sp1), getAreaFromSPID(sp1), getIndexFromSPID(sp1));

    sp1 = getSuperPixel(0, 0, -200);
    printf("layer=%d, area=%d, index=%d\n", getLayerFromSPID(sp1), getAreaFromSPID(sp1), getIndexFromSPID(sp1));
    sp1 = getSuperPixel(20, 0, 0);
    printf("layer=%d, area=%d, index=%d\n", getLayerFromSPID(sp1), getAreaFromSPID(sp1), getIndexFromSPID(sp1));
    sp1 = getSuperPixel(45, 0, 100);
    printf("layer=%d, area=%d, index=%d\n", getLayerFromSPID(sp1), getAreaFromSPID(sp1), getIndexFromSPID(sp1));
    sp1 = getSuperPixel(55, 0, 200);
    printf("layer=%d, area=%d, index=%d\n", getLayerFromSPID(sp1), getAreaFromSPID(sp1), getIndexFromSPID(sp1));
    sp1 = getSuperPixel(65, 0, 200);
    printf("layer=%d, area=%d, index=%d\n", getLayerFromSPID(sp1), getAreaFromSPID(sp1), getIndexFromSPID(sp1));

    float x=55.0;
    float y=3.0;
    float z=100;

    unsigned int sp2=getSuperPixel(x,y,z);
    int phiSP = getPhiSP(getHitPhi(x, y));
    int zSP = getZSP(z);
    int index = computeIndex(zSP, phiSP, this->wBinCount);
    int layer = getLayer(getRadius(x,y));
    int area=0;
    printf("-- custom evaluation: x=%f, y=%f, z=%f, phiSP=%d, zSP=%d, SPindex=%d, layer=%d, area=%d\n",x,y,z,phiSP, zSP, index, layer, area);
    unsigned int SPID = computeSPID(layer, area, index);
    int reslayer = getLayerFromSPID(SPID);
    int resarea = getAreaFromSPID(SPID);
    int resindex = getIndexFromSPID(SPID);
    printf("-- result for SPID dec=%u hex=%#X    ->   layer=%d, area=%d, index=%d\n", SPID, SPID, reslayer, resarea, resindex);

    assert(reslayer == layer);
    assert(resindex == index);
    assert(resarea == area);

//    getIndexFromSPID(0x1234F000); //assertion must fail
}
//// <------------  TESTS (ONLY WORK ON CENTRAL DETECTOR SO FAR)



//// PRINT GRAPHS AND PLOTS ----------->
void PatternEngineSingle::displayBinBoundaries() {
    /**
     * Makes plot that shows the super pixel mapping. Also writes this to the root file, where gROOT is pointing to.
     */

    auto *canvas = new TCanvas(("bin_boundaries" + this->getAreaTag()).c_str(), "bin_boundaries", 1200, 900);
    canvas->SetLeftMargin(0.15);
    canvas->SetGrid(1,1);
    canvas->SetTicks(1, 1);

    TH2F * grid = new TH2F(("h_bins" + this->getAreaTag()).c_str(), ("Superpixel bin grid boundaries "+ this->getAreaTag()).c_str() , zBinCount, centralDetectorZmin, centralDetectorZmax, wBinCount, widthMin, widthMax);
    grid->SetStats(false);
    grid->GetYaxis()->SetBinLabel(1, "-#pi");
    grid->GetYaxis()->SetBinLabel((int) (grid->GetNbinsY() / 2), "0");
    grid->GetYaxis()->SetBinLabel(grid->GetNbinsY(), "#pi");
    labelAxis(grid, "z [mm]", "#phi", 0.07);
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

    for(int i = 0; i< boundaries.size(); i++){
        boundaries[i]->Draw();
    }

    canvas->Write();
    canvas->Print((this->plottingpath + this->runspecs + pltf()).c_str(), "pdf");

}

void PatternEngineSingle::displayBinWeightDistributionLayer(int layer) {
    /**
     * Makes 2D hist, that shows the frequency if how often each super pixel was hit.
     * Also writes this to the root file, where gROOT is pointing to.
     */

    std::string layerno = get_string(layer);
    auto *canvas = new TCanvas(("bin_weights"+layerno).c_str(), ("bin_weights"+layerno).c_str(), 1200, 900);
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.15);
    canvas->SetGrid(1,1);
    canvas->SetTicks(1, 1);

    TH2F * h_binweights = new TH2F(("h_binweights" + layerno + getAreaTag()).c_str(), ("Superpixel bin weights layer " + layerno + " " + getAreaTag()).c_str(), zBinCount, 0, zBinCount, wBinCount, 0, wBinCount);
    h_binweights->SetStats(true);
    labelAxis(h_binweights, "z bins", "#phi bins");

    //Fill 2d hist with weights
    for(int i=0; i<totalBinCount; i++) {
        for(int j=0; j<spWeights[layer].at(i); j++) {
            h_binweights->Fill(getSPZcoord(i), getSPWcoord(i));
        }
    }

    // Norm
    h_binweights -> Scale(1.0 / h_binweights->Integral());

//    std::cout << "Bin z :" << h_binweights->GetXaxis()->GetBinWidth(1) << std::endl;

    h_binweights->Draw("colz");
    h_binweights->Write();

    canvas->Print((this->plottingpath + "PatternEngine_" + this->runspecs + pltf()).c_str(), "pdf");

//    TFile * tF = new TFile("output.root", "RECREATE");
//    h_binweights->Write();
//    TH1F * tP = (TH1F*)h_binweights->ProjectionY();
//    tP->Write();
//    tF->Close();

}

void PatternEngineSingle::displayBinWeightDistribution() {
    for(int i=(this->area == 0 ? 0 : 2); i<4; i++) {
        this->displayBinWeightDistributionLayer(i);
    }
}

std::string PatternEngineSingle::pltf() {
    //opens plotting pdf, if necessary.
    if(this->firstplot == 1) {
        printf("first file\n");
        this->firstplot = 0;
        return plottingfile + "(";
    } else {
        return plottingfile;
    }
}

std::string PatternEngineSingle::getRunSpecs() {
    return this->runspecs;
}

void PatternEngineSingle::closePlot() {
    auto *c_final = new TCanvas("c_final", "c_final");
    c_final->Print((this->plottingpath + this->runspecs + pltf() + ")").c_str(), "pdf");
}

std::string PatternEngineSingle::getAreaTag() {
    if(0 <= this->area && this->area <= 2) {
        return std::string(areaDescript[this->area]);
    } else {
        return "area" + get_string(this->area);
    }
}

//// <------------ PRINT GRAPHS AND PLOTS






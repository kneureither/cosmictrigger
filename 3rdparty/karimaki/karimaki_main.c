#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>
#include <string>
#include <TH1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <fstream>
#include <filesystem>
#include <iostream>

#include "karimaki.h"
#include "../../util/utility_functions.h"
#include "../../util/custom_types.h"

using std::endl;
using std::cout;

int main(int argc, char *argv[]) {

//    TApplication* app = new TApplication("app",&argc, argv);

    int karimaki_hit(KariFit&, int , double *, double *, double *, double *, double *, double *);

    int npoints = 4;
    double xp[4] = {85.8653, 71.3723,-36.469,-49.4697};
    double yp[4] = {-2.68301,-9.67875,-62.8534, -69.3326};
    double zp[4] = {-493.95, -496.46, -514.88, -517.09};
    double tres[4] = {1.0,1.0,1.0,1.0};
    double zres[4] = {1.0,1.0,1.0,1.0};
    double rres[4] = {1.0,1.0,1.0,1.0};

    KARITRACK karires;
    karimaki_hit(karires, npoints, xp, yp, zp, tres, zres, rres);

    printf("success!");
    printf("kari.phi=%4.2f", karires.phi);

}
//#include <string>
//#include <fstream>
//#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "../3rdparty/karimaki/karimaki_hit.c"
#include "../3rdparty/karimaki/karimaki.h"
#include "../util/utility_functions.h"
#include "../util/custom_types.h"

using std::cout;
using std::endl;

int main(int argc, char *argv[]) {

    int npoints = 4;
    double xp[85.8653, 71.3723,-36.469,-49.4697];
    double yp[-2.68301,-9.67875,-62.8534, -69.3326];
    double zp[-493.95, -496.46, -514.88, -517.09];
    double tres[1.0,1.0,1.0,1.0];
    double zres[1.0,1.0,1.0,1.0];
    double rres[1.0,1.0,1.0,1.0];

    KARITRACK karidata;

    karimaki_hit(npoints, xp, yp, zp, tres, zres, rres);

    cout << "success!" << endl;
    cout << "kari.phi=" << kari.phi << endl;

}
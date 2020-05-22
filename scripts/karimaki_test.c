#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "../3rdparty/karimaki/karimaki_hit.c"

int main(int argc, char *argv[]) {

    int npoints = 4;
    double xp[4] = {85.8653, 71.3723,-36.469,-49.4697};
    double yp[4] = {-2.68301,-9.67875,-62.8534, -69.3326};
    double zp[4] = {-493.95, -496.46, -514.88, -517.09};
    double tres[4] = {1.0,1.0,1.0,1.0};
    double zres[4] = {1.0,1.0,1.0,1.0};
    double rres[4] = {1.0,1.0,1.0,1.0};

    KARITRACK karidata;

    karimaki_hit(npoints, xp, yp, zp, tres, zres, rres);

    printf("success!");
    printf("kari.phi=%4.2f", kari.phi);

}
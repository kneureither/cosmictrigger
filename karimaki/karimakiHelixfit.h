#ifndef COSMICTRIGGER_KARIMAKIHELIXFIT_H
#define COSMICTRIGGER_KARIMAKIHELIXFIT_H

#define MAXLAYER 16


typedef struct KariFit
{
  float r3d;
  float rad;
  float dca;
  float phi;
  float tchi2n;
  float z0;
  float theta;
  float zchi2n;
} KARITRACK;

int karimakiHelixfit(KariFit &karires, int npoints, double xp[MAXLAYER], double yp[MAXLAYER], double zp[MAXLAYER],
        double phip[MAXLAYER], double thetas[MAXLAYER], double tres[MAXLAYER], double zres[MAXLAYER], double rres[MAXLAYER]);

template <typename T>
static int sign(const T val) {
    return (T(0) < val) - (val < T(0));
}

static void swapKariMomentum(KariFit &kari, float rec_r, int &count) {
    if(sign(kari.rad) != sign(rec_r)) {
        kari.rad = -kari.rad;
        kari.dca = -kari.dca;
        kari.r3d = -kari.r3d;
        count++;
    }
}

#endif //COSMICTRIGGER_KARIMAKIHELIXFIT_H

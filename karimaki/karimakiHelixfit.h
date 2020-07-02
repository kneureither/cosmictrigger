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

int karimakiHelixfit(KariFit &karires, int npoints, double xp[16], double yp[16], double zp[16], double phip[16], double thetas[16], double tres[16], double zres[16], double rres[16]);

template <typename T>
static int sign(const T val) {
    return (T(0) < val) - (val < T(0));
}

void swapKariMomentum(KariFit &kari, float rec_r, int &count) {
    if(sign(kari.rad) != sign(rec_r)) {
        kari.rad = -kari.rad;
        kari.dca = -kari.dca;
        kari.r3d = -kari.r3d;
        count++;
    }
}

#endif //COSMICTRIGGER_KARIMAKIHELIXFIT_H

// basics.h needed!!

// technical parameters
#define MAXSEED  (20*MAXTRACK)
#define MAXCAND  (60*MAXTRACK)
#define MAXFINAL (4*MAXTRACK)

#define MAXLAYER 12

#define NSECTORS 1000

#define NHITCUT  7  // minimum number of hits to be reconstrcuted


// reconstruction parameters
#define CHI2LOCAL 20.0

#define CHI2GLOBAL4 30.0
#define CHI2GLOBALSTEP 10.0
/*
#define CHI2GLOBAL4 1e10
#define CHI2GLOBALSTEP 0.0
*/
#define CHI2MAXLOCAL 25.0
#define R3DCHI2PASS 1e10

// acceptance cuts
#define ETAMIN 0.0
#define ETAMAX 1.0

//#define  R3DMIN3 0.3  // triplet cut
//#define  RMIN3 0.3
#define  R3DMIN3 0.6  // triplet cut
#define  RMIN3 0.6

// R=0.6 corresponds to pt~400 MeV for B=0.6 (2 Tesla)
//#define  R3DMIN 0.6   // track cut
//#define  RMIN 0.6
#define  R3DMIN 0.8   // track cut
#define  RMIN 0.8
#define  DZCUT  (0.5*ZVERTEX)
#define  PRESEL_COTTHETA 0.2  // cut to find hit pairs


#define RESRING0 (PIXRING0/sqrt(12.))
#define RESRING1 (PIXRING1/sqrt(12.))
#define RESRING2 (PIXRING2/sqrt(12.))

#define RESENDC1 (PIXENDC1/sqrt(12.))
#define RESENDC2 (PIXENDC2/sqrt(12.))

#define SEARCHWINDOWSIZE 5.0 // extrapolation search window size in units of scattering RMS


typedef struct Triplet
{
  int nlog;
  int hit1;
  int hit2;
  int hit3;
  double sms1;
  double sms2;
  double sms3;
  double chi2;
  double r3d;
  double sr3d;
  double th1,th2;
  double ph1,ph2;
  double sth2;
  double sph2;
  double sig_thetarot;
  double sig_phirot;
} TRIPLET;

#define TYPEBARR 0
#define TYPERING 1
#define TYPEENDC 2

typedef struct Track
{
  int nlog;
  int nfit;
  enum detgroups linked_group;
  int dettypes[MAXLAYER];            // 0=barrel, 1=ring, 2=endcap
  int layers[MAXLAYER];              // layer number intgrated over stattions!
  int hits[MAXLAYER];
  double sms[MAXLAYER];
  double chi2;
  double r3d;
  double sr3d;
  double th1,th2;
  double ph1,ph2;
  double sth2;
  double sph2;
  double sig_thetarot;
  double sig_phirot;
  //  double r3d_alt;
  //  double sr3d_alt;
} TRACK;

typedef struct RecTrack
{
  int nfit;             // number hits fitted
  int hits[MAXLAYER];
  int dettypes[MAXLAYER];
  int layers[MAXLAYER];
  int track_type;        // see list below

  double chi2;          // chi2 of global fit
  double r3d;
  double sr3d;
  double th1,th2;       // polar angle at beginning and end of track
  double ph1,ph2;       // azimuth at beginning and end of track
  double vtx_phi,dca,z0;


  int flag_correct;    // all linked hits belong to same true track
  int flag_ambig;      // ambiguity flag: 1=internal ambiguity, 2=external ambiguity
  int flag_noise;      // noise flag: 1=noise fitted, 2=all noise
  int gen_itr;         // linked track number (generated)
  double gen_r3d;       // radius of generated track number
  double gen_theta;     // theta of generated track number
  double gen_phi;       // phi of generated track number
  double gen_z0;        // z0 of generated track number (on beamline)
  double kari_r3d;
  double kari_dca;
  double kari_phi;
  double kari_tchi2n;
  double kari_z0;
  double kari_theta;
  double kari_zchi2n;
  double gbl_r3d;
  double gbl_dca;
  double gbl_phi;
  double gbl_z0;
  double gbl_theta;
  double gbl_chi2n;
} RECTRACK;



#define MIN_AMBIG_REJECT 1
typedef struct Track_List
{
int itr;
int nambig;
double target;
} TRACK_LIST;


// track parameter resolution
typedef struct Resol_Track
{
  int n;
  double ptmean;
  double etamean;
  double d0;
  double z0;
  double relpt;
  double relp;
  double theta;
  double phi;
} RESOLTRACK;


// global variable
#ifdef MAIN
// Book keeping variables why here?
int nrecfound,ncorrect,ngenfound,ntracc;
int nrecfound0,ncorrect0,ngenfound0,ntracc0;

int nrecfound_x,ncorrect_x,ngenfound_x,ntracc_x;
int nrecfound0_x,ncorrect0_x,ngenfound0_x,ntracc0_x;
int nrectrack;

RECTRACK track_rec[MAXFINAL];
TRACK_LIST track_ambig[MAXCAND];
TRACK_LIST track_work[MAXCAND];
TRACK_LIST track_list[MAXCAND];
RESOLTRACK track_rms_sum[3];    // 0=MS-fit; 1=Karimaki-fit; 2=GBL


#else
// Book keeping variables why here?
extern int nrecfound,ncorrect,ngenfound,ntracc;
extern int nrecfound0,ncorrect0,ngenfound0,ntracc0;

extern int nrecfound_x,ncorrect_x,ngenfound_x,ntracc_x;
extern int nrecfound0_x,ncorrect0_x,ngenfound0_x,ntracc0_x;
extern int nrectrack;

extern RECTRACK track_rec[MAXFINAL];
extern TRACK_LIST track_ambig[MAXCAND];
extern TRACK_LIST track_list[MAXCAND];
extern RESOLTRACK track_rms_sum[3];    // 0=MS-fit; 1=Karimaki-fit; 2=GBL
#endif

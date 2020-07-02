
#define SIMULATEPHOTONS
#define BREMSSTRAHLUNG
#define NUCLEARIA
#define SIM_SECONDARIES


typedef struct SIMVECTOR
{
  double x;
  double y;
  double z;
  double theta;
  double phi;
  double r3d;
  int code;
  int ivtx;         // vertex number: 1=primary, 2=gamma conv, 3=secondary (external input), 4=nuclear IA
  int idaughter;    // end of life:   0=live, 1=Brems,  2=gamma conv, 4=nuclear IA
  int imother;     // track number of mother
  int gennr;        // generated track number
  int recnr;        // reconstructed track number
  int trignr;        // reconstructed track-trigger number
  int trignr2;        // reconstructed track-trigger number
  int endc3nr;        // reconstructed track-trigger number
  int endc2nr;        // reconstructed track-trigger number
  double lrad;
  int nbarr;
  int nring;
  int nendc;
  int flag_reco;
  int flag_trigger;
  int flag_trigger2;
  int flag_endcap3;
  int flag_endcap2;
  double kap_tri;
  double kap_doub;
  double kap_endc3;
  double kap_endc2;
} SIMVEC;

#define NDETGROUPS 9
// oder should not be changed
enum detgroups {barr0, barr1, ring0, ring1, endc1, barr2, ring2, endc2, undef };

// global variable
#ifdef MAIN
int nlayers[NDETGROUPS];
int layer0[NDETGROUPS];

double bmsgroup[NDETGROUPS]; // scattering parameter
double xgroup[NDETGROUPS];   // rad-length parameter

//double tsmgroup[NDETGROUPS];  // transverse smearing pixel size (barrel)
//double zsmgroup[NDETGROUPS];  // longitudinal smearing pixel size (barrel)

double tsm_barr[NBARRTOT];  // transverse smearing pixel size (barrel)
double zsm_barr[NBARRTOT];  // longitudinal smearing pixel size (barrel)

double xsmgroup[NDETGROUPS];  // x-smearing pixel size (discs)
double ysmgroup[NDETGROUPS];  // y-smearing pixel size (discs)

SIMVEC track_arr[MAXTRACK];
SIMVEC phot_arr[MAXTRACK];
#else
extern int nlayers[NDETGROUPS];
extern int layer0[NDETGROUPS];

extern double bmsgroup[NDETGROUPS]; // scattering parameter
extern double xgroup[NDETGROUPS];   // rad-length parameter

//double tsmgroup[NDETGROUPS];  // transverse smearing pixel size (barrel)
//double zsmgroup[NDETGROUPS];  // longitudinal smearing pixel size (barrel)

extern double tsm_barr[NBARRTOT];  // transverse smearing pixel size (barrel)
extern double zsm_barr[NBARRTOT];  // longitudinal smearing pixel size (barrel)

extern double xsmgroup[NDETGROUPS];  // x-smearing pixel size (discs)
extern double ysmgroup[NDETGROUPS];  // y-smearing pixel size (discs)

extern SIMVEC track_arr[MAXTRACK];
extern SIMVEC phot_arr[MAXTRACK];
#endif

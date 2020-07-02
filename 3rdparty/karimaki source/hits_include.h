// basics.h has to be included!
// geometry.h has to be included!

/* foresee recurlers x2; x4 from clusters */
#define MAXHITS (10*MAXTRACK)

typedef struct HITTYPE_struct
{
  double x;
  double y;
  double z;
  int ivalid;   // can be used for reconstruction (no veto, not inefficient, etc. )
  /* 0 = inefficient
     1 = hit in cluster
     2 = single hit or cluster
  */
  int status;   // reconstruction status: number of track matches
  int itr;      // simulated track number+1!!! (negative means recurling), 0=noise
  // itr+1 should be fixed!!!

} HITTYPE;

#ifdef MAIN
// Global Memory
/* memory for hit management: outgoing track */
HITTYPE hit_mem_barr[MAXHITS][NBARRTOT]; // hit record barrel
HITTYPE hit_mem_ring[MAXHITS][NRINGTOT]; // hit record ring
HITTYPE hit_mem_endc[MAXHITS][NENDCTOT]; // hit record endcap

int nhit_barr[NBARRTOT]; // number of hits in barrel layers
int nhit_ring[NRINGTOT]; // number of hits in ring layers
int nhit_endc[NENDCTOT]; // number of hits in endcap layers


/* fake hit statistics */
int hit_barr_corr[NBARRTOT][2];
int hit_ring_corr[NRINGTOT][2];
int hit_endc_corr[NENDCTOT][2];
#else
extern HITTYPE hit_mem_barr[MAXHITS][NBARRTOT]; // hit record barrel
extern HITTYPE hit_mem_ring[MAXHITS][NRINGTOT]; // hit record ring
extern HITTYPE hit_mem_endc[MAXHITS][NENDCTOT]; // hit record endcap

extern int nhit_barr[NBARRTOT]; // number of hits in barrel layers
extern int nhit_ring[NRINGTOT]; // number of hits in ring layers
extern int nhit_endc[NENDCTOT]; // number of hits in endcap layers


/* fake hit statistics */
extern int hit_barr_corr[NBARRTOT][2];
extern int hit_ring_corr[NRINGTOT][2];
extern int hit_endc_corr[NENDCTOT][2];
#endif

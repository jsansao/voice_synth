/*
  Adapted in 02.10.95 to output jitter values
  Expanded 26 july 96 for more egg measures
  Modified in 01 may 2001: per-phase measures were exclude; spline options
			   excluded; measures based on a single cycle were
			   time-alined with the cycles;

  Modified in 22 May 2001: teste in "checar aqui" allows larger cycle-to-cycle
			  variations; (F0Min, F0Max, RatioMin, RatioMax inlcuded)

  Modified 10 abr 2006: upper F0 limit (for countertenors)

  Modified 24 mai 2010:
       crossing level for speed index via comand line
       crossing level for contact quotient (%) via comand line

		     2) PosPeak
		 ----------                              ---------
	       /            \                           /
	      /               \                        /
	 3)  / UpLevelx         \ DownLevelx          /
=====================================================o============CrossLevel
	   /                        \               /
     -----                            --------------
	1) NegPeak

 loop:
   1) Find NegPeak
   2) Find Next PosPeak
   3) Find single UpLevelx betwen NegPeak and PosPeak
      if(there is a previous single UpLevelx)
	4) Chaining
	5) measures
	   Measures are taken when a second up-going zero crossing is detected,
	   following other conditions in the chaining; the values are put into
	   the structure NewCycle but the output to the file is aligned with
	   the previous level crossing ("LastOutTime")

*/

#include "stdio.h"
//#include "io.h"
#include "fcntl.h"
//#include "conio.h"
#include "stdlib.h"
//#include "sys\stat.h"
#include "math.h"
#include "string.h"
#include "ctype.h"
#include "errno.h"
#include "time.h"

#define NOISE_THR	 15    	/* percentage of (dynamic) peakpct 
*/
#define CROSS_LEVEL	 0
#define MIN_PEAK	 400	/* absolute minumum allowed peak value */
#define TRUE		 1
#define FALSE		 !TRUE
#define ON		 1
#define NR_END		 1	/* used in spline */
#define FREE_ARG	 char*  /*  idem   */
#define CQLEVELX	0.25	/* passed in -t option */
#define INITCONTACT     0.10  /* default crossing level for Speed Index measurements */

#define PITCH		 1
#define JITTER		 2	/* based on up going zero crossings */
#define SHIMMER_PK	 3	/* based on peaks */
#define SHIMMER_RMS	 4	/* based on rms values */
#define SNR		 5       /* based on comb filters */
#define CLOSED_QUOT	 6
#define CLOSED_QUOT_AREA 7
#define SPEED_INDEX	 8
#define FREQUENCY	11

#define F0Max	1000.0
#define F0Min	40.0
#define RatioMax 5.0
#define RatioMin .20


/* FUNCTION PROTOTYPES */
void usage(void),
     msg(void),
     initialization(void),
     overlap(int),
     shift_structures(void),
     postprocessor(char *argv[]),
     shimmer(void),
     closed_quotient(void),
     speed_index(void),
     findphases(void),
     snr(void),
     statistics(int, int);

float *vector(long nl, long nh);

char* strrev (char *str);

int NegSignificantPoint(int),
    PosSignificantPoint(int);


/* HEADER OF .WAV FILE */
struct {
  char riff[4];               	/* "RIFF" */
  long int filesize;		/* size of file - 8 bytes */
  char wave[4];              	/* "WAVE" */
  char fmt[4];                	/* "fmt " */
  long int fmtsize;             /* 16, in general */
  short wFormatTag;		/* 1 for PCM */
  unsigned short nChannels;		/* 1 for mono, 2 for stereo */
  unsigned long nSamplesPerSec;	/* 44100, 22050, or 11025 */
  unsigned long nAvgBytesPerSec;	
  unsigned short nBlockAlign;		  
  unsigned short wBitsPerSample;	
  char data[4];			/* "data" */
  unsigned long datasize;		/* size of speech data */
} header;

/* FILE RELATED */
FILE *fdr = NULL,		/* data file header */
     *fpw = NULL,               /* header for output file (*.wav) */
     *tmp = NULL, 		/* 20 cycles of "spline" */
     *cum = NULL;		/* cummulative statistics */

char outfile[60], infile[60], cumfile[60];
signed short *x = NULL,				/* speech samples buffer */
    ni;					/* nb. of samples actually read */

/* TIME RELATED */
int milisec1,			/* integer corresponding to 1 ms */
    Tmax,			/* 1/Tmax = f0 min */
    Lframe,			/* length of frame (buffer) */
    Nframe = 0,			/* frame counter */
    N,				/* extra buffer around Lframe */
    S;				/* Sig point search interval */

/* GLOTTAL CYCLES ("humps") */
struct {
  int BreakFlag,
      PosPeakFrame,		/* frame number for Positive Peak; => out */
      PosPeak_i,		/* Pos peak position inside the frame */
      PosPeakValue,
      NegPeakFrame,

      NegPeak_i,
      NegPeakValue,
      UpLevelx_i,
      DownLevelx_i;

  float UpLevelxTime,		/* abs time (sec) of up going (interp.) */
	DownLevelxTime,         /* (not being used) */
	PosPeakTime,		/* abs time of pos peak (non-interpolated) */
	NegPeakTime,

	pitch,
	F0,
	jitter,
	shimmerRMS,
	shimmerPeak,
	snr,
	ClosedQuot,
	ClosedQuotArea,
	SpeedIndex,

	Closing,
	Closed,
	Opening,
	Open,

	rms;

  int InitClosing, EndClosing,
      InitOpening, EndOpening;

} LastCycle = {FALSE, 0, 0, 0, 0,
		   0, 0, 0, 0,
		 0.0, 0.0, 0.0, 0.0, 0.0,
		 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		 0.0, 0.0, 0.0, 0.0,
		 0.0, -1, -1, -1, -1},
  NewCycle =  {FALSE, 0, 0, 0, 0,
		   0, 0, 0, 0,
		 0.0, 0.0, 0.0, 0.0, 0.0,
		 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		 0.0, 0.0, 0.0, 0.0,
		 0.0, -1, -1, -1, -1};


/* UPGOING LEVEL CROSSING SEARCH */
struct {
  int NegPeakFlag,	/* Neg peak was found=> searched for Pos peak */
      Only1Levelx,	/* a single zero crossing was found */
      Levelx_count;	/* # of Levelx between pos and neg peak */
   float    CQCrossLevel;	/* Threshold pct for level crossing for Contact Quotient */
   int   m_points;		/* number of points for slope calculation */
   float PosDynLim,	/* Positive Dynamic Limiar for peak picking */
	 NegDynLim,
	 peakpct,	/* limiar for peak-picking */
	 InitContact;   /* assumed percentage for init of closing phase */
} search = {FALSE, FALSE, 0, CQLEVELX, 10, MIN_PEAK, -MIN_PEAK, NOISE_THR, INITCONTACT};
/* CHAINING OF GLOTTAL CYCLES */
struct {
  float LastOutTime,
	last_output,	/* last measure written to file */
	AvgPitch,	/* average pitch */
	AvgNegPeak,
	AvgPosPeak,
	UnvoicedTime,
	InitOfTracking,
	EndOfTracking;

  int PitchBreaks,
      f0flag;      	/* TRUE: previous output was a true pitch value
			   FALSE:   "       "     "  a chain break (-1 -1) */
} chain = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, FALSE};

/* CMD LINE ARGUMENTS */
struct {
  int input,
      output,
      CQCrossLevel,
      pitch,
      jitter,
      shimmer,
      snr,
      ClosedQuot,
      SpeedIndex,
      PerPhase,
      cummulative,
      InitContact;
} arg = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

struct STATISTICS {
  int n;
  double SumX, SumX2;
  float mean, StdDev;
};

int algorithm;		/* algorithm selected in the cmd line */

int main(int argc, char *argv[])

{
  int i, j;

  char str[60], str2[60], *p = NULL;
  float PreviousOut[3] = {0.0, 0.0, 0.0};

  float aux, aux2,	/* general use */
	ratio,		/* gereral use */
	pitch,		/* candidate value to NewCycle.pitch */
	deltat;         /* size of time correction after zx interpolation */

  /* CHECK PROGRAM CALL */
  if(argc<4)  usage();

  for(i=1; i<argc && *argv[i] == '-'; i++) {
    switch(argv[i][1]) {
      default:
	if (argc <= i+1) usage();

	switch(argv[i++][1]) {
	  case 'i':
	  case 'I':
	    arg.input = i;
	    strcpy(infile, argv[i]);
	  break;

	  case 'o':
	  case 'O':
	    arg.output = i;
	    strcpy(outfile, argv[i]);
	  break;

	  case 'c':
	  case 'C':
	    arg.cummulative = i;
	    strcpy(cumfile, argv[i]);
	  break;

	  case 'p':
	  case 'P':
	    arg.pitch = i;
	    j = atoi(argv[i]);
	    if(j==1) algorithm = PITCH; /* seconds */
	    else if(j==2) algorithm = FREQUENCY; /* hertz */
	    else usage();
	  break;

	  case 'j':
	  case 'J':
	    arg.jitter = i;
	    j = atoi(argv[i]);
	    if(j==1) algorithm = JITTER;
	    else usage();
	  break;

	  case 's':
	  case 'S':
	    arg.shimmer = i;
	    j = atoi(argv[i]);
	    if(j==1) algorithm = SHIMMER_PK;
	    else if(j==2) algorithm = SHIMMER_RMS;
	    else usage();
	  break;

	  case 'n':
	  case 'N':
	    arg.snr = i;
	    j = atoi(argv[i]);
	    if(j==1) algorithm = SNR;
	    else usage();
	  break;

	  case 'f':
	  case 'F':
	    printf("Option not available in this version\n");
	    exit(0);
	  break;

	  case 'q':
	  case 'Q':
	    arg.ClosedQuot = i;
	    j = atoi(argv[i]);
	    if(j==1) algorithm = CLOSED_QUOT;
	    else if(j==2) algorithm = CLOSED_QUOT_AREA;
	    else usage();
	  break;

	  case 'x':
	  case 'X':
	    arg.SpeedIndex = i;
	    j = atoi(argv[i]);
	    if(j==1) algorithm = SPEED_INDEX;
	    else usage();
	  break;

	  case 'T':
	  case 't':
	    arg.CQCrossLevel = i;
	    search.CQCrossLevel = atof(argv[arg.CQCrossLevel]);
	    if(search.CQCrossLevel < 10 ||
	       search.CQCrossLevel >  50) usage();
	    search.CQCrossLevel = search.CQCrossLevel/100.0;
	  break;

	  case 'M':
	  case 'm':
	    arg.InitContact = i;
	    search.InitContact = atof(argv[arg.InitContact]);

	    if(search.InitContact < 10 ||
	       search.InitContact > 50) usage();
	    search.InitContact = search.InitContact/100;
	    printf("InitContact = %f \n", search.InitContact);
	  break;

	  default:
	    usage();
	  break;
	}
	break;
    }
  }

  if((i != argc && *argv[i] != 'i') || arg.input == -1 || arg.output==-1) {
    usage();
  }

  printf("CQCrossLevel = %f\n", search.CQCrossLevel);

  if(arg.CQCrossLevel == -1 && arg.cummulative == -1 && argc != 9 ||
     arg.CQCrossLevel != -1 && arg.cummulative == -1 && argc != 11 ||
     arg.CQCrossLevel == -1 && arg.cummulative != -1 && argc != 11 ||
     arg.CQCrossLevel != -1 && arg.cummulative != -1 && argc != 13) usage();


  /* OPEN SPEECH FILE */
  if((fdr = fopen(argv[arg.input], "rb")) !=  NULL) {
    fread(&header, sizeof(header), 1, fdr);

    if(header.wFormatTag != 1) {
      printf(".wav file is not PCM");
      exit(0);
    }

    if(header.wBitsPerSample != 16) {
      printf(".wav file is not 16 bits per sample!");
    }
  }
  else {
    printf(".wav file not found\n");
    exit(0);
  }

  /* OPEN OUTPUT  FILE */
  if((fpw = fopen(argv[arg.output], "wt"))==NULL) {
    printf("Error while creating output file\n");
    exit(1);
  }

  if(arg.cummulative != -1) {
    if((cum = fopen(argv[arg.cummulative], "a+t"))==NULL) {
      printf("Error while creating cummulative file\n");
      exit(1);
    }

    /*-- take the 1st 6 characters of file name (PID) --*/
    strcpy(str2, infile); /* path...db0964la.wav */
    strrev(str2);  /* str2=> vaw.al4790bd...thap */
    strcpy(str, "12345678");

    p = strchr(str2, '.');
    if(p != NULL) strncpy(str, p+1, 8);
    else strncpy(str, str2, 8);

    strrev(str); /* str=> db0964la */
    strcpy(str2, "123456");
    strncpy(str2, str, 6); /* str2 = db0964 */
    fprintf(cum, "%s ", str2);
    /*printf("Cummulative file: %s\n", cumfile);*/
  }

  /* INITIALIZE VARIABLES, BUFFERS, ETC */
  fprintf(fpw,"0.0 ");
  gcvt((float) 8.0*header.datasize/
		   (header.wBitsPerSample*header.nSamplesPerSec), 12, str);
  fprintf(fpw, str);
  fprintf(fpw, "\n");

  /* teste aqui */
  fprintf(fpw, "-1 -1\n");
  chain.last_output = -1;



  initialization();


  /* PUT MESSAGES ON THE SCREEN */
  msg();

  /*
  printf("Noise Threshold: %5.2f%%\n", search.peakpct);
  */

  /*---------------------- MAIN LOOP ---------------------------------*/



  while((ni = (fread(x + 2*N, sizeof(int), Lframe, fdr))) != NULL) {
	       /* fill buffer for 2N <= x < Lframe+N */

   /*
   Buffers:                                             (ni)
   i=  0                      ...                 2N + Lframe - 1
       |                                           |
      [nnnnnnnnnnnn...nnnnnnnnnnnnnnnnn] = x[] (int)
       |-N-|-N-|

       |-N-|-------------Lframe------------|-N-|
	20                 100               20  [ms]

      o=old samples (overlap), n=new samples,
      Lframe >= ni (ni < Lframe only when eof is reached)
   */


    for(i = N; i < N+ni; i++) {
      /* 1) find NegPeak */
      if( (x[i]<search.NegDynLim) || (x[i]>search.PosDynLim) ) {
	if(NegSignificantPoint(i)) {
	  NewCycle.NegPeakFrame = Nframe;
	  NewCycle.NegPeakValue = x[i];
	  NewCycle.NegPeak_i = i;

	  search.NegPeakFlag = TRUE;
	}
	else if(search.NegPeakFlag) {
	  /* 2) Find PosPeak */
	  if(PosSignificantPoint(i)) {
	    NewCycle.PosPeakFrame = Nframe;
	    NewCycle.PosPeakValue = x[i];
	    NewCycle.PosPeak_i = i;
	    NewCycle.PosPeakTime = (i +
			((float) Nframe)*Lframe)/header.nSamplesPerSec;

	    /* test for min f0 (Tmax) */

	    /*
	    if(NewCycle.PosPeakValue < 0) {
	      printf("NewPosPeak < 0");
	      getch();
	    }
	    */

	    aux = NewCycle.PosPeak_i - NewCycle.NegPeak_i;
	    if(aux < Lframe/2) {
	      search.Only1Levelx = FALSE;
	      search.Levelx_count = 0;

	      /* 3) Find single up going Level crossing between peaks */
	      for(j = NewCycle.NegPeak_i; j<i; j++) {
		if((x[j]<=search.CQCrossLevel) &&
		   (x[j+1]>=search.CQCrossLevel)) {
		  NewCycle.UpLevelx_i = j;
		  search.Only1Levelx = TRUE;

		  if(++search.Levelx_count > 1) {
		    /* more than one up zero crossing was found */
		    search.Only1Levelx = FALSE;
		    break;
		  }
		}
	      }
	      if(search.Only1Levelx) { /* only one zero crossing was found */
		/* calculate increment in seconds (linear interpolation) */
		deltat = (float) -x[NewCycle.UpLevelx_i]/
			((x[NewCycle.UpLevelx_i+1] - x[NewCycle.UpLevelx_i])*
						    header.nSamplesPerSec);

		/* calculate absolute time of zero crossing in seconds */
		NewCycle.UpLevelxTime = (NewCycle.UpLevelx_i  +
		    ((float) Nframe)*Lframe)/header.nSamplesPerSec + deltat;

		/* 4) chaining ... */
		if(LastCycle.UpLevelxTime != 0) {
		  pitch = NewCycle.UpLevelxTime - LastCycle.UpLevelxTime;

		  /* compare with positive peak distance */
		  if( (NewCycle.PosPeakTime - LastCycle.PosPeakTime)==0)
		    ratio = 1;
		  else
		    ratio = pitch/(NewCycle.PosPeakTime - LastCycle.PosPeakTime);

		  /* checar aqui */
		  /*
		  if( (pitch > 1.0/500.0) && (pitch < 1.0/50.0) &&
		      (ratio>0.8 && ratio < 1.2) ) {
		  */
		  if( (pitch > 1.0/F0Max) && (pitch < 1.0/F0Min) &&
		      (ratio>RatioMin && ratio < RatioMax) ) {

		    NewCycle.pitch = pitch;

		    /* ... chaining */
		    if(chain.LastOutTime != 0.0) { /* nao e' inicio de cadeia */
		      if( (NewCycle.UpLevelxTime - chain.LastOutTime)/
						 LastCycle.pitch > 1.8) {
			/* trap against period doubling */
			if(chain.last_output != -1) {
			  fprintf(fpw, "-1 -1\n");
			  chain.PitchBreaks++;

			  PreviousOut[0] = PreviousOut[1];
			  PreviousOut[1] = PreviousOut[2];
			  PreviousOut[2] = LastCycle.UpLevelxTime;
			}

			chain.f0flag = FALSE;
			chain.LastOutTime = 0;
			chain.last_output = -1;
		      }
		      else {

			/*:::::::::::::::::::::::::::::::::*/
			NewCycle.F0 = 1.0/pitch;
			NewCycle.jitter =
			  200.0*fabs(pitch - LastCycle.pitch)/
				    (pitch + LastCycle.pitch);

			findphases();
			shimmer();
			snr();
			speed_index();
			closed_quotient(); /* time and area */

			statistics(0, arg.cummulative);
			/*::::::::::::::::::::::::::::::::::*/


			/* OUTPUT TO FILE (Time series) */
			switch(algorithm) {
			  case PITCH:
			    aux = 1000.0*NewCycle.pitch;
			  break;

			  case FREQUENCY:
			    aux = NewCycle.F0;
			  break;

			  case JITTER:
			    aux = NewCycle.jitter;
			  break;

			  case SHIMMER_PK:
			    aux = NewCycle.shimmerPeak;
			  break;

			  case SHIMMER_RMS:
			    aux = NewCycle.shimmerRMS;
			  break;

			  case SNR:
			    aux = NewCycle.snr;
			  break;

			  case CLOSED_QUOT:
			    aux = NewCycle.ClosedQuot;
			  break;

			  case CLOSED_QUOT_AREA:
			    aux = NewCycle.ClosedQuotArea;
			  break;

			  case SPEED_INDEX:
			    aux = NewCycle.SpeedIndex;
			  break;

			  default:
			  break;
			}

			/* teste aqui 01/05/2001 ***** */
			if(algorithm == JITTER || algorithm == SHIMMER_PK
			   || algorithm == SHIMMER_RMS || algorithm == SNR) {
			  aux2 = chain.LastOutTime;
			}
			else {
			  aux2 = chain.LastOutTime - LastCycle.pitch;
			}
			/******************************/

			fprintf(fpw, "%f %5.2f\n", aux2, aux);

			if(chain.InitOfTracking == 0.0) {
			 /* used in the estimation of UnvoicedTime */
			  chain.InitOfTracking = NewCycle.UpLevelxTime;
			}

			chain.last_output = NewCycle.jitter;
			chain.LastOutTime = NewCycle.UpLevelxTime;
			chain.f0flag = TRUE;

                        PreviousOut[0] = PreviousOut[1];
			PreviousOut[1] = PreviousOut[2];
			PreviousOut[2] = LastCycle.UpLevelxTime;
			chain.EndOfTracking = LastCycle.UpLevelxTime;

			if(PreviousOut[1]==-1 &&
			   PreviousOut[2] > PreviousOut[0] &&
			   PreviousOut[0] != 0.0) {
			  /* add to unvoiced intervals */
			  chain.UnvoicedTime += PreviousOut[2] -
						PreviousOut[0];
			}
		      }
		    }
		    else {
		      chain.LastOutTime = NewCycle.UpLevelxTime;
		      chain.f0flag = TRUE;
		    }

		    if(chain.AvgPitch==0.0) {
		      /* new chain is starting */
		      chain.AvgPitch = pitch;
		      chain.AvgPosPeak = NewCycle.PosPeakValue;
		      chain.AvgNegPeak = NewCycle.NegPeakValue;
		    }
		    else {
		      chain.AvgPitch = 0.8*chain.AvgPitch + 0.2*pitch;
		      chain.AvgPosPeak = 0.8*chain.AvgPosPeak +
					 0.2*NewCycle.PosPeakValue;
		      chain.AvgNegPeak = 0.8*chain.AvgNegPeak +
					 0.2*NewCycle.NegPeakValue;

		      search.PosDynLim = chain.AvgPosPeak*search.peakpct/100.0;
		      if(search.PosDynLim < MIN_PEAK) {
			search.PosDynLim = MIN_PEAK;
		      }
		      search.NegDynLim = chain.AvgNegPeak*search.peakpct/100.0;
		      if(search.NegDynLim > -MIN_PEAK) {
			search.NegDynLim = -MIN_PEAK;
		      }
		    }
		  }
		  else if(chain.f0flag) {
		    if(chain.last_output != -1) {
		      fprintf(fpw, "-1 -1\n");
		      chain.PitchBreaks ++;

		      PreviousOut[0] = PreviousOut[1];
		      PreviousOut[1] = PreviousOut[2];
		      PreviousOut[2] = -1;

		      if(chain.last_output == 0.0) {
			chain.PitchBreaks --; /* discount for first -1 -1 */
		      }
		    }
		    chain.f0flag = FALSE;
		    chain.LastOutTime = 0.0;
		    chain.AvgPitch=0.0;
		    chain.last_output = -1;
		  }
		}
	      }
	    }
	    else { /* aux >= Lframe+2N */
	       /* break in the chain */
	      if(LastCycle.BreakFlag) {
		NewCycle.BreakFlag = FALSE;
		if(chain.f0flag) {
                  if(chain.last_output != -1) {
		    fprintf(fpw, "-1 -1\n");
		    chain.PitchBreaks++;

		    PreviousOut[0] = PreviousOut[1];
		    PreviousOut[1] = PreviousOut[2];
		    PreviousOut[2] = -1.0;
		  }

		  chain.f0flag = FALSE;
		  chain.AvgPitch = chain.AvgPosPeak = chain.AvgNegPeak = 0.0;
		  chain.LastOutTime = 0.0;
		  chain.last_output = -1;
		}
	      }
	    }
	    shift_structures();
	  }
	}
      }
    }
    overlap(ni);
    Nframe++;
  }

  statistics(1, arg.cummulative);

  free(x);
  fclose(fdr);
  fclose(fpw);


  if(arg.cummulative) fclose(cum);

/*  postprocessor(argv);*/

  exit(0);
}

/*---------------------------- FUNCTIONS ---------------------------*/

void usage()
{
  printf("\negg2 -i inputfile -o seriesfile [-c cummfile] [-options {default}]\n");
  printf("Last modified 24 mai 2010\n");
  printf("   Lower crossing level for speed index via comand line\n");
  printf("   Crossing level for contact quotient (time) via comand line\n");

  printf("-i   .wav file, mono\n");
  printf("-o   ascii file (time series)\n");
  printf("-c   cummulative statistical file\n\n");

  printf("Time series output:\n");
  printf("-p # (fundamental freq/period, 1=ms, 2=Hz)\n");
  printf("-j # (jitter [PF1%], 1=interpolated zero crossings)\n");
  printf("-s # (shimmer [PF1%], 1=peaks, 2=rms \n");
  printf("-n # (signal-to-noise ratio [HNR dB], 1=Comb\n");
  printf("-q # (Closed Quotient [%] 1 = based on time only, 2 = based on areas\n");
  printf("-x # (Speed Index, 1=(Opening-Closing)/(Opening+closing)\n\n");

  printf("F0 detection (common to all measurements):\n");
  printf("-t # (Crossing Level for contact quotient [between 10-50]) {25}\n");
  printf("-m # (Lower Crossing level for Speed Index (pct), [between 10-50 pct]) {10}\n\n");

/*
  printf("Postprocessor settings:\n");
  printf("-m   (median non linear filter)\n");
  printf("-l   (linear smoothing)\n");
  printf("-g x (maximum filled gap [milliseconds]) {x=25}\n");
*/
  exit(0);
}

void initialization()
{
  milisec1  = (int) (header.nSamplesPerSec * 0.001/2.0)*2;

  Lframe = 50*milisec1;

  S = 0.5*milisec1; /* used in the search of sig peaks */
  N = 20*milisec1;

  Tmax = 20*milisec1;

  /* CREATE INPUT BUFFER */
  if((x = malloc(sizeof(signed short) * 2*(N + Lframe))) == NULL) {
    printf("out of memory in call to malloc(x).\n");
    exit(1);
    }

  /* INITIALIZE BUFFERS for 0 <= x < N */
  fread(x, sizeof(signed short), 2*N, fdr);
}


void overlap(int ni)
/* ni <= Lframe */
{
  register int i;

  for(i=0; i < 2*N; i++) {
    /*
    x[i] = x[Lframe+i];
    */
    x[i] = x[ni + i];

  }

  LastCycle.PosPeak_i -= ni;    NewCycle.PosPeak_i  -= ni;
  LastCycle.NegPeak_i -= ni;    NewCycle.NegPeak_i  -= ni;
  LastCycle.UpLevelx_i -= ni;   NewCycle.UpLevelx_i  -= ni;
  LastCycle.DownLevelx_i -= ni; NewCycle.DownLevelx_i  -= ni;
  LastCycle.InitClosing -= ni;  NewCycle.InitClosing -= ni;
  LastCycle.EndClosing -= ni;   NewCycle.EndClosing -= ni;
  LastCycle.InitOpening -= ni;  NewCycle.InitOpening -= ni;
  LastCycle.EndOpening -= ni;   NewCycle.EndOpening -= ni;
}


void msg(void)
{
  printf(" \n (c) Maurilio N. Vieira, 26 July 96. \n");
  printf(" EGG measures\n");

  printf("Time series: ");
  if(algorithm == JITTER) printf("jitter\n");
  else if(algorithm == PITCH) printf("pitch (ms)\n");
  else if(algorithm == FREQUENCY) printf("F0 (Hz)\n");
  else if(algorithm == SHIMMER_PK) printf("shimmer (peaks)\n");
  else if(algorithm == SHIMMER_RMS) printf("shimmer (rms)\n");
  else if(algorithm == CLOSED_QUOT) printf("closed quot (time)\n");
  else if(algorithm == CLOSED_QUOT_AREA) printf("closed quot (area)\n");
  else if(algorithm == SNR) printf("signal-to-noise ratio (comb)\n");
  else if(algorithm == SPEED_INDEX) printf("speed index\n");
}


int NegSignificantPoint(int i)
{
  int k,
      type = -1;

      /* aqui */
  for(k=1; k <= S; k++) {     /* value of S modified in 10 abr 2006 */
/*
    if(!((x[i]<0)&&(x[i]<=x[i+k])&&(x[i]<x[i-k]))) {
*/
    if(!((x[i]<search.CQCrossLevel)&&(x[i]<=x[i+k])&&(x[i]<x[i-k]))) {

      type=1;
      break;
    }
  }

  if(type==-1) return TRUE;
  else return FALSE;
}


int PosSignificantPoint(int i)
{
  int k,
      type = 1;

  for(k=1; k <= S; k++) {
    if(!((x[i]>0)&&(x[i]>=x[i+k])&&(x[i]>x[i-k]))) {
      type=-1;
      break;
    }
  }

  if(type==1) return TRUE;
  else return FALSE;
}

void shift_structures(void)
{
  /* everything is being shift in case a difference is needed */
  LastCycle.BreakFlag     = NewCycle.BreakFlag;
  LastCycle.PosPeakFrame  = NewCycle.PosPeakFrame;
  LastCycle.PosPeak_i     = NewCycle.PosPeak_i;
  LastCycle.PosPeakValue  = NewCycle.PosPeakValue;
  LastCycle.NegPeakFrame  = NewCycle.NegPeakFrame;
  LastCycle.NegPeak_i     = NewCycle.NegPeak_i;
  LastCycle.NegPeakValue  = NewCycle.NegPeakValue;
  LastCycle.UpLevelx_i    = NewCycle.UpLevelx_i;
  LastCycle.DownLevelx_i  = NewCycle.DownLevelx_i;

  LastCycle.UpLevelxTime   = NewCycle.UpLevelxTime;
  LastCycle.DownLevelxTime = NewCycle.DownLevelxTime;
  LastCycle.PosPeakTime   = NewCycle.PosPeakTime;
  LastCycle.NegPeakTime   = NewCycle.NegPeakTime;

  LastCycle.pitch         = NewCycle.pitch;
  LastCycle.F0		  = NewCycle.F0;
  LastCycle.jitter	  = NewCycle.jitter;
  LastCycle.shimmerRMS	  = NewCycle.shimmerRMS;
  LastCycle.shimmerPeak   = LastCycle.shimmerPeak;
  LastCycle.snr    	  = NewCycle.snr;
  LastCycle.ClosedQuot    = NewCycle.ClosedQuot;
  LastCycle.ClosedQuotArea = NewCycle.ClosedQuotArea;

  LastCycle.SpeedIndex    = NewCycle.SpeedIndex;
  LastCycle.Closing       = NewCycle.Closing;
  LastCycle.Closed	  = NewCycle.Closed;
  LastCycle.Opening       = NewCycle.Opening;
  LastCycle.Open	  = NewCycle.Open;
  LastCycle.rms		  = NewCycle.rms;
  LastCycle.SpeedIndex    = NewCycle.SpeedIndex;
  LastCycle.InitClosing   = NewCycle.InitClosing;
  LastCycle.EndClosing    = NewCycle.EndClosing;
  LastCycle.InitOpening   = NewCycle.InitOpening;
  LastCycle.EndOpening    = NewCycle.EndOpening;
}


void shimmer(void)
{
  int j;

  NewCycle.rms = 0.0;
  for(j = NewCycle.UpLevelx_i; j>LastCycle.UpLevelx_i; j--) {
    NewCycle.rms += (float) x[j]*x[j];
  }
  NewCycle.rms = sqrt(NewCycle.rms/NewCycle.pitch);

  if(NewCycle.rms <= 0.0) {
    printf("rms<0.0\n");
    NewCycle.shimmerRMS = -1;
  }
  else if (LastCycle.rms != -1) {
    NewCycle.shimmerRMS = 200.0*fabs(NewCycle.rms - LastCycle.rms)/
			       (NewCycle.rms + LastCycle.rms);
    if(NewCycle.shimmerRMS < 0) {
      printf("shimmerRMS < 0\n");
      NewCycle.shimmerRMS = -1;
    }
  }

  if(NewCycle.PosPeakValue<0 || LastCycle.PosPeakValue<0) {
    printf("Pospeak < 0\n");
    NewCycle.shimmerPeak = -1;
  }
  else {
    NewCycle.shimmerPeak = 200.0*
	   fabs((float) NewCycle.PosPeakValue - LastCycle.PosPeakValue)/
	       ((float) NewCycle.PosPeakValue + LastCycle.PosPeakValue);
    if(NewCycle.shimmerPeak < 0) {
      printf("shimmerPeak < 0\n");
      NewCycle.shimmerPeak = -1;
    }
  }
}

void closed_quotient(void)
/*
	    PosPeak
	    -+ -
	  /      \
	 /         \             /
  ======/===========|===========/=========  CrossLevel
  -----+ta----------+tb--------+tc--------  CqLeveLx
   -+-               - - - - -
    NegPeak

    ClosedQuot(time) = (tb - ta)/(tc-ta)
*/

{
  int j, ta, tb, tc;
  float aux, ClosedPhaseArea, OpenPhaseArea, CqLeveLx;

  if(LastCycle.PosPeak_i == -1 || LastCycle.NegPeak_i == -1) {
    NewCycle.ClosedQuot = NewCycle.ClosedQuotArea = -1;
    return;
  }

  /*  pega o cruzamento em 25% a partir do pico negativo
  CqLeveLx =   (float) LastCycle.PosPeakValue -
	   ( (float) LastCycle.PosPeakValue - LastCycle.NegPeakValue)*CQLEVELX;
  */

  /* pega o cruzamento a partir do percentual passado em -t da cmd line */
  aux = ( (float) LastCycle.PosPeakValue - LastCycle.NegPeakValue)*CQLEVELX;
  CqLeveLx = (float) LastCycle.NegPeakValue + aux;

  ClosedPhaseArea = OpenPhaseArea = 0.0;

  /* find ta */
  for(j = LastCycle.PosPeak_i; j>LastCycle.NegPeak_i; j--) {
    if(x[j] >= CqLeveLx) {
      ClosedPhaseArea += fabs(x[j]);
    }
    else {
      ta = j;
      break;
    }
  }

  /* find tb; tc = ta + pitch */
  for(j = LastCycle.PosPeak_i + 1; j < NewCycle.NegPeak_i; j++) {
    if(x[j] >= CqLeveLx) {
      ClosedPhaseArea += fabs(x[j]);
    }
    else {
      tb = j;
      tc = ta + (int) ceil(NewCycle.pitch*header.nSamplesPerSec);
		/* avoid humps in the open phase */
      break;
    }
  }

  for(j = tb+1; j < tc; j++) {
    OpenPhaseArea += fabs(x[j]);
  }

  NewCycle.ClosedQuot = 100.0*((float) tb-ta)/((float) tc - ta);
  if(NewCycle.ClosedQuot < 0) {
    printf("CQ(time) < 0\n");
    NewCycle.ClosedQuot = -1;
  }

  NewCycle.ClosedQuotArea = 100.0*ClosedPhaseArea/
			       (ClosedPhaseArea+OpenPhaseArea);
  if(NewCycle.ClosedQuotArea < 0) {
    printf("CQ(area) < 0\n");
    NewCycle.ClosedQuotArea = -1;
  }
}

void speed_index(void)
{
  float num, den;

  num = ( (float) LastCycle.EndClosing - LastCycle.InitClosing) -
	( (float) LastCycle.EndOpening - LastCycle.InitOpening);
  den = ( (float) LastCycle.EndClosing - LastCycle.InitClosing) +
	( (float) LastCycle.EndOpening - LastCycle.InitOpening);

  if(LastCycle.InitClosing != -1 && LastCycle.EndClosing != -1 &&
    LastCycle.InitOpening != -1 && LastCycle.EndOpening != -1 &&
    den != 0.0 && num != 0.0) {
    NewCycle.SpeedIndex = num/den;
  }
  else NewCycle.SpeedIndex = 2.0; /* invalid number */
}

void findphases(void)
{
  int k;
  float pctLow, pctHigh, aux;

  aux = search.InitContact*( (float) NewCycle.PosPeakValue - NewCycle.NegPeakValue);
  pctLow = (float)  NewCycle.NegPeakValue + aux;

  /* colocar 0.10 como parametro de linha de cmd */
  aux = 0.10*( (float) NewCycle.PosPeakValue - NewCycle.NegPeakValue);
  pctHigh = (float) NewCycle.PosPeakValue - aux;


  NewCycle.InitClosing = NewCycle.EndClosing =
  NewCycle.InitOpening = NewCycle.EndOpening = -1;

  for(k=NewCycle.PosPeak_i; k>LastCycle.PosPeak_i; k--) {
    if( (x[k] >= pctHigh) && (x[k-1] <= pctHigh) ) {
      NewCycle.EndClosing = k;
    }
    if( (x[k] > pctLow) && (x[k-1] <= pctLow) ) {
      NewCycle.InitClosing = k;
      break;
    }
  }

  aux = ceil(NewCycle.pitch*header.nSamplesPerSec) +
	(float) NewCycle.PosPeak_i;
  for(k=NewCycle.PosPeak_i; k < aux; k++) {
    if(x[k] >= pctHigh && x[k+1] <= pctHigh) {
      NewCycle.InitOpening = k;
    }
    if(x[k] >= pctLow && x[k+1] <= pctLow) {
      NewCycle.EndOpening = k;
      break;
    }
  }
 /* checar aqui 24 mai 2010; usado para CQ!! */
  /* take first down going crossing between levels of SearchCrossLevel */
  NewCycle.DownLevelx_i = -1;
  if(NewCycle.InitOpening != -1 && NewCycle.EndOpening != -1) {
    for(k=NewCycle.InitOpening; k<NewCycle.EndOpening; k++) {
      if( NewCycle.DownLevelx_i == -1 &&
	  x[k] >= search.CQCrossLevel  && x[k+1] <= search.CQCrossLevel) {
	NewCycle.DownLevelx_i = k;
	break;
      }
    }
  }
}


void snr(void)
{
  int i, j, T, aux;
  float noise, harmonic;

  noise = harmonic = 0.0;
  T = (int) ceil(NewCycle.pitch*header.nSamplesPerSec);

  aux = LastCycle.UpLevelx_i + T;
  for(i=LastCycle.UpLevelx_i, j = NewCycle.UpLevelx_i; i < aux; i++, j++) {
    harmonic += ((float) x[i] + x[j]) *( (float) x[i] + x[j]);
    noise    += ((float) x[i] - x[j]) *( (float) x[i] - x[j]);
  }

  if(noise != 0.0) {
    NewCycle.snr = 10*log10(harmonic/noise);
  }
  else {
    NewCycle.snr = 1e20; /* large value (= 1/0 ) */
  }
}


void statistics(int code, int cummulative)
{
  char str[60];
  FILE *stat = NULL;
  float dur, UnvoicedPct, aux;

  static struct STATISTICS
	      F0 = {0, 0.0, 0.0, 0.0, 0.0},
	      jitter = {0, 0.0, 0.0, 0.0, 0.0},
	      shimmerRMS = {0, 0.0, 0.0, 0.0, 0.0},
	      shimmerPeak = {0, 0.0, 0.0, 0.0, 0.0},
	      snr = {0, 0.0, 0.0, 0.0, 0.0},
	      ClosedQuot = {0, 0.0, 0.0, 0.0, 0.0},
	      ClosedQuotArea = {0, 0.0, 0.0, 0.0, 0.0},
	      SpeedIndex = {0, 0.0, 0.0, 0.0, 0.0},
	      Closing = {0, 0.0, 0.0, 0.0, 0.0},
	      Closed = {0, 0.0, 0.0, 0.0, 0.0},
	      Opening = {0, 0.0, 0.0, 0.0, 0.0},
	      Open = {0, 0.0, 0.0, 0.0, 0.0},
	      ClosingPct = {0, 0.0, 0.0, 0.0, 0.0},
	      OpeningPct = {0, 0.0, 0.0, 0.0, 0.0};

  if(code == 1) { /* end statistics */
    if(F0.n) {
      F0.mean = F0.SumX/F0.n;
      F0.StdDev = sqrt(
		   (F0.SumX2 - F0.n*F0.mean*F0.mean)/F0.n );
    }
    if(jitter.n) {
      jitter.mean = jitter.SumX/jitter.n;
      jitter.StdDev = sqrt(
	(jitter.SumX2 - jitter.n*jitter.mean*jitter.mean)/jitter.n);
    }

    if(shimmerRMS.n) {
      shimmerRMS.mean = shimmerRMS.SumX/shimmerRMS.n;
      shimmerRMS.StdDev = sqrt(
	(shimmerRMS.SumX2 - shimmerRMS.n*shimmerRMS.mean*shimmerRMS.mean)/
	shimmerRMS.n);
    }

    if(shimmerPeak.n) {
      shimmerPeak.mean = shimmerPeak.SumX/shimmerPeak.n;
      shimmerPeak.StdDev = sqrt(
	(shimmerPeak.SumX2 - shimmerPeak.n*shimmerPeak.mean*shimmerPeak.mean)/
	shimmerPeak.n);
    }

    if(snr.n) {
      snr.mean = (float) (snr.SumX/snr.n);
      snr.StdDev = (float) sqrt((snr.SumX2 - snr.n*snr.mean*snr.mean)/snr.n);
    }

    if(ClosedQuot.n) {
      ClosedQuot.mean = (float) (ClosedQuot.SumX/ClosedQuot.n);
      ClosedQuot.StdDev = (float) sqrt(
	(ClosedQuot.SumX2 - ClosedQuot.n*ClosedQuot.mean*ClosedQuot.mean)/
	  ClosedQuot.n);
    }

    if(ClosedQuotArea.n) {
      ClosedQuotArea.mean = (float) (ClosedQuotArea.SumX/ClosedQuotArea.n);
      ClosedQuotArea.StdDev = (float) sqrt(
	(ClosedQuotArea.SumX2 - ClosedQuotArea.n*ClosedQuotArea.mean*ClosedQuotArea.mean)/
	  ClosedQuotArea.n);
    }

    if(SpeedIndex.n) {
      SpeedIndex.mean = (float) (SpeedIndex.SumX/SpeedIndex.n);
      SpeedIndex.StdDev = (float) sqrt(
	(SpeedIndex.SumX2 - SpeedIndex.n*SpeedIndex.mean*SpeedIndex.mean)/
	  SpeedIndex.n);
    }

    if(Closing.n) {
      Closing.mean = (float) (Closing.SumX/Closing.n);
      Closing.StdDev = (float) sqrt(
	(Closing.SumX2 - (double) Closing.n*Closing.mean*Closing.mean)/Closing.n);
    }

    if(Closed.n) {
      Closed.mean = (float) (Closed.SumX/Closed.n);
      Closed.StdDev = (float) sqrt(
	(Closed.SumX2 - (double) Closed.n*Closed.mean*Closed.mean)/Closed.n);
    }

    if(Opening.n) {
      Opening.mean = (float) (Opening.SumX/Opening.n);
      Opening.StdDev = (float) sqrt((Opening.SumX2 -
		  (double) Opening.n*Opening.mean*Opening.mean)/Opening.n);
    }

    if(Open.n) {
      Open.mean = (float) (Open.SumX/Open.n);
      Open.StdDev = (float) sqrt((Open.SumX2 -
			      (double) Open.n*Open.mean*Open.mean)/Open.n);
    }

    if(ClosingPct.n) {
      ClosingPct.mean = (float) (ClosingPct.SumX/ClosingPct.n);
      ClosingPct.StdDev = (float) sqrt((ClosingPct.SumX2 -
			      (double) ClosingPct.n*ClosingPct.mean*ClosingPct.mean)/ClosingPct.n);
      ClosingPct.mean = 100.0*ClosingPct.mean/(float) header.nSamplesPerSec;
      ClosingPct.StdDev = 100.0*ClosingPct.StdDev/(float) header.nSamplesPerSec;
    }

    if(OpeningPct.n) {
      OpeningPct.mean = (float) (OpeningPct.SumX/OpeningPct.n);
      OpeningPct.StdDev = (float) sqrt((OpeningPct.SumX2 -
			      (double) OpeningPct.n*OpeningPct.mean*OpeningPct.mean)/OpeningPct.n);
      OpeningPct.mean = 100.0*OpeningPct.mean/ (float) header.nSamplesPerSec;
      OpeningPct.StdDev = 100.0*OpeningPct.StdDev/ (float) header.nSamplesPerSec;
    }

    /* OPEN *.stt FILE */
    strcpy(str, outfile);
    if(!strchr(str, '.')) {
      strcat(str, ".stt");
    }
    else {
      strcpy(strchr(str, '.'), ".stt");
    }




    if((stat = fopen(str, "wt"))==NULL) {
      printf("Error while creating statitistical file\n");
      exit(1);
    }

    /*   Unvoiced/silent intervals
       |------------------dur----------------------------+
       |-InitOfTrack-+                               EndOfTrack.
       +.............+------------------+..+----+.+--+...+
       ....=
       0 to InitOfTracking  => assumed as being silence;
       EndOfTracking to Dur => assumed as being unvoiced.
     */

    dur = (float) header.datasize/(2.0 * header.nSamplesPerSec);
    if(dur > chain.InitOfTracking) {
      UnvoicedPct = 100.0*chain.UnvoicedTime/
		     (dur - chain.InitOfTracking);
    }
    else {
      UnvoicedPct = 0.0;
    }

    printf("\n Statistics in file %s\n", str);
    printf("------------------------------------------\n");
    printf(" MEASURE               MEAN       STD DEV\n");
    printf(" Dur          (sec)   %5.2f\n", dur);
    printf(" Unvoiced     (%%)     %5.2f\n", UnvoicedPct);
    printf(" Pitch Breaks (#)      %i\n", chain.PitchBreaks);
    printf(" F0           (Hz)   %5.2f        %5.2f\n", F0.mean, F0.StdDev);
    printf(" jitter       (%%)     %5.2f        %5.2f\n", jitter.mean, jitter.StdDev);
    printf(" shimmerRMS   (%%)     %5.2f        %5.2f\n", shimmerRMS.mean, shimmerRMS.StdDev);
    printf(" shimmerPeak  (%%)     %5.2f        %5.2f\n", shimmerPeak.mean, shimmerPeak.StdDev);
    printf(" snr          (dB)    %5.2f        %5.2f\n", snr.mean, snr.StdDev);
    printf(" CQ (Eq 3.8)  (%%)     %5.2f        %5.2f\n", ClosedQuot.mean, ClosedQuot.StdDev);
    printf(" CQ Area      (%%)     %5.2f        %5.2f\n", ClosedQuotArea.mean, ClosedQuotArea.StdDev);
    printf(" SI (Eq 3.7) (none)   %5.2f        %5.2f\n", SpeedIndex.mean, SpeedIndex.StdDev);

    printf("------------------------------------------\n");
    printf(" Cp (Eq 3.5)  (%%)     %5.2f        %5.2f\n", ClosingPct.mean, ClosingPct.StdDev);
    printf(" Op (Eq 3.6)  (%%)     %5.2f        %5.2f\n", OpeningPct.mean, OpeningPct.StdDev);
    printf("Cummulative file: %s\n", cumfile);

    fprintf(stat, " MEASURE                MEAN       STD DEV\n");
    fprintf(stat, " Dur          (sec)    %5.2f\n", dur);
    fprintf(stat, " Unvoiced     (%%)      %5.2f\n", UnvoicedPct);
    fprintf(stat, " Pitch Breaks (#)       %i\n", chain.PitchBreaks);
    fprintf(stat, " F0           (Hz)    %5.2f       %5.2f\n", F0.mean, F0.StdDev);
    fprintf(stat, " jitter       (%%)      %5.2f       %5.2f\n", jitter.mean, jitter.StdDev);
    fprintf(stat, " shimmerRMS   (%%)      %5.2f       %5.2f\n", shimmerRMS.mean, shimmerRMS.StdDev);
    fprintf(stat, " shimmerPeak  (%%)      %5.2f       %5.2f\n", shimmerPeak.mean, shimmerPeak.StdDev);
    fprintf(stat, " snr          (dB)     %5.2f       %5.2f\n", snr.mean, snr.StdDev);
    fprintf(stat, " CQ Time      (%%)      %5.2f       %5.2f\n", ClosedQuot.mean, ClosedQuot.StdDev);
    fprintf(stat, " CQ Area      (%%)      %5.2f       %5.2f\n", ClosedQuotArea.mean, ClosedQuotArea.StdDev);
    fprintf(stat, " SpeedIndex  (none)    %5.2f       %5.2f\n", SpeedIndex.mean, SpeedIndex.StdDev);
    fprintf(stat, " ClosingPct   (%%)     %5.2f        %5.2f\n", ClosingPct.mean, ClosingPct.StdDev);
    fprintf(stat, " OpeningPct   (%%)     %5.2f        %5.2f\n", OpeningPct.mean, OpeningPct.StdDev);

    fclose (stat);

    if(cummulative != -1) {

      fprintf(cum, "%5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",

	dur,   	               UnvoicedPct,     (float) chain.PitchBreaks,
	F0.mean,               F0.StdDev,
	jitter.mean,           jitter.StdDev,
	shimmerRMS.mean,       shimmerRMS.StdDev,
	shimmerPeak.mean,      shimmerPeak.StdDev,
	snr.mean,              snr.StdDev,
	ClosedQuot.mean,       ClosedQuot.StdDev,
	ClosedQuotArea.mean,   ClosedQuotArea.StdDev,
	SpeedIndex.mean,       SpeedIndex.StdDev,
	ClosingPct.mean,       ClosingPct.StdDev,
	OpeningPct.mean,       OpeningPct.StdDev

      );
    }

    return;
  }

  /* ================= Add values ======================= */

  F0.n ++;
  F0.SumX += NewCycle.F0;
  F0.SumX2 += (NewCycle.F0*NewCycle.F0);

  if(NewCycle.jitter <= 10.0) {
    jitter.n++;
    jitter.SumX += NewCycle.jitter;
    jitter.SumX2 += (NewCycle.jitter*NewCycle.jitter);
  }

  if(NewCycle.shimmerRMS != -1) {
    shimmerRMS.n++;
    shimmerRMS.SumX += NewCycle.shimmerRMS;
    shimmerRMS.SumX2 += (NewCycle.shimmerRMS*NewCycle.shimmerRMS);
  }

  if(NewCycle.shimmerPeak != -1) {
    shimmerPeak.n++;
    shimmerPeak.SumX += NewCycle.shimmerPeak;
    shimmerPeak.SumX2 += (NewCycle.shimmerPeak*NewCycle.shimmerPeak);
  }

  if(NewCycle.snr != 1e20) {
    snr.n++;
    snr.SumX += NewCycle.snr;
    snr.SumX2 += (NewCycle.snr*NewCycle.snr);
  }

  if(NewCycle.ClosedQuot != -1) {
    ClosedQuot.n++;
    ClosedQuot.SumX += NewCycle.ClosedQuot;
    ClosedQuot.SumX2 += (NewCycle.ClosedQuot*NewCycle.ClosedQuot);
  }

  if(NewCycle.ClosedQuotArea != -1) {
    ClosedQuotArea.n++;
    ClosedQuotArea.SumX += NewCycle.ClosedQuotArea;
    ClosedQuotArea.SumX2 += (NewCycle.ClosedQuotArea*NewCycle.ClosedQuotArea);
  }

  if(NewCycle.SpeedIndex != 2) {
    SpeedIndex.n++;
    SpeedIndex.SumX += NewCycle.SpeedIndex;
    SpeedIndex.SumX2 += (NewCycle.SpeedIndex*NewCycle.SpeedIndex);
  }

  if(NewCycle.Closing != -1 && NewCycle.Closed != -1 &&
     NewCycle.Opening != -1 && NewCycle.Open != -1 &&
     NewCycle.Closing != 1e20 &&  NewCycle.Closed != 1e20 &&
     NewCycle.Opening != 1e20 && NewCycle.Open != 1e20) {

    Closing.n++;
    Closing.SumX += (double) NewCycle.Closing;
    Closing.SumX2 += ((double) NewCycle.Closing*NewCycle.Closing);

    Closed.n++;
    Closed.SumX += (double) NewCycle.Closed;
    Closed.SumX2 += ((double) NewCycle.Closed*NewCycle.Closed);

    Opening.n++;
    Opening.SumX += (double) NewCycle.Opening;
    Opening.SumX2 += ( (double) NewCycle.Opening*NewCycle.Opening);

    Open.n++;
    Open.SumX += (double) NewCycle.Open;
    Open.SumX2 += ((double) NewCycle.Open*NewCycle.Open);

    if(LastCycle.InitClosing != -1 && LastCycle.EndClosing != -1 &&
       LastCycle.InitOpening != -1 && LastCycle.EndOpening != -1 &&
       (LastCycle.InitClosing < LastCycle.EndClosing) &&
       (LastCycle.InitOpening < LastCycle.EndOpening)) {
      OpeningPct.n++;
      aux = ((float) LastCycle.EndOpening - LastCycle.InitOpening)/
			LastCycle.pitch;
      /* pitch is in seconds; normalised by the sampling frequency at
	 the end of statistics */

      OpeningPct.SumX +=  (double) aux;
      OpeningPct.SumX2 += (double) aux*aux;

      ClosingPct.n++;
      aux = ((float) LastCycle.EndClosing - LastCycle.InitClosing)/
			LastCycle.pitch;
      ClosingPct.SumX +=  (double) aux;
      ClosingPct.SumX2 += (double) aux*aux;
    }
  }
}


void postprocessor(char *argv[])
{
  FILE *CorrectedFile;
  FILE *OldFile;

  int i,
      fim = FALSE;

  float aux,
	aux2;

  char x[40], y[20];
  struct {
    float t,
	  f0;
  }p[3];

  if((CorrectedFile = fopen(argv[arg.output], "wt"))==NULL) {
    printf("Error while opening tmp.tmp file in Postprocessor\n");
    exit(1);
  }

  if((OldFile = fopen("tmp.tmp", "rt"))==NULL) {
    printf("Error while opening OldFile in Postprocessor \n");
    exit(1);
  }

  /* read 1st line (.wav duration) */
  fscanf(OldFile, "%s%s", x, y);

  /* write to new file */
  fprintf(CorrectedFile, "0.0 "); fprintf(CorrectedFile, y);
  fprintf(CorrectedFile, "\n");

  /* initialize first 2 positions of buffer p[] */
  for(i=0; i<2; i++) {
    if(fscanf(OldFile, "%s%s\n", x, y)==2) {
      p[i].t = (float) atof(x);
      p[i].f0 = (float) atof(y);
    }
    else {
      fim = TRUE;
      break;
    }
  }

  if(!fim) {
    for(;;) {
      /* get new line from OldFile */
      if(fscanf(OldFile, "%s%s\n", x, y)==2){
	p[2].t = (float) atof(x);
	p[2].f0 = (float) atof(y);
      }
      else {
	break;
      }

      /* TIME CORRECTION */
      if(p[1].t == -1.0) {
	/* "spurious" -1 -1 ? */
	aux = (p[2].t - p[0].t)*p[0].f0;

	if(aux > 0.80 && aux < 1.20) {
	  /* overwrite -1 -1 in p[1] */
	  p[1].f0 = 1.0/(p[2].t - p[0].t);
	  p[1].t = p[2].t;
	}
	else if(p[0].t == p[2].t) {
	  if(fscanf(OldFile, "%s%s\n", x, y)==2){
	    p[1].t = (float) atof(x);
	    p[1].f0 = (float) atof(y);
	  }
	  else {
	    break;
	  }
	}
	else {
	  /* WRITE 1st buffer position */
	  fprintf(CorrectedFile, "%f %f\n", p[0].t, p[0].f0);

	  /* SHIFT p[] */
	  p[0].t = p[1].t; p[0].f0 = p[1].f0;
	  p[1].t = p[2].t; p[1].f0 = p[2].f0;
	}
      }
      else if( (p[0].t == -1) && (p[2].t == -1) ) {
	/* remove "isolated" f0
	   -1   x   -1
	    t0  t1  t2
	*/

	if(fscanf(OldFile, "%s%s\n", x, y)==2){
	  p[1].t = (float) atof(x);
	  p[1].f0 = (float) atof(y);
	}
	else {
	  break;
	}
      }
      else {
	/* sporadic halvig (zero crossing was missed) */
	aux  = 2.0*p[1].f0/(p[0].f0 + p[2].f0);
	aux2 = (p[2].t - p[1].t)/(p[1].t - p[0].t);

	if(aux>0.45 && aux < .55 && aux2>0.45 && aux2 < .55) {
	  p[1].f0 = (p[0].f0 + p[2].f0)/2;
	}

	/* WRITE 1st buffer position */
	fprintf(CorrectedFile, "%f %f\n", p[0].t, p[0].f0);

	/* SHIFT p[] */
	p[0].t = p[1].t; p[0].f0 = p[1].f0;
	p[1].t = p[2].t; p[1].f0 = p[2].f0;
      }
    }
  }

  /* WRITE REMAINING p[] */
  fprintf(CorrectedFile, "%f %f\n", p[0].t, p[0].f0);
  fprintf(CorrectedFile, "%f %f\n", p[1].t, p[1].f0);

  fclose(CorrectedFile);
  fclose(OldFile);

  /* del OldFile */
  system("del tmp.tmp");
}



char* strrev (char *str)
{
        char *left = str;
        char *right = str + strlen(str) -1;
        char ch;

       while (left < right)
        {
                ch = *left;
                *left++ = *right;
                *right-- = ch;
        }

        return(str);
}


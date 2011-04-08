
/*
    Acoustic.c
    (c) Maurilio Nunes Vieira 1996

    F0 tracker based on
      Kurt Schafer Vincent, "Pitch Period Detection and Chaining: Method
      and Evaluation," Phonetica 40: 177-202 (1983).

    compile in the compact memory model
*/


/* 

Modified 9 Feb. 2007 - Joao SANSAO 

* Ported to gcc  

Main modifications: 
* excluded borland libraries io.h, conio.h and stat.h
* wave header format 
* changed int to signed short type for the data array (x) 
* implemented strrev not present in gcc library

To do:
* correct time function callings to gcc library

* seems to work fine under GNU/Linux
* tested with gcc 4.1.2 20060928

* Known bugs: segmentation fault when option -f not specified, 
workaround: always add -f [m|f|t] 

* gcc compilation procedure :

$gcc -lm -Wall acoustic.c -o acoustic

*/








#define	TRUE		1
#define	FALSE 		!TRUE

#define	YES		0
#define	NO 		!YES

#define	MAX		1
#define	MIN	       -1
#define	NOT_MAX_MIN 	0

#define	DEFAULT_NOISE   600     /* default value */

#define	SIG_LEN		20      /* significant points (SP) buffer length */
#define	TWIN_LEN	50     /* twin period (TP) buffer length */

#define NO_CHAIN	100     /* flag for twin[].link */
#define FPW_BUF_LEN	64      /* length of the output-file write buffer */
#define MAXBINS		21

#define PEAKS		1
#define IntPEAKS	2
#define NZP	        3
#define PZN	        4
#define RMS		2
#define MAXLAG 		10

#define JITTER		1
#define SHIMMER		2
#define HNR		3
#define COMB		1

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


/* for statistical analysis */
//#include "dos.h"
#include "time.h"


/* FUNCTION PROTOTYPES */
void usage(void),
     initialization(char *argv[]),
     link(void),
     update(void),
     install(void),
     overlap(void),
     postprocessor(void),
     shimmerms(int, int),
     GlottalNoise(int, int, int),
     comb2(int, int, int),
     statistics(int, int, int, float);

int significant_point(void),
    test1(void),
    test2(void),
    test4(void),
    test5(void),
    test6(void),
    test6a(void),
    output(void),
    sign(float value),
    fit(int, int, int),
    outofmemory(void);

/*float timedif(struct tm, struct tm),*/
float FindZx(int, int) ; /*, fmin(float, float)*/

char* strrev (char *str);

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
;

/* BUFFER OF TWIN PERIODS (TP) */
struct {
  long int dist3,	/* distance of the TP right end to the present time */
	   ntotal,	/* number of TP in the chain */
	   dtotal;	/* duration of the chain (in # samples) */


  int frame,            /* frame number */
      t3,               /* time of 3rd SP of the TP */
      x1,               /* amplitude of x1 */
      x2,               /*     "      " x2 */
      x3,               /*     "      " x3 */
      dur12,         	/* duration of 1st period */
      dur23,		/* duration of 2nd period */
      type,          	/* type of SP: MAX (+)  or MIN (-) */
      link,          	/* pointer to previous PG */
      rightend,      	/* boolean: is('nt) the right end of the chain */
      empty;         	/* boolean: buffer position empty or full */
  double int_dur23,	/* interpolated  dur23, float isn't enough */
       int_dur12,
	t3_int;		/* ", t3 */
  float zxjitterN_,	/* NZP */
	zxjitter_P,

	zxjitterP_,     /* PZN */
	zxjitter_N,
	shimmer_peak,	/* shimmer */
	shimmer_rms,
	hnr;		/* signal- or harmonic-to-noise ratio */
} twin[TWIN_LEN];


/* BUFFER OF SIGNIFICANT POINTS (SP) */
struct {
  int x,             	/* amplitude */
      type;          	/* MAXimum  or  MINimum */
  long int time;     	/* position */
  double interp_time;	/* float isn't precise enough */
  float zxP_, zx_N,	/* interpolated zero crossing instant */
	zxN_, zx_P;
} sig[SIG_LEN];


/* GENERAL USE */
int i, k, j, m, n, q, r,	/* counters */
    d;				/* distances */
float razao,            	/* ratios */
      aux;

/* SKELETIZATION */
int tipo;  			/* temporary: MAX or MIN */


/* SEARCH OF TP */
int i1, i2, i3;       		/* pointers */
int x1, x2, x3;			/* magnitudes */
int c1, c2;			/* counters */
long int t1, t2, t3; 		/* time intervals */
double t1_int = 0.0,            /* used in parabolic interpolation */
       t2_int = 0.0;
float E12, E23;			/* energy between t12 and t23 (for shimmer) */
double hnr;			/* harmonic or signal to noise ratio */

int stop = FALSE;		/* program control flag */

float slope23,        	/* slope of the line between x2 and x3 */
slope12;        	/*   "   "   "    "     "    x1 and  x2 */

long tin;		/* used in test6 */
int w12, w23,           /*  "    "   "    (time slots) */
    L12, L23;           /*  "    "   "    (length of time slots) */

float am12, am23,       	/* in test6: av. magnitude of periods of TP */
      amv12[8], amv23[8],	/* "   "   : av. mag. of time slots */
      var12, var23,		/* "   "   : abs. mag. var. of each period */
      dif;

float zx1P_, zx1_N, zx1N_, zx1_P,
      zx2P_, zx2_N, zx2N_, zx2_P,
      zx3P_, zx3_N, zx3N_, zx3_P; /* interpolated zero crossings */


/* CHAINING */
long int delta,	      	/* auxiliar: t3 - t3old */
	 t3old,         /* t3 of previous TP */
	 gap;         	/* dist. between end of 1st and beg. of 2nd period */
int alivept,          	/* rigth end of and 'alive' chain */
    precede,        	/* temporay: pointer to the previous chain element */
    append,         	/* boolean: present TP is to be appended */
    done,           	/* boolean: loop control */
    perfect_match,  	/* boolean: perfect-match of TP */
    atual,          	/* time of present SP */
    inhibit,        	/* boolean: inhibit chaining if TP is the double or
			   triple of (assumed) true period */
    birth_ready;       	/* boolean: chain gained birth condition */


/* FILE RELATED */
FILE *fdr;				/* data file's header */
FILE *fpcum;				/* cumulative file's header*/
FILE *fpwt;                             /* time series file's header*/
FILE *TmpNZP, *TmpPZN;			/* used in maxratio method */

char InFile[50];			/* wav input file */
char OutFile[50];			/* time series */
char CumFile[50];			/* cummulative file */

int cum_flag = 0;

signed short *x = NULL,				/* speech samples buffer */
    Nframe,                             /* frame counter */
    ni;					/* nb. of samples actually read */

char *fpwt_buf = NULL;			/* buffer for file output */
int fpwt_buf_cnt;			/* buffer counter */


/* TIME INTERVALS */
int milisec1,
    Lframe;


/* PARAMETERISATION */
struct {
  int   limiar,		/* noise threshold level (in ADC units) */
	N,		/* 1/N = max f0 detectable */

	t1_Pmax;        /* 1/t1_Pmax = min f0 detectable */

  float t2_Rmin,        /* boundary for max jitter 'within' twin period */
	t2_Rmax,	/* idem */
	t5_Rmin,	/* boundary for max shimmer within twin period */
	t5_Rmax,	/* idem */

	t6_K,         	/* allowance for similarity test */

	link_Rmin,      /* boundary for max jitter 'between' twin periods */
	link_Rmax;	/* idem */

  int  	outp_BirthL,	/* min chain length required for birth */
	outp_BirthN;	/* min required # elements in chain for birth */
  float outp_K23;	/* how many times a chain must be greater than dur23
			   to get birth condition */
  float gap;		/* used in postprocessor() */
  int pitchdetection_alg, /* NZP or PZN */
      jitter_alg,	/* from peaks (1) or rms values (2) */
      shimmer_alg,
      hnr_alg;
} par;


/* STATISTICS */
struct STATISTICS {
  int n;
  double SumX, SumX2,
	 mean, StdDev;
};

time_t  prog_start,  prog_end,
	     loop_start,  loop_end;

int pitchbreaks = 0;

float UnvoicedTime = 0.0, InitOfTrack = 0.0;

/* HISTOGRAM */
/*
 to do later: distribution of f0: 50-500 Hz at 5 Hz step -> 90 bins
 this may be useful in the "siren" -> gross f0 range and intermediate values
 -> paper for ICLSP with selected patients of eg paralysis?
*/


int infile_arg = -1,
    outfile_arg = -1,
    threshold_arg = -1,
    /*f0range_arg = -1,*/ 
    f0range_arg = 0,
    pitchdetection_arg = -1,
    jitter_arg = -1,
    shimmer_arg = -1,
    hnr_arg = -1,

    linsmooth_arg = -1,
    gap_arg = -1,
    medfilter_arg = -1;

struct {
  float f0_peak,
	f0_intpeak,
	shimmer_peak,
	shimmer_rms,
	jitterP_, jitter_P,
	jitterN_, jitter_N,
	hnr_comb;
} outvalue;



void main(argc,argv)


int argc;
char *argv[];

{
  char str[50], str2[50], *p;

  /*gettime(&prog_start);*/
  
 time (&prog_start);
  

  /* CHECK PROGRAM CALL */
  if(argc<2)  usage();

  for(i=1; i<argc && *argv[i] == '-'; i++) {
    j = i+1;

    switch(argv[i][1]) {
      case 'M':
      case 'm':
      case 'L':
      case 'l':
	j--;

      default:
	if (argc <= j) usage();

	switch(argv[i++][1]) {
	  case 'c':
	  case 'C':
	    cum_flag = 1;
	    strcpy(CumFile, argv[i]);
		
	    break;
	  case 'i':
	  case 'I':
	    infile_arg = i;		
	    strcpy(InFile, argv[i]);
	    break;
	  case 'o':
	  case 'O':
	    outfile_arg = i;
	    strcpy(OutFile, argv[outfile_arg]);
	    break;
	  case 't':
	  case 'T':
	    threshold_arg = i;
	    break;
	  case 'f':
	  case 'F':
	    f0range_arg = i;
	    break;
	  case 'p':
	  case 'P':
	    pitchdetection_arg = i;
	    break;
	  case 'j':
	  case 'J':
	    jitter_arg = i;
	    break;
	  case 's':
	  case 'S':
	    shimmer_arg = i;
	    break;
	  case 'n':
	  case 'N':
	    hnr_arg = i;
	    break;
	  case 'M':	/* median filter */
	  case 'm':
	    medfilter_arg = 1;
	    i--;
	    break;
	  case 'L':	/* linear smoothing */
	  case 'l':
	    linsmooth_arg = 1;
	    i--;
	    break;
	  case 'g':
	  case 'G':
	    gap_arg = i;
	    break;
	  default:
	    usage();
	    break;
	}
	break;
    }
  }
  if ((i != argc && *argv[i] != 'i')
       || outfile_arg == -1 || infile_arg == -1) usage();

  if ( ((pitchdetection_arg != -1) && (jitter_arg != -1)) ||
      ((pitchdetection_arg != -1) && (shimmer_arg != -1)) ||
      ((shimmer_arg != -1) && (jitter_arg != -1)) ||
      ((pitchdetection_arg != -1) && (hnr_arg != -1)) ||
      ((jitter_arg != -1) && (hnr_arg != -1)) ||
      ((shimmer_arg != -1) && (hnr_arg != -1)) ||
      (((pitchdetection_arg == -1) && (jitter_arg == -1)) &&
      ((shimmer_arg == -1) && (hnr_arg == -1)))) {
    printf("Choose one of the (mutually exclusive) [-p], [-j], [-s], or [-n]\n");
    exit(0);
  }



  /* OPEN SPEECH FILE */
  if((fdr = fopen(argv[infile_arg], "rb")) !=  NULL) {
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

  /* OPEN tmp.tmp TMP OUTPUT FILE */
  if((fpwt = fopen("tmp.tmp", "w"))==NULL) {
    printf("Error while opening output file (%s)\n", str);
    exit(1);
  }



  /* OPEN cummulative FILE */
  if(cum_flag) {
    if((fpcum = fopen(CumFile, "a+t"))==NULL) {
      printf("Error while creating cummulative file\n");
      exit(1);
    }

    /*
    fprintf(fpcum, "PID Dur Unv PB Fav FSd FSmAv FSmSd Jav ratio NzPzN JSd Sav Ssd SNRav SNRsd\n");
    */

    /*-- take the 1st 6 characters of file name (PID) --*/
    strcpy(str2, InFile); /* path...db0964la.wav */
    strrev(str2);  /* str2=> vaw.al4790bd...thap */
    strcpy(str, "12345678");

    p = strchr(str2, '.');
    if(p != NULL) strncpy(str, p+1, 8);
    else strncpy(str, str2, 8);

    strrev(str); /* str=> db0964la */
    strcpy(str2, "123456");
    strncpy(str2, str, 6); /* str2 = db0964 */
    fprintf(fpcum, "%s ", str2);
  }

   


  /* INITIALIZE VARIABLES, BUFFERS, ETC */
  initialization(argv);

  fprintf(fpwt,"0.0 ");
  gcvt((float) 8.0*header.datasize/
	       (header.wBitsPerSample*header.nSamplesPerSec), 10, str);
  fprintf(fpwt, str);
  fprintf(fpwt, "\n");


  /* PUT MESSAGES ON THE SCREEN */
  printf(" (c) Maurilio N. Vieira, Measures of Acoustic Perturbation\n");
  printf(" Time series: ");
  if(pitchdetection_arg != -1)
    printf("Pitch estimation\n");
  else if(jitter_arg != -1)
    printf("Jitter estimation\n");
  else if(par.shimmer_alg == PEAKS)
    printf("Shimmer based on peaks\n");
  else if(par.shimmer_alg == RMS)
    printf("Shimmer based on RMS values\n");
  else if(par.hnr_alg == COMB)
    printf("HNR based on Comb Filtering\n");


  if(par.pitchdetection_alg == PEAKS || par.jitter_alg == PEAKS)
    printf("Parabolic Interpolation = OFF\n");
  else if(par.pitchdetection_alg == IntPEAKS || par.jitter_alg == IntPEAKS)
    printf("Parabolic Interpolation = ON\n");

  if(par.jitter_alg == PZN)  printf("PZN zero crossing pattern\n");
  else if(par.jitter_alg == NZP) printf("NZP zero crossing pattern\n");

  if(linsmooth_arg == 1) printf("Linear Smoothing = ON\n");
  else printf("3 point Hamming Window linear smoothing = OFF\n");

  if(medfilter_arg == 1) printf("Median Filter = ON\n");
  else printf("3 point Median Filter non-linear smoothing = OFF\n");

  if(gap_arg != -1) printf("Gap in post processor = %5.2f ms\n", par.gap*1000.0);


  /*---------------------- ANALYSIS LOOP ------------------------------*/
  printf("\nWait...");

 /* gettime(&loop_start);*/
  time (&loop_start);

  /* 20 jul 95 */
  while((ni =
	(fread(x+ Lframe + 2*par.N, sizeof(int), Lframe, fdr))) != NULL) {
   if(ni <= par.N)
     break;

   /* 20 jul 95 */
   for(i = Lframe + par.N; i<Lframe+par.N+ni; i++) {
     if(abs(x[i]) > par.limiar) {
       if(significant_point()) {
	 stop = FALSE; /* my change */
	 for(c2=1; c2<=i3; c2++) {
	   if(test1()) { /* f0 max */
	     for(c1=1; c1<=i2; c1++) {
	       if(test2()) { /* aprox. duration: jitter<10% */
		 if(test4()) { /* envelope */
		   if(test5()) { /* low shimmer within twin-period */
		     if(test6a()) { /* similiarity of period (my change) */
		       /* a twin-period was detected */
		       update();
		       link();
		       install();
		       if(output()) {
			 stop = TRUE;  /* my change */
			 break;
		       }
		     }
		   }
		 }
	       }
	     }
	   }
	  if(stop) break;  /* my change */
	 }
       }
     }
   }
   overlap();
  }

  if(fpwt_buf_cnt > 0) {
    fprintf(fpwt, fpwt_buf);
    fprintf(fpwt, "-1 -1\n");
  }

  free(fpwt_buf);
  free(x);

  fclose(fdr);
  fclose(fpwt);
  /*gettime(&loop_end);*/
  //  time (&loop_end);
  //  printf("Correcting time...\n");

  if(pitchdetection_arg != -1 ) postprocessor();

  statistics(1, 0, cum_flag, 0.0);
  if(cum_flag) fclose(fpcum);

  /* pzn.pzn and nzp.nzp are closed in statistics() */
 
  //  time(&prog_end);




 //printf("\n-------------------PERFORMANCE-----------------\n");
 // printf("Total Time: %f\n", difftime(prog_start, prog_end));
 //  printf("Initialization Time: %f\n", difftime(prog_start, loop_start));
 //  printf("Loop Time: %f\n", difftime(loop_start, loop_end));
 //  printf("Postproc Time: %f \n", difftime(loop_end, prog_end));

  printf("\ndone.");

  return(0);
}

/*---------------------------- FUNCTIONS ---------------------------*/

void usage()
{
  printf("acoustic -i infile -o outfile [-c CumFile] [-options {default}]\n");
  printf("-i  .wav file, mono\n");
  printf("-o  ascii file (time series)\n");
  printf("-c  cummulative file (ascii)\n\n");

  printf("F0 detection (common to all measurements):\n");
  printf("-t x (noise Threshold in ADC units) {x=600}\n");
  printf("-f x (f0: (m)ale 77-250 Hz, (f)emale 100-500 Hz, (t)otal 50-500 Hz) {x=t}\n\n");

  printf("Output measurement:\n");
  printf("-p x (pitch [F0 Hz], 1=peaks, 2=interp. peaks )\n");
  printf("-j x (jitter [PF1%], 1=peaks, 2=interp. peaks, 3=NZP, 4=PZN)\n");
  printf("-s x (shimmer [PF1%], 1=peaks, 2=rms values\n");
  printf("-n x (glottal noise [SNR dB], 1: Comb Filter\n\n");

  printf("Postprocessor settings:\n");
  printf("-m   (median non linear filter)\n");
  printf("-l   (linear smoothing)\n");
  printf("-g x (maximum gap to be filled [milliseconds]) {x=25}\n");

  exit(0);
}

void initialization(char *argv[])
{
  milisec1  = (int) (header.nSamplesPerSec * 0.001/2.0)*2;

  Lframe =50*milisec1; /* should be > 2*par.t1_Pmax*/

  par.pitchdetection_alg = -1;
  par.jitter_alg = -1;
  par.shimmer_alg = -1;
  par.hnr_alg = -1;

  /* set control parameters of the algorithm */
  if(pitchdetection_arg != -1) {
    par.pitchdetection_alg = atoi(argv[pitchdetection_arg]);
    if(par.pitchdetection_alg != IntPEAKS && par.pitchdetection_alg != PEAKS)
      usage();
  }

  if(shimmer_arg != -1) {
    par.shimmer_alg = atoi(argv[shimmer_arg]);
    if(par.shimmer_alg != 1 && par.shimmer_alg != 2)
      usage();
  }

  if(jitter_arg != -1) {
    par.jitter_alg = atoi(argv[jitter_arg]);

    if(par.jitter_alg < 1 || par.jitter_alg > 5)
     usage();
  }

  if(hnr_arg != -1) {
    par.hnr_alg = atoi(argv[hnr_arg]);

    if( abs(par.hnr_alg) != 1 ) usage();
  }

  if(threshold_arg == -1)
    par.limiar = DEFAULT_NOISE;
  else
    par.limiar = atoi(argv[threshold_arg]);

  /* maximum gap allowed in postprocessor() */
  if(gap_arg == -1)
    par.gap = 0.025; /* milliseconds */
  else  {
    par.gap = atof(argv[gap_arg])/1000.0;
    if(par.gap < 0 ) usage();
  }




  if (tolower(argv[f0range_arg][0]) == 'm') {
    /* my settings for male voice */

    par.N = 4*milisec1;			/* -> f0 max = 250 Hz */
    par.t1_Pmax = 13*milisec1;		/* -> f0 min =  77 Hz */
    par.t2_Rmin = 0.9;
    par.t2_Rmax = 1.1;
    par.t5_Rmin = 0.3;
    par.t5_Rmax = 3.3;
    par.t6_K = 1.0;
    par.link_Rmin = 0.7;
    par.link_Rmax = 1.3;
    par.outp_BirthL = 30*milisec1;
    par.outp_BirthN = 3;
    par.outp_K23 = 2.5;
  }
  else if (tolower(argv[f0range_arg][0]) == 'f') {
    /* my settings for female voice */
    par.N = 2*milisec1;			/* f0 max = 500 Hz */
    par.t1_Pmax = 10*milisec1;		/* f0 min = 100 Hz */
    par.t2_Rmin = 0.9;
    par.t2_Rmax = 1.1;
    par.t5_Rmin = 0.3;
    par.t5_Rmax = 3.3;
    par.t6_K = 1.0;
    par.link_Rmin = 0.7;
    par.link_Rmax = 1.3;
    par.outp_BirthL = 30*milisec1;
    par.outp_BirthN = 3;
    par.outp_K23 = 2.5;
  }
  else if (tolower(argv[f0range_arg][0]) == 't') {
    /* my all range settings -> test  */
    par.N = 2*milisec1;			/* f0 max = 500 Hz */
    par.t1_Pmax = 20*milisec1;		/* f0 min = 50 Hz */
    par.t2_Rmin = 0.9;
    par.t2_Rmax = 1.1;
    par.t5_Rmin = 0.3;
    par.t5_Rmax = 3.3;
    par.t6_K = 1.0;
    par.link_Rmin = 0.7;
    par.link_Rmax = 1.3;
    par.outp_BirthL = 30*milisec1;
    par.outp_BirthN = 3;
    par.outp_K23 = 2.5;
  }
  else {
    /* default settings from original paper */

    par.N = 2*milisec1;
    par.t1_Pmax = 20*milisec1;
    par.t2_Rmin = 0.9;
    par.t2_Rmax = 1.1;
    par.t5_Rmin = 0.5;
    par.t5_Rmax = 2.0;
    par.t6_K = 0.5;
    par.link_Rmin = 0.86;
    par.link_Rmax = 1.16;
    par.outp_BirthL = 30*milisec1;
    par.outp_BirthN = 3;
    par.outp_K23 = 2.8;

  }

  if((x = malloc( sizeof(int) * 2*(par.N + Lframe))) == NULL) {
    outofmemory();
  }


  for(i=0; i<Lframe; i++) {
    x[i] = 0;
  }

  fread(x + Lframe, sizeof(int), 2*par.N, fdr);

  if((fpwt_buf = malloc( 15*2*500  )) == NULL) {
    outofmemory();
  }

  fpwt_buf[0] = '\0';
  fpwt_buf_cnt = 0;

  /* initialize buffers */
  for(i=0; i<SIG_LEN; i++) {
    sig[i].x = 0;
    sig[i].time = 0L;
    sig[i].interp_time = 0.0;
    sig[i].type = 0;
    sig[i].zxP_ = 0.0; sig[i].zx_N = 0.0;
    sig[i].zxN_ = 0.0; sig[i].zx_P = 0.0;
  }

  for(i=0; i<TWIN_LEN; i++) {
    twin[i].empty = YES;
    twin[i].t3 = 0;
    twin[i].t3_int = 0.0;
    twin[i].int_dur12 = 0.0;
    twin[i].int_dur23 = 0.0;
    twin[i].zxjitterP_ = 0.0; twin[i].zxjitter_N = 0.0;
    twin[i].zxjitterN_ = 0.0; twin[i].zxjitter_P = 0.0;
  }

  Nframe = 0;
  alivept = NO_CHAIN;
  t3old = 0;
}



int significant_point()
{

  /* teste 3: pontos significativos devem ser notavelmente maior
		    que o ruido de quantizacao */


  /* ESQUELETIZACAO */
  /* procura ponto significativo positivo */
  tipo = MAX;
  for(k=1; k <= par.N; k++) {
    if(!((x[i]>0)&&(x[i]>=x[i+k])&&(x[i]>x[i-k]))) {
      tipo=NOT_MAX_MIN;
      break;
    }
  }

  /* procura ponto significativo negativo */
  if(tipo==NOT_MAX_MIN) {
    tipo = MIN;
    for(k=1; k <=  par.N; k++) {
      if(!((x[i]<0)&&(x[i]<=x[i+k])&&(x[i]<x[i-k]))) {
	tipo = NOT_MAX_MIN;
	break;
      }
    }
  }

  if(tipo!=NOT_MAX_MIN) {

    /* x[i] E' UM PONTO SIGNIFICATIVO:
       atualiza buffer de pontos significativos */
    for(k=0; k <= SIG_LEN - 2; k++) {
      sig[k].time = sig[k+1].time;
      sig[k].interp_time = sig[k+1].interp_time;
      sig[k].x    = sig[k+1].x;
      sig[k].type = sig[k+1].type;
      sig[k].zxP_   = sig[k+1].zxP_; sig[k].zx_N   = sig[k+1].zx_N;
      sig[k].zxN_   = sig[k+1].zxN_; sig[k].zx_P   = sig[k+1].zx_P;
    }

    sig[SIG_LEN-1].time = i;

    /* possibility of ovwerflow! */
    sig[SIG_LEN-1].interp_time = (float) i +
	  (( (float) x[i-1] - x[i+1])/(2.0*x[i+1] - 4.0*x[i] + 2.0*x[i-1]))
	  /header.nSamplesPerSec;

    sig[SIG_LEN-1].x = x[i];
    sig[SIG_LEN-1].type = tipo;

    i3 = SIG_LEN-1;
    t3 = sig[i3].time;
    x3 = sig[i3].x;


    /* interpolated zero crossing */
    sig[i3].zxP_ = 0.0; sig[i3].zx_N = 0.0;
    sig[i3].zxN_ = 0.0; sig[i3].zx_P = 0.0;

    zx3P_ = 0.0; zx3_N = 0.0;
    zx3N_ = 0.0; zx3_P = 0.0;

    /*
     create variables; open two files and decide after statistics;
     do statistics for zx_pzn and _nzp
    */
    if(tipo == MAX) {
      zx3P_ = sig[i3].zxP_ = FindZx(tipo, PZN);
      zx3_P = sig[i3].zx_P = FindZx(tipo, NZP);
    }
    else if (tipo == MIN) {

      zx3N_ = sig[i3].zxN_ = FindZx(tipo, NZP);
      zx3_N = sig[i3].zx_N = FindZx(tipo, PZN);
    }

    return TRUE;
  }
  else {
    return FALSE;
  }
}


int test1()
{
  /* procura t2 */
  i2 = SIG_LEN-c2-1;

  /* teste #1 */

  if((sig[i2].type == tipo) && ((t3 - sig[i2].time) <= par.t1_Pmax)) {
    t2 = sig[i2].time;
    t2_int = sig[i2].interp_time;
    x2 = sig[i2].x;

    if(tipo == MAX) {
      zx2P_ = sig[i2].zxP_;
      zx2_P = sig[i2].zx_P;
    }
    else {
      zx2N_ = sig[i2].zxN_;
      zx2_N = sig[i2].zx_N;
    }

    return TRUE;
  }
  else {
    return FALSE;
  }
}


int test2()
{

  /* procura t1 */
  i1=i2-c1;

  if(sig[i1].type == tipo) {
    /* extremos de mesma polaridade */
    t1 = (long) sig[i1].time;
    t1_int = sig[i1].interp_time;
    x1 = sig[i1].x;
    if(tipo==MAX) {
      zx1P_ = sig[i1].zxP_;
      zx1_P = sig[i1].zx_P;
    }
    else if (tipo == MIN) {
      zx1N_ = sig[i1].zxN_;
      zx1_N = sig[i1].zx_N;
    }

    /* calcula relacao entre duracao dos periodos */
    razao = (float) (t2-t1)/(t3-t2);

    /* teste #2: a diferenca na variacao da duracao entre
       os dois periodos do PG esta' limitada a 10% */

    if((razao >= par.t2_Rmin)&&(razao <= par.t2_Rmax)) {
      slope23 = (float) (x3-x2)/(t3-t2);

      return TRUE;
    }
    else {
      return FALSE;
    }
  }
  else {
    return FALSE;
  }
}


int test4()
{
  /* teste #4: linhas retas entre (x1,t1)-(x2,t2) e (x2,t2)-(x3,t3)
     devem envolver quaisquer outros pontos intermediarios */

  /* testa envoltoria do periodo 2-3 */
  for(m=i2+1; m<i3; m++) {
    if((sig[m].type==tipo)&&(abs(sig[m].x) >
		      fabs(x2+(sig[m].time-t2)*slope23))) {

      return FALSE;
    }
  }

  /* testa envoltoria do periodo 1-2 */
  slope12 = (float)(x2-x1)/(t2-t1);
  for(m=i1+1; m<i2; m++) {
    if((sig[m].type==tipo)&&(abs(sig[m].x)>
				fabs(x1+(sig[m].time-t1)*slope12))) {
      return FALSE;
    }
  }

  return TRUE;
}


int test5()
{
  /* calcula x2/x_medio */
  razao =  (float) x2 * 2.0/( (float) x1 + (float) x3 );

  /* teste #5: a envoltoria determinada por
     x1, x2 e x3 nao varia mais que ... */
  if((razao >= par.t5_Rmin) && (razao <= par.t5_Rmax)) {
    return TRUE;
  }
  else {
    return FALSE;
  }
}


/*
  Waveform matching by linear time warping between t12 and t23
*/

int test6a()
{
  register int z, w, index1, index2;
  float e, E, ratio, alpha;

  /*
     => t1 and t2 are subtracted by L frame in overlap(),
	which is called before a new frame is read from the speech file;
     => par.N is an offset;
     => t3 is always in the current frame.
  */

  if(t1 < 0) index1 = (int) t1 + Lframe + par.N; /* SP do frame anterior*/
  else index1 = (int) t1; /* SP no frame corrente */

  if(t2 < 0) index2 = (int) t2 + Lframe + par.N;
  else index2 = (int) t2;

  e = 0.0;
  E = (t3 - index1)*( (float) abs(x1) + abs(x2) + abs(x3))/3.0;
  alpha = (float) (t3-index2-1.0)/(index2-index1-1.0);

  for(z=index1; z<index2; z+=2) {
    w = (int) ((z - index1)*alpha + index2);
    e +=  fabs(x[z] - x[w]);
  }

  ratio = e/E;

  if(ratio < 0.07) { /* empirical */
    shimmerms(index1, index2);
    GlottalNoise(index1, index2, par.hnr_alg);
    return TRUE;
  }
  else {
   return FALSE;
  }
}




void update()
{
  delta=t3-t3old;
  t3old=t3;
  for(q=0; q<TWIN_LEN; q++) {
    if(twin[q].empty == NO) {

      /* desloca buffer */
      twin[q].dist3 += delta;

      if(twin[q].rightend) {
	/* otimizacao: libera cadeias que nao passariam no link() */

	if(!((twin[q].dist3 <= 3*twin[q].dur23) ||
	       ((twin[q].dist3 <= 5*twin[q].dur23) && (alivept==q)))) {

	  /* libera extremidade da cadeia */
	  twin[q].empty = YES;

	  if(alivept != q) {
	    /* o elemento atual nao e' extremo de cadeia viva */

	    r = twin[q].link; /* aponta para elemento anterior */

	    while(r != NO_CHAIN) {
	      twin[r].empty = YES;  /* libera */
	      r=twin[r].link;  /* aponta para elemento anterior */
	    }
	  }
	  else alivept = NO_CHAIN;
	}
      }
    }
  }
}


void link()
{
  append = FALSE;

  if(alivept != NO_CHAIN) {
    /* se ha' uma cadeia viva, tenta encaixar o novo
       PG nesta cadeia */

    razao=(float)twin[alivept].dur23/(t2-t1);

    if((razao >= par.link_Rmin) && (razao <= par.link_Rmax)) {

      gap=twin[alivept].dist3-(t3-t1);

      if((float) gap <= 1.1*(t2-t1)) { /* parameter? */
	append = TRUE;
	precede = alivept;
      }
    }
  }
  if(!append) {
    /* o novo PG nao pode ser apendado
       a uma cadeia viva */

    for(q=0; q < TWIN_LEN; q++) {
        /* procura cadeias formadas e nao vivas (pre'-natais) */
      if((twin[q].empty==NO) && (twin[q].rightend)) {

	/* testa de forma mais rigorosa se o novo
	   PG encaixa numa cadeia na condicao pre'-natal */

	razao=(float)twin[q].dur23/(t2-t1);

        /* faz um pre' teste */
	if((razao >= par.link_Rmin) && (razao <= par.link_Rmax)) {
	  gap=twin[q].dist3 - (t3-t1);
	  done = FALSE;

	  /* faz o teste mais rigoroso para evitar que novas cadeias
	     nao sejam iniciadas em situacoes complicadas */
	  perfect_match = FALSE;
	  r=q;

	  while((gap<=0)&&(!done)) {
	    if(twin[r].type == tipo) {
	      if(((t3-t1)==twin[r].dist3) ||
		 ((t3-t1)==(twin[r].dist3+twin[r].dur23))||
		 ((t3-t1)==(twin[r].dist3+twin[r].dur23 + twin[r].dur12))||
		 ((t3-t2)==twin[r].dist3)) {

		perfect_match = TRUE;
	      }
	    }
	    if((perfect_match) || (twin[r].type!=tipo)) {
	      append=TRUE;
	      precede=q;
	    }
	    if(append || twin[r].link == NO_CHAIN)
	      done=TRUE;

	    if(!done) {
	      r=twin[r].link;
	      gap=twin[r].dist3-(t3-t1);
	    }
	  }
	}
      }
      if(append) break;
    }
  }
}


void install()
{
  for(q=0; q<TWIN_LEN; q++) {
    if(twin[q].empty==YES) { /* procura um elemento vazio no buffer */
      atual=q;
      twin[atual].dist3 = 0;
      twin[atual].dur12 = (int) (t2-t1);
      twin[atual].dur23 = (int) (t3-t2);
      twin[atual].type = tipo;
      twin[atual].rightend = TRUE;
      twin[atual].empty = NO;

      /* parabolic interpolation */
      /* possibility of ovwerflow */
      twin[atual].t3_int = (float) i +
	  (( (float) x[i-1] - x[i+1])/(2.0*x[i+1] - 4.0*x[i] + 2.0*x[i-1]))
	  /header.nSamplesPerSec;

      /* elementos introduzidos nesta modificacao do algoritmo */
      twin[atual].t3 = i;

      twin[atual].int_dur23 = twin[atual].t3_int - t2_int;
      twin[atual].int_dur12 = t2_int - t1_int;

      twin[atual].frame = Nframe;
      twin[atual].x1 = x1;
      twin[atual].x2 = x2;
      twin[atual].x3 = x3;

      /*----------------------MEASUREMENTS--------------------*/
      switch(twin[atual].type) {
	case MAX:
	  if( (zx1P_ != 0.0) && (zx2P_ != 0.0) && (zx3P_ != 0.0) ) {

	  /*                |(zx3 - zx2) - (zx2 - zx1)|
	   jitter = 100 * --------------------------------
			  0.5*[(zx3 - zx2) + (zx2 - zx1)]
	  */

	  twin[atual].zxjitterP_ =
	   200.0*fabs( (zx3P_ - zx2P_) - (zx2P_ - zx1P_) )/(zx3P_ - zx1P_);
	  }
	  else {
	    twin[atual].zxjitterP_ = 0.0000001;
	  }

	  if( (zx1_P != 0.0) && (zx2_P != 0.0) && (zx3_P != 0.0) ) {
	  twin[atual].zxjitter_P =
	   200.0*fabs( (zx3_P - zx2_P) - (zx2_P - zx1_P) )/(zx3_P - zx1_P);
	  }
	  else {
	    twin[atual].zxjitter_P = 0.0000001;
	  }
	break;

	case MIN:
	  if( (zx1N_ != 0.0) && (zx2N_ != 0.0) && (zx3N_ != 0.0) ) {
	  twin[atual].zxjitterN_ =
	   200.0*fabs( (zx3N_ - zx2N_) - (zx2N_ - zx1N_) )/(zx3N_ - zx1N_);
	  }
	  else {
	    twin[atual].zxjitterN_ = 0.0000001;
	  }

	  if( (zx1_N != 0.0) && (zx2_N != 0.0) && (zx3_N != 0.0) ) {
	  twin[atual].zxjitter_N =
	   200.0*fabs( (zx3_N - zx2_N) - (zx2_N - zx1_N) )/(zx3_N - zx1_N);
	  }
	  else {
	    twin[atual].zxjitter_N = 0.0000001;
	  }
	break;
      }

      twin[atual].shimmer_peak =  200.0*fabs(x1 - x2)/fabs(x1 + x2);

      twin[atual].shimmer_rms =
	    200.0*fabs( sqrt(E12)/(t2-t1) - sqrt(E23)/(t3-t2) )/
		      ( sqrt(E12)/(t2-t1) + sqrt(E23)/(t3-t2) );

      /*if(hnr>0) twin[atual].hnr = (float) 10.0*log10(hnr);
      else twin[atual].hnr = 1e20;
*/
	if(hnr>0) twin[atual].hnr = (float) (hnr);
      	else twin[atual].hnr = 1e20;


      /*-----------------------------------------------------*/

      if(append) { /* vai apendar a uma cadeia formada (que nao passou,
		      necessariamente, da condicao pre'natal */
	twin[atual].link = precede; /* aponta p/ elem. anterior na cadeia */
	twin[atual].ntotal = twin[precede].ntotal+1; /* incr. contador */
	twin[atual].dtotal = twin[precede].dtotal + twin[precede].dist3;
			     /* incrementa medida do comprimento da cadeia */
	twin[precede].rightend = FALSE; /* elem. anterior nao e' mais a
				extremidade da cadeia */

	if(alivept==precede) /* se a cadeia ja' nasceu, desloca o apontador */
	  alivept=atual;
      }
      else { /* senao, inicia nova cadeia (ainda em condicao pre'natal) */
	twin[atual].link = NO_CHAIN;
	twin[atual].ntotal = 1;
	twin[atual].dtotal = t3-t1;
      }
      break;
    }
  }
}


int output()
{
  unsigned long int output_buf[100];  /* usado para ordenar a saida de PG */
  int cnt = 0;         /* contador para output_buf[] */
  char tmp[200];
  float measure;  /* for time series */
  long int time_offset = 0L;	/* to align with zero crossings */

  /* sera' dada a saida de uma cadeia, sempre que um novo PG fizer c/ que
     ele atinja a condicao de nascimento (birth_ready = TRUE ) */

  if(append) { /* o PG foi apendado a uma cadeia */
    if(alivept==atual) { /* a cadeia ja' esta' viva (alive) */
      twin[atual].link = NO_CHAIN;
      twin[precede].empty = YES;
      output_buf[cnt++] = atual;
    }
    else { /* a cadeia `a qual o PG foi apendado nao esta' viva */
      /* verifica se a cadeia pode ser recem-nascida */
      if((twin[atual].dtotal >= par.outp_BirthL) &&
	 (twin[atual].ntotal >= par.outp_BirthN) &&
	 (twin[atual].dtotal >= par.outp_K23*(t3-t2)))
	birth_ready=TRUE;
      else birth_ready = FALSE;

      if(birth_ready) {

	/* a cadeia e' recem-nascida e e' candidata a alive */

	if(alivept == NO_CHAIN) { /* nao ha' cadeia viva competindo */

	  /* verifica se ha' qualquer outra cadeia com periodo menor */
	  inhibit = FALSE;

	  for(q=0; q<TWIN_LEN; q++) {

/* obs: fabio used ... (twin[atual].dur23 > 1.5*twin[q].dur23) */
/***** original test **********/
	    if((twin[q].empty == NO) && (twin[q].rightend) &&
		(twin[atual].dur23 > twin[q].dur23) ){
	      inhibit = TRUE;
	      break;
	    }
/******************************/

/****** my test **********
	    if((twin[q].empty == NO) && (twin[q].rightend) &&
		(twin[q].ntotal > twin[atual].ntotal) ){
	      razao = (float) twin[atual].dur23/twin[q].dur23;

	      if( (razao<=0.6 && razao >= 0.4) || (razao<=2.2 && razao>=1.8) ) {
		inhibit=TRUE;
		break;
	      }
	    }
*******************************/

	  }
	  if(!inhibit) { /* nao ha' cadeia com periodo menor */
	    /* a cadeia do novo PG e' considerada alive */
	    alivept=atual;

	    strcat(fpwt_buf, "-1 -1\n");

	    /* SAIDA DO PG ATUAL (da cadeia alive e recem-nascida) */
	    output_buf[cnt++]=atual;

	    /* SAIDA DOS PG ANTERIORES (da cadeia recem-nascida) */
	    r=twin[atual].link;
	    while(r!=NO_CHAIN) {
	      output_buf[cnt++] = r;
	      twin[r].empty = YES;
	      m=r;
	      r=twin[r].link;
	    }
	    twin[atual].link =  NO_CHAIN;
	  }
	}
	else { /* ha' uma cadeia alive competindo */

	  /* verifica se o novo PG pode ser apendado `a cadeia  */

/* obs: fabio used ... (twin[atual].dur23 > 1.5*twin[q].dur23) */
/********* original test *******/
	  if(twin[atual].dur23 > twin[alivept].dur23) {
/*******************************/

/*******  my test *********
          razao = (float) twin[atual].dur23/twin[alivept].dur23;

	  if( (razao<=0.6 && razao >= 0.4) || (razao<=2.2 && razao>=1.8) ) {
**************************/

	    /* nao pode: libera elemento atual e sua cadeia */
	    twin[atual].empty = YES;
	    r = twin[atual].link;

	    /* libera precedentes */
	    while(r != NO_CHAIN) {
	      twin[r].empty=YES;
	      r=twin[r].link;
	    }
	  }
	  else {
	    /* inclui PG atual na cadeia viva ja' existente */
	    twin[alivept].empty = YES;
	    alivept=atual;

	    /* SAIDA DO PG ATUAL (apendado a cadeia viva existente) */

	    output_buf[cnt++] = atual;

	    /* libera PG atual no buffer */
	    r=twin[atual].link;

	    /* SAIDA DE PG's ANTERIORES */
	    while(r!=NO_CHAIN) {
	      output_buf[cnt++] = r;
	      twin[r].empty=YES; /* libera */
	      r=twin[r].link; /* aponta para anterior */
	    }

	    twin[atual].link = NO_CHAIN;
	  }
	}
      }
    }

    /* descarrega output_buf[] (ordenador) */
    if(cnt) {
      if( (fpwt_buf_cnt += (30*cnt)) > 400) { /* write buffer to file */
	fprintf(fpwt, fpwt_buf);
	fpwt_buf_cnt = 2*cnt;
	fpwt_buf[0] = '\0';
      }

      /* 20 jul 95 */
      for(r = cnt-1; r>=0; r--) { /* fill buffer */

	/* OUTPUT VALUES */
	outvalue.f0_intpeak = (float) header.nSamplesPerSec/
			    (float) twin[output_buf[r]].int_dur23;
	outvalue.f0_peak = (float) header.nSamplesPerSec/
			    (float) twin[output_buf[r]].dur23;
	outvalue.shimmer_peak = (float) twin[output_buf[r]].shimmer_peak;
	outvalue.shimmer_rms = (float) twin[output_buf[r]].shimmer_rms;

	if(twin[output_buf[r]].type == MAX) {
	  outvalue.jitterP_ = (float) twin[output_buf[r]].zxjitterP_;
	  outvalue.jitter_P = (float) twin[output_buf[r]].zxjitter_P;
	  outvalue.jitterN_ = outvalue.jitter_N = 0.0000001;
	}
	else if(twin[output_buf[r]].type == MIN) {
	  outvalue.jitterN_ = (float) twin[output_buf[r]].zxjitterN_;
	  outvalue.jitter_N = (float) twin[output_buf[r]].zxjitter_N;
	  outvalue.jitterP_ = outvalue.jitter_P = 0.0000001;
	}

	outvalue.hnr_comb = (float) twin[output_buf[r]].hnr;


	if(par.pitchdetection_alg == IntPEAKS) {
	  measure = outvalue.f0_intpeak;
	  time_offset = 0L;
	}

	else if(par.pitchdetection_alg == PEAKS) {
	  measure = outvalue.f0_peak;
	  time_offset = 0L;
	}

	else if(par.jitter_alg == NZP) {
	  if(twin[output_buf[r]].type == MIN)
	    measure = outvalue.jitterN_;
	  else
	    measure = outvalue.jitter_P;

	  time_offset = (long int) twin[output_buf[r]].int_dur23;
	}

	else if(par.jitter_alg == PZN) {
	  if(twin[output_buf[r]].type == MAX)
	    measure = outvalue.jitterP_;
	  else
	    measure = outvalue.jitter_N;

	  time_offset = (long int) twin[output_buf[r]].int_dur23;
	}

	else if(par.jitter_alg == PEAKS) {
	  measure = 200.0*
	    fabs(twin[output_buf[r]].dur12 - twin[output_buf[r]].dur23)/
		   ((float) twin[output_buf[r]].dur12 +
					  (float) twin[output_buf[r]].dur23);
	  time_offset = 0L;
	}

	else if(par.jitter_alg == IntPEAKS) {
	  measure = 200.0*
	    fabs(twin[output_buf[r]].int_dur12 - twin[output_buf[r]].int_dur23)/
	     (twin[output_buf[r]].int_dur12 + twin[output_buf[r]].int_dur23);
	  time_offset = 0L;
	}

	else if(par.shimmer_alg == 2) {
	  measure = twin[output_buf[r]].shimmer_rms;
	  time_offset = 0L;
	}

	else if(par.shimmer_alg == 1) {
	  measure = twin[output_buf[r]].shimmer_peak;
	  time_offset = 0L;
	}

	else if(hnr_arg != -1) {
	  measure = twin[output_buf[r]].hnr;
	  time_offset = 0L;
	}

	if(measure == 0.0) {
	  measure = 0.0000001;
	}

	sprintf(tmp, "%11.7f %11.7f\n",
	  (float) ((long int) twin[output_buf[r]].frame*Lframe -
	      Lframe + twin[output_buf[r]].t3 -
		twin[output_buf[r]].dur12 + time_offset)/header.nSamplesPerSec,
	  (float) twin[output_buf[r]].type*measure);

	strcat(fpwt_buf, tmp);

	statistics(0, twin[output_buf[r]].type, cum_flag, 0.0);

      }

      return TRUE; /* a twin period was found and linked; stop search */

    }
  }

  return FALSE;
}


void overlap(void)
{
  register int count;

  /* 20 jul 95 */
  for(count=0; count<Lframe+2*par.N; count++) {
    x[count] = x[Lframe+count];
  }

  /* desloca tempos do buffer de ptos. signficativos */
  for(i=0; i<SIG_LEN; i++) {
     sig[i].time -= Lframe;
     sig[i].interp_time -= Lframe;
     if(sig[i].zxP_ != 0.0 ) {
       sig[i].zxP_ -= Lframe;
     }
     if(sig[i].zx_P != 0.0 ) {
       sig[i].zx_P -= Lframe;
     }
     if(sig[i].zxN_ != 0.0 ) {
       sig[i].zxN_ -= Lframe;
     }
     if(sig[i].zx_N != 0.0 ) {
       sig[i].zx_N -= Lframe;
     }
  }

  t3old-=Lframe;

  Nframe++;

}

void postprocessor(void)
/*
  - make time correction on the .tmp file;
  - fill small gaps in pitch tracking;
  - depending on the parameters from cmd line:
    . run a non-linear filtering (3 point median filter) and/or
    . run a linear smoothing (3 point hanning window).
*/
{
  FILE *CorrectedFile;
  FILE *OldFile;

  int i, j, to, from,		/* counters */
      last_max,                 /* dynamic pointer to the end of p[] */
      terminate = FALSE,
      bufok;			/* informe whether buf[] contains -1 values */

  float f0ratio,
	buf[3],
	tmp;

  char x[15], y[15];

  struct {
    float t,
	  f0;
  }p[10];

  if((CorrectedFile = fopen(OutFile, "w"))==NULL) {
    printf("Error while opening time-corrected file\n");
    exit(1);
  }


  if((OldFile = fopen("tmp.tmp", "r"))==NULL) {
    printf("Error while opening .tmp file\n");
    exit(1);
  }

  for(i=0; i<10; i++) {
    p[i].f0 = 0.0;
  }

  /* read 1st line (.wav duration) */
  fscanf(OldFile, "%s%s", x, y);

  /* write to new file */
  fprintf(CorrectedFile, "0.0 "); fprintf(CorrectedFile, y);
  fprintf(CorrectedFile, "\n");

  /* initialize first 9 positions of buffer p[] */
  for(i=0; i<9; i++) {
    if(fscanf(OldFile, "%s%s\n", x, y)==2) {
      p[i].t = (float) atof(x);
      p[i].f0 = (float) atof(y);
      last_max = i+1;
    }
    else {
      terminate = TRUE;
      break;
    }
  }

  do {
    /* fill the remaining positions of p[] from last_max to 9 */
    for(i=last_max; i<10; i++){
      if(fscanf(OldFile, "%s%s\n", x, y)==2){
	p[i].t = (float) atof(x);
	p[i].f0 = (float) atof(y);
	last_max = i+1;
      }
      else {
	terminate = TRUE;
	break;
      }
    }

    /* TIME CORRECTION */
    for(i=1; i<last_max; i++) {
      if(p[i].t == -1)  continue;
      for(j=0; j<i; j++) {
	if(p[j].t == -1) continue;
	if(p[i].t <= p[j].t) {
	  printf("shift (%f)\n", p[i].t);
	  for(to=j, from=i; from<last_max; to++, from++) {
	    p[to].t  = p[from].t;
	    p[to].f0 = p[from].f0;
	  }
	  last_max = last_max - (i-j);
	  i = j;
	}
      }
    }

    /* FILL "BREAKS" UP TO par.gap: create rules for other time series */
    if(last_max>2) {
      for(i=1; i<last_max-1; i++) {

	/* test for buffer extremities and size of the time gap */
	if( (p[i].t == -1) && ((p[i+1].t - p[i-1].t) < par.gap) ) {

	  /* test similarity of f0 values */
	  f0ratio = fabs(p[i+1].f0/p[i-1].f0);

	  if( (0.8<f0ratio) && (f0ratio<1.2) ) {
	    printf("Fill gap (%3.0f ms) at %f seconds\n",
			    1000*(p[i+1].t - p[i-1].t), p[i+1].t);
	    /* shift buffer */
	    for(to=i, from=i+1; from<last_max; to++, from++) {
	      p[to].t  = p[from].t;
	      p[to].f0 = p[from].f0;
	    }
	    last_max--;
	  }
	}
      }
    }


    /* small buffer for median filtering */
    bufok = TRUE;
    for(i=0; i<3; i++) {
      buf[i] = p[i].f0;
      if(buf[i]==-1) {
	bufok = FALSE;
	break;
      }
    }

    /* NON LINEAR MEDIAN FILTER (length = 3) */
    /* note: this does not run for all of the *last* buffer before eof */
    if((medfilter_arg==1) && bufok) {
      for(i=0; i<2; i++) {
	for(j=i+1; j<3; j++) {
	  if(fabs(buf[j]) < fabs(buf[i])) {
	    tmp = buf[i];
	    buf[i] = buf[j];
	    buf[j] = tmp;
	  }
	}
      }
      /* takes the median value */
      p[1].f0 = buf[1];
    }

    /* LINEAR SMOOTHING (Hanning window, length = 3) */
    if((linsmooth_arg==1) && bufok) {
      tmp = 0.25*fabs(p[0].f0) + 0.50*fabs(p[1].f0) + 0.25*fabs(p[2].f0);
      if(p[2].f0 < 0) {
	p[2].f0 = -tmp;
      }
      else {
	p[2].f0 = tmp;
      }
    }

    fprintf(CorrectedFile, "%11.7f %11.7f\n", p[0].t, p[0].f0);
    statistics(2, 0, 0, tmp);

    /* voicing detection */
    if(InitOfTrack == 0.0 && p[0].t == 0.0) {
      if(p[1].t == -1 && p[2].t != -1) InitOfTrack = p[2].t;
      if(p[1].t != -1) InitOfTrack = p[1].t;
    }

    if(p[0].t != -1 && p[1].t == -1 && p[2].t != -1) {
      UnvoicedTime += p[2].t - p[0].t;
      pitchbreaks ++;
    }

    /* shift buffer */
    for(i=0; i<last_max-1; i++) {
      p[i].t = p[i+1].t;
      p[i].f0 = p[i+1].f0;
    }
    last_max--;

  } while(!terminate);

  /* write remaining values of the buffer into the time-corrected file */
  /* note: 1st position is not written in the last run of do{}while */
  for(i=0; i<last_max; i++) {
    fprintf(CorrectedFile, "%11.7f %11.7f\n", p[i].t, p[i].f0);
  }

  fclose(CorrectedFile);
  fclose(OldFile);

}

/*
float timedif(struct tm &t1, struct tm &t2)
{
  return (float)
    ( t2->tm_hour*360000L + t2->tm_min*6000L + t2->tm_sec*100L 
    - t1->tm_hour*360000L - t1->tm_min*6000L - t1->tm_sec*100L)/100.0;

}
*/

int sign (float value)
{
  if(value<0) return -1;
  else if(value>0) return 1;
  else return 0;
}


float FindZx(int polarity, int pattern)

/*
   "int i" is a global variable with the current sample

   polarity = 1, pattern=PZN  -> positive peak, zx ahead
	      1          NZP  -> positive peak, zx behind
	     -1          PZN  -> negative peak, zx behind
	     -1          NZP  -> negative peak, zx ahead

		       1                     1
	              *                     *
	            *   *                 *   *
	          *       *            *
	     �������-*�����*������-
	                  /  *      *  \
	                /     *   *      \
	      (pattern = PZN    *        (pattern = NZP
		search this      -1        search this zero crossing)
		zero crossing

*/

{
  int k, found;
  float delta;

  found = 0;
  if(((polarity>0) && (pattern == PZN)) ||
     ((polarity<0) && (pattern == NZP))) { /* look ahead */
    for(k=i+1; k < (int) (i+milisec1*1.3); k++) {
      if(polarity==1) { /* down going zero crossing */
	if(x[k]>=0.0 && x[k+1]<0.0) {
	  found = 1;
	  break;
	}
      }
      else { /* up going zero crossing */
	if(x[k]<=0.0 && x[k+1]>0.0) {
	  found = 1;
	  break;
	}
      }
    }
    if(found) {
      delta = (x[k] / ( (float) x[k] - x[k+1]));
      return (float) k + delta;
    }
    else {
      return 0.0;
    }
  }
  else if (((polarity<0) && (pattern == PZN)) ||
     ((polarity>0) && (pattern == NZP))) { /* look behind */
    for(k=i-1; k > (int) (i-milisec1*1.3); k--) {
      if(polarity == -1) { /* up going zero crossing */
	if(x[k]<=0.0 && x[k-1]>0.0) {
	  found = 1;
	  break;
	}
      }
      else { /* down going zero crossing */
	if(x[k]>=0.0 && x[k-1]<0.0) {
	  found = 1;
	  break;
	}
      }
    }

    if(found) {
      delta = (float)  (x[k] / ( (float) x[k] - x[k-1]));
      return k - delta;
    }
    else {
      return 0.0;
    }
  }
}


void GlottalNoise(int index1, int index2, int algorithm)
/*
  Global variables: t3, X[], Y[]
  algorithm was used to call dtw, wiener, etc (excluded here)
*/
{
  int N;

  /* take N as the smallest period */
  if( (index2-index1) <= (t3-index2) ) {
    N = (int) (index2 - index1 + 1);
  }
  else {
    N = (int) (t3 - index2 + 1);
  }
   /* original 0.8*N */
  comb2(index1, index2, (int) (0.7*N) );
}



/*
  Maximises the SNR
  between x[t1]...x[t2] and x[t2+offset]...x[t2+N+offset]
*/
#define MAXLAG 10
void comb2(index1,index2, N)
{
  int i, z, offset;
  float SNR[2*MAXLAG + 1],
	noise,
	signal,
	max;

  for(z = -MAXLAG; z <= MAXLAG; z++) {
    noise = signal = 0.0;
    for(i=0; i < N + z; i++) {
      noise += ( (float) x[index1+i] - (float) x[index2+i+z] )*
	       ( (float) x[index1+i] - (float) x[index2+i+z] );
      signal += ( (float) x[index1+i] + (float) x[index2+i+z] )*
		( (float) x[index1+i] + (float) x[index2+i+z] );
    }

    if(noise != 0.0) SNR[z+MAXLAG] = signal/noise;
    else SNR[z+MAXLAG] = 1e20;
  }

  max = SNR[0]; offset = 0;
  for(i=1; i < 2*MAXLAG+1; i++) {
    if(SNR[i] > max) {
      max = SNR[i];
      offset = i;
    }
  }

  hnr = SNR[offset];
}



/*
  calculate energy in t12 and t23 for shimmer based on rms values
  shimmer is calculated during the "installation" of the twin-period
  float E12, E23; int x[], t2, t3 = global variables
*/

void shimmerms(index1, index2)
/*
 t3 is a global variable
*/
int index1, index2;
{
  int z;

  E12 = 0.0;
  for(z=index1; z< index2; z++)
    E12 += ( (float) x[z]*x[z]);

  E23 = 0.0;
  for(z=index2; z < (int) t3; z++)
    E23 += ( (float) x[z]*x[z]);
}


int outofmemory()
{
  printf("malloc failure\n");
  exit(1);
}

void statistics(int code, int polarity, int cummulative, float FiltF0)
/*
  code: 0 = incremental; 1 = end; 2 = special for smoothed f0 (FiltF0)
  polarity of twin period: MAX (1, pos) or MIN (-1, neg)
  cummulative: 1/0 = write/don't write do cummulative file
*/
{
  char str[60];
  FILE *stat = NULL;
  float dur, UnvoicedPct, ratioPZN, ratioNZP, ratio;

  static struct STATISTICS
	      F0 = {0, 0.0, 0.0, 0.0, 0.0},
	      SmoothF0 = {0, 0.0, 0.0, 0.0, 0.0},

	      /* PZN */
	      jitterP_ = {0, 0.0, 0.0, 0.0, 0.0},
	      jitter_N = {0, 0.0, 0.0, 0.0, 0.0},

	      /*NZP*/
	      jitterN_ = {0, 0.0, 0.0, 0.0, 0.0},
	      jitter_P = {0, 0.0, 0.0, 0.0, 0.0},

	      jitter = {0, 0.0, 0.0, 0.0, 0.0},
	      shimmerRMS = {0, 0.0, 0.0, 0.0, 0.0},
	      shimmerPeak = {0, 0.0, 0.0, 0.0, 0.0},
	      snr = {0, 0.0, 0.0, 0.0, 0.0};

  switch(code) {
    case 1: /* end statistics */
      if(F0.n) {
	F0.mean = F0.SumX / (double) F0.n;
	F0.StdDev = sqrt(
		   (F0.SumX2 - F0.n*F0.mean*F0.mean)/F0.n );
      }

      if(SmoothF0.n) {
	SmoothF0.mean = SmoothF0.SumX / (double) SmoothF0.n;
	SmoothF0.StdDev =
	  sqrt((SmoothF0.SumX2 - SmoothF0.n*SmoothF0.mean*SmoothF0.mean)/SmoothF0.n);
      }

      if(jitterP_.n) {
	jitterP_.mean = jitterP_.SumX/ (double) jitterP_.n;
	jitterP_.StdDev = sqrt(
	  (jitterP_.SumX2 - jitterP_.n*jitterP_.mean*jitterP_.mean)/jitterP_.n);
      }

      if(jitterN_.n) {
	jitterN_.mean = jitterN_.SumX/ (double) jitterN_.n;
	jitterN_.StdDev = sqrt(
	  (jitterN_.SumX2 - jitterN_.n*jitterN_.mean*jitterN_.mean)/jitterN_.n);
      }

      if(jitter_N.n) {
	jitter_N.mean = jitter_N.SumX/ (double) jitter_N.n;
	jitter_N.StdDev = sqrt(
	  (jitter_N.SumX2 - jitter_N.n*jitter_N.mean*jitter_N.mean)/jitter_N.n);
      }

      if(jitter_P.n) {
	jitter_P.mean = jitter_P.SumX/ (double) jitter_P.n;
	jitter_P.StdDev = sqrt(
	  (jitter_P.SumX2 - jitter_P.n*jitter_P.mean*jitter_P.mean)/jitter_P.n);
      }

      if(shimmerRMS.n) {
	shimmerRMS.mean = shimmerRMS.SumX/ (double) shimmerRMS.n;
	shimmerRMS.StdDev = sqrt(
	  (shimmerRMS.SumX2 - shimmerRMS.n*shimmerRMS.mean*shimmerRMS.mean)/
	  shimmerRMS.n);
      }

      if(shimmerPeak.n) {
	shimmerPeak.mean = shimmerPeak.SumX/ (double) shimmerPeak.n;
	shimmerPeak.StdDev = sqrt(
	  (shimmerPeak.SumX2 - shimmerPeak.n*shimmerPeak.mean*shimmerPeak.mean)/
	  shimmerPeak.n);
      }

      if(snr.n) {
	snr.mean = 10*log10(snr.SumX/ (double) snr.n);
	snr.StdDev =  sqrt((snr.SumX2 - snr.n*snr.mean*snr.mean)/snr.n);
      }


      /* OPEN *.stt FILE */
      strcpy(str, OutFile);
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

      dur = (float) header.datasize/(2.0 * header.nSamplesPerSec);
      if(dur > InitOfTrack) {
	UnvoicedPct = 100.0*UnvoicedTime/
		     (dur -InitOfTrack);
      }
      else {
	UnvoicedPct = 0.0;
      }

      /* choose between PZN and NZP jitter */
      ratioPZN = ratioNZP = 0.0;
      if(jitterP_.mean != 0.0 && jitter_N.mean != 0.0) {
	ratioPZN = (float)fminf( (float)(jitterP_.mean/jitter_N.mean),
		       (float)(jitter_N.mean/jitterP_.mean));
      }

      if(jitterN_.mean != 0.0 && jitter_P.mean != 0.0) {
	ratioNZP = (float)fminf( (float)(jitterN_.mean/jitter_P.mean),
		       (float)(jitter_P.mean/jitterN_.mean));
      }

      if(ratioPZN >= 0.8 && ratioPZN >= ratioNZP) {
	jitter.mean = (jitterP_.mean + jitter_N.mean)/2.0;
	jitter.StdDev = (jitterP_.StdDev + jitter_N.StdDev)/2.0;
	ratio = ratioPZN;
	strcpy(str, "PZN ");
      }
      else if(ratioNZP >= 0.8 && ratioNZP > ratioPZN) {
	jitter.mean = (jitterN_.mean + jitter_P.mean)/2.0;
	jitter.StdDev = (jitterN_.StdDev + jitter_P.StdDev)/2.0;
	ratio = ratioNZP;
	strcpy(str, "NZP ");
      }
      else {
	ratio = 0.0;
	strcpy(str, "NONE");
      }

      printf("\n Statistics in file %s\n", str);
      printf("------------------------------------------\n");
      printf(" MEASURE               MEAN       STD DEV\n");
      printf(" Dur          (sec)   %5.2f\n", dur);
      printf(" Unvoiced     (%%)     %5.2f\n", UnvoicedPct);
      printf(" Pitch Breaks (#)      %i\n", pitchbreaks);
      printf(" F0           (Hz)   %5.2f        %5.2f\n", F0.mean, F0.StdDev);
      printf(" SmoothF0     (Hz)   %5.2f        %5.2f\n", SmoothF0.mean, SmoothF0.StdDev);

      printf(" jitterP_     (%%)     %5.2f        %5.2f\n", jitterP_.mean, jitterP_.StdDev);
      printf(" jitter_N     (%%)     %5.2f        %5.2f\n", jitter_N.mean, jitter_N.StdDev);
      printf(" jitterN_     (%%)     %5.2f        %5.2f\n", jitterN_.mean, jitterN_.StdDev);
      printf(" jitter_P     (%%)     %5.2f        %5.2f\n", jitter_P.mean, jitter_P.StdDev);

      printf(" JITTER       (%%)     %5.2f        %5.2f\n", jitter.mean, jitter.StdDev);
      printf(" ratio        (%%)     %5.2f         %s\n", ratio, str);
      printf(" shimmerRMS   (%%)     %5.2f        %5.2f\n", shimmerRMS.mean, shimmerRMS.StdDev);
      printf(" shimmerPeak  (%%)     %5.2f        %5.2f\n", shimmerPeak.mean, shimmerPeak.StdDev);
      printf(" snr          (dB)    %5.2f        %5.2f\n", snr.mean, snr.StdDev);
      printf("------------------------------------------\n");

      printf("Cummulative file: %s\n", CumFile);

      fprintf(stat, " MEASURE                MEAN       STD DEV\n");
      fprintf(stat, " Dur          (sec)    %5.2f\n", dur);
      fprintf(stat, " Unvoiced     (%%)      %5.2f\n", UnvoicedPct);
      fprintf(stat, " Pitch Breaks (#)       %i\n", pitchbreaks);
      fprintf(stat, " F0           (Hz)    %5.2f       %5.2f\n", F0.mean, F0.StdDev);
      fprintf(stat, " SmoothF0     (Hz)    %5.2f       %5.2f\n", SmoothF0.mean, SmoothF0.StdDev);

      fprintf(stat, " jitterP_     (%%)      %5.2f       %5.2f\n", jitterP_.mean, jitterP_.StdDev);
      fprintf(stat, " jitter_N     (%%)      %5.2f       %5.2f\n", jitter_N.mean, jitter_N.StdDev);
      fprintf(stat, " jitterN_     (%%)      %5.2f       %5.2f\n", jitterN_.mean, jitterN_.StdDev);
      fprintf(stat, " jitter_P     (%%)      %5.2f       %5.2f\n", jitter_P.mean, jitter_P.StdDev);

      fprintf(stat, " jitter       (%%)      %5.2f       %5.2f\n", jitter.mean, jitter.StdDev);
      fprintf(stat, " ratio        (%%)      %5.2f         %s\n", ratio, str);

      fprintf(stat, " shimmerRMS   (%%)      %5.2f       %5.2f\n", shimmerRMS.mean, shimmerRMS.StdDev);
      fprintf(stat, " shimmerPeak  (%%)      %5.2f       %5.2f\n", shimmerPeak.mean, shimmerPeak.StdDev);
      fprintf(stat, " snr          (dB)     %5.2f       %5.2f\n", snr.mean, snr.StdDev);
      fclose (stat);

      if(cummulative != -1) {
	fprintf(fpcum, "%5.2f %5.2f %i %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %s %5.2f %5.2f %5.2f %5.2f %5.2f\n",
	  dur,
	  UnvoicedPct,           pitchbreaks,
	  F0.mean,               F0.StdDev,
	  SmoothF0.mean,	 SmoothF0.StdDev,
	  jitter.mean,        	 ratio,		 str,
	  jitter.StdDev,         shimmerRMS.mean, shimmerRMS.StdDev,
	  snr.mean,              snr.StdDev
	);
      }
    break;

    case 0: /* increment statistics */
      F0.n ++;
      F0.SumX += outvalue.f0_peak;
      F0.SumX2 += (outvalue.f0_peak*outvalue.f0_peak);

      switch(polarity) {
	case MAX:
	  if(outvalue.jitterP_ <= 10.0) {
	    jitterP_.n++;
	    jitterP_.SumX += outvalue.jitterP_;
	    jitterP_.SumX2 += (outvalue.jitterP_*outvalue.jitterP_);
	  }

	  if(outvalue.jitter_P <= 10.0) {
	    jitter_P.n++;
	    jitter_P.SumX += outvalue.jitter_P;
	    jitter_P.SumX2 += (outvalue.jitter_P*outvalue.jitter_P);
	  }

	break;

	case MIN:
	  if(fabs(outvalue.jitterN_) <= 10.0) {
	    jitterN_.n++;
	    jitterN_.SumX += (float) fabs(outvalue.jitterN_);
	    jitterN_.SumX2 += (outvalue.jitterN_*outvalue.jitterN_);
	  }

	  if(fabs(outvalue.jitter_N) <= 10.0) {
	    jitter_N.n++;
	    jitter_N.SumX += (float) fabs(outvalue.jitter_N);
	    jitter_N.SumX2 += (outvalue.jitter_N*outvalue.jitter_N);
	  }
	break;

	default:
	break;
      }

      shimmerRMS.n++;
      shimmerRMS.SumX += outvalue.shimmer_rms;
      shimmerRMS.SumX2 += (outvalue.shimmer_rms*outvalue.shimmer_rms);

      shimmerPeak.n++;
      shimmerPeak.SumX += outvalue.shimmer_peak;
      shimmerPeak.SumX2 += (outvalue.shimmer_peak*outvalue.shimmer_peak);

      if(outvalue.hnr_comb != 1e20) {
	snr.n++;
	snr.SumX += outvalue.hnr_comb;
	snr.SumX2 += (outvalue.hnr_comb*outvalue.hnr_comb);
      }
    break;

    case 2: /* special increment for Smoothed F0 */
      SmoothF0.n ++;
      SmoothF0.SumX += FiltF0;
      SmoothF0.SumX2 += (FiltF0*FiltF0);
    break;

    default:
    break;
  }

  return;
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



/*
float fmin(float a, float b)
{
  if( a <= b ) return a;
  else return b;
}
*/


/* Boxes
 ALT + ...  sigma = 228

	 218   196   191
  (180)  179         179  (195)
	 192   196   217

*/

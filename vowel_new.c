/*
  Cascade Formant Synthesiser
  Formant frequencies/bandwidths from Rabiner & Schafer (1978) (book)
*/

/* 
Original Implementation: Maurilio Nunes Vieira (1997)

Modified 9 Feb. 2007 - Joao SANSAO 

* Ported to gcc  

Main modifications: 
* random() calling --> Number between 0 and RAND_MAX, takes no parameter
* excluded borland libraries io.h, conio.h and stat.h
* wave header format 
* changed int to signed short type for the data array (x) 

* seems to work fine under GNU/Linux
* tested with gcc 4.1.2 20060928

* gcc compilation procedure :

$gcc -lm -Wall vowel.c -o vowel

*/
#include "stdio.h"
//#include "conio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

#define MAX_ORDER	40

/* FUNCTION PROTOTYPES */
void usage(),
     msg(void),
     initialization(char *argv[]),
     overlap(void),
     coefficients(double *B, double *A, int, int);
signed short round2int(double);


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

/* FILE RELATED */
FILE *fdr;				/* data file header */
signed short *x = NULL,				/* speech samples buffer */
    ni;					/* nb. of samples actually read */

FILE *fpwb;                             /* header for output file (*.wav) */
signed short *y = NULL;				/* buffer for output file */


/* TIME INTERVALS */
int milisec1,
    Lframe,			/* length of buffer */
    N;				/* overlap */

float pre_emphasis = 1.0,
      gain = 10.0; /* to accommodate range; for airflow in the 0-1000 range:
		   0.0 for oral flow, ~100 for oral pressure */

/* noise generation */
float snr, NoiseDistWidth, sig_power, noise_power, aux, noiseval;

void main(argc,argv)
int argc;
char *argv[];

{
  int i, j; /* used in main loop */

  double y_double[MAX_ORDER+1];		/* keep the old y values with max precision */

  int input_arg = -1,
      output_arg = -1,
      algorithm_arg = -1,
      noise_arg = -1,
      alg; /*tmp*/


  /* filter coefficients, H(z) = B(z)/A(z);*/

/**  BP5_8.CH2
      | | |  |
      | | |  +---> Chebyshev type II
      | | +------> Order = 8
      | +--------> Attenuation in Stop band = -5 dB
      +----------> Band Pass
***/

  double B[MAX_ORDER+1],
	 A[MAX_ORDER+1];

  int Order;


  /* CHECK PROGRAM CALL */
  if(argc<2)  usage();

  for(i=1; i<argc && *argv[i] == '-'; i++) {
    switch(argv[i][1]) {
      default:
	if (argc <= i+1) usage();

	switch(argv[i++][1]) {
	  case 'p':
	  case 'P':
	    pre_emphasis = atof(argv[i]);
	    if(pre_emphasis < 0.0 || pre_emphasis > 1.0) usage();
	    break;
	  case 'g':
	  case 'G':
	    gain = atof(argv[i]);
	    if(gain < 1) usage();
	    break;
	  case 'i':
	  case 'I':
	    input_arg = i;
	    break;
	  case 'n':
	  case 'N':
	    noise_arg = i;
	    snr = atof(argv[i]);
	    if(snr <= 0) usage();
	    else snr = pow(10, snr/10);
	    break;
	  case 'o':
	  case 'O':
	    output_arg = i;
	    break;
	  case 'v':
	  case 'V':
	    algorithm_arg = i;
	    j = (int) argv[i][0];
	    if(j == 'i' || j == 'a' || j == 'u' ||
	       j == 'I' || j == 'A' || j == 'U' || 
	       j== '1'|| j== '2'|| j== '3' || j== '4'
	       || j== '5'|| j== '6'|| j== '7') 
	      {
	      alg = j;

	      switch(alg) 
		{
		case 'a':
		case 'i':
		case 'u':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		  Order = 22;
		coefficients(B, A, Order, alg);
		break;
	      
		default:
		  break;
		}
	      }
	    else usage();

	    break;

	  default:
	    usage();
	    break;
	}
	break;
    }
  }
  if ((i != argc && *argv[i] != 'i') || input_arg == -1 ||
       algorithm_arg == -1) usage();

  /* OPEN SPEECH FILE */
  if((fdr = fopen(argv[input_arg], "rb")) !=  NULL) {
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
  if((fpwb = fopen(argv[output_arg], "wb"))==NULL) {
    printf("Error while creating file (%s)\n", argv[output_arg]);
    exit(1);
  }

  /* copy the header from the input file */
  fwrite(&header, sizeof(header), 1, fpwb);

  /* INITIALIZE VARIABLES, BUFFERS, ETC */
  for(i=0; i<Order+1; i++) {
    y_double[i] = 0.0;
  }
  initialization(argv);

  /* PUT MESSAGES ON THE SCREEN */
  msg();

  /*---------------------- MAIN LOOP ---------------------------------*/
  printf("pre_emphasis=%5.2f, gain=%5.2f, snr=%5.2f\n", pre_emphasis, gain, snr);
  printf("Wait...");

srandom(time(NULL)); /* iniciando o gerador de números aleatórios */


  while((ni = (fread(x+N, sizeof(signed short), Lframe, fdr))) != NULL) {
	       /* fill buffer for N <= x < Lframe+N */

   /*
   Buffers:                                        (ni)
   i=  0         N             ...                Lframe+N-1
      [oooooooooonnnnnnnnnnnnnn...nnnnnnnnnnnnnnnnn] = x[] (int)
      [ooooooooooffffffffffffff...fffffffffffffffff] = y[] (int)

      o=old samples (overlap), n=new samples, f=filtered samples
      Lframe >= ni (ni < Lframe only if eof was reached)
   */

   /*.......................... filters ...............................*/

   for(i=N; i<N+ni; i++) {

     switch(alg) {
       case 'a':
       case 'u':
       case 'i':
       case '1':
       case '2':
       case '3':
       case '4':
       case '5':
       case '6':
       case '7':
	 /* zeros */
	 y_double[0] = 0.0;
	 for(j=0; j<Order+1; j++) {
	   y_double[0] = y_double[0] + B[j]*x[i-j]*gain;
	 }

	 /* poles */
	 /* y_double[]:
	    +---+---+---+---+---+---+---+
	    | 6 | 5 | 4 | 3 | 2 | 1 | 0 |
	    +---+---+---+---+---+---+---+
	     n-6	 ...     n-1  n
	 */

	 for(j=1; j<Order+1; j++) {
	   y_double[0] = y_double[0] - A[j]*y_double[j];
	 }

	 /* pre-emphasis and rounding of output value */
	 y[i] = round2int(y_double[0] - pre_emphasis*y_double[1]);

	 /* shift y_double[] */
	 for(j=Order; j>0; j--) {
	   y_double[j] = y_double[j-1];
	 }
	 break;

       default:
	 break;
     }

   }
   /*.................................................................*/

   /*
     Noise generation
   */
   if(noise_arg != -1) {
     aux = 0.0;
     for(i=N; i<N+ni; i++) {
       aux += (float) y[i]*y[i];
     }
     sig_power = aux/(float) ni;

     NoiseDistWidth = sqrt(12*sig_power/snr);

     noise_power = 0.0;
     //randomize();
     
     for(i=N; i < N+ni; i++) {
       noiseval = (1.0*random())/RAND_MAX;
       aux =  NoiseDistWidth * (noiseval - 0.5);
       y[i] = round2int( 1.0 * y[i] + 1.0*aux);
       noise_power += aux*aux;
     }
     noise_power = noise_power/( (float) ni);
    
    // printf("SNRdb = %5.2f\n", 10.0*log10(sig_power/noise_power));
    
   }

   /* write processed samples into file */
   fwrite(y+N, sizeof(signed short), ni, fpwb);

   if(ni <= N)  break; /* eof */
   overlap();
  }

  free(x);
  free(y);
  fclose(fdr);
  fclose(fpwb);

  printf("done\n");

  exit(0);
}

/*---------------------------- FUNCTIONS ---------------------------*/

void usage()
{
  printf("lxfilter -i file1.wav  -o file2.wav [- arg {descr (defaults)}]\n");
  printf("-i   {input file  (.wav, pcm, mono, 16 bits/sample, 22050 Hz)}\n");
  printf("-o   {output file (.wav, pcm, mono, 16 bits/sample, 22050 Hz)}\n");
  printf("-v x {vowel: n = x, i, or u}\n");
  printf("-p x {pre_emphasis (0.0 <= x <= 1.0, 1.0}\n");
  printf("-g x {gain (x > 0.0, 10.0)}\n");
  printf("-n x {SNR ratio (white noise added to oral pressure), x > 0}\n");
  exit(0);
}

void initialization(char *argv[])
{
  int i;

  milisec1  = (int) (header.nSamplesPerSec * 0.001/2.0)*2;

  Lframe = 50*milisec1;

  N = 20;

  /* CREATE INPUT BUFFER */
  if((x = malloc( sizeof(signed short) * (N + Lframe))) == NULL) {
    printf("out of memory in call to malloc(x).\n");
    exit(1);
  }

  /* CREATE OUTPUT BUFFER */
  /* this may be of float type in other applications */
  if((y = malloc( sizeof(signed short) * (N + Lframe))) == NULL) {
    printf("out of memory in call to malloc(y).\n");
    exit(1);
  }

  /* INITIALIZE BUFFERS for 0 <= x < N */
  /* non zero initial conditions: */
  /* fread(x, sizeof(int), N, fdr); */

  /* zero initial conditions */
  for(i=0; i<N; i++) {
    x[i] = 0;
    y[i] = 0;
  }
}


void overlap(void)
{
  register int i;

  for(i=0; i<N; i++) {
    x[i] = x[Lframe+i];
    y[i] = y[Lframe+i];
  }
}



void msg(void)
{
  printf(" \nMaurilio N. Vieira, 28 mar 97. \n");
  printf(" Cascade Formant Synthesiser\n");
  printf(" Formant frequencies/bandwithds from Rabiner & Schafer (1978),\n");
  printf(" Digital Processing of Speech Signals, Prentice Hall, pp. 74-77\n");
}


signed short round2int(double x)
{
  double dec;


  dec = x - floor(x);
  if(dec > 0.5) {
    x = x + 1;
  }

  if(x>32767) x = 32767;
  else if(x<-32767) x = -32767;

  return (signed short) floor(x);
}


void coefficients(double *B, double *A, int order, int vowel)
{
  int i;

 /* /a/ vowel (Rabiner & Schafer) */
 double B_a[22+1] = {
    1.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,  
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0, 
    0.0, 0.0 
  };

 double A_a[22+1] = {  1.0000,   -1.0834,   -0.0633,   -0.0748,    0.9822,   -0.9253,    0.0508, -0.0586,    1.2470,   -0.6341,   -0.2670,   -0.3795,    0.6916,   -0.3385, 0.1811,   -0.1897,    0.4624,   -0.3500,    0.3861,   -0.3103,   -0.0433,  -0.1641,    0.3319};


 /* /i/ vowel (Rabiner & Schafer) */
  double B_i[22+1] = {
    1.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,  
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0, 
    0.0, 0.0 
  };

  double A_i[22+1] = {  
    1.0000,   -0.5402,    0.5386,   -0.8206,    1.8397,   -1.6048,    0.9437,
    -1.9194,    2.3815,   -1.9621,    1.1651,   -2.0608,    2.0610,   -1.5520,
     0.9860,   -1.1767,    1.1736,   -0.5621,    0.5424,   -0.3210,    0.4121,
    -0.0800,    0.1109
 };

 /* /u/ vowel (Rabiner & Schafer) */
  double B_u[22+1] = {
    1.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,  
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0, 
    0.0, 0.0 
  };


  double A_u[22+1] = {
    1.0000,   -0.7191,   -0.0635,   -0.7478,    0.9152,   -0.8705,    0.6543,
   -0.7978,    1.1919,   -0.7942,    0.9960,   -0.8106,    0.7699,   -0.8660,
    0.4908,   -0.5704,    0.4488,   -0.2506,    0.2409,   -0.2666,    0.2400,
   -0.1351,    0.0990
  };


  /* double A_zz1[22+1] = { */
  /*   1,-2.8612,3.1841,-2.0855,1.4903,-1.0136,0.20456,-0.18686,0.98001,-1.2302,1.041,-1.1055,0.72005,0.17129,-0.33227,-0.16355,0.41607,-0.68266,0.8933,-0.5762,0.20119,-0.06272,0.0037288}; */

  /* double A_zz2[22+1] = { */
  /*   1,-2.5745,2.343,-0.94085,0.70452,-0.97762,0.45166,-0.25279,0.69516,-0.70511,0.53763,-0.72364,0.55629,0.10976,-0.17062,-0.04305,-0.1669,0.14387,0.30413,-0.4663,0.33867,-0.28592,0.12997}; */

  /* double A_zz3[22+1] = { */
  /*   1,-2.3339,2.0297,-0.69428,0.073279,-0.44429,0.48562,-0.2367,0.21391,-0.44347,0.62499,-0.51608,0.2139,0.23154,-0.057067,-0.36845,0.38443,-0.094472,-0.084909,0.16286,-0.18416,-0.048406,0.090362}; */

  /* double A_zz4[22+1] = { */
  /*   1,-1.9594,1.3269,-0.44981,0.10627,-0.11209,-0.063338,-0.045392,0.26162,-0.10648,-0.024688,0.031364,-0.080296,0.27955,0.080201,-0.30531,0.13132,0.0017249,-0.11225,0.1998,-0.19607,-0.10478,0.14557}; */


  /* double A_zz5[22+1] = { */
  /*   1,-2.3861,1.9264,-0.73108,0.51213,-0.64305,-0.065166,0.40315,0.59226,-0.53892,-0.24842,-0.21442,1.1287,-0.81991,0.12439,-0.39315,0.19686,0.48666,-0.37824,-0.17531,0.21037,0.19977,-0.18278}; */

  /* double A_zz6[22+1] = { */
  /*   1,-2.6091,2.2411,-0.31572,-0.55109,-0.034761,0.32006,0.076746,0.027201,-0.15949,-0.058595,0.088923,0.017173,0.064056,-0.12183,-0.10195,0.023969,0.24124,-0.13137,-0.15734,0.13959,0.14095,-0.13489}; */

  /* double A_zz7[22+1] = { */
  /*   1,-1.6893,0.49945,0.1387,0.025117,-0.18234,0.11312,-0.0070647,0.22668,0.11767,-0.10212,-0.12436,-0.01977,0.15953,-0.099512,-0.089738,-0.054842,-0.042349,-0.013237,0.11297,0.1187,0.0043178,-0.087704}; */



double A_zz1[22+1] = {
  1,-1.9209,1.3822,-0.79498,0.73304,-0.28773,-0.088081,-0.26752,0.73022,-0.57144,0.53328,-0.57742,0.11902,0.33303,-0.078531,-0.18882,0.24966,-0.46435,0.40522,-0.095661,-0.031051,0.13504,-0.08825};

double A_zz2[22+1] = {
1,-1.6188,0.78654,-0.16682,0.5267,-0.48262,0.018241,-0.25822,0.46733,-0.29664,0.29636,-0.44338,0.1071,0.22995,0.038005,0.0086797,-0.17205,0.0010315,0.27459,-0.18409,0.10782,-0.046035,-0.054049};

double A_zz3[22+1] = {
1,-1.3779,0.72802,-0.036978,0.063248,-0.39709,0.093791,-0.088007,0.066124,-0.37397,0.28998,-0.20462,-0.062905,0.27393,0.12884,-0.20283,0.14209,0.11943,-0.024631,0.1333,-0.14007,0.085528,-0.13335};

double A_zz4[22+1] = {
1,-0.98772,0.37199,-0.12661,0.024009,-0.1072,-0.17165,-0.1928,0.034513,-0.058549,-0.045453,-0.024005,-0.10138,0.18711,0.24479,-0.035571,0.088531,0.081769,-0.047952,0.17438,-0.10188,-0.00059236,-0.14424};

double A_zz5[22+1] = {
1,-1.4418,0.57161,-0.20901,0.32433,-0.32691,-0.3872,0.033461,0.63796,0.051528,-0.1644,-0.41906,0.75986,-0.11283,0.039019,-0.37524,-0.17132,0.33358,-0.045904,-0.23838,0.011798,0.15687,0.034466};

double A_zz6[22+1] = {
1,-1.6708,0.66938,0.32891,-0.25474,-0.28354,0.071261,0.14345,0.15291,-0.020042,-0.072386,0.020538,0.042754,0.10146,-0.036888,-0.13553,-0.096515,0.17165,0.026696,-0.17227,-0.041513,0.25834,-0.079834};

double A_zz7[22+1] = {
1,-0.74371,-0.20562,-0.057142,-0.027027,-0.20498,-0.079405,-0.081376,0.15025,0.25524,0.1427,0.0083873,-0.012021,0.1418,0.03303,-0.057638,-0.10839,-0.14086,-0.14557,-0.025487,0.10014,0.10982,0.0027533};
  


  switch(vowel) {
    case 'a':
      printf("vowel /a/ JPHS");
      for(i=0; i<order+1; i++) {
	B[i] = B_a[i];
	A[i] = A_a[i];
      }
      break;

     case 'i':
      printf("vowel /i/ JPHS");
      for(i=0; i<order+1; i++) {
	B[i] = B_i[i];
	A[i] = A_i[i];
      }
      break;

    case 'u':
      printf("vowel /u/ JPHS");
      for(i=0; i<order+1; i++) {
	B[i] = B_u[i];
	A[i] = A_u[i];
      }
      break;

    case '1':
      printf("vowel /a/ MNV");
      for(i=0; i<order+1; i++) {
	B[i] = B_a[i];
	A[i] = A_zz1[i];
      }
      break;

    case '2':
      printf("vowel /e/ MNV");
      for(i=0; i<order+1; i++) {
	B[i] = B_a[i];
	A[i] = A_zz2[i];
      }
      break;

    case '3':
      printf("vowel /e_/ MNV");
      for(i=0; i<order+1; i++) {
	B[i] = B_a[i];
	A[i] = A_zz3[i];
      }
      break;

    case '4':
      printf("vowel /i/ MNV");
      for(i=0; i<order+1; i++) {
	B[i] = B_a[i];
	A[i] = A_zz4[i];
      }
      break;

    case '5':
      printf("vowel /o_/ MNV");
      for(i=0; i<order+1; i++) {
	B[i] = B_a[i];
	A[i] = A_zz5[i];
      }
      break;

    case '6':
      printf("vowel /o/ MNV");
      for(i=0; i<order+1; i++) {
	B[i] = B_a[i];
	A[i] = A_zz6[i];
      }
      break;

    case '7':
      printf("vowel /u/ MNV");
      for(i=0; i<order+1; i++) {
	B[i] = B_a[i];
	A[i] = A_zz7[i];
      }
      break;

    default:
      usage();
      break;
  }
}




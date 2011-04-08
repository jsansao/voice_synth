/*
  Simulated airflow based on Fant (1979),
  Glottal Source and Excitation Analysis,
  STL-QPSR 1/1979, pp. 85-107.
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

$gcc -lm -Wall flowgen.c -o flowgen

*/


#include "stdio.h"
#include "fcntl.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "ctype.h"
#include "time.h"

#define PI 4.0*atan(1.0)

/* FUNCTION PROTOTYPES */
void usage(void),
     msg(char *argv[]),
     initialization(char *argv[]);
signed short truncate(float);


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
char WaveFile[30];
FILE *outfile;

signed short *x = NULL;				/* speech samples buffer */

struct PAR {
  float dur,		/* duration of the created file */
	jitter,
	cq,
	K,		/* speed of closure */
	Fg,		/* glottal formant; Fg = 1/(2.T2) */
	F0,		/* Fg > F0 */
	DC,		/* DC flow (% of max amplitude) */
	noise;		/* pow(10, arg.noise/10) */

  long fs;		/* sampling rate */
   int amp,             /* maximum amplitude */
       NoiseDistWidth;	/* noise = uniform(0,...,NoiseDistWidth) */
  float Kvar,Shimmer;           /* speed closure variation, Shimmer */
} par = {1.0, .0, 0.55, 0.65, 125, 120, 0.0, 0.0, 22050L, 12000, 0.0, 0.0,0.0};


struct ARG {
    int wav,
	dur,
	jitter,
	cq,
	K,
	Fg,
	F0,
	DC,
	noise,
	fs,
	amp, Kvar, Shimmer;
} arg = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,-1,-1};


void main(argc,argv)
int argc;
char *argv[];
{
  int i, j, k,			/* general use */
      P,			/* nominal Pitch period */
      T,			/* real Pitch period (jittered) */
      T2,			/* = 1/(w.Fg) */
      T3,			/* point in the closing phase where flow = 0 */
      T4,			/* instant where the rising flow = DCflow */
      w[500];			/* vector for additive white noise */
  unsigned long nSamples, CountSamples;
  float aux,
	x_pow,
	w_pow,
	J,                      /* random jitter */
	DeltaPer[2] = {0, 0};
  float S, DeltaShimmer[2] = {0, 0};

  char *p;			/* used to extract WaveFile from cmd line */


  /* CHECK PROGRAM CALL */
  if(argc<2)  usage();

  for(i=1; i<argc && *argv[i] == '-'; i++) {
    j = i+1;

    switch(argv[i][1]) {
      default:
	if (argc <= j) usage();

	switch(argv[i++][1]) {
	  case 'o':
	  case 'O':
	    arg.wav = i;
	    /* p points to the last occurrence of '\' on the path */
	    if( (p = strrchr(argv[arg.wav], '\\')) != NULL) {
	      strcpy(WaveFile, p+1);
	    }
	    else {
	      strcpy(WaveFile, argv[arg.wav]);
	    }
	    break;

	  case 'g':
	  case 'G':
	    arg.Fg = i;
	    break;

	  case 'f':
	  case 'F':
	    arg.F0 = i;
	    break;

	  case 'd':
	  case 'D':
	    arg.dur = i;
	    break;

	  case 'c':
	  case 'C':
	    arg.cq = i;
	    break;

	  case 'j':
	  case 'J':
	    arg.jitter = i;
	    break;

	  case 'k':
	  case 'K':
	    arg.K = i;
	    break;

	  case 'n':
	  case 'N':
	    par.DC = .25;
	    arg.noise = i;
	    break;

	  case 'r':
	  case 'R':
	    arg.fs = i;
	    break;

	  case 'a':
	  case 'A':
	    arg.amp = i;
	    break;

	  case 'l':
	  case 'L':
	    arg.DC = i;
	    break;

	  case 'z':
	  case 'Z':
	    arg.Kvar = i;
	    break;
	
	  case 's':
	  case 'S':
	    arg.Shimmer = i;
	    break;


	  default:
	    usage();
	    break;
	}
	break;
    }
  }
  if ((i != argc && *argv[i] != 'i') || arg.wav == -1) usage();

  /* INITIALIZE VARIABLES, BUFFERS, ETC */
  initialization(argv);

  /* OPEN SPEECH FILE */
  if((outfile = fopen(argv[arg.wav], "wb")) ==  NULL) {
    printf("Error while creating %s\n", argv[arg.wav]);
  }

  /* WRITE HEADER */
  if(fwrite(&header, sizeof(header), 1, outfile) != 1) {
    printf("Error while writing header to %f\n", argv[arg.wav]);
    exit(0);
  }

  /* PUT MESSAGES ON THE SCREEN */
  msg(argv);
  printf("Wait...");


  /* ---------------------- MAIN LOOP ------------------------------ */
  srandom(time(NULL));
  nSamples = (unsigned long) par.fs*par.dur; /* # of samples in the file */
  CountSamples = 0L;
  P = T = (int) ( (float) par.fs/par.F0);

  do {

    if(arg.jitter != -1 && par.jitter != 0.0) {

      /*
	par.jitter = specified MEAN jitter / 100   [=> initialization()]
	J = RANDOM jitter, uniformely distributed,
	    with average value = par.jitter;

				  p(J)
				   
		 ������������������
		                                   
		                                   
	     �����������������������J
	-2*par.jitter -par.jitter  0   par.jitter  2*par.jitter

		 argument for random():
		 �������������������		 0                                 4(*10000)*jitter/[10000)]


	The recursion for DeltaPer[0] (below) comes from:

				 (n)	            (n-1)
		   {(P + DeltaPer[0]) - (P + DeltaPer[1])}
	   J  = ����������������������		(1/2)*[(P + DeltaPer[0]) + (P + DeltaPer[1])]
				    (n)                 (n-1)
      */

      /* shift */
      DeltaPer[1] = DeltaPer[0];
      do {
	/* for the old borland compiler
	J = (float) random( (int) ceil(40000.0*par.jitter) )/10000.0 -
	    2.0*par.jitter;*/

         /* making gcc happy */
	J = (  random() / (RAND_MAX * 10000.0) ) * 40000.0 * par.jitter -
	    2.0*par.jitter;

	DeltaPer[0] = DeltaPer[1]*(2.0 + J)/(2.0 - J) +
		      2.0*P*J/(2.0 - J);

	T = (signed short) ceil((float) P + DeltaPer[0]);
      } while ( (float) T > (float) 1.2*P || (float) T < (float) 0.8*P);
    }

    float Amplitude;

   if(arg.Shimmer != -1 && par.Shimmer != 0.0) {
   DeltaShimmer[1] = DeltaShimmer[0];
      do {
	float epsilon =  ((float)random()) / RAND_MAX;
	
	S =  epsilon  * 4.0 * par.Shimmer - 2.0*par.Shimmer;

	DeltaShimmer[0] = DeltaShimmer[1]*(2.0 + S)/(2.0 - S) +
		      2.0*par.amp*S/(2.0 - S);

	Amplitude = ((float) par.amp + DeltaShimmer[0]);
      } while (  (Amplitude > (float) 1.8*par.amp) || (Amplitude < (float) 0.2*par.amp));
	printf("%5.2f \n", S);
	
    }
    else
    {
	Amplitude = (float) par.amp;
    }


    /* generate glottal flow */
    T2 = ceil(0.5*par.cq*P);
    for(i=0; i<T2; i++) {
      x[i] = ceil( Amplitude * 0.5*(1.0 - cos(PI*i/T2) ));
      if(x[i] < par.DC) {
	x[i] = par.DC;
	T4 = i;
      }
    }
    float Knew = par.K * ( 1 + 2 * par.Kvar *  ( ( (1.0 * random())/RAND_MAX )  - 0.5));
    //printf("valor de K: %f \n",Knew); 
    for(i=T2; i<2*T2; i++) {
      x[i] = ceil( (float) Amplitude*(Knew*cos(PI*(i-T2)/T2) - Knew + 1.0));
      if(x[i] < par.DC) break;
    }

    T3 = i;

    for(i=T3; i<T; i++) {
      x[i] = par.DC;
    }


    /* NOISE GENERATION

		       1
		      ��SUM[x(i)]
		       T   i=T4...T3 (see figure below)
      SNR_db = 10*log ��������= specified SNR (cmd line)
		       1
		      ��SUM[w(i)]
		       T   i=0...T4 and T3...T (figure see below)

			   ��> Denominator = <w(n)>

      w(i) = random variable, uniform(0, ... M) => E[w(i)] = M/12

      Assuming ergordicity,  <w(n)> = E[w(i)], whence

						      
      M = par.NoiseDistWidth = [12*<x(n)>/par.noise)]

					  par.noise = pow(10, arg.noise/10)
			.                    .
		      .   .                .   .
	  |	     .     .              .      .
	  | noise   .       .  noise     .        .  noise
      DC  +xxxxxxxx+         +xxxxxxxxxx+          +xxxxxxx
	  |       .
	0 +------+-+---+-----+----------+---------------------> time
		 0     T2   T3          T
		 +-+
		  T4
    */



    if(arg.noise != -1) {
      aux = 0.0;
      for(i=T4; i<T3; i++) {
	aux += (float) x[i]*x[i];
      }
      x_pow = aux/((float) T3 - T4);


      aux = 1.0 + ((float) T3 - T4)/((float) T);
      par.NoiseDistWidth = sqrt(12*aux*(x_pow)/par.noise);

      aux = 0.0;
      for(i=0; i<T4; i++) {

	w[i] = (signed short) ceil( ((  1.0 * random()) / RAND_MAX )* par.NoiseDistWidth - par.NoiseDistWidth/2.0);

	/*
	w[i] = (int) ceil(random(par.NoiseDistWidth));
	*/

	aux += (float) w[i]*w[i];
	x[i] = truncate( (float) x[i] + w[i]);
      }
      for(i=T3; i<T; i++) {

	w[i] = (signed short) ceil( ( ( 1.0*random() )/RAND_MAX )* par.NoiseDistWidth - par.NoiseDistWidth/2.0);

	/*
	w[i] = (int) ceil(random(par.NoiseDistWidth));
	*/

	aux += (float) w[i]*w[i];
	x[i] = truncate( (float) x[i] + w[i]);
      }
      w_pow = aux/( (float) T);

      printf("SNRdb = %5.2f\n", 10.0*log10(x_pow/w_pow));

    }

    CountSamples += T;
    if(CountSamples > nSamples) k = T - (CountSamples - nSamples);
    else k = T;


    if(fwrite(x, sizeof(signed short), k, outfile) != k) {
      printf("Error while writing samples to %f\n", argv[arg.wav]);
      exit(0);
    }

  } while(CountSamples < nSamples);


  free(x);

  fclose(outfile);

  printf("done\n");
  exit(0);
}

/*---------------------------- FUNCTIONS ---------------------------*/

void usage(void)
{
  printf("(c) Maurilio N. Vieira, 28 mar 1997\n");
  printf(" Simulated airflow based on Fant (1979),\n");
  printf(" Glottal Source and Excitation Analysis,\n");
  printf(" STL-QPSR 1/1979, pp. 85-107\n\n");
  printf("usage:\n");
  printf("%s -o file [-args {description (defaults <range>)}]\n\n",
							   "voicegen");
  printf("-o x {Output file (.wav, pcm, 16 bits/sample)}\n");
  printf("-r x {sampling Rate (22050 Hz <44100, 22050, or 11025>)}\n");
  printf("-d x {Duration (> 0.5 seconds)}\n");
  printf("-j x {jitter (0% <0-10%>)}\n");
  printf("-c x {closed quotient (.55 <0-1>)}\n");
  printf("-f x {Fundamental frequency, F0, 120 Hz }\n");
  printf("-g x {Glottal formant, Fg > F0, in Fant's (1979) model, 125 Hz }\n");
  printf("-k x {Speed of closure, K, in Fant's (1979) model 0.65 <0.55-1.00>}\n");
  printf("-z x {Variation of speed of closure (0.0 <0-1.0>))}\n");
  printf("-s x {shimmer (0.0 <0-10%>))}\n");
  printf("-n x {cycle-to-cycle SNR (0 dB  <0-50>) \n");
  printf("      aditive noise, uniforme distribution, closed phase}\n");
  printf("-a x {maximum amplitude (12000 <0-32767>)}\n");
  printf("-l x {DC flow, proportion of max amplitude (0.0 <0-0.30>))}\n");

  exit(0);
}

void initialization(char *argv[])
{
  float f;
  int i;
  long l;

  /* COMMAND LINE PARAMETERS */
  if(arg.dur != -1) {
    f = atof(argv[arg.dur]);
    if(f >= 0.5) par.dur = f;
    else usage();
  }

  if(arg.jitter != -1) {
    f = atof(argv[arg.jitter])/100.0;
    if(f >= 0.0 && f <= 10.0) par.jitter = f;
    else usage();
  }

  if(arg.K != -1) {
    f = atof(argv[arg.K]);
    if(f >= 0.50) par.K = f;
    else usage();
  }

  if(arg.cq != -1) {
    f = atof(argv[arg.cq]);
    if(f >= 0.0 && f<=1.0) par.cq = f;
    else usage();
  }

  if(arg.Fg != -1) {
    f = atof(argv[arg.Fg]);
    if(f >= 50) par.Fg = f;
    else usage();
  }



  if(arg.F0 != -1) {
    f = atof(argv[arg.F0]);
    if((f >= 50) && (f<par.Fg)) par.F0 = f;
    else usage();
  }

  if(arg.noise != -1) {
    f = atof(argv[arg.noise]);
    if(f >= 0.0 && f <= 50) {
      par.noise = pow(10, f/10);
    }
    else usage();
  }

  if(arg.amp != -1) {
    i = atoi(argv[arg.amp]);
    if(i >= 0 && i < 32767) par.amp = i;
    else usage();
  }

  if(arg.DC != -1) {
    f = atof(argv[arg.DC]);
    if(f >= 0 && f<= 0.3) par.DC = f*par.amp;
    else usage();
  }

 if(arg.Kvar != -1) {
    f = atof(argv[arg.Kvar]);
    if(f >= 0 && f <= 1) par.Kvar = f;
    else usage();
  }


  if(arg.fs != -1) {
    l = atol(argv[arg.fs]);
    if( (l == 44100L) || (l != 22050L) || (l == 11025L) )
       par.fs = l;
    else usage();
  }

  if(arg.Shimmer != -1) {
    f = atof(argv[arg.Shimmer]);
    if(f >= 0 && f <= 100) par.Shimmer = f/100;
    else usage();
  }


  /* FILE HEADER */
  memcpy(header.riff, "RIFF", 4);
  memcpy(header.wave, "WAVE", 4);
  memcpy(header.fmt, "fmt ", 4);
  memcpy(header.data, "data", 4);

  header.datasize = (long int) (par.dur*par.fs*2);
  header.filesize = (long int) header.datasize + 44L - 8L;
  header.fmtsize = 16L;
  header.wFormatTag = (int) 1;
  header.nChannels = (int) 1;
  header.nSamplesPerSec = par.fs;
  header.wBitsPerSample = (int) 16;
  header.nBlockAlign =
    (int) (header.wBitsPerSample/8 * header.nChannels);
  header.nAvgBytesPerSec =
     (long int) (header.nBlockAlign*header.nSamplesPerSec);


  /* alloc memory buffers */
  if((x = malloc(sizeof(signed short) * (par.fs/par.Fg) * 2)) == NULL) {
    printf("out of memory in call to malloc(x).\n");
    exit(1);
  }

}

void msg(char *argv[])
{
  printf("(c) Maurilio N. Vieira, 1996\nSynthetic vowel generator\n");
  printf("ported to gcc - Joao SANSAO, Feb. 2007");
  printf("Output file = %s\n", argv[arg.wav]);

  if(arg.noise != -1) {
    printf("SNR: %5.2f dB, ", 10.0*log10(par.noise) );
  }

  printf("Fs=%ld Hz, Dur=%5.2f s, Fg=%d Hz, Amp = %d, DCflow=%5.2f\n",
	  par.fs, par.dur, (int) par.Fg, par.amp, par.DC);
}


signed short truncate(float aux)
{
  signed short i;
  if(aux > 32767) i = 32767;
  else if(aux < -32767) i = -32767;
  else i = (signed short) ceil(aux);

  return i;

}


/* Boxes
 ALT + ...  sigma = 228
	    2     = 253
	    infinite = 236

	 218   196   191
  (180)  179         179  (195)
	 192   196   217

*/

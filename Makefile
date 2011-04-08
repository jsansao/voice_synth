all : vowel flowgen_shimmer

vowel : vowel_new.o
	cc -lm -o  vowel vowel_new.o

flowgen_shimmer : flowgen_shimmer.o
	cc -lm -o flowgen_shimmer flowgen_shimmer.o 

vowel_new.o : vowel_new.c
	cc -c -lm vowel_new.c

flowgen_shimmer.o : flowgen_shimmer.c
	cc -c -lm flowgen_shimmer.c

clean: 
	rm *.o flowgen_shimmer vowel



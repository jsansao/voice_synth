all : vowel flowgen_shimmer acoustic egg 

vowel : vowel_new.o
	cc -lm -o  vowel vowel_new.o

acoustic : acoustic.o
	cc -lm -o  acoustic acoustic.o

flowgen_shimmer : flowgen_shimmer.o
	cc -lm -o flowgen_shimmer flowgen_shimmer.o 

egg : egg.o
	cc -lm -o egg egg.o

vowel_new.o : vowel_new.c
	cc -c -lm vowel_new.c

acoustic.o : acoustic.c
	cc -c -lm acoustic.c

flowgen_shimmer.o : flowgen_shimmer.c
	cc -c -lm flowgen_shimmer.c

egg.o : egg.c
	cc -c -lm -Wall egg.c


clean: 
	rm *.o flowgen_shimmer vowel acoustic egg



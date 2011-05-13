APPS= cycleprunner

all: $(APPS)

cycleprunner: cycleprunner.c
	mpicc -o cycleprunner -O3 -fopenmp -Wall -lm cycleprunner.c cmaes.c forwardmodel.c dkhdebug.o

cmaes.o: cmaes.h cmaes.c
	gcc -o cmaes.o -c cmaes.c

dkhdebug.o: dkhdebug.c dkhdebug.h
	gcc -o dkhdebug.o -c dkhdebug.c
	
clean:
	rm -f *.dat *.o $(APPS)

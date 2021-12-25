#include <stdio.h>
#include <mpi.h>
#include <unistd.h> //Header file optarg.
#include <stdlib.h> //For atoi 


/* This is the struct inputs will be loaded into. Rather than using global variables,
    the struct will be initialised in main and its reference will be passed to the functions
    that need this information. */
struct inputStruct {
    int dimension;
    int numOfProcessors;
    double precision;
};


void handleInput(int argc, char *argv[], struct inputStruct *input);

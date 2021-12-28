#include <stdio.h>
#include <mpi.h>
#include <unistd.h> //Header file optarg.
#include <stdlib.h> //For atoi.
#include <stdbool.h> //For bool.
#include <math.h> //For fabs().

/* This is the struct inputs will be loaded into. Rather than using global variables,
    the struct will be initialised in main and its reference will be passed to the functions
    that need this information. */
struct inputStruct {
    int dimension;
    int totalProcessorCOUNT;
    double precision;
};

struct rankRowInfo {
    int startRow;
    int endRow;
    int rank;
};


void handleInput(int argc, char *argv[], struct inputStruct *input);
struct rankRowInfo computeResponsibleRows(int rank, int dimension, int numOfProcs);
double** allocateMemory(bool *errorCheck, double **bufferAddress, int dimension);
void initMatrix(double **matrix, int dimension);
void printMatrix(double **matrix, int dimension);
double** allocateSubMatrixMemory(bool *errorCheck, double **bufferAddress, int width, int height);
void printSubMatrix(double **matrix, int width, int height);
double computeNeighborAverage(double **matrix, int row, int col);
void copyMatrix(double **originalMatrix, double **copyMatrix, int width, int height);

#include "main.h" // Include the header files I defined, because the top of this file was too cluttered.

int main(int argc, char** argv) {
    int rc, myrank, nproc, namelen;
    bool firstIteration = true;
    char name[MPI_MAX_PROCESSOR_NAME];

    rc = MPI_Init(&argc,  &argv);
    if (rc != MPI_SUCCESS) {
        printf ("Error starting MPI program\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Get_processor_name(name, &namelen);
    
    struct inputStruct inputs;
    inputs.numOfProcessors = nproc;
    handleInput(argc, argv, &inputs);



    /* rank 0 will be the main thread, running all the workers, and checking the precision. */
    if(myrank == 0) {
        
        // This is the first iteration. Node 0 should generate the random square matrix, cut it up accordingly, and distribute it across the procs. 
        if(firstIteration == true) {
            firstIteration = false;
            
            bool memoryError = false;
            double *originalMatrixBuffer;   // Buffer that will be freed at the end. 

            double **originalMatrix = allocateMemory(&memoryError, &originalMatrixBuffer, inputs.dimension);


            if(memoryError == true) {
                printf("ERROR OCCURRED DURING MEMORY ALLOCATION TO MAIN MATRIX\n");
            }

            initMatrix(originalMatrix, inputs.dimension); // Initialise the matrix. Seeded random values, and non symmetric boundary conditions.

            // Now need to allocate the right rows to the right machine. 
            for(int remoteRank = 1; remoteRank < inputs.numOfProcessors; remoteRank++) {
                struct rankRowInfo currentInfo = computeResponsibleRows(remoteRank, inputs.dimension, nproc);
                int sendCOUNT = (currentInfo.endRow-currentInfo.startRow+1+2) * inputs.dimension;

                //printf("My rank is: %i, start row: %i, end row: %i, send count: %i\n", remoteRank, currentInfo.startRow, currentInfo.endRow, sendCOUNT);
                MPI_Send(originalMatrix[currentInfo.startRow-1], sendCOUNT, MPI_DOUBLE, currentInfo.rank, 0, MPI_COMM_WORLD);
            }
        }

        // Check to see if all other ranks (which are workers) report that they have converged.
        // If so, collect the pieces of the matrix, and stitch together.
        else if(firstIteration == false) {

        }
    } 
    else {
        struct rankRowInfo myRankRowInfo = computeResponsibleRows(myrank, inputs.dimension, nproc);
        
        if(firstIteration == true) {
            firstIteration = false;
            bool memoryError = false;

            // Get the difference, and then add two, because the boundary condition rows also need to be inserted, and add one because its an inclusive range.
            int height = myRankRowInfo.endRow - myRankRowInfo.startRow + 2 + 1; 
            int receiveCOUNT = height * inputs.dimension;

            double *subMatrixBuffer;
            double **subMatrix = allocateSubMatrixMemory(&memoryError, &subMatrixBuffer, inputs.dimension, height);

            printf("This is rank: %i, receive count: %i, \n", myRankRowInfo.rank, receiveCOUNT);
            MPI_Recv(subMatrix[0], receiveCOUNT, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            for(int i = 1; i < inputs.numOfProcessors; i++) {
                if(i == myrank) break;
                sleep(0.5*i);
            }
            printSubMatrix(subMatrix, inputs.dimension, height);
        }

        else if(firstIteration == false) {

        }

        //printf("My rank is %i, Start Row: %i, End Row: %i \n", myRankRowInfo.rank, myRankRowInfo.startRow, myRankRowInfo.endRow);
    }

    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}



void handleInput(int argc, char *argv[], struct inputStruct *input) {
    int c;
    // The while loop below parses inputs, using getopt. 
    while((c = getopt(argc, argv, "d:p:" )) != -1) 
    {
        switch(c)
        {
            case 'd':
                input->dimension = atoi(optarg);
                break;

            case 'p':
                input->precision = atof(optarg);
                break;

            case '?': 
                printf("Unknown option: %c\n", optopt);
                break; 
        }
    }

    for(; optind < argc; optind++){     
        printf("extra arguments: %s\n", argv[optind]); 
    }

    // Check for input errors.
    
    if(input->dimension <= 3) {
        printf("INPUT ERROR: The dimension you have entered is too little. Setting dimension to 4 (a 4 by 4 square will be used).\n");
        input->dimension = 4;
    }

    if(input->precision <= 0.0005) {
        printf("INPUT ERROR: Desired precision too low. Setting it to 0.0005.\n");
        input->precision = 0.0005;
    }

    if(input->numOfProcessors-1 > input->dimension - 2) {
        printf("INPUT ERROR: There are too many cores for this job, will only be using %i cores actively.\n", input->dimension - 2);
    }

    //printf("Size is: %i\n", input->dimension);
    //printf("Desired precision is: %f\n", input->precision);
}

struct rankRowInfo computeResponsibleRows(int rank, int dimension, int numOfProcs) {
    int workerProcCOUNT = numOfProcs-1;      // Count of worker processes. DO I NEED TO REDUCE THIS IF NPROC IS BIGGER THAN DIMENSION?
    int rowsToRelaxCOUNT = dimension - 2;    // Rows to relax, meaning it excludes the boundary rows. 
    int extraRowCOUNT = rowsToRelaxCOUNT % workerProcCOUNT;     // This many procs need to take care of an additional row.
    int startRow, endRow, height = rowsToRelaxCOUNT / workerProcCOUNT;

    if(rank <= extraRowCOUNT + 1) {
        startRow = ((rank-1)*(height+1)) + 1;
        if(rank <= extraRowCOUNT) {
            height++;
        }
        endRow = startRow + height - 1;
    }
    else {
        startRow = ((extraRowCOUNT) * (height+1)) + ((rank - extraRowCOUNT - 1) * height) + 1;
        endRow = startRow + height - 1;
    }
    struct rankRowInfo result;
    result.startRow = startRow;
    result.endRow = endRow;
    result.rank = rank;

    return result;
}

double** allocateMemory(bool *errorCheck, double **bufferAddress, int dimension){

    // Allocate memory for 2-D array of size dimension * dimension.
    double **a = calloc((long unsigned int) dimension, (long unsigned int) sizeof(double*));
    if(a == NULL) {
        printf("ERROR: Could not malloc original list of pointers to memory!\n");
        *errorCheck = true;
    }

    double *aBuf = calloc((long unsigned int)dimension, (long unsigned int)dimension*sizeof(double));
    *bufferAddress = aBuf;
    if(aBuf == NULL) {
        printf("ERROR: Could not malloc buffer!\n");
        *errorCheck = true;
    }

    for(int i = 0; i < dimension; i++) {
        a[i] = aBuf + dimension*i;
    }

    return a;
}

double** allocateSubMatrixMemory(bool *errorCheck, double **bufferAddress, int width, int height) {
    double **a = calloc((long unsigned int) height, (long unsigned int) sizeof(double*));
    if(a == NULL) {
        printf("ERROR: Could not malloc original list of pointers to memory!\n");
        *errorCheck = true;
    }

    double *aBuf = calloc((long unsigned int)height, (long unsigned int)width*sizeof(double));
    *bufferAddress = aBuf;
    if(aBuf == NULL) {
        printf("ERROR: Could not malloc buffer!\n");
        *errorCheck = true;
    }

    for(int i = 0; i < height; i++) {
        a[i] = aBuf + width*i;
    }

    return a;

}

void initMatrix(double **matrix, int dimension) {

    srand((unsigned)66);    // Seed the random number generator with any random number.
    // What matters is consistency, and any seed achieves this.

    for(int row = 0; row < dimension; row++) {
        for(int col = 0; col < dimension; col++) {
            matrix[row][col] = rand() / ((double) RAND_MAX);
        }
    }

    for(int i = 0; i < dimension; i++) {
        matrix[0][i] = 1;
        matrix[i][0] = 1;

        // Make the right and bottom border conditions 0.5, so that the matrix isn't symmetric. 
        if(i < dimension - 1) {
            matrix[dimension - 1][i+1] = 0.5;
            matrix[i+1][dimension - 1] = 0.5;
        }

    }
}

void printMatrix(double **matrix, int dimension) {
    for(int row = 0; row < dimension; row++) {
        for(int col = 0; col < dimension; col++) {
            printf("%f ", matrix[row][col]);
        }
        printf("\n");
    }
    printf("\n");
    printf("\n");
    printf("\n");

}

void printSubMatrix(double **matrix, int width, int height) {
    for(int row = 0; row < height; row++) {
        for(int col = 0; col < width; col++) {
            printf("%f ", matrix[row][col]);
        }
        printf("\n");
    }
    printf("\n");
    printf("\n");

}
#include "main.h" // Include the header files I defined, because the top of this file was too cluttered.

int main(int argc, char** argv) {
    int rc, myrank, nproc, namelen;
    char name[MPI_MAX_PROCESSOR_NAME];

    bool allSubMatrixSolved = false;
    bool localSubMatrixSolved = false;

    rc = MPI_Init(&argc,  &argv);
    if (rc != MPI_SUCCESS) {
        printf ("Error starting MPI program\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Get_processor_name(name, &namelen);
    
    struct inputStruct inputs;
    inputs.totalProcessorCOUNT = nproc;
    handleInput(argc, argv, &inputs);

    /* rank 0 will be the main thread, running all the workers, and checking the precision. */
    if(myrank == 0) {
        // Start recording time for benchmarking.
        double beginTime = MPI_Wtime();

        // This is the first iteration. Node 0 should generate the random square matrix, cut it up accordingly, and distribute it across the procs. 
        bool memoryError = false;
        double *originalMatrixBuffer;   // Buffer that will be freed at the end. 
        double **originalMatrix = allocateMemory(&memoryError, &originalMatrixBuffer, inputs.dimension);

        if(memoryError == true) {
            printf("ERROR OCCURRED DURING MEMORY ALLOCATION TO MAIN MATRIX\n");
        }
        
        // Initialise the matrix. Seeded random values, and non symmetric boundary conditions.
        initMatrix(originalMatrix, inputs.dimension);

        // Now need to allocate the right rows to the right machine. 
        for(int remoteRank = 1; remoteRank < inputs.totalProcessorCOUNT; remoteRank++) {
            struct rankRowInfo currentInfo = computeResponsibleRows(remoteRank, inputs.dimension, inputs.totalProcessorCOUNT);
            int sendCOUNT = (currentInfo.endRow-currentInfo.startRow+1+2) * inputs.dimension;

            MPI_Send(originalMatrix[currentInfo.startRow-1], sendCOUNT, MPI_DOUBLE, currentInfo.rank, 0, MPI_COMM_WORLD);
        }

        // Check to see if all other ranks (which are workers) report that they have converged.
        // If so, collect the pieces of the matrix, and stitch together.
        while(!allSubMatrixSolved) {
            localSubMatrixSolved = true;
            MPI_Allreduce(&localSubMatrixSolved, &allSubMatrixSolved, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
        }

        // The computations have been completed, and all sub matrices are solved.
        // Stitch them back together. 
        for(int remoteRank = 1; remoteRank < inputs.totalProcessorCOUNT; remoteRank++) {
            struct rankRowInfo currentInfo = computeResponsibleRows(remoteRank, inputs.dimension, inputs.totalProcessorCOUNT);
            int recvCOUNT = (currentInfo.endRow-currentInfo.startRow+1+2) * inputs.dimension;

            MPI_Recv(originalMatrix[currentInfo.startRow], recvCOUNT, MPI_DOUBLE, currentInfo.rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        double endTime = MPI_Wtime();
        //printMatrix(originalMatrix, inputs.dimension);
        printf("Time taken is %f seconds!\n", endTime-beginTime);
    }
    else if (myrank >= 1 && myrank < inputs.totalProcessorCOUNT) 
    {
        struct rankRowInfo myRankRowInfo = computeResponsibleRows(myrank, inputs.dimension, inputs.totalProcessorCOUNT);

        // Initialise the subMatrix memory. 
        bool memoryError = false;
        int absoluteHeight = myRankRowInfo.endRow - myRankRowInfo.startRow + 2 + 1; // NOT an index, literal height. 
        double *subMatrixBuffer;
        double **subMatrix = allocateSubMatrixMemory(&memoryError, &subMatrixBuffer, inputs.dimension, absoluteHeight);
        
        // Initialise the read matrix as well.
        double *writeMatrixBuffer;
        double **writeMatrix = allocateSubMatrixMemory(&memoryError, &writeMatrixBuffer, inputs.dimension, absoluteHeight);

        // Get the difference, and then add two, because the boundary condition rows also need to be inserted, and add one because its an inclusive range.
        int receiveCOUNT = absoluteHeight * inputs.dimension;

        // Receive subMatrix from the main proc. 
        MPI_Recv(subMatrix[0], receiveCOUNT, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        copyMatrix(subMatrix, writeMatrix, inputs.dimension, absoluteHeight);

        while(!allSubMatrixSolved) {
            localSubMatrixSolved = true;
            for(int row = 1; row < absoluteHeight - 1; row++) {
                for(int col = 1; col < inputs.dimension - 1; col++) {
                    writeMatrix[row][col] = computeNeighborAverage(subMatrix, row, col);
                    if(fabs(writeMatrix[row][col] - subMatrix[row][col]) > inputs.precision) localSubMatrixSolved = false;
                }
            }
            copyMatrix(writeMatrix, subMatrix, inputs.dimension, absoluteHeight);
            
            // Send out this ranks updated boundary conditions
            // Dont send the first row anywhere, because there is nowhere logical for it to be sent. 
            if(myrank > 1) 
            {
                MPI_Send(subMatrix[1], inputs.dimension, MPI_DOUBLE, myrank-1, 0, MPI_COMM_WORLD);
            }
            // Similarly, send the final rows on each rank, except the one containing the final row for the entire matrix.
            if(myrank < inputs.totalProcessorCOUNT-1) 
            {
                MPI_Send(subMatrix[absoluteHeight-2], inputs.dimension, MPI_DOUBLE, myrank+1, 0, MPI_COMM_WORLD);
            }

            // Receive the updated boundary conditions.
            // Rank 1's boundary conditions are set in stone - it is the uppermost row of the entire matrix.
            if(myrank > 1)
            {
                //printf("rank %i in recv %i\n", myrank, cnt);
                MPI_Recv(subMatrix[0], inputs.dimension, MPI_DOUBLE, myrank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            // Similarly, the final rank's lower row does not need to be updated. 
            if(myrank < inputs.totalProcessorCOUNT - 1)
            {
                //printf("rank %i in recv %i\n", myrank, cnt);
                MPI_Recv(subMatrix[absoluteHeight-1], inputs.dimension, MPI_DOUBLE, myrank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            // Reduce the status of all seperate nodes and broadcast it back, to see if we still need to compute more. 
            MPI_Allreduce(&localSubMatrixSolved, &allSubMatrixSolved, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
        }

        int sendCOUNT = (myRankRowInfo.endRow-myRankRowInfo.startRow+1) * inputs.dimension;
        MPI_Send(subMatrix[1], sendCOUNT, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    // This final else statement is added because the communicator COMM_WORLD is being used.
    // Since all ranks in the communicator must reply to the reduce, in the case that there are 
    // too many procs created initially (too many cores for the matrix size), these procs simply 
    // have to contribute to the Allreduce which expects input from them. 
    else {
        while(!allSubMatrixSolved) {
            localSubMatrixSolved = true;
            MPI_Allreduce(&localSubMatrixSolved, &allSubMatrixSolved, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
        }
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

    if(input->precision < 0.0001) {
        printf("INPUT ERROR: Desired precision too low. Setting it to 0.0001.\n");
        input->precision = 0.0001;
    }

    if(input->precision > 0.2) {
        printf("INPUT ERROR: Desired precision too high, no meaningful computations will be performed. Setting it to 0.01.\n");
        input->precision = 0.01;
    }
    
    // Check if there are too many processor to compute this matrix, using a row slice method. 
    if((input->totalProcessorCOUNT - 1) > (input->dimension-2)) {
        printf("INPUT ERROR: There are too many cores for this job, will only be using %i procs to compute the matrix.\n", input->dimension - 2);
        input->totalProcessorCOUNT = input->dimension - 1;
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

/* Deep copies originalMatrix into copyMatrix */
void copyMatrix(double **originalMatrix, double **copyMatrix, int width, int height) {
    for(int row = 0; row < height; row++) {
        for(int col = 0; col < width; col++) {
            copyMatrix[row][col] = originalMatrix[row][col];
        }
    }
    return;
}

double computeNeighborAverage(double **matrix, int row, int col) {
    return (matrix[row-1][col] + matrix[row+1][col] + matrix[row][col+1] + matrix[row][col-1]) / 4;    
}
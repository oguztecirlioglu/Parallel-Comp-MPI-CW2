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
            
        }

        // Check to see if all other ranks (which are workers) report that they have converged.
        // If so, collect the pieces of the matrix, and stitch together.
        else if(firstIteration == false) {

        }
    } 
    else {
        struct rankRowInfo myRankRowInfo = computeResponsibleRows(myrank, inputs.dimension, nproc);
        
        if(firstIteration == true) {

        }

        else if(firstIteration == false) {
            
        }

        printf("My rank is %i, Start Row: %i, End Row: %i \n", myRankRowInfo.rank, myRankRowInfo.startRow, myRankRowInfo.endRow);
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

    if(input->numOfProcessors > input->dimension - 2) {
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

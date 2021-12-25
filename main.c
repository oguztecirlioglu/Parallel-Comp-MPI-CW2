#include "main.h" // Include the header files I defined, because the top of this file was too cluttered.

int main(int argc, char** argv) {
    int rc, myrank, nproc, namelen;
    char name[MPI_MAX_PROCESSOR_NAME];

    rc = MPI_Init(&argc,  &argv);
    if (rc != MPI_SUCCESS) {
        printf ("Error starting MPI program\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Get_processor_name(name, &namelen);
    
    if(myrank == 0) {
        printf("I am the head node! \n");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    printf("I am processor %s, I am rank %i out of %i processors.\n", name, myrank, nproc);

    struct inputStruct inputs;
    inputs.numOfProcessors = nproc;
    handleInput(argc, argv, &inputs);


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

    printf("Size is: %i\n", input->dimension);
    printf("Desired precision is: %f\n", input->precision);
}

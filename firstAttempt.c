#include <stdio.h>
#include <mpi.h>
#include <unistd.h> //Header file optarg.

void handleInput(int argc, char *argv[], struct inputStruct input);


/* This is the struct inputs will be loaded into. Rather than using global variables,
    the struct will be initialised in main and its reference will be passed to the functions
    that need this information. */
struct inputStruct {
    int dimension;
    double precision;
};  

int main(int argc, char** argv) {
    
    MPI_Init(&argc,  &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        
    struct inputStruct userInputs;


    
    
    char proc_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(proc_name, &name_len);
    printf("Hello world from processor %s, rank %d out of %d processors\n", proc_name, world_rank, world_size);

    // Finalize the MPI environment.
    MPI_Finalize();

}



void handleInput(int argc, char *argv[], struct inputStruct input) {
    int c;
    // The while loop below parses inputs, using getopt. 
    while((c = getopt(argc, argv, "n:t:p:" )) != -1) 
    {
        switch(c)
        {
            case 'd':
                input.dimension = atoi(optarg);
                break;

            case 'p':
                input.precision = atof(optarg);
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
    
    if(input.dimension <= 3) {
        printf("INPUT ERROR: The dimension you have entered is too little. Setting dimension to 4 (a 4 by 4 square will be used).\n");
        input.dimension = 4;
    }

    if(input.threads < 1) {
        printf("INPUT ERROR: You have entered too few threads. Setting thread count to 1 \n");
        input.threads = 1;
    }

    if(input.threads > 44) {
        printf("INPUT ERROR: The machine this program was written for only has 44 cores. Setting thread count to 44. \n");
        input.threads = 44;
    }

    if(input.threads > input.dimension - 2) {
        printf("INPUT ERROR: Too many threads for this dimension. Setting thread count to dimension or 44 (minimum of the two).\n");
        input.threads = input.dimension > 44 ? 44 : input.dimension - 2; 
    }

    if(input.precision > 0.2) {
        printf("INPUT ERROR: Precision required is too high. Reducing the precision to produce more refined results. Setting precision to 0.001\n");
        input.precision = 0.2;
    }

    if(input.precision <= 0.0005 ) {
        printf("INPUT ERROR: Precision is either too low or 0. Setting precision to 0.0005\n");
        input.precision = 0.0005;
    }

    printf("Thread count is %i\n", input.threads);
    printf("Size is: %i\n", input.dimension);
    printf("Desired precision is: %f\n", input.precision);
}

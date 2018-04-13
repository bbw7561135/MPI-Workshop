#include <stdio.h>
#include <mpi.h>

#define MPI_DEFAULT_TAG 0
#define MPI_DEFAULT_COMM MPI_COMM_WORLD
#define DEFAULT_REORDER 0

#define NDIMS 1

int
main(int argc, char *argv[])
{

/* -------------------------------------------------------------------------- */

    int n_ranks, current_rank;
    int data_to_send, received_data;
    int rank_sum = 0;
    int right_rank, left_rank;

    int dims[NDIMS]; 
    int periods[NDIMS]; 

    MPI_Status recv_status;
    MPI_Request rank_request_right, rank_request_left;
    MPI_Comm cart_comm;

/* -------------------------------------------------------------------------- */

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &current_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);

    if (current_rank == 0)
        printf("Number of processes: %d.\n", n_ranks);
    
    // calculate the dimensions to use for the cart grid
    // set all of dims to be 0
    for (int i = 0; i < NDIMS; i++)
    {
        dims[i] = 0;
    }
    
    periods[0] = 1;  // periodic boundary

    MPI_Dims_create(n_ranks, NDIMS, dims);
    MPI_Cart_create(MPI_DEFAULT_COMM, NDIMS, dims, periods, DEFAULT_REORDER,
                    &cart_comm);
    
    // initialise the start of the sum
    data_to_send = current_rank; 
    for (int rank = 0; rank < n_ranks; rank++)
    {
        MPI_Cart_shift(cart_comm, 0, 1, &current_rank, &right_rank);
        MPI_Cart_shift(cart_comm, 0, -1, &current_rank, &left_rank);
        
        // send (to right) and receive (to left) -- non-blocking 
        MPI_Issend(&data_to_send, 1, MPI_INT, right_rank, MPI_DEFAULT_TAG,
                   cart_comm, &rank_request_right); 
        MPI_Irecv(&received_data, 1, MPI_INT, left_rank, MPI_DEFAULT_TAG,
                  cart_comm, &rank_request_left);

        // stop the non-blocking send and receives
        MPI_Wait(&rank_request_right, &recv_status);
        MPI_Wait(&rank_request_left, &recv_status);

        // sum everything up
        rank_sum += received_data;
        data_to_send = received_data;
    }

    printf("Rank %d: sum = %d.\n", current_rank, rank_sum);

    MPI_Finalize();
    return 0;
}
        

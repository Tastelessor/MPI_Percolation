#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>

#include "percolate.h"
#include "arralloc.h"

/*
 * Simple parallel program to test for percolation of a cluster.
 */

int main(int argc, char *argv[]) {
    /*
     * Define and initialize variables
     * Noted that this program can produce correct answer with any number of process
     * The first 2 arguments are accepted in this sequence:
     * mpirun -n <process number> ./percolate <seed> <rho>
     * If there are only 4 arguments, the command can be interpreted as
     * mpirun -n <process number> ./percolate <seed> <rho> <L>
     * If there are 5 arguments, the command will be
     * mpirun -n <process number> ./percolate <seed> <rho> <width> <height>
     * where width is the length of the first dimension, height is for the second dimension
     */
    // Basic variables
    int width = 480;
    int height = 480;
    int seed = 8759;
    double rho = 0.4064;
    int size, rank;

    // Check arguments
    if (argc >= 2) {
        seed = atoi(argv[1]);
    }
    if (argc >= 3) {
        rho = atof(argv[2]);
    }
    if (argc >= 4) {
        height = atoi(argv[3]);
        width = height;
    }
    if (argc >= 5) {
        height = atoi(argv[4]);
    }
    if (argc > 5) {
        printf("Too many arguments!\n");
        return -1;
    }

    // Init MPI and check arguments
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm new_comm;
    MPI_Request request[4];
    MPI_Status status[4];

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if(size > width * height){
        if(rank == 0){
            printf("Too many processes! \n"
                   "This is a %d * %d map but assigned with %d processes!\n", width, height, size);
        }
        MPI_Finalize();
        return -1;
    }

    // Local variables for iterations 
    int i, j, nhole, step, maxstep, oldval, newval;
    int nchangelocal, nchange, printfreq;
    int boundary_start, boundary_stop, perc;
    bool is_percolate = false;
    int displacement = width / size;
    double random_value; // Temporary random variable for map initialization
    int tag = 0;
    double average_value = 0; // Average value of the map array for correctness testing
    double average = 0; // The reduced target

    // Variables for 2D decomposition
    int dims[ndims] = {0, 0};
    int period[ndims] = {1, 0};
    int coords[ndims];
    int reorder = 1;

    // Variables for identifying process position
    int P = width / 6;
    int top, down, left, right;
    int rstart, rend, cstart, cend;
    int row_elements, col_elements, number_of_elements;

    // Decomposition of map
    MPI_Dims_create(size, ndims, dims);
    MPI_Cart_create(comm, ndims, dims, period, reorder, &new_comm);
    MPI_Cart_coords(new_comm, rank, ndims, coords);

    // Temporary variables for identifying the position of each block
    int M = width / dims[0];
    int N = height / dims[1];

    // Assign neighbors
    MPI_Cart_shift(new_comm, 0, 1, &top, &down);
    MPI_Cart_shift(new_comm, 1, 1, &left, &right);

    /*
     * Declaration for start index and end index for each dimension is not redundant here.
     * Once the map size is not divisible by the number of processes, then according to the
     * Pigeon-Hole Principle, there must be at least 1 process that has different workload than
     * others in a dimension. These indices are used for identifying Datatype length, identification
     * of position of communication related variables.
     */

    // Identify the exact start and end index for 2 dimensions
    rstart = coords[0] * M;
    cstart = coords[1] * N;
    if (coords[0] == dims[0] - 1)
        rend = width;
    else
        rend = rstart + M;
    if (coords[1] == dims[1] - 1)
        cend = height;
    else
        cend = cstart + N;

    col_elements = cend - cstart;
    row_elements = rend - rstart;
    number_of_elements = row_elements * col_elements;

    // Define map
    int **old;
    int **new;
    int **map;
    old = (int **) arralloc(sizeof(int), 2, row_elements + 2, col_elements + 2);
    new = (int **) arralloc(sizeof(int), 2, row_elements + 2, col_elements + 2);
    map = (int **) arralloc(sizeof(int), 2, width, height);

    // Define MPI datatype for
    MPI_Datatype row_boundary;
    MPI_Datatype col_boundary;

    MPI_Type_vector(row_elements, 1, col_elements + 2, MPI_INT, &row_boundary);
    MPI_Type_contiguous(col_elements, MPI_INT, &col_boundary);

    MPI_Type_commit(&row_boundary);
    MPI_Type_commit(&col_boundary);

    // Variables for performance test
    double total_start, total_end, par_start, par_end;

    /*
     *  Update for a fixed number of steps, periodically report progress
     */

    maxstep = width > height ? 5 * width : 5 * height;
    printfreq = 100;

    if (rank == 0) {
        total_start = MPI_Wtime();
        printf("percolate: running on %d process(es)\n", size);

        // Set the randum number seed and initialise the generator
        printf("percolate: width = %d, height = %d, rho = %f, seed = %d, maxstep = %d\n", width, height, rho, seed,
               maxstep);
    }
    rinit(seed);
    /*
     * Initialize map
     * This the parallelized version of map initialization. This is not applied because it generates
     * a different map with the same seed. To satisfy the requirements of the coursework, I have to
     * make it an additional version.
     */

    nhole = rank * width * height / size; // start value of nhole
    int initial_nholt = nhole; // record the start value to calculate how many holes are initialized
    int sum_hole; // The variable for reduction
    int row_start = rank * width / size; // start row
    int row_end = row_start + width / size; // end row

    // In case of the width is not divisible by the number of processes
    if(rank == size - 1)
        row_end = width;

    for (i = row_start; i < row_end; i++) {
        for (j = 0; j < height; j++) {
            random_value = uni();
            if (random_value < rho) {
                map[i][j] = 0;
            } else {
                nhole++;
                map[i][j] = nhole;
            }
        }
    }
    // Calculate how many squares are assigned with non-zero values
    nhole -= initial_nholt;

//    MPI_Reduce(&average_value, &average, 1, MPI_DOUBLE, MPI_SUM, 0, new_comm);
    MPI_Reduce(&nhole, &sum_hole, 1, MPI_INT, MPI_SUM, 0, comm);
    if(rank == 0)
        printf("percolate: rho = %f, actual density = %f\n", rho, 1.0 - ((double) sum_hole) / ((double) width * height));

    // Reduce map to all processes
    MPI_Allreduce(&map[rstart][0], &map[0][0], (row_end-row_start)*height, MPI_INT, MPI_MAX, comm);

    /*
     * Initialise the old map: assign meaningful values to its centre
     * leave boundaries to be halos
     */

    // Zero the top and bottom halos
    for (i = 0; i < row_elements + 2; i++) {
        old[i][0] = 0;
        old[i][col_elements + 1] = 0;
    }

    // Zero the left and right halos
    for (j = 0; j < col_elements + 2; j++) {
        old[0][j] = 0;
        old[row_elements + 1][j] = 0;
    }

    // Initialize the old array
    for (i = 1; i < row_elements + 1; i++) {
        for (j = 1; j < col_elements + 1; j++) {
            old[i][j] = map[rstart + i - 1][cstart + j - 1];
        }
    }

    // Prepare to loop for update maps
    step = 1;
    nchange = 1;

    // Record parallel updating time cost
    if (rank == 0) {
        par_start = MPI_Wtime();
    }
    while (step <= maxstep) {
        // Exchange boundaries
        MPI_Issend(&old[1][1], 1, row_boundary, left, 1, new_comm, &request[0]);
        MPI_Issend(&old[1][col_elements], 1, row_boundary, right, 2, new_comm, &request[1]);
        MPI_Issend(&old[1][1], 1, col_boundary, top, 3, new_comm, &request[2]);
        MPI_Issend(&old[row_elements][1], 1, col_boundary, down, 4, new_comm, &request[3]);

        MPI_Recv(&old[1][col_elements + 1], 1, row_boundary, right, 1, new_comm, &status[0]);
        MPI_Recv(&old[1][0], 1, row_boundary, left, 2, new_comm, &status[1]);
        MPI_Recv(&old[row_elements + 1][1], 1, col_boundary, down, 3, new_comm, &status[2]);
        MPI_Recv(&old[0][1], 1, col_boundary, top, 4, new_comm, &status[3]);

        MPI_Waitall(4, request, status);

        /*
         * This is the additional implementation for new boundary condition
         * The program judge if a rank is on extreme top or bottom, if it is, then start filtering
         * However, if the number of process is not perfect, or in other words, cannot divide map
         * perfectly, and the number of process is not enough. Then, a process may on top and bottom
         * at the same time. Additional judgement is added to deal with such conditions.
         * */

        // New boundary condition
        if (coords[0] == 0 || coords[0] == dims[0] - 1) {
            // Judge if it's top layer
            int crd = coords[0] == 0 ? 0 : row_elements + 1;
            for (j = 0; j < col_elements + 2 && coords[1] * col_elements + j < P + 1; j++) {
                old[crd][j] = 0;
            }
            for (j = col_elements + 1; j >= 0 && coords[1] * col_elements + j > 5 * P; j--) {
                old[crd][j] = 0;
            }
            // Additional judgement in case of a block is both top and bottom
            if (coords[0] == dims[0] - 1) {
                crd = row_elements + 1;
                for (j = 0; j < col_elements + 2 && coords[1] * col_elements + j < P + 1; j++) {
                    old[crd][j] = 0;
                }
                for (j = col_elements + 1; j >= 0 && coords[1] * col_elements + j > 5 * P; j--) {
                    old[crd][j] = 0;
                }
            }
        }

        // Update maps and check if there is any change
        nchangelocal = 0;

        for (i = 1; i < row_elements + 1; i++) {
            for (j = 1; j < col_elements + 1; j++) {
                oldval = old[i][j];
                newval = oldval;

                // Assign the element with the greatest value from its 4 neighbours
                if (oldval != 0) {
                    if (old[i][j - 1] > newval) newval = old[i][j - 1];
                    if (old[i][j + 1] > newval) newval = old[i][j + 1];
                    if (old[i - 1][j] > newval) newval = old[i - 1][j];
                    if (old[i + 1][j] > newval) newval = old[i + 1][j];

                    if (newval != oldval) {
                        ++nchangelocal;
                    }
                }
                new[i][j] = newval;
            }
        }

        // Compute global number of changes on rank 0
        MPI_Allreduce(&nchangelocal, &nchange, 1, MPI_INT, MPI_SUM, new_comm);

        // Report progress every print frequency number of times
        if (step % printfreq == 0) {
            /*
             * This part is for calculating average value of the whole map, thus verifying the
             * correctness of this program by comparing results with the serial program.
             * Invoking this section every printfreq times is to reduce the influence to performance
             */

            average_value = 0;

            for (i = 1; i < row_elements + 1; i++) {
                for (j = 1; j < col_elements + 1; j++) {
                    average_value += (double) new[i][j];
                }
            }

            MPI_Reduce(&average_value, &average, 1, MPI_DOUBLE, MPI_SUM, 0, new_comm);

            if (rank == 0) {
                average /= (double) (width * height);
                printf("percolate: changes on step %d is %d, average map value: %f\n", step, nchange, average);
            }
    }

    // Break if all maps stop changing
        if (nchange == 0)
            break;

        // Copy back in preparation for next step, omitting halos
        for (i = 1; i < row_elements + 1; i++) {
            for (j = 1; j < col_elements + 1; j++) {
                old[i][j] = new[i][j];
            }
        }
        step++;
    }

    // Record end time of updating
    if (rank == 0) {
        par_end = MPI_Wtime();
    }

    /*
     * For each process, put the final results into map according to their coordinates.
     * Then, apply MAX operation to reduce all maps to every process. The reason why this
     * works is that empty squares can only be increased. Thus, after reducing, squares will
     * be updated as the greatest value, but 0 is still 0.
     */

    // Synchronize values to map
    for (i = rstart; i < rend; i++) {
        for (j = cstart; j < cend; j++) {
            map[i][j] = old[i - rstart + 1][j - cstart + 1];
        }
    }

    // Check if step limit is too small
    if (rank == 0 && nchange != 0) {
        printf("percolate: WARNING max steps = %d reached but nchange != 0\n", maxstep);
    }

    // Reduce updated maps from each process together
    MPI_Allreduce(&map[0][0], &map[0][0], width * height, MPI_INT, MPI_MAX, new_comm);

    // Distribute verification work into processes
    boundary_start = rank * displacement;
    if (rank == size - 1)
        boundary_stop = width;
    else
        boundary_stop = boundary_start + displacement;

    // parallel percolate verification
    for (i = boundary_start; i < boundary_stop; i++) {
        for (j = 0; j < width; j++) {
            if (map[i][0] == map[j][height - 1]) {
                is_percolate = true;
                break;
            }
        }
    }

    // Reduce verification result
    MPI_Allreduce(&is_percolate, &is_percolate, 1, MPI_C_BOOL, MPI_LOR, comm);

    // Deal with output
    if (rank == 0) {
        // Print time cost results
        total_end = MPI_Wtime();
        printf("Total time cost: %f\nMap updating time cost: %f\n", total_end - total_start,
               (par_end - par_start) / step);

        // Draw cluster into pgm file
        if (is_percolate) {
            printf("percolate: cluster DOES percolate\n");
        } else {
            printf("percolate: cluster DOES NOT percolate\n");
        }

        // Write results with specific number of clusters
        mapwrite("map.pgm", map, width, height, 2);
    }

    MPI_Finalize();
    return 0;
}

/******************************************************************************	
*                  Cpt S 411, Introduction to Parallel Computing
*	                          Programming Project #2
* 
* Programmer: Jessica Cuevas
* 
* Goal: this project implements the Conway’s Game of Life in MPI. 
* The game is a cellular automaton containing m*n cells arranged 
* as m rows and n columns.
*******************************************************************************/

/*--------------How to compile------------------//
Login to Leidaes:
ssh jessica.cuevas2@pleiades.tricity.wsu.edu
myWSU password

Write on terminal:
nano golp.c
mpicc -o tst golp.c
ls -t
squeue
sbatch -N "num" -n "num" sub.sh 
nano sub.sh ---> change 
cat slurm-"number".out

scancel ** if want to cancel process

sun.sh
mpirun ~/test 
/----------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <mpi.h>
#include <sys/time.h>
#include <assert.h>
//#include <unistd.h>

//------------------------------------------------------------------------------
#define RAND_MAX 4294967295 
#define MAX_SIZE 32768
//bit manipulation marco functions
#define SetBit(A,k)    ( A[(k/32)] |= (1 << (k%32)) )
#define ClearBit(A,k)  ( A[(k/32)] &= ~(1 << (k%32)) )
#define GetBit(A,k)    ( A[(k/32)] & (1 << (k%32)) )  
#define Index(r,c,n)     ( c + (n*r) )
//------------------------------------------------------------------------------
// eauclidean division algorithm for modular.
int mod(int x, int y)
{
   int t = x - ((x / y) * y);
   if (t < 0) t += y;
   return t;
}

//------------------------------------------------------------------------------
//Checks if cell (bit) is alive (or is one).
bool alive(int *A, int row, int col, int n)
{
    return GetBit(A, Index(row, col, n)) ? 1 : 0;
} 

//------------------------------------------------------------------------------
//Print the whole matrix.
// ASSUMES ROW MAJOR FORMAT
int print_bitmat(int *A, int m, int n)
{
    if ((m * n) / 32 > MAX_SIZE) return -1;
 
    for (int r = 0 ; r <= m - 1 ; ++r)
    {
        for (int c = 0 ; c <= n - 1 ; ++c)
        {
            printf("%d ", alive(A, r, c, n));
        }
        printf("\n");
    }
    
    return 0;
}

//------------------------------------------------------------------------------
//prints a part (or block that is assigned to each processor) of the matrix.
int print_block(int *A, int m, int n)
{
    if ((m * n) / 32 > MAX_SIZE) return -1;
 
    for (int r = 1 ; r < m - 1 ; ++r)
    {
        //printf("ROW %d\n", r);
        for (int c = 0 ; c <= n - 1 ; ++c)
        {
            if (alive(A, r, c, n)) 
            {
                printf("@ ");
            }
            else 
            { 
                printf("- "); 
            }   
        }
        printf("\n");
    }
    
    return 0;
}

//------------------------------------------------------------------------------
// Psuedo-random bitmatix generator with zeros and ones.
// This is an function equavalent to GenerateInitialGoL()
int init_bitmat(int *A, int m, int n, int rank)
{
    if ((m * n) / 32 > MAX_SIZE) return -1;
    
    // For Parallel, include Rank in calculation.
    unsigned int seed = (time(NULL) % MAX_SIZE) * (rank + 1); 

    for (int k = 0 ; k <= m * n - 1 ; ++k)
    {
        if (rand_r(&seed) % 2) SetBit(A, k);
    }

    return 0;
}

//------------------------------------------------------------------------------
// Perform deep copy of B into A (A <-- B). (Both are matrixes))
int deep_copy_bitmat(int *A, int *B, int m, int n)
{
    if ((m * n) / 32 > MAX_SIZE) return -1;

    int i_max = (m * n) % 32 ? (m * n) / 32 : (m * n) / 32 - 1;

    for (int i = 0 ; i <= i_max ; ++i)
    {
        A[i] = B[i];
    }

    return 0;
}

//------------------------------------------------------------------------------
// B is changed, A is referenced. Updates each block with the game of life rules.
// Equivalent to Simulate() function.
int gol_update(int *A, int *B, int m, int n)
{
    assert(m > 2);
    
    int nb_count = 0;
    for (int row = 1; row < m - 1; ++row)
    {
        for (int col = 0; col < n; ++col)
        {
            // Reset for each new cell
            nb_count = 0;
            
            // Count neighbors
            nb_count += alive(A, (row - 1), col, n);             // Up
            nb_count += alive(A, (row - 1), mod(col + 1, n), n); // UpRight
            nb_count += alive(A, row, mod(col + 1, n), n);       // Right
            nb_count += alive(A, (row + 1), mod(col + 1, n), n); // DownRight
            nb_count += alive(A, (row + 1), col, n);             // Down
            nb_count += alive(A, (row + 1), mod(col - 1, n), n); // DownLeft
            nb_count += alive(A, row, mod(col - 1, n), n);       // Left
            nb_count += alive(A, (row - 1), mod(col - 1, n), n); // UpLeft

            if (alive(A, row, col, n))
            {
                if (nb_count < 3 || nb_count > 5)
                {
                    ClearBit(B, Index(row, col, n));
                }
            }
            else
            {
                if (2 < nb_count && nb_count < 6)
                {
                    SetBit(B, Index(row, col, n));
                } 
            }
        }
    }
}

//------------------------------------------------------------------------------
//Basic max function to find the largest int time
int max_func(int *A, int size)
{
    int max = A[0];
    //printf("Time p0: %d\n", max);//DEBUG
    for(int i = 1; i < size; ++i)
    {
        if (A[i] > max) max = A[i];
        //printf("Time p%d: %d\n", i, A[i]);// DEBUG
    }
    
    return max; 
}
//================================= MAIN =======================================
int main(int argc, char* argv[])
{

//------------------------------------------------------------------------------
    //create a struct to store times
    struct timeval t1, t2, ct1, ct2;
//------------------------------------------------------------------------------
    assert(argc == 3);
    // n == number of rows and cols, g == total number of generations, 
    // p == number of processes, rank == id of each processor. 
    int n, g, p, rank, dest, delta, block_height, max_size;
    int tag, res, ct_time;

    int block[MAX_SIZE] = {0};
    int old_block[MAX_SIZE] = {0};
    
    //is neighbor alive? T or F
    bool nb_alive = 0; 
    bool nb_info = 0;

    //initialising all MPI variables
    MPI_Init(&argc,&argv); 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
    MPI_Comm_size(MPI_COMM_WORLD,&p);
    MPI_Status status;

    //taking values from terminals
    n = atoi(argv[1]);
    g = atoi(argv[2]);
    //gen_interval = atoi(argv[3]);

    //checking for assumptions
    assert((n > p) && !(n % p)); // n have to be greater than p and 
                                 // n is divisible by p

    //common eq used that divides procs into rows
    delta = n / p;
    
    // each processor initalizes their own matrix sizes.
    // why delta + 2? we add two for neighboring rows.
    block_height = delta + 2;
    max_size = (block_height * n) % 32 ? (block_height * n) / 32 + 1 : (block_height * n) / 32;


//============================== GOL LOOP ======================================
    
//----------------------- start total timer ------------------------------------
    gettimeofday(&t1,NULL);
//------------------------------------------------------------------------------
    //randomize general matrix
    init_bitmat(block, block_height, n, rank);

    //simulate game of life.
    for(int gen = 0; gen < g; ++gen)
    {
// STEP 0 ----------------------------------------------------------------------
        // 0) Receive / send correct information of neighboring rows.
        // Inner (sending) neighbors are on row indices 1 & delta.
        // Outer (receiving) neighbors are on row indices 0 & delta + 1.
//-----------start comm time ---------------------------------------------------
        gettimeofday(&ct1, NULL);
//------------------------------------------------------------------------------
        for (int col = 0; col < n; ++col)
        {  
            //have to initialize delta again, because macro is been called.
            int bottom = delta + 1;
            
            // Sending (top) [1 --> delta + 1]:
            nb_info = alive(block, 1, col, n);
            dest = mod(rank - 1,  p);
            tag = col + n * dest;
            MPI_Send(&nb_info, 1, MPI_C_BOOL, dest, tag, MPI_COMM_WORLD);
            
            
            // Receiving (bottom):
            MPI_Recv(
                &nb_alive, 
                1, 
                MPI_C_BOOL, 
                MPI_ANY_SOURCE, 
                MPI_ANY_TAG, 
                MPI_COMM_WORLD, &status);
            
            if (nb_alive)
            {
                SetBit(block, Index(bottom, col, n));
            }
            else
            {
                ClearBit(block, Index(bottom, col, n));
            }
 
            // Sending (bottom) [delta -+> 0]:
            nb_info = alive(block, delta, col, n);
            dest = mod(rank + 1,  p);
            tag = col + n * dest;
            MPI_Send(&nb_info, 1, MPI_C_BOOL, dest, tag, MPI_COMM_WORLD);
            

            // Receiving (top):
            MPI_Recv(
                &nb_alive, 
                1, 
                MPI_C_BOOL, 
                MPI_ANY_SOURCE, 
                MPI_ANY_TAG, 
                MPI_COMM_WORLD, &status);
            
            if (nb_alive)
            {
                SetBit(block, Index(0, col, n));
            }
            else
            {
                ClearBit(block, Index(0, col, n));
            }
        }
//-------------------end comm time ---------------------------------------------
        gettimeofday(&ct2, NULL);
        ct_time += ((ct2.tv_sec - ct1.tv_sec) * 1000 + (ct2.tv_usec - ct1.tv_usec) / 1000);
//STEP 1------------------------------------------------------------------------
        // 1) Store memory of current block.
        deep_copy_bitmat(old_block, block, block_height, n);

//STEP 2------------------------------------------------------------------------
        // 2) Update status of inner block cells based on GoL rules.
        // aka simulates the game of life
        gol_update(old_block, block, block_height, n); 

        // BARRIER: Wait until all gol updates are complete
        MPI_Barrier(MPI_COMM_WORLD);
//STEP 3------------------------------------------------------------------------   
        // 3) Print the combined matrix in order 
        // Equivalent to function DisplayGoL() 
       /* if (gen % gen_interval == 0 || gen == g - 1)
        {
            //Rank 0 receives all updated blocks from other ranks 
            if (rank == 0)
            {
                printf("Generation: %d\n", gen);
                print_block(block, block_height, n);
                
                int* block_buf;
                
                for (int proc = 1; proc < p; ++proc)
                {
                    //allocating the memory based on the size max_size times
                    // the total number of processors            
                    block_buf = (int*)malloc(sizeof(int) * max_size);

                    MPI_Recv(
                    block_buf,
                    max_size,
                    MPI_INT, 
                    proc, 
                    MPI_ANY_TAG, 
                    MPI_COMM_WORLD,&status);

                    print_block(block_buf, block_height, n);

                    free(block_buf);
                }
            }
            else
            {
                // processors send to rank 0
                MPI_Send(
                    &block, 
                    max_size, 
                    MPI_INT, 
                    0, 
                    0, 
                    MPI_COMM_WORLD);
            }
        }*/
    }
//----------------------- end total timer --------------------------------------
    gettimeofday(&t2,NULL);

//-------------------------calculate total runtime -----------------------------
    int totalRuntime = ((t2.tv_sec-t1.tv_sec) * 1000 + (t2.tv_usec-t1.tv_usec) / 1000);

//Rank 0 receives all total, and comm runtimes----------------------------------
    if (rank == 0)
    {

        int* time_buf = (int*)malloc(sizeof(int) * p);
        int* comm_buf = (int*)malloc(sizeof(int) * p);

        time_buf[0] = totalRuntime;
        comm_buf[0] = ct_time;

        int time_box;
        
        for (int proc = 1; proc < p; ++proc)
        {
            MPI_Recv(
            &time_box,
            1,
            MPI_INT, 
            proc, 
            0, 
            MPI_COMM_WORLD,&status);

            time_buf[proc] = time_box;
            //printf("proc time: %d\n",time_box); //DEBUG             

            MPI_Recv(
            &time_box,
            1,
            MPI_INT, 
            proc, 
            1, 
            MPI_COMM_WORLD,&status);
            
            comm_buf[proc] = time_box;
        }

        int max_time = max_func(time_buf, p);
        int max_ct = max_func(comm_buf,p); 
        float ave_time = (float)max_time / (float)g;

        printf("\nTotal runtime: %d milliseconds\n", max_time);
        printf("\nAverage time per generation: %f milliseconds\n", ave_time);
        printf("\nTotal communication time: %d milliseconds\n", max_ct);
        printf("\nTotal computation time: %d milliseconds\n", max_time - max_ct);
   
        free(time_buf);
        free(comm_buf);
    }
    else // other ranks send their total and comm runtimes to rank 0
    {
        
        MPI_Send(
            &totalRuntime, 
            1, 
            MPI_INT, 
            0, 
            0, 
            MPI_COMM_WORLD);

        MPI_Send(
            &ct_time, 
            1, 
            MPI_INT, 
            0, 
            1, 
            MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}

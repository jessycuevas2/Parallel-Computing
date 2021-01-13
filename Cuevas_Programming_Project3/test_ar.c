/******************************************************************************	
*                  Cpt S 411, Introduction to Parallel Computing
*	                          Programming Project #3
* 
* Programmer: Jessica Cuevas
* 
* Goal: Implement and compare the performace of three different versions 
* of the parallel reduction operation.
* (myAllReduce(), naiveAllReduce(), and MPI_Allreduce(). )
* 
*******************************************************************************/

            /*--------------How to compile------------------//
            Login to Leidaes:
            ssh jessica.cuevas2@pleiades.tricity.wsu.edu
            myWSU password

            Write on terminal:
            nano golp.c
            mpicc -o out all_reduce.c
            ls -t
            squeue
            sbatch -N "num" -n "num" sub.sh 
            nano sub.sh ---> change 
            cat slurm-"number".out

            scancel ** if want to cancel process

            sun.sh
            mpirun ~/PA3/out
            /----------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
//#include <mpi.h>
#include <sys/time.h>
#include <assert.h>

#define MAX_SIZE 32768
#define GetBit(n,k)    ( n & (1 << (k%32)) )


int generate_block(int *A, int size, int rank)
{
    unsigned int seed = (time(NULL) % MAX_SIZE) * (rank + 1); 

    for (int i = 0; i < size; ++i)
    {
        A[i] = rand_r(&seed) % 101;
    }

    return 0;
}

// 0 <= k <= 31 (indexed bit from right)
bool check_bit(int n, int k)
{
    return GetBit(n, k) ? 1 : 0;
}

bool is_pow_two(int n)
{
    int bit_count = 0;

    for (int k = 0; k < 32; ++k)
    {
        bit_count += check_bit(n, k);
    }

    if (bit_count == 1) return true;
    
    return false;
}


// Log base 2 of power of 2 number
int log_pow_two(int n)
{
    if (n == 0) return -1;

    int k = 0;
    while (!check_bit(n, k))
    {
        ++k;
    }
    
    return k;
}

int max_func(int x, int y)
{
    return x ^ ((x ^ y) & -(x < y));
}

/*int max_fo(int x, int y)
{
    return (x > y) ? x : y; 
}*/

//------------------------------------------------------------------------------
// the binary associative operator function. 
int bop(int x, int y)
{
    int sum;
    //int max;

    sum = x + y;
    //max = max_func(x ,y);

    return sum;
    //return max;   
}
//------------------------------------------------------------------------------
// Applies operation on each block and returns output value (x_i). 
int calc_output_x(int *A, int size)
{
    if (size <= 0) return 0;
    if (size == 1) return A[0];

    int output = bop(A[0], A[1]);

    for (int i = 2; i < size; ++i)
    {
        output = bop(output, A[i]);        
    }

    return output;
}

void print_array(int *A, int size)
{
    printf("[");

    for (int i = 0; i < (size - 1); ++i)
    {
        printf("%d ", A[i]);
    }
    printf("%d]\n", A[size - 1]);
}

int main(int argc, char* argv[])
{
    //n == input array size, p == total number of processes
    // rank == processor ID, delta == n / p
    int n, p, rank, delta, log_p, x;

    assert(argc == 2); 
    //taking array size from terminal
    n = atoi(argv[1]);

    //initialising all MPI variables
    //MPI_Init(&argc,&argv); 
    //MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
    //MPI_Comm_size(MPI_COMM_WORLD,&p);
    //MPI_Status status;
    rank = 0;
    p = 1;

    //Check all assumptions apply
    assert((n > p) &&  ((n % p) == 0) && is_pow_two(p));
    // n have to be greater than p and 
    // n is a multiple of p

    //common eq used that divides procs into blocks
    delta = n / p;
    log_p = log_pow_two(p);

    //STEP 0: generate array A[0..n-1]
    //array
    int *block = (int*)malloc(sizeof(int) * delta);
    generate_block(block, delta, rank);
    
    print_array(block, delta);
    x = calc_output_x(block, delta);
    printf("%d\n", x);
    //int m1 = max_f(1, 5);
    //int m2 = max_fo(1, 5);
    //printf(" %d    %d\n", m1, m2);
    //printf(" %d\n", log_p);

    //STEP 1: Implement MyAllReduce()

    //MPI_Finalize();

    return 0;
}

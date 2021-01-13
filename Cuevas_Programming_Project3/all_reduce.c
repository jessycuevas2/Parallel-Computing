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
#include <mpi.h>
#include <sys/time.h>
#include <assert.h>

#define MAX_SIZE 32768
#define GetBit(n,k)    ( n & (1 << (k%32)) )


//------------------------------------------------------------------------------
//generates random blocks of an array.
// Equivalent to GenerateArray()

int generate_block(int *A, int size, int rank)
{
    unsigned int seed = (time(NULL) % MAX_SIZE) * (rank + 1); 

    for (int i = 0; i < size; ++i)
    {
        A[i] = rand_r(&seed) % 11; // generate ran numbers, range 1-100
    }

    return 0;
}

/*
int generate_block(int *A, int size, int rank) //DEBUG
{
    for (int i = 0; i < size; ++i)
    {
        A[i] = rank * size + i; // generate ran numbers, range 1-100
    }

    return 0;
}*/

//------------------------------------------------------------------------------
// 0 <= k <= 31 (indexed bit from right)
bool check_bit(int n, int k)
{
    return GetBit(n, k) ? 1 : 0;
}

//------------------------------------------------------------------------------
//checks if p is a power of 2. 
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

//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
// eauclidean division algorithm for modular.
int mod(int x, int y)
{
   int t = x - ((x / y) * y);
   if (t < 0) t += y;
   return t;
}

//------------------------------------------------------------------------------
// basic max function to compare two integers.
int max_func(int x, int y)
{
    return x ^ ((x ^ y) & -(x < y));
}

//------------------------------------------------------------------------------
// the binary associative operator function. 
int bop(int x, int y)
{
    int val;

    val = x + y;
    //val = max_func(x, y);
    //val = y;

    return val;
       
}

//------------------------------------------------------------------------------
// Applies operation on each block and returns output value (x_i). 
int output_x(int *A, int size)
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

//------------------------------------------------------------------------------
int get_partner(int t, int rank)
{
    int value;
    value = 1 << t;
    return rank ^ value;
}

//------------------------------------------------------------------------------
void print_array(int *A, int size)
{
    printf("[");

    for (int i = 0; i < (size - 1); ++i)
    {
        printf("%d ", A[i]);
    }
    printf("%d]\n", A[size - 1]);
}

//------------------------------------------------------------------------------
// myAllReduce() implementation
int my_all_reduce(int rank, int x_i, int k, MPI_Status *status_ptr)
{
    int x_t, partner_rank, y;
    x_t = x_i;
    //printf("rank: %d, x_i = %d\n",rank, x_t); //DEBUG

    for (int t = 0; t < k; ++t)
    {
        partner_rank = get_partner(t, rank);
        //printf("t = %d; rank = %d; partner = %d\n",t,rank,partner_rank); //DEBUG

        MPI_Send(
                &x_t, 
                1, 
                MPI_INT, 
                partner_rank, 
                0, 
                MPI_COMM_WORLD);

        MPI_Recv(
                &y, 
                1, 
                MPI_INT, 
                partner_rank, 
                MPI_ANY_TAG, 
                MPI_COMM_WORLD, 
                status_ptr);

        if (rank < partner_rank) 
        {
            x_t = bop(x_t, y);
        }
        else
        {
            x_t = bop(y, x_t);
        }
        //printf("rank: %d, x_%d = %d\n",rank, t, x_t); //DEBUG
    }
    
    return x_t;
    
}

//------------------------------------------------------------------------------
// NaiveAllReduce() implementation: Array/Bus permutation all reduce.  
int naive_all_reduce(int rank, int x_i, int p, MPI_Status *status_ptr)
{
    int y, new_x, recv_rank, send_rank;

    new_x = x_i;

    if (p == 1)
    {
        return new_x;
    }

    if (rank == 0)
    {
        //SEND TO RANK 1
        MPI_Send(
            &new_x, 
            1, 
            MPI_INT, 
            1, 
            0, 
            MPI_COMM_WORLD);
        
        recv_rank = p - 1;

        MPI_Recv(
                &y, 
                1, 
                MPI_INT, 
                recv_rank, 
                MPI_ANY_TAG, 
                MPI_COMM_WORLD, 
                status_ptr);

        new_x = y;     
        //printf("sending value %d to rank 1\n", new_x); //DEBUG

    }
    else // if rank > 0
    {
        //which rank send it? 
        recv_rank = rank - 1;
        
        //receive x_prime(product of operation from previous rank) 
        MPI_Recv(
                &y, 
                1, 
                MPI_INT, 
                recv_rank, 
                MPI_ANY_TAG, 
                MPI_COMM_WORLD, 
                status_ptr);

        //printf("receiving %d from rank %d\n", y, recv_rank); //DEBUG
        
        // then computes new x_prime
        new_x = bop(y, x_i);

        // then sends new x_prime to the next rank (rank + 1) % p 
        send_rank = mod(rank + 1, p);

        MPI_Send(
            &new_x, 
            1, 
            MPI_INT, 
            send_rank, 
            0, 
            MPI_COMM_WORLD);

        //printf("sending %d to rank %d\n", new_x, send_rank); //DEBUG
        // Note: the mod p will make the last rank send the final x to rank 0 
    }
    
    return new_x;
}

//------------------------------------------------------------------------------
int mpi_library_all_reduce(int x_i)
{
    int output;

    MPI_Allreduce(
        &x_i,
        &output,
        1,
        MPI_INT,
        MPI_SUM, //change to MPI_MAX
        MPI_COMM_WORLD);

    return output;

}

//------------------------------------------------------------------------------
//Basic max function to find the largest int time
int max(int *A, int size)
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
    struct timeval t1, t2, t3, t4, t5, t6;
//------------------------------------------------------------------------------
    // n == input array size,    p == total number of processes
    // rank == processor ID,     delta == n / p
    int n, p, rank, delta, log_p, x;
    int all_red, n_red, mpi_red, all_time, n_time, mpi_time;

// check size n is written on terminal.-----------------------------------------
    assert(argc == 2); 
    n = atoi(argv[1]);

// initialising all MPI variables.----------------------------------------------
    MPI_Init(&argc, &argv); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Status status;

// Check all conditions before before running parallel program.-----------------
    assert((n > p) &&  ((n % p) == 0) && is_pow_two(p));

//common eq used that divides procs into blocks & log of p ---------------------
    delta = n / p;

    log_p = log_pow_two(p);

// STEP 0: generate array (A[0..n-1]) block for each processor.-----------------
    int *block = (int*)malloc(sizeof(int) * delta);

    generate_block(block, delta, rank);
    
// STEP 1: Implement local operations.------------------------------------------

    x = output_x(block, delta);   

// STEP 2: Implement MyAllReduce()----------------------------------------------

    //gettimeofday(&t1,NULL);
    //all_red = my_all_reduce(rank, x, log_p, &status);
    //gettimeofday(&t2,NULL);
    //all_time = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
    
// STEP 3: Implement NaiveAllReduce()-------------------------------------------

    gettimeofday(&t3,NULL);
    n_red = naive_all_reduce(rank, x, p, &status);
    gettimeofday(&t4,NULL);
    n_time = (t4.tv_sec-t3.tv_sec)*1000000 + (t4.tv_usec-t3.tv_usec);

// STEP 4: Implement  MPILibraryAllReduce()-------------------------------------
    //gettimeofday(&t5,NULL);
    //mpi_red = mpi_library_all_reduce(x);
    //gettimeofday(&t6,NULL);
    //mpi_time = (t6.tv_sec-t5.tv_sec)*1000000 + (t6.tv_usec-t5.tv_usec);
// PRINT -----------------------------------------------------------------------  
    if (rank == 0)
    {
        int *A = (int*)malloc(sizeof(int) * n);
        int *block_buf;
        int index = 0;

        //place rank 0 block into array A[]
        for (int i = 0; i < delta; ++i)
        {
            A[index] = block[i];
            ++index;
        }
        
        for (int proc = 1; proc < p; ++proc)
        {
            block_buf = (int*)malloc(sizeof(int) * delta);

            MPI_Recv(
                    block_buf,
                    delta,
                    MPI_INT,
                    proc, 
                    MPI_ANY_TAG, 
                    MPI_COMM_WORLD, 
                    &status);

            //print_array(block_buf,delta); //DEBUG

            for (int i = 0; i < delta; ++i)
            {
                A[index] = block_buf[i];
                ++index;
            }

            free(block_buf);   
        }
        //printf("Array: \n");
        //print_array(A, n);       
        
        //RESULTS
        //printf("My All Reduce final output: %d\n", all_red);
        //printf("Naive all reduce final output: %d\n", n_red);
        //printf("MPI all reduce final output: %d\n", mpi_red); 

        free(A);
    }
    else
    {
        MPI_Send(
                block, 
                delta, 
                MPI_INT, 
                0, 
                0, 
                MPI_COMM_WORLD);
    }
//------------------------------------------------------------------------------
    if (rank == 0)
    {
        int *time_buf = (int*)malloc(sizeof(int) * p);
        
        //time_buf[0] = all_time;
        time_buf[0] = mpi_time;
        int time_box;

        //printf("proc time: %d\n",time_buf[0]); //DEBUG

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
        }

        int max_time = max(time_buf, p);
        //=================== Timing Testing =============================
        //printf("\nTotal runtime: %.2f microseconds\n", (float)max_time);
        printf("\nTotal runtime: %.0f microseconds\n", (float)max_time);
        //printf("\nTotal runtime: %.0f microseconds\n", (float)max_time);

        free(time_buf);
    }
   else // other ranks send their total runtimes to rank 0
    {
        
        MPI_Send(
            &n_time, 
            1, 
            MPI_INT, 
            0, 
            0, 
            MPI_COMM_WORLD);
    }

    free(block);
    MPI_Finalize();

    return 0;
}

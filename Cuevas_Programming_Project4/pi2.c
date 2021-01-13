/*******************************************************************************	
*                  Cpt S 411, Introduction to Parallel Computing
*	                          Programming Project #4
* 
* Programmer: Jessica Cuevas
* 
* Goal: implement an OpenMP multithreaded PI value estimator using the PI
* algorithm discussed in class.
*
* This algorithm essentially throws a dart n times into a unit square 
* and computes the fraction of times that dart falls into the
* embedded unit circle. This fraction multiplied by 4 gives an estimate for PI.
*
* Note: the precision of the PI estimate could potentially vary with the number
*of threads (assuming n is fixed).
*******************************************************************************/

            /*--------------How to compile------------------//
            Login to Leidaes:
            ssh jessica.cuevas2@pleiades.tricity.wsu.edu
            myWSU password

            Write on terminal:
            nano pi.c
            (mpicc) gcc -o out -fopenmp pi.c
            ls -t
            squeue
            sbatch -N 1 sub.sh
            nano sub.sh ---> change 
            cat slurm-"number".out

            scancel ** if want to cancel process

            sun.sh
            mpirun ~/PA4/out
            /----------------------------------------------*/

// input: iterations "n" and number of threads " th"
// output: the PI value and total time 
  
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

//------------------------------------------------------------------------------
//toss dart simulator at a squared board with a target circle.
long int count_in_circle(long int t_tosses)
{

    long int count = 0; //counter, and iterator
    double x, y; // coordinates
    unsigned int seed = ((unsigned) time(NULL));
    srand(seed);
    
    for(long int i = 0; i < t_tosses; ++i)
    {
        //set random numbers within squareboard
        x = rand_r(&seed)/ (double) RAND_MAX;
        y = rand_r(&seed)/ (double) RAND_MAX;
    
        //check if coordinates is within circle, if it is keep count.        
        if((x * x + y * y) < 1.0) count++;
    }

    //return 4.0 * count / (double) t_tosses;
    return count;
}

double run_test(long int n, long int t_block, int t_threads)
{
    long int in_count;
    in_count = 0;

    //get start timer-----------------------------------------------------------
    double dif;
    double start = omp_get_wtime();

    #pragma omp parallel for shared(t_block, t_threads) reduction(+:in_count)
    for (long int i = 0; i < t_threads; ++i)
    {
        in_count += count_in_circle(t_block);
    }
    
    
    double pi = (double) in_count / n * 4.0;

    //End of timer--------------------------------------------------------------
    double end = omp_get_wtime();

    dif = end - start;

    printf("Pi estimate is: %.16f\n", pi);
    printf("\n %f seconds \n ", dif);

    return dif;
}
//================================= MAIN =====================================//
int main(int argc, char* argv[])
{
    // n == number of dart tosses/throws, how many darts you throw?
    int t_threads;

    long int n, t_block;

    double test_time;
    //double PI25DT = 3.141592653589793238462643;    

    //Get input
/*
    assert(argc == 3); 
    n = atoi(argv[1]);
    t_threads = atoi(argv[2]);

    //check t_thread is greater than 0, and iternations is creater than 1
    //and check t_threads divides n
    assert(t_threads > 0 && n > 0 && n % t_threads == 0);

    //Initialize threads
    omp_set_num_threads(t_threads);

    //distribute number of tosses pre thread.
    t_block = n / t_threads;
*/
    // FILE OPENED OUTSIDE
    FILE *fp;
    fp = fopen("results.csv", "w");

    // WRITE HEADER "n,time,threads\n"
    fputs("threads,n,time\n", fp);

    // RUN EXPERIMENT (long int n, int t, FILE *outfile)
    // DOUBLE FOR LOOP, t_threads: 1 -> 8 (pow 2 or linear),
    // n: 840 -> 840 * (1<<21) (linear) or (1<<10) to (1<<30)
    for (int i = 0; i < 4; ++i)
    {
        t_threads = (1 << i);
        omp_set_num_threads(t_threads);
        
        for (int j = 10; j < 31; ++j)
        {
            n = (1 << j);
            t_block = n / t_threads;
            test_time = run_test(n, t_block, t_threads);
            fprintf(fp, "%d,%ld,%lf\n", t_threads, n, test_time);
        }
    } 
    
/*
    in_count = 0;
    
    //get start timer-----------------------------------------------------------
    double dif;
    double start = omp_get_wtime();

    #pragma omp parallel for shared(t_block, t_threads) reduction(+:in_count)
    for (long int i = 0; i < t_threads; ++i)
    {
        in_count += count_in_circle(t_block);
    }
    
    pi = (double) in_count / n * 4.0;

    //End of timer--------------------------------------------------------------
    double end = omp_get_wtime();

    dif = end - start;

    printf("Pi estimate is: %.16f\n", pi);
    printf("\n %f seconds \n ", dif);
*/
    // WRITE fprintf "%ld,%lf,%d\n", n, time, t_threads
    
    fclose(fp);
    return 0;
}

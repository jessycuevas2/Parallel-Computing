/******************************************************************************	
*                  Cpt S 411, Introduction to Parallel Computing
*	                          Programming Project #1
* 
* Programmer: Jessica Cuevas
* 
* Goal: The goal of this project is to empirically estimate the network 
* parameters latency and bandwidth, and the network buffer size, 
* for the network connecting the nodes of the compute cluster. 
* 
* Send receive test:
* rank 1 sends to rank 0 (all other ranks sit idle)
* For timing use of C gettimeofday() is recommended.
*******************************************************************************/

/*--------------How to compile------------------//
Login to Leidaes:
ssh jessica.cuevas2@pleiades.tricity.wsu.edu
myWSU password

Write on terminal:
nano send_recv_test.c
mpicc -o send_recv_test send_recv_test.c
ls -t
squeue
sbatch -N 2 -n 2 sub.sh 
squeue
cat slurm-"number".out

scancel ** if want to cancel process

sun.sh
mpirun ~/test 
/----------------------------------------------*/

#include <stdio.h>
#include <mpi.h> 
#include <assert.h>
#include <sys/time.h>
#include <stdlib.h>

/******************* MPI Blockig Send/Receive program *********************/
int main(int argc,char *argv[])
{
    // rank represents the id of a processor, p represents number of processes
    int rank,p; 
    
    //create a struct ot store times
    struct timeval t1,t2;
    
    //initialize MPI system, rank and size
    MPI_Init(&argc,&argv); // creates ids
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
    MPI_Comm_size(MPI_COMM_WORLD,&p);

    assert(p>=2);

    //int size=1;
    int loops=10; // 10 sends 10 receives.
    int i = 0;
    int j = 0;

//------------------------------- SENDER -------------------------------------//
    if(rank==1)
    {
        //destination tag
        int dest = 0;

        char test_x = 0;
        
        //----------------------LATENCY TEST----------------------------------//
        //start timer
        gettimeofday(&t1,NULL);

        MPI_Send(&test_x,1,MPI_CHAR,dest,0,MPI_COMM_WORLD);

        //End of timer.
        gettimeofday(&t2,NULL);

        int tSend = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
        printf("%d,%d,%d\n", rank,1,tSend);
        //----------------------LATENCY TEST----------------------------------//

        for (j= 0; j <= 20; ++j)
        {
            //size will be sifted to be power of two.
            size = 1 << j; 
            //allocate memory for the message that will be sent
            char *x = (char*) malloc(sizeof(char)*size);
            
            //start timer
            gettimeofday(&t1,NULL);
            //for loop to repeat sending 10 times to calcualte average time. 
            for(i = 0; i < loops; ++i)
            {
                MPI_Send(x,size,MPI_CHAR,dest,0,MPI_COMM_WORLD);
                /*MPI Send has the following parameters:
                (1) initial address of send buffer (choice)
                (2) number of elements in send buffer,size of message(noneg int)
                (3) datatype of each send buffer element (handle)
                (4) rank or id of destination (int)
                (5) processor the message is going to be sent to (int)
                (6) communicator (handle)
                */
            }
            //End of timer
            gettimeofday(&t2,NULL);
            //Total time in millisecs
            tSend = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
            printf("%d,%d,%d\n", rank,size,tSend/loops);
            //deallocate memory so we can assign a new size for message x
            free(x)
        }
    }
    
//-----------------------------RECEIVER --------------------------------------//
    else if (rank==0)
    {
        //value for MPI system 
        MPI_Status status;

        char test_y = 0;
        
        // LATENCY TEST
        //start timer
        gettimeofday(&t1,NULL);

        MPI_Recv(&test_y,1,MPI_CHAR,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);

        //End timer
        gettimeofday(&t2,NULL);

        int tRecv = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
        printf("%d,%d,%d\n", rank,1,tRecv);

        for (j= 0; j <= 20; ++j)
        {
            //size will be sifted to be power of two.
            size = 1 << j;
            //allocate memory for the message that will be sent
            char *y = (char*) malloc(sizeof(char)*size);
            
            //start timer
            gettimeofday(&t1,NULL);
            //for loop to repeat recieved 10 times to calcualte average time. 
            for(i=0; i < loops; ++i)
            {
                MPI_Recv(y,size,MPI_CHAR,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
            }
            //End timer
            gettimeofday(&t2,NULL);
            //Total time in millisecs
            tRecv = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
            // printf("ZERO\n\n");
            printf("%d,%d,%d\n", rank,size,tRecv/loops);
            //deallocate memory so we can assign a new size for message y
            free(y)
        }
    }

    //this terminates a paralled programm---------------------------------------
    MPI_Finalize();
}

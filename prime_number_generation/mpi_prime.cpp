#include <iostream>
#include <algorithm>
#include <cmath>
#include "mpi.h"

#define LIMIT     2500000+11    /* Increase this to find more primes */
#define PRINT     100000      /* Print a line after this many numbers */
#define MASTER 0

bool isprime(const int& n){
    int sqr_root;
    if(n>10){
        sqr_root = (int) std::sqrt(n);
        for(int i=3; i<=sqr_root; i+=2){
            if((n%i)==0){return false;}
        }
        return true;
    }
    else{
        return false;
    }
}



int main(int argc, char* argv[]){

    MPI_Init(&argc, &argv);
    int comm_size, task_rank, local_prime_counter{0}, local_recent_prime, local_chunk_size, local_start_point, local_end_point;
    int prime_counter, recent_prime;

    MPI_Comm_size( MPI_COMM_WORLD , &comm_size);
    MPI_Comm_rank( MPI_COMM_WORLD , &task_rank);

    if(task_rank == MASTER){
        std::cout<<"Starting! Numbers to be scanned= "<<LIMIT<<std::endl;
    }

    local_chunk_size = (LIMIT-11)/comm_size;

    if(local_chunk_size%2==1){
        if(task_rank%2==1){
            local_start_point = 11+(task_rank*local_chunk_size)+1;
        }
        else{
            local_start_point = 11+(task_rank*local_chunk_size);
        }
    }
    else{
        local_start_point = 11+(task_rank*local_chunk_size);
    }


    if(task_rank==comm_size-1){
        local_end_point = LIMIT+1;
    }
    else{
        local_end_point = local_start_point+local_chunk_size;
    }

    for(int i=local_start_point; i<local_end_point; i+=2){
        if(isprime(i)){
            ++local_prime_counter;
            local_recent_prime=i;
        }
    }

    MPI_Reduce( &local_prime_counter, &prime_counter, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
    MPI_Reduce( &local_recent_prime, &recent_prime, 1, MPI_INT, MPI_MAX, MASTER, MPI_COMM_WORLD);

    if(task_rank==MASTER){
        std::cout<<"Done. Largest prime is "<<recent_prime<<". Total number of primes are "<<prime_counter+4<<"."<<std::endl;
    }
    MPI_Finalize();
    return 0;
}
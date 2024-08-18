#include <iostream>
#include "mpi.h"
#define MASTER 0

int main(int argc, char* argv[]){
    int num_tasks, task_rank, message;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_rank);

    MPI_Request request;

    if(num_tasks%2!=0){
        if(task_rank == MASTER){
            std::cout<<"Not doing any communicatioin as there are Odd number of tasks in the communicator"<<std::endl;
        }
    }
    else{
        if(task_rank<num_tasks/2){
            int partner = (num_tasks/2)+task_rank;
            MPI_Isend( &task_rank, 1, MPI_INT, partner, 000, MPI_COMM_WORLD, &request);
            MPI_Irecv( &message, 1, MPI_INT, partner, 000, MPI_COMM_WORLD, &request);
            std::cout<<"Message: \""<<message<<"\", is sent to "<<partner<<" by task: "<<task_rank<<std::endl;
        }
        else{
            int partner =task_rank - num_tasks/2;
            MPI_Isend( &task_rank, 1, MPI_INT, partner, 000, MPI_COMM_WORLD, &request);
            MPI_Irecv( &message, 1, MPI_INT, partner, 000, MPI_COMM_WORLD, &request);
            std::cout<<"Message: \""<<message<<"\", is sent to "<<partner<<" by task: "<<task_rank<<std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}
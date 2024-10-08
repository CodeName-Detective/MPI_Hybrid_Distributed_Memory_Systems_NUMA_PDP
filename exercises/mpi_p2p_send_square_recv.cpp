#include <iostream>
#include "mpi.h"
#include <cmath>

int main(int argc, char* argv[]){
    MPI_Init(&argc , &argv);

    int task_id;
    MPI_Comm_rank(MPI_COMM_WORLD , &task_id);

    if(task_id==0){
        int num_tasks;
        MPI_Comm_size(MPI_COMM_WORLD , &num_tasks);

        for(int i=1; i<num_tasks; ++i){
            MPI_Send( &i, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

        for(int i=1; i<num_tasks; ++i){
            int buffer;
            MPI_Recv( &buffer, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::cout<<"Recieved: "<<buffer<<" from process "<<i<<std::endl;
        }
    }
    else{
        int buffer;
        MPI_Recv( &buffer, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        buffer = std::pow(buffer,2);
        MPI_Send(&buffer, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
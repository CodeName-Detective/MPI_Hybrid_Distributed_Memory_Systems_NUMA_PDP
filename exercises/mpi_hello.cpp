#include <iostream>
#include "mpi.h"

#define MASTER 0

int main(int argc, char* argv[]){
    int num_tasks, task_rank, len_process_name;
    char host_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_rank);

    MPI_Get_processor_name(host_name, &len_process_name);
    
    if (task_rank == MASTER)
    {
        std::cout<<"MASTER: Number of MPI tasks is: "<<num_tasks<<std::endl;
    }

    std::cout<<"Printing from host: "<<host_name<< " and task: "<<task_rank<<"/"<<num_tasks<<std::endl;

    MPI_Finalize();
    return 0;
}

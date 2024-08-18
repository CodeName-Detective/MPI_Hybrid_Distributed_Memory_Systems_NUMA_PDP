#include <iostream>
#include <memory>
#include "mpi.h"

#define MASTER 0
#define ARRAY_SIZE 16

void const print_array(const float* const array, const int dim){
    for(int i=0; i<dim; ++i){
        std::cout<<array[i]<<" ";
    }
    std::cout<<std::endl;
}

int main(int argc, char* argv[]){
    int num_tasks, task_rank;

    MPI_Init( &argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_rank);

    int buffer_size = ARRAY_SIZE/num_tasks;

    float* array;

    float* buffer = new float[buffer_size];

    if (task_rank == MASTER){
        array = new float[ARRAY_SIZE];
        for(int i=0; i<ARRAY_SIZE; ++i){
            array[i] = i*1.0;
        }
        print_array(array, ARRAY_SIZE);
    }

    MPI_Scatter( array, buffer_size, MPI_FLOAT, buffer, buffer_size, MPI_FLOAT, MASTER, MPI_COMM_WORLD);

    for(int local_idx=0; local_idx<buffer_size; ++local_idx){
        int global_idx = task_rank*buffer_size+local_idx;
        buffer[local_idx] = buffer[local_idx] + global_idx*1.0;
    }

    MPI_Gather( buffer, buffer_size, MPI_FLOAT, array, buffer_size, MPI_FLOAT, MASTER, MPI_COMM_WORLD);

    delete[] buffer;

    if (task_rank == MASTER){
        print_array(array, ARRAY_SIZE);
        delete[] array;
    }

    MPI_Finalize();
    return 0;
}

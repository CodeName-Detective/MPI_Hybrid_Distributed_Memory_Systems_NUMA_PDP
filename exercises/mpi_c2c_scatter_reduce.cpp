#include <iostream>
#include <memory>
#include <random>

#include "mpi.h"

int main(int argc, char* argv[]){
    MPI_Init(&argc, &argv);

    int num_size;
    MPI_Comm_size( MPI_COMM_WORLD, &num_size);

    const int num_elements = 1 << 10; // Binary 1 left shift 1 by 10 positions.
    const int chunk_size = num_elements/num_size;

    int task_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);

    std::unique_ptr<int[]>send_ptr;

    if(task_id==0){
        send_ptr = std::make_unique<int[]>(num_elements);

        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_int_distribution dist(1,1);

        std::generate(send_ptr.get(), send_ptr.get()+num_elements, [&] {return dist(mt);});
    }

    auto recv_buffer = std::make_unique<int[]>(chunk_size);

    // Scatter the data
    MPI_Scatter(send_ptr.get(), chunk_size, MPI_INT, recv_buffer.get(), chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Local Reduction
    auto local_sum = std::reduce(recv_buffer.get(), recv_buffer.get()+chunk_size);

    int global_sum;

    MPI_Reduce( &local_sum, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(task_id==0){
        std::cout<<"At process: "<<task_id<<" the global sum is "<<global_sum<<std::endl;
    }
    else if(task_id==1){
        std::cout<<"At process: "<<task_id<<" the global sum is "<<global_sum<<std::endl;
    }

    MPI_Finalize();
    return 0;
}
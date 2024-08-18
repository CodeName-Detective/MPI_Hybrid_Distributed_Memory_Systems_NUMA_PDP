#include <iostream>
#include "mpi.h"

#define MASTER 0


int main(int argc, char* argv[]){
    MPI_Init(&argc, &argv);
    int comm_world_size, task_world_rank;

    MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);
    MPI_Comm_rank( MPI_COMM_WORLD, &task_world_rank);

    int x;
    double y;

    int int_size, double_size, total_size;
    MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &int_size);
    MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &double_size);
    total_size = int_size+double_size+MPI_BSEND_OVERHEAD; // MPI_BSEND_OVERHEAD -> Additional Overhead.
    char* pack_buffer = new char[total_size]; //4 Bytes for integer and 8 bytes for double

    int pack_position{0};

    if(task_world_rank==MASTER){
        x = 9;
        y = 3.14325;

        MPI_Pack( &x, 1, MPI_INT, pack_buffer,total_size, &pack_position, MPI_COMM_WORLD);
        MPI_Pack( &y, 1, MPI_DOUBLE, pack_buffer,total_size, &pack_position, MPI_COMM_WORLD);
        MPI_Bcast( pack_buffer, total_size, MPI_PACKED, MASTER, MPI_COMM_WORLD);
    }
    else{
        MPI_Bcast( pack_buffer, total_size, MPI_PACKED, MASTER, MPI_COMM_WORLD);
        MPI_Unpack( pack_buffer, total_size, &pack_position, &x, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack( pack_buffer, total_size, &pack_position, &y, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    }
    std::cout<<"For rank: "<<task_world_rank<<" the value of integer: "<<x<<" the value of the double: "<<y<<std::endl;

    delete[] pack_buffer;
    MPI_Finalize();
    return 0;
}
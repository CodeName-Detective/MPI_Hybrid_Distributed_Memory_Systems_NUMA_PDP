#include <iostream>
#include "mpi.h"

#define MASTER 0

struct Variable{
        int x;
        double y;
    };


void create_variable_struct_datatype(MPI_Datatype& mpi_var_struct_type){

    const int number_of_blocks{2}; // use const other wise will get error.
    // Defining the memory layout of the new datatype for the MPI to understand
    int block_lengths[number_of_blocks] = {1,1}; // Length of each block.

    /* MPI_Aint -> MPI Address Integer - datatype defined in the MPI (Message Passing Interface) standard to represent addresses and address differences.
    Portability: Different systems may have different sizes for pointers and integers. MPI_Aint ensures that address arithmetic works correctly across all platforms.
    Consistency: Using MPI_Aint guarantees that you can correctly handle memory addresses and displacements in MPI operations, which is critical for defining custom datatypes.*/
    MPI_Aint offsets[number_of_blocks], block1_address, block2_address; // Offsets of each block

    MPI_Datatype primitive_datatypes[number_of_blocks] = {MPI_INT, MPI_DOUBLE}; // Primitive data types of each block

    MPI_Get_address(&((struct Variable*)0)->x, &block1_address); //  saying "imagine there's a Variable struct at memory address 0". 
    //Since we're starting from address 0, this gives us the offset of 1st block within the struct.
    MPI_Get_address(&((struct Variable*)0)->y, &block2_address);
    offsets[0] = 0;
    offsets[1] = block2_address-block1_address;

    MPI_Type_create_struct(number_of_blocks, block_lengths, offsets, primitive_datatypes, &mpi_var_struct_type); //Creating MPI Datatype

    MPI_Type_commit(&mpi_var_struct_type);// Commit the datatype to MPI.
}


int main(int argc, char* argv[]){
    MPI_Init(&argc, &argv);
    int comm_world_size, task_world_rank;

    MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);
    MPI_Comm_rank( MPI_COMM_WORLD, &task_world_rank);

    Variable var;

    if(task_world_rank==MASTER){
        var.x = 9;
        var.y = 3.14325;
    }

    MPI_Datatype mpi_var_struct_type;
    create_variable_struct_datatype(mpi_var_struct_type); // Create Variable Struct Datatype

    MPI_Bcast( &var, 1, mpi_var_struct_type, MASTER, MPI_COMM_WORLD);
    std::cout<<"For rank: "<<task_world_rank<<" the value of integer: "<<var.x<<" the value of the double: "<<var.y<<std::endl;
    MPI_Type_free(&mpi_var_struct_type);
    MPI_Finalize();
    return 0;
}
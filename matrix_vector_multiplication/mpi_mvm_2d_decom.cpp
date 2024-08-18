#include <chrono>
#include <iostream>
#include <vector>
#include <random>
#include <mpi.h>


#define MASTER 0

int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    int comm_size, task_comm_world_rank;

    MPI_Comm_size( MPI_COMM_WORLD , &comm_size);
    MPI_Comm_rank( MPI_COMM_WORLD , &task_comm_world_rank);

    int matrix_size, vector_size;
    int  process_buffer_rank,recv_buffer_rank = task_comm_world_rank; // To get rank of the input vector buffers of which we are processing

    std::vector<double> input_matrix, input_vector, output_vector;

    if(task_comm_world_rank == 0){
        if(argc != 3){
            std::cout << "Need input 2 parameters. Input should be mpiexec -np 8 ./mpiapp matrix_size(#rows) vector_size"<<std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        matrix_size = std::atoi(argv[1]);
        vector_size = std::atoi(argv[2]);

        std::mt19937 random_number_algorithm(11113);
        std::uniform_real_distribution<double> dist(-1,1);
        auto random_number_generator = std::bind(dist, random_number_algorithm);

        input_matrix.resize(matrix_size*vector_size);
        input_vector.resize(vector_size);
        output_vector.resize(matrix_size);

        std::generate(std::begin(input_matrix), std::end(input_matrix), random_number_generator);
        std::generate(std::begin(input_vector), std::end(input_vector), random_number_generator);
    }

    // * Sending the input data to every matrix
    MPI_Bcast( &matrix_size, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast( &vector_size, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}
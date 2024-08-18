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

    int row_chunk, matrix_buffer_size, vector_buffer_size, output_buffer_size;

    // 1D Decomposition of Matrix across rows
    row_chunk = (matrix_size+(comm_size-1))/comm_size;

    matrix_buffer_size = row_chunk*vector_size;
    output_buffer_size = row_chunk;

    vector_buffer_size = (vector_size+(comm_size-1))/comm_size;



    std::vector<double> input_matrix_buffer(matrix_buffer_size), input_vector_buffer(vector_buffer_size), recv_input_vector_buffer(vector_buffer_size),output_vector_buffer(output_buffer_size, 0);


    int matrix_local_size = ((task_comm_world_rank == comm_size-1) ? matrix_size - (row_chunk*(comm_size-1)): row_chunk);
    int vector_local_size = ((task_comm_world_rank == comm_size - 1) ? vector_size - (vector_buffer_size * (comm_size - 1)) : vector_buffer_size);

    int send_counts[comm_size], displacements[comm_size];

    // Scattering Matrix in 1D Decomposition
    for(int i =0; i<comm_size; ++i){
        send_counts[i] = ((i!=comm_size-1) ? row_chunk : matrix_size - (row_chunk*(comm_size-1)))*vector_size;
        displacements[i] = i*row_chunk*vector_size;
    }

    MPI_Scatterv(input_matrix.data(), send_counts, displacements, MPI_DOUBLE, input_matrix_buffer.data(), send_counts[task_comm_world_rank], MPI_DOUBLE, MASTER, MPI_COMM_WORLD);


    // Scattering Vector
    for(int i =0; i<comm_size; ++i){
        send_counts[i] = (i!=comm_size-1) ? vector_buffer_size : vector_size - (vector_buffer_size*(comm_size-1));
        displacements[i] = i*vector_buffer_size;
    }

    MPI_Scatterv(input_vector.data(), send_counts, displacements, MPI_DOUBLE, input_vector_buffer.data(), send_counts[task_comm_world_rank], MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    int process_col_start, process_col_end;

    int src = (task_comm_world_rank+1)%comm_size;
    int  dest= (comm_size+(task_comm_world_rank-1))%comm_size;

    MPI_Request send_recv_request[4];

    for(int communication=0; communication<comm_size; ++communication){
        process_buffer_rank = recv_buffer_rank;
        for(int row = 0; row<matrix_local_size; ++row){
            process_col_start = process_buffer_rank*vector_buffer_size;
            process_col_end = ((process_buffer_rank!=comm_size-1) ? process_col_start+vector_buffer_size : vector_size);
            for(int col=process_col_start, local_col=0; col<process_col_end; ++col, ++local_col){
                output_vector_buffer[row] +=  input_matrix_buffer[row*vector_size+col]*input_vector_buffer[local_col];
            }
        }

        if(communication!=comm_size-1){

            // * communicating the vector buffer
            MPI_Isend(input_vector_buffer.data(), vector_buffer_size , MPI_DOUBLE , dest , 111, MPI_COMM_WORLD, &send_recv_request[0]);
            MPI_Irecv(recv_input_vector_buffer.data(), vector_buffer_size, MPI_DOUBLE, src, 111, MPI_COMM_WORLD, &send_recv_request[1]);            //communicating the vector buffer and its rank

            // * Communicating rank
            MPI_Isend( &process_buffer_rank, 1, MPI_INT, dest , 222, MPI_COMM_WORLD, &send_recv_request[2]);
            MPI_Irecv( &recv_buffer_rank, 1, MPI_INT, src, 222, MPI_COMM_WORLD, &send_recv_request[3]);

            MPI_Waitall(4, send_recv_request, MPI_STATUS_IGNORE);

            std::copy(recv_input_vector_buffer.begin(), recv_input_vector_buffer.end(), input_vector_buffer.begin());
        }
    }

    for(int i =0; i<comm_size; i++){
        send_counts[i] = (i!=comm_size-1) ? output_buffer_size : matrix_size - (output_buffer_size*(comm_size-1));
        displacements[i] = i*output_buffer_size;
    }

    MPI_Gatherv(output_vector_buffer.data(), send_counts[task_comm_world_rank], MPI_DOUBLE, output_vector.data(), send_counts, displacements, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}
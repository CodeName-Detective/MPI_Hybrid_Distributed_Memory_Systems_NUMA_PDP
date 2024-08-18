#include <algorithm>
#include <iostream>
#include <random>
#include <vector>
#include <memory>
#include "mpi.h"


void print_matrix(const float* matrix, int dim){
    for(int i=0; i<dim; ++i){
        for(int j=0; j<dim; ++j){
            std::cout<<matrix[i*dim+j]<<" ";
        }
        std::cout<<std::endl;
    }
}


int main(int argc, char* argv[]){

    MPI_Init(&argc, &argv);


    int num_size;
    MPI_Comm_size( MPI_COMM_WORLD, &num_size);

    const int dim = 4096;
    const int chunk_rows = dim/num_size;

    int task_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);

    std::unique_ptr<float[]> m;
    auto m_chunk = std::make_unique<float[]>(dim*chunk_rows);
    auto pivot_chunk = std::make_unique<float[]>(dim);

    if (task_id==0){

        std::mt19937 mt(123);
        std::uniform_real_distribution dist(1.0f, 2.0f);

        m = std::make_unique<float[]>(dim*dim);
        std::generate_n(m.get(), dim*dim, [&] {return dist(mt);});

    }


    MPI_Scatter(m.get(), dim*chunk_rows, MPI_FLOAT, m_chunk.get(), dim*chunk_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Store requests that for non-blocking sends
    std::vector<MPI_Request> requests(num_size);

    const int start_row_chunk = task_id*chunk_rows;
    const int end_row_chunk = start_row_chunk + chunk_rows;


    // Gaussian Elimination
    for(int row=0; row<dim; ++row){

        int mapped_rank = row/chunk_rows;


        if (task_id==mapped_rank){
            
            int local_row = row - start_row_chunk;

            float pivot=m_chunk[local_row*dim+row];

            for(int col=0; col<dim; ++col){
                m_chunk[local_row*dim+col] = m_chunk[local_row*dim+col]/pivot;
            }

            for(int dest_task = mapped_rank+1; dest_task<num_size; ++dest_task){
                MPI_Isend( m_chunk.get()+local_row*dim, dim, MPI_FLOAT, dest_task, 0, MPI_COMM_WORLD, &requests[dest_task]);
            }


            for(int eliminate_row = row+1; eliminate_row<end_row_chunk; ++eliminate_row){
                float scale = m_chunk[(eliminate_row - start_row_chunk)*dim+row];

                for(int eliminate_col = 0; eliminate_col<dim; ++eliminate_col){
                    m_chunk[(eliminate_row - start_row_chunk)*dim+eliminate_col] = m_chunk[(eliminate_row - start_row_chunk)*dim+eliminate_col] - (m_chunk[(row - start_row_chunk)*dim+eliminate_col]*scale);
                }
            }

            for(int dest_task = mapped_rank+1; dest_task<num_size; ++dest_task){
                MPI_Wait( &requests[dest_task], MPI_STATUS_IGNORE);
            }

        }

        else if (task_id > mapped_rank) {

            MPI_Recv(pivot_chunk.get(), dim, MPI_FLOAT, mapped_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for(int eliminate_row = 0; eliminate_row<chunk_rows; ++eliminate_row){
                float scale = m_chunk[eliminate_row*dim+row];

                for(int eliminate_col = 0; eliminate_col<dim; ++eliminate_col){
                    m_chunk[eliminate_row*dim+eliminate_col] = m_chunk[eliminate_row*dim+eliminate_col] - (pivot_chunk[eliminate_col]*scale);
                }
            }

        }

    }


    MPI_Gather(m_chunk.get(), dim*chunk_rows, MPI_FLOAT, m.get(), dim*chunk_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}
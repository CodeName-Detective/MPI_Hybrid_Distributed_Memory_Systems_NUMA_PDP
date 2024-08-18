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

    const int dim = 1<<3;
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

    // Scattering in Roundrobin Fashion
    for(int i=0; i<chunk_rows; ++i){
        MPI_Scatter( m.get()+i*dim*num_size , dim , MPI_FLOAT,m_chunk.get()+i*dim, dim, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }
 

    for(int row=0; row<dim; ++row){
        int mapped_rank = row % num_size;

        if(task_id==mapped_rank){

            int local_row = row / num_size;

            float pivot = m_chunk[local_row*dim+row];

            for(int col=0; col<dim; ++col){
                m_chunk[local_row*dim+col] = m_chunk[local_row*dim+col]/pivot;
            }

            MPI_Bcast( m_chunk.get()+local_row*dim, dim, MPI_FLOAT, task_id, MPI_COMM_WORLD);

            for(int eliminate_row=local_row+1; eliminate_row<chunk_rows; ++eliminate_row){
                float scale = m_chunk[eliminate_row*dim+row];
                for(int col=0; col<dim; ++col){
                    m_chunk[eliminate_row*dim+col] = m_chunk[eliminate_row*dim+col] - (scale*m_chunk[local_row*dim+col]);
                }
            }
        }

        else{
            MPI_Bcast( pivot_chunk.get(), dim, MPI_FLOAT, mapped_rank, MPI_COMM_WORLD);

            for(int eliminate_row=0; eliminate_row<chunk_rows; ++eliminate_row){
                int eliminate_row_global = num_size*eliminate_row+task_id;
                if(eliminate_row_global>row){
                    float scale = m_chunk[eliminate_row*dim+row];
                    for(int col=0; col<dim; ++col){
                        m_chunk[eliminate_row*dim+col] = m_chunk[eliminate_row*dim+col] - (scale*pivot_chunk[col]);
                    }
                }
            }
        }

    }

    // Gathering in Roundrobin Fashion
    for(int i=0; i<chunk_rows; ++i){
        MPI_Gather(m_chunk.get()+i*dim, dim, MPI_FLOAT, m.get()+i*dim*num_size, dim, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }

    if(task_id==0){
        print_matrix(m.get(), dim);
    }

    MPI_Finalize();
    return 0;
}
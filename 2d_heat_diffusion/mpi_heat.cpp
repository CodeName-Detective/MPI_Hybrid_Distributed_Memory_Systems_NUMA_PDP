#include <iostream>
#include <fstream>
#include <iomanip>
#include "mpi.h"

#define NXPROB 100
#define NYPROB 100
#define MASTER 0


struct Params {
    float x;
    float y;
    int nts;
};


void intidat(const int nx, const int ny, float* u1, const int& start_idx_global, const int& end_idx_global){

    for(int ix=start_idx_global; ix<end_idx_global; ++ix){
        for(int iy=0; iy<ny; ++iy){
            u1[((ix-start_idx_global)*NYPROB)+iy] = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
        }
    }
}


void prtdat(const int nx, const int ny, float* u1, const char* fnam){
    std::ofstream fout;
    fout << std::fixed << std::setprecision(3);
    fout.open(fnam, std::ios::out);

    for(int iy=ny-1; iy>=0; --iy){
        for(int ix=0; ix<nx; ++ix){
            fout<<u1[(ix*NYPROB)+iy];
            if(ix != nx-1){
                fout<<" ";
            }
            else{
                fout<<"\n";
            }
        }
    }
    fout.close();
}


void update(const int nx, const int ny, float* u1, float* u2, const Params& parms, 
            const int& max_buffer_size, const int& current_buffer_size, const int& task_comm_world_rank, const int& comm_size){

    float buffer[NYPROB];

    MPI_Request send_recv_req[2];

    for(int ix=0; ix<max_buffer_size; ++ix){
        
        if (ix==0){
            int dest = (task_comm_world_rank+1)%comm_size;
            int src = (comm_size+(task_comm_world_rank-1))%comm_size;

            if (task_comm_world_rank != comm_size-1)
            {
                MPI_Isend( &u1[(max_buffer_size-1)*NYPROB+0], NYPROB, MPI_FLOAT, dest, 0, MPI_COMM_WORLD, &send_recv_req[0]);
            }
            else{
                MPI_Isend( &u1[(current_buffer_size-1)*NYPROB+0], NYPROB, MPI_FLOAT, dest, 0, MPI_COMM_WORLD, &send_recv_req[0]);
            }

            MPI_Irecv( &buffer, NYPROB, MPI_FLOAT, src, 0, MPI_COMM_WORLD, &send_recv_req[1]);

            MPI_Waitall( 2, send_recv_req, MPI_STATUS_IGNORE);
        }
        else if(ix==max_buffer_size-1){
            int src = (task_comm_world_rank+1)%comm_size;
            int dest = (comm_size+(task_comm_world_rank-1))%comm_size;

            MPI_Isend( &u1[0*NYPROB+0], NYPROB, MPI_FLOAT, dest, 0, MPI_COMM_WORLD, &send_recv_req[0]);
            MPI_Irecv( &buffer, NYPROB, MPI_FLOAT, src, 0, MPI_COMM_WORLD, &send_recv_req[1]);
            MPI_Waitall( 2, send_recv_req, MPI_STATUS_IGNORE);
        }
        
        for(int iy=1; iy<ny-1; ++iy){
            if (ix==0 && task_comm_world_rank==MASTER){
                continue; //Intintially skipping the iteration
            }
            else if(ix>=current_buffer_size-1 && task_comm_world_rank==comm_size-1){
                continue; //Intintially skipping the iteration
            }
            else if(ix==0){
                u2[ix*NYPROB+iy] = u1[ix*NYPROB+iy]+ (parms.x*(u1[(ix+1)*NYPROB+iy]+buffer[iy]-2*u1[ix*NYPROB+iy])) + (parms.y*(u1[(ix*NYPROB)+(iy+1)]+u1[(ix*NYPROB)+(iy-1)]-2*u1[ix*NYPROB+iy]));
            }
            else if(ix==max_buffer_size-1){
                u2[ix*NYPROB+iy] = u1[ix*NYPROB+iy]+ (parms.x*(buffer[iy]+u1[(ix-1)*NYPROB+iy]-2*u1[ix*NYPROB+iy])) + (parms.y*(u1[(ix*NYPROB)+(iy+1)]+u1[(ix*NYPROB)+(iy-1)]-2*u1[ix*NYPROB+iy]));
            }
            else{
                u2[ix*NYPROB+iy] = u1[ix*NYPROB+iy]+ (parms.x*(u1[(ix+1)*NYPROB+iy]+u1[(ix-1)*NYPROB+iy]-2*u1[ix*NYPROB+iy])) + (parms.y*(u1[(ix*NYPROB)+(iy+1)]+u1[(ix*NYPROB)+(iy-1)]-2*u1[ix*NYPROB+iy]));   
            }     
            
        }

        
    }
}

int main(int argc, char* argv[]){
    
    MPI_Init( &argc, &argv);
    Params parms = {0.1, 0.1, 100};
    int comm_size, task_comm_world_rank;

    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_comm_world_rank);

    float* input, * output, *input_buffer, *output_buffer;


    if(task_comm_world_rank == MASTER){
        input = new float [NXPROB*NYPROB];
        output = new float [NXPROB*NYPROB];


        std::cout<<"Starting serial version of 2D heat example..."<<std::endl;
        std::cout<<"Using "<<NXPROB<<"x"<<NYPROB<<" grid"<<std::endl;
        std::cout<<"Initializing grid and creating input file:"<<std::endl;
    }

    int max_buffer_size = NXPROB/comm_size;
    int reminder = NXPROB%comm_size;

    if(reminder>0){++max_buffer_size;}
    
    int start_index_global = max_buffer_size*task_comm_world_rank, end_idx_global;

    if(task_comm_world_rank == comm_size-1){
        end_idx_global= NXPROB;
    }
    else{
        end_idx_global = start_index_global+max_buffer_size;
    }

    int current_buffer_size = end_idx_global - start_index_global;


    input_buffer = new float [max_buffer_size*NYPROB];
    output_buffer = new float [max_buffer_size*NYPROB];

    
    intidat(NXPROB, NYPROB, input_buffer, start_index_global, end_idx_global); // In this context u[0] will decay to &u[0][0] automatically in function call. Just like int arr[5];

    MPI_Request input_request;
    MPI_Igather( input_buffer, max_buffer_size*NYPROB, MPI_FLOAT, task_comm_world_rank == MASTER ? input:nullptr, max_buffer_size*NYPROB, MPI_FLOAT, MASTER, MPI_COMM_WORLD, &input_request);

    std::cout<<"Waiting for data to gather";

    for(int ix=0; ix<current_buffer_size; ++ix){
        output_buffer[(ix*NYPROB)+0] = input_buffer[(ix*NYPROB)+0];
        output_buffer[(ix*NYPROB)+(NYPROB-1)] = input_buffer[(ix*NYPROB)+(NYPROB-1)];
    }

    for(int iy=0; iy<NYPROB; ++iy){
        output_buffer[(0*NYPROB)+iy] = input_buffer[(0*NYPROB)+iy];
        output_buffer[((NXPROB-1)*NYPROB)+iy] = input_buffer[((NXPROB-1)*NYPROB)+iy];
    }

    std::cout<<"Iterating over "<<parms.nts<<" time steps..."<<std::endl;

    for(int it=1; it<=parms.nts; ++it){

        if(it%2==1){
            update(NXPROB, NYPROB, input_buffer, output_buffer, parms, max_buffer_size, current_buffer_size, task_comm_world_rank, comm_size);
        }
        else{
            update(NXPROB, NYPROB, output_buffer, input_buffer, parms, max_buffer_size, current_buffer_size, task_comm_world_rank, comm_size);
        }
    }

    MPI_Wait( &input_request , MPI_STATUS_IGNORE);
    MPI_Gather( output_buffer, max_buffer_size*NYPROB, MPI_FLOAT, task_comm_world_rank == MASTER ? output:nullptr, max_buffer_size*NYPROB, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
    if(task_comm_world_rank == MASTER){
        prtdat(NXPROB, NYPROB, input, "mpi_initial.dat");
        std::cout<<"input file created!";
        prtdat(NXPROB, NYPROB, output, "mpi_final.dat");
        std::cout<<"output file created!";

        delete[] input;
        delete[] output;
    }

    delete[] input_buffer;
    delete[] output_buffer;

    MPI_Finalize();
    return 0;
}
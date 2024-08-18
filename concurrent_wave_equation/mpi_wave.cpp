#include <iostream>
#include <cmath>
#include <cstring>

#include "mpi.h"

#define PI 3.14159265
#define MASTER 0


void init_param(int& total_steps, int& total_points){
    std::cout<<"Enter total number of points vibrating along the string between [20-1000] such that value+2 should be multiple of 8:";
    std::cin>>total_points;
    std::cout<<"Enter total number of time steps [1-1000]:";
    std::cin>>total_steps;
}

void init_line(const int& total_points, double* const values, double* const old_values){
    double x, fac{2.0*PI}, k{0.0}, tmp{(double)total_points-1};

    for(int i=0; i<total_points; ++i,++k){
        x=k/tmp;
        values[i] = std::sin(fac*x);
    }

    std::memcpy(old_values, values, (total_points)*sizeof(double));
}


void do_math(int idx, double* const new_values, const double* const values, const double* const old_values, const double& prev_value, const double& nxt_value){
    double dtime{0.3}, c{1.0}, dx{1.0}, tau, sqtau;
    tau = (c * dtime / dx);
    sqtau = tau * tau;
    new_values[idx] = (2.0*values[idx]-old_values[idx]+ (sqtau*(prev_value + nxt_value - (2.0*values[idx]))));
}


void update(const int& total_steps, const int&  total_points, double* const values, double* const old_values){

    double *values_buffer, *old_values_buffer, *new_values_buffer;
    double prev_value_buffer, nxt_value_buffer;

    int num_tasks, task_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_rank);

    int local_points = total_points/num_tasks;
    int reminder = total_points%num_tasks;

    if(reminder>0){++local_points;}

    values_buffer = new double[local_points];
    old_values_buffer = new double[local_points];
    new_values_buffer = new double[local_points];

    int start_index_global = local_points*num_tasks, end_idx_global;

    if(task_rank == num_tasks-1){
        end_idx_global= total_points-1;
    }
    else{
        end_idx_global = start_index_global+local_points;
    }

    MPI_Scatter(values, local_points, MPI_DOUBLE, values_buffer, local_points, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Scatter(old_values, local_points, MPI_DOUBLE, old_values_buffer, local_points, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    MPI_Request send_recv_request[2]; 

    for(int i=0; i<total_steps; ++i){
        for(int local_idx=0; local_idx<local_points; ++local_idx){
            int global_idx = task_rank*local_points+local_idx;
 

            if(global_idx==0 || global_idx==total_points-1){
                new_values_buffer[local_idx] = 0.0;
            }
            if (local_idx==0){
                int dest = (task_rank+1)%num_tasks;
                int src = (num_tasks+(task_rank-1))%num_tasks;
                MPI_Isend( &values_buffer[local_points-1], 1, MPI_DOUBLE, dest, 111, MPI_COMM_WORLD, &send_recv_request[0]);
                MPI_Irecv( &prev_value_buffer, 1, MPI_DOUBLE, src, 111, MPI_COMM_WORLD, &send_recv_request[1]);
                
                nxt_value_buffer = values_buffer[local_idx+1];
                MPI_Waitall(2, send_recv_request, MPI_STATUS_IGNORE);
            }
            else if(local_idx==local_points-1){
                int dest = (num_tasks+(task_rank-1))%num_tasks ;
                int src = (task_rank+1)%num_tasks;
                MPI_Isend( &values_buffer[0], 1, MPI_DOUBLE, dest, 222, MPI_COMM_WORLD, &send_recv_request[0]);
                MPI_Irecv( &nxt_value_buffer, 1, MPI_DOUBLE, src, 222, MPI_COMM_WORLD, &send_recv_request[1]);

                prev_value_buffer = values_buffer[local_idx-1];
                MPI_Waitall(2, send_recv_request, MPI_STATUS_IGNORE);
            }
            else{
                prev_value_buffer = values_buffer[local_idx-1];
                nxt_value_buffer = values_buffer[local_idx+1];
            }

            if(global_idx==0 || global_idx>=total_points-1){
                new_values_buffer[local_idx] = 0.0;
            }

            else{
                do_math(local_idx, new_values_buffer, values_buffer, old_values_buffer, prev_value_buffer, nxt_value_buffer);
            }

        }
        std::memcpy(old_values_buffer, values_buffer, (local_points)*sizeof(double));
        std::memcpy(values_buffer, new_values_buffer, (local_points)*sizeof(double));
    }

    MPI_Gather(values_buffer, local_points, MPI_DOUBLE, values, local_points, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    delete[] values_buffer;
    delete[] old_values_buffer;
    delete[] new_values_buffer;
}


void print_final( const int& total_points,const double* const values){
    for(int i=0; i<total_points; ++i){
        std::cout<<values[i]<<", ";
        if((i+1)%10 == 0){std::cout<<std::endl;}
    }
}

int main(int argc, char* argv[]){
    int total_steps{1}, total_points{35};
    double *values, *old_values;

    MPI_Init(&argc, &argv);
    int num_tasks, task_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_rank);

    if(task_rank==MASTER){

        //init_param(total_steps, total_points);

        values = new double[total_points];
        old_values = new double[total_points];

        init_line(total_points, values, old_values);
    }

    update(total_points, total_points, values, old_values);


    if(task_rank==MASTER){
        print_final(total_points, values);
        delete[] values;
        delete[] old_values;
    }

    MPI_Finalize();
    return 0;
}
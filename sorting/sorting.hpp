#ifndef A2_HPP
#define A2_HPP

#include <vector>
#include <mpi.h>

int get_start_ptr(const int& input_range_size, std::vector<int>& cumulative_sum, const int& start_index_global){
    for(int i=0; i<input_range_size; ++i){
        if(start_index_global+i < cumulative_sum[i]){
            return i;
        }
    }
    return -1;
}

void isort(std::vector<short int>& Xi, MPI_Comm comm) {
    int comm_size, task_comm_world_rank;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &task_comm_world_rank);

    const int input_range_size = 2*comm_size-1;
    std::vector<int> local_counts(input_range_size, 0), global_counts(input_range_size), cumulative_sum(input_range_size);

    for (short int a : Xi) {
        if (a + (comm_size-1) >= 0 && a + (comm_size-1) < input_range_size) {
            ++local_counts[a + (comm_size-1)];
        }
    }

    MPI_Allreduce(local_counts.data(), global_counts.data(), input_range_size, MPI_INT, MPI_SUM, comm);

    for (int i = 0, temp_sum = 0; i < input_range_size; ++i) {
        temp_sum += global_counts[i];
        cumulative_sum[i] = temp_sum;
    }

    int local_vector_size = Xi.size();
    int start_index_global = 0;

    MPI_Exscan(&local_vector_size, &start_index_global, 1, MPI_INT, MPI_SUM, comm);

    int start_ptr = get_start_ptr(input_range_size, cumulative_sum, start_index_global);

    int local_ptr = 0;

    while (start_ptr < input_range_size && (start_index_global + local_vector_size) > cumulative_sum[start_ptr]) {
        int end_idx = std::min(local_vector_size, cumulative_sum[start_ptr] - start_index_global);
        for (int idx = local_ptr; idx < end_idx; ++idx) {
            Xi[idx] = start_ptr - comm_size + 1;
        }
        local_ptr = end_idx;
        ++start_ptr;
    }

    for (int idx = local_ptr; idx < local_vector_size; ++idx) {
        Xi[idx] = start_ptr - comm_size + 1;
    }
}

#endif // A2_HPP

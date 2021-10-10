#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>
#include "auxiliary.c"
#include "auxiliary.h"


void mergeSort(int *arr, int l, int r, int *mirror);
void BMM(int world_rank, int* chunk_offset, int* A_block_IDs, int* A_locations, int* A_nz_ptr, int* B_block_IDs, int* B_locations, 
            int* B_nz_ptr, int Blocks_per_row, int* A_nz_blocks_ptr, int* B_nz_blocks_ptr, int* C_i, int* C_j, int* C_nz);
void BMMfiltered(int world_rank, int* chunk_offset, int* A_block_IDs, int* A_locations, int* A_nz_ptr, int* B_block_IDs, int* B_locations, int* B_nz_ptr, int Blocks_per_row, 
            int* A_nz_blocks_ptr, int* B_nz_blocks_ptr, unsigned int* F_ptr, unsigned int* F_row, int* C_i, int* C_j, int* C_nz);
void COO_to_CSC(unsigned int* C_ptr, unsigned int* C_row, int* C_i, int* C_j, int C_nz);
void preprocess_A(unsigned int* A_ptr, unsigned int* A_row, int* A_block_IDs, int* A_locations, 
            int* A_nz_ptr, int* A_nz_blocks_per_col, int* A_nz_blocks_ptr, int* nz_count_ptr, int* ID_count_ptr);
void preprocess_B(unsigned int* B_ptr, unsigned int* B_row, int* B_block_IDs, int* B_locations, 
            int* B_nz_ptr, int* B_nz_blocks_per_col, int* B_nz_blocks_ptr, int* nz_count_ptr, int* ID_count_ptr);


int main(){
    srand(42);

	/* ---------------------------------
		Sparse Boolean Matrix Creation
	   --------------------------------- */
	   
	   
    MPI_Init(NULL, NULL);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int processes_count;
    MPI_Comm_size(MPI_COMM_WORLD, &processes_count);

    unsigned int* A_row = (unsigned int *)malloc(d*n*sizeof(unsigned int));
    if (A_row==NULL) exit(-1);
    unsigned int* B_row = (unsigned int *)malloc(d*n*sizeof(unsigned int));
    if (B_row==NULL) exit(-1);
	
	unsigned int F_nz_per_col = 10 * d;
    unsigned int* F_row = (unsigned int *)malloc(F_nz_per_col*n*sizeof(unsigned int));
    if (F_row==NULL) exit(-1);

    if(world_rank==0){
        sparseBoolean(A_row, d);
        sparseBoolean(B_row, d);
        sparseBoolean(F_row, F_nz_per_col);
    }
    MPI_Bcast(A_row, d*n, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(B_row, d*n, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(F_row, F_nz_per_col*n, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    unsigned int* A_ptr = (unsigned int *)malloc((n+1)*sizeof(unsigned int));
    if (A_ptr==NULL) exit(-1);
    A_ptr[0] = 1;
    for (int i=1; i<n+1; i++)
        A_ptr[i] = A_ptr[i-1] + d;

    unsigned int* B_ptr = (unsigned int *)malloc((n+1)*sizeof(unsigned int));
    if (B_ptr==NULL) exit(-1);
    B_ptr[0] = 1;
    for (int i=1; i<n+1; i++)
        B_ptr[i] = B_ptr[i-1] + d;

    unsigned int* F_ptr = (unsigned int *)malloc((n+1)*sizeof(unsigned int));
    if (F_ptr==NULL) exit(-1);
    F_ptr[0] = 1;
    for (int i=1; i<n+1; i++)
        F_ptr[i] = F_ptr[i-1] + F_nz_per_col;

    int total_blocks = (n/b)*(n/b);
	
	/* ---------------------------------
		      Preprocessing step
	   --------------------------------- */
	   
	   
    int* A_block_IDs = (int *)malloc(total_blocks*sizeof(int));
    if (A_block_IDs==NULL) exit(-1);

    int length_A_row = d*n;
    int* A_locations = (int *)malloc(3*length_A_row*sizeof(int));
    if (A_locations==NULL) exit(-1);

    int* A_nz_ptr = (int *)malloc(total_blocks*sizeof(int));
    if (A_nz_ptr==NULL) exit(-1);

    int* A_nz_blocks_per_row = (int *)malloc((n/b)*sizeof(int));
    if (A_nz_blocks_per_row==NULL) exit(-1);

    int* A_nz_blocks_ptr = (int *)malloc((n/b)*sizeof(int));
    if (A_nz_blocks_ptr==NULL) exit(-1);

    int* B_block_IDs = (int *)malloc(total_blocks*sizeof(int));
    if (B_block_IDs==NULL) exit(-1);

    int length_B_row = d*n;
    int* B_locations = (int *)malloc(3*length_B_row*sizeof(int));
    if (B_locations==NULL) exit(-1);

    int* B_nz_ptr = (int *)malloc(total_blocks*sizeof(int));
    if (B_nz_ptr==NULL) exit(-1);

    int* B_nz_blocks_per_col = (int *)malloc((n/b)*sizeof(int));
    if (B_nz_blocks_per_col==NULL) exit(-1);

    int* B_nz_blocks_ptr = (int *)malloc((n/b)*sizeof(int));
    if (B_nz_blocks_ptr==NULL) exit(-1);

    MPI_Barrier(MPI_COMM_WORLD);
    struct timespec start;
	struct timespec total_start;   
	
    if(world_rank==0){
		clock_gettime(CLOCK_MONOTONIC, &total_start);
        clock_gettime(CLOCK_MONOTONIC, &start);
    }

    int nz_count_ptr_A, ID_count_ptr_A;
    int nz_count_ptr_B, ID_count_ptr_B;
	
    if(world_rank==0)
        preprocess_A(A_ptr, A_row, A_block_IDs, A_locations, A_nz_ptr, A_nz_blocks_per_row, A_nz_blocks_ptr, &nz_count_ptr_A, &ID_count_ptr_A);
    if(world_rank==1 || processes_count == 1)
        preprocess_B(B_ptr, B_row, B_block_IDs, B_locations, B_nz_ptr, B_nz_blocks_per_col, B_nz_blocks_ptr, &nz_count_ptr_B, &ID_count_ptr_B);

    MPI_Bcast(&nz_count_ptr_A, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ID_count_ptr_A, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(A_block_IDs, ID_count_ptr_A, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(A_locations, 3*nz_count_ptr_A, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(A_nz_blocks_ptr, n/b, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(A_nz_ptr, ID_count_ptr_A, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    MPI_Bcast(&nz_count_ptr_B, 1, MPI_INT, 1, MPI_COMM_WORLD);
    MPI_Bcast(&ID_count_ptr_B, 1, MPI_INT, 1, MPI_COMM_WORLD);
    MPI_Bcast(B_block_IDs, ID_count_ptr_B, MPI_UNSIGNED, 1, MPI_COMM_WORLD);
    MPI_Bcast(B_locations, 3*nz_count_ptr_B, MPI_UNSIGNED, 1, MPI_COMM_WORLD);
    MPI_Bcast(B_nz_blocks_ptr, n/b, MPI_UNSIGNED, 1, MPI_COMM_WORLD);
    MPI_Bcast(B_nz_ptr, ID_count_ptr_B, MPI_UNSIGNED, 1, MPI_COMM_WORLD);

    
    struct timespec end;
    double time;
    if(world_rank==0){
        clock_gettime(CLOCK_MONOTONIC, &end);
        time = time_spent(start, end);
        printf("Preprocessing: %lf\n", time);
    }
	MPI_Barrier(MPI_COMM_WORLD);
    free(B_ptr); free(B_row);
    free(A_ptr); free(A_row);
	
	/* ---------------------------------
		         BMM algorithm
	   --------------------------------- */
	   
	   
    int* C_i;
    int* C_j;

    int C_nz = 0;
    int Blocks_per_row = n / b;

    if(world_rank==0){
        clock_gettime(CLOCK_MONOTONIC, &start);

        C_i = (int *)malloc(d*d*n*sizeof(int));
        if(C_i==NULL) exit(-1);

        C_j = (int *)malloc(d*d*n*sizeof(int));
        if(C_j==NULL) exit(-1);
    }

    int *chunk_offset = (int *)calloc(processes_count+1, sizeof(int));
    if(chunk_offset==NULL) exit(-1);

    if(world_rank==0){
        int chunk = total_blocks / processes_count;
        int remaining_blocks = total_blocks - chunk * processes_count;

        int *chunk_size = (int *)malloc(processes_count*sizeof(int));
        if(chunk_size==NULL) exit(-1);

        for(int process_iter=0; process_iter<processes_count; process_iter++){
            chunk_size[process_iter] = chunk;
        }

        while(remaining_blocks != 0){
            for(int process_iter=0; process_iter<processes_count; process_iter++){
                chunk_size[process_iter] += 1;
                remaining_blocks --;
                if(remaining_blocks == 0){
                    break;
                }
            }
        }
        for(int process_iter=1; process_iter<processes_count+1; process_iter++){
            chunk_offset[process_iter] = chunk_offset[process_iter-1] + chunk_size[process_iter-1];
        }
        free(chunk_size);
    }
    MPI_Bcast(chunk_offset, processes_count+1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
	
	if(FILTERED == 0){
		BMM(world_rank, chunk_offset, A_block_IDs, A_locations, A_nz_ptr, B_block_IDs, B_locations, B_nz_ptr, 
            Blocks_per_row, A_nz_blocks_ptr, B_nz_blocks_ptr, C_i, C_j, &C_nz);
    }
	else{
		BMMfiltered(world_rank, chunk_offset, A_block_IDs, A_locations, A_nz_ptr, B_block_IDs, B_locations, B_nz_ptr, 
            Blocks_per_row, A_nz_blocks_ptr, B_nz_blocks_ptr, F_ptr, F_row, C_i, C_j, &C_nz);
	}
    MPI_Barrier(MPI_COMM_WORLD);
    if(world_rank==0){
        clock_gettime(CLOCK_MONOTONIC, &end);
        time = time_spent(start, end)-1;
        printf("BMM: %lf\n", time);
    }

    free(A_block_IDs);
    free(B_block_IDs);
    free(A_locations);
    free(B_locations);
    free(A_nz_ptr);
    free(B_nz_ptr);
    free(A_nz_blocks_per_row);
    free(B_nz_blocks_per_col);
    free(A_nz_blocks_ptr);
    free(B_nz_blocks_ptr);
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	/* ---------------------------------
		      COO to CSC conversion
	   --------------------------------- */
	   
	   
    if(world_rank==0){
		struct timespec start;
		
        unsigned int* C_row = (unsigned int*)malloc(C_nz*sizeof(unsigned int));
        if(C_row==NULL) exit(-1);

        unsigned int* C_ptr = (unsigned int*)malloc((n+1)*sizeof(unsigned int));
        if(C_ptr==NULL) exit(-1);
    
        clock_gettime(CLOCK_MONOTONIC, &start);

        COO_to_CSC(C_ptr, C_row, C_i, C_j, C_nz);
        
		struct timespec end;
        clock_gettime(CLOCK_MONOTONIC, &end);
        double time = time_spent(start, end);
        printf("Conversion to CSC: %lf\n", time);
		struct timespec total_end;
		clock_gettime(CLOCK_MONOTONIC, &total_end);
		time = time_spent(total_start, total_end)-1;
		
		
		printf("---\nTotal time: %lf seconds\n---\n", time);
        free(C_row);
        free(C_ptr);
        free(C_i);
        free(C_j);
    }

    free(F_row);
    free(F_ptr);
    MPI_Finalize();
    return 0;
}


void preprocess_A(unsigned int* A_ptr, unsigned int* A_row, int* A_block_IDs, int* A_locations, 
    int* A_nz_ptr, int* A_nz_blocks_per_row, int* A_nz_blocks_ptr, int *nz_count_ptr, int *ID_count_ptr){


    A_nz_ptr[0] = 0;
    
    int nz_count = 0;
    int ID_count = 0;
    int nz_blocks_per_row = 0;
    int current_row = 0;

    for(int p=0; p<(n/b); p++){ // iterate rows of blocks
        for(int q=0; q<(n/b); q++){ //iterate blocks of each row
            bool block_is_nz = false;
            for(int j=q*b; j<(q*b)+b; j++){ //iterate each column inside block
                for(int k=(A_ptr[j]-1); k<(A_ptr[j+1]-1); k++){ //iterate col_non_zeros
                    int current_row = A_row[k];
                    if(current_row >= p*b && current_row < (p+1)*b){
                        A_locations[3*nz_count] = current_row;
                        A_locations[3*nz_count + 1] = j;
                        A_locations[3*nz_count + 2] = j % b;
                        nz_count++;
                        if(block_is_nz == false){
                            A_block_IDs[ID_count] = p*(n/b) + q;
                            ID_count++;
                            nz_blocks_per_row++;
                            block_is_nz = true;
                        }
                    }
                    else if (A_row[k] >= (p+1)*b)
                        break;
                }
            }
            A_nz_ptr[ID_count] = nz_count;
        }
        A_nz_blocks_per_row[p] = nz_blocks_per_row;
        nz_blocks_per_row = 0;
    }
    A_nz_blocks_ptr[0] = 0;
    for(int i=1; i<(n/b); i++){
        A_nz_blocks_ptr[i] = A_nz_blocks_ptr[i-1] + A_nz_blocks_per_row[i-1];
    }

    
    A_block_IDs = (int *)realloc(A_block_IDs, ID_count*sizeof(int));
    A_nz_ptr = (int *)realloc(A_nz_ptr, ID_count*sizeof(int));
    A_locations = (int *)realloc(A_locations, 3*nz_count*sizeof(int));
    *ID_count_ptr = ID_count;
    *nz_count_ptr = nz_count;
}


void preprocess_B(unsigned int* B_ptr, unsigned int* B_row, int* B_block_IDs, int* B_locations, 
    int* B_nz_ptr, int* B_nz_blocks_per_col, int* B_nz_blocks_ptr, int *nz_count_ptr, int *ID_count_ptr){


    B_nz_ptr[0] = 0;

    int nz_count = 0;
    int ID_count = 0;
    int nz_blocks_per_col = 0;
    int current_row = 0;

    for(int q=0; q<(n/b); q++){ // iterate columns of blocks
        for(int p=0; p<(n/b); p++){ //iterate blocks of each column
            bool block_is_nz = false;
            for(int j=q*b; j<(q*b)+b; j++){ //iterate each column inside block
                for(int k=(B_ptr[j]-1); k<(B_ptr[j+1]-1); k++){ //iterate col_non_zeros
                    current_row = B_row[k];
                    if(current_row >= p*b && current_row < (p+1)*b){
                        B_locations[3*nz_count] = current_row;
                        B_locations[3*nz_count + 1] = j;
                        B_locations[3*nz_count + 2] = current_row % b; 
                        nz_count++;
                        if(block_is_nz == false){
                            B_block_IDs[ID_count] = q*(n/b) + p;
                            ID_count++;
                            nz_blocks_per_col++;
                            block_is_nz = true;
                        }
                    }
                    else if (B_row[k] >= (p+1)*b)
                        break;
                }
            }
            B_nz_ptr[ID_count] = nz_count;
        }
        B_nz_blocks_per_col[q] = nz_blocks_per_col;
        nz_blocks_per_col = 0;
    }
    B_nz_blocks_ptr[0] = 0;
    for(int i=1; i<(n/b); i++){
        B_nz_blocks_ptr[i] = B_nz_blocks_ptr[i-1] + B_nz_blocks_per_col[i-1];
    }

    
    B_block_IDs = (int *)realloc(B_block_IDs, ID_count*sizeof(int));
    B_nz_ptr = (int *)realloc(B_nz_ptr, ID_count*sizeof(int));
    B_locations = (int *)realloc(B_locations, 3*nz_count*sizeof(int));
    *ID_count_ptr = ID_count;
    *nz_count_ptr = nz_count;
}


void COO_to_CSC(unsigned int* C_ptr, unsigned int* C_row, int* C_i, int* C_j, int C_nz){
    
	
	mergeSort(C_j, 0, C_nz-1, C_i);
	
    C_ptr[0] = 1;
    int ptr = 0;
	
	// for each column find how many non zeros exist
    for(int i=0; i<n; i++){
        C_ptr[i+1] = C_ptr[i];
        for(int j=ptr; j<C_nz; j++){
            if (C_j[j] != i)
                break;
            else{
                C_row[C_ptr[i+1] - 1] = C_i[j];
                C_ptr[i+1]++;
                ptr++;
            }
        }
    }
}


void BMMfiltered(int world_rank, int* chunk_offset, int* A_block_IDs, int* A_locations, int* A_nz_ptr, int* B_block_IDs, int* B_locations, int* B_nz_ptr, 
            int Blocks_per_row, int* A_nz_blocks_ptr, int* B_nz_blocks_ptr, unsigned int* F_ptr, unsigned int* F_row, int* C_i, int* C_j, int* C_nz){
    
	/* 
	Iterate every block of A and multiply it with every block of B. Depending on
	the "matched" bits coordinates activate the corresponding bit of output matrix C.
	The multiplication is done by iterating every non-zero of block_A and checking
	every non-zero of the corresponding block_B. Finally, the coordinates are checked.
	*/
	
    int total_blocks = (n/b)*(n/b);
    int blocks_per_row = n / b;
    int nz_counter = 0;
    
    unsigned int* C_i_process = (unsigned int *)malloc(b*b*(chunk_offset[world_rank+1]-chunk_offset[world_rank])*sizeof(unsigned int));
    if(C_i_process==NULL) exit(-1);
    unsigned int* C_j_process = (unsigned int *)malloc(b*b*(chunk_offset[world_rank+1]-chunk_offset[world_rank])*sizeof(unsigned int));
    if(C_j_process==NULL) exit(-1);

    #pragma omp parallel for
    for (int block = chunk_offset[world_rank]; block<chunk_offset[world_rank+1]; block++){
        int blocks_row = block / blocks_per_row;
        int blocks_col = block % blocks_per_row;
        int block_nz = 0;

        unsigned int* C_i_block = (unsigned int *)malloc(b*b*sizeof(unsigned int));
        if(C_i_block==NULL) exit(-1);
        unsigned int* C_j_block = (unsigned int *)malloc(b*b*sizeof(unsigned int));
        if(C_j_block==NULL) exit(-1);

        // for each NZ block of A in "blocks_row" row
        for(int p = A_nz_blocks_ptr[blocks_row]; p<A_nz_blocks_ptr[blocks_row+1]; p++){

            int current_A_block = A_block_IDs[p];

            // for each NZ block of b in "blocks_col" column
            for(int q = B_nz_blocks_ptr[blocks_col]; q<B_nz_blocks_ptr[blocks_col+1]; q++){

                int current_B_block = B_block_IDs[q];

                if ((current_A_block%blocks_per_row) != (current_B_block%blocks_per_row)){
                    if ((current_A_block%blocks_per_row) < (current_B_block%blocks_per_row))
                        break;
                    continue;
                }
   
                // for each NZ bit in NZ block of A
                for(int nz_A = A_nz_ptr[p]; nz_A<A_nz_ptr[p+1]; nz_A++){
                    
                    int offset_A = A_locations[3*nz_A + 2];

                    // for each NZ bit in NZ block of B
                    for(int nz_B = B_nz_ptr[q]; nz_B<B_nz_ptr[q+1]; nz_B++){
                        int offset_B = B_locations[3*nz_B + 2];
                        
                        if(offset_A == offset_B){
                            // check if bit exists in F
                            for(int k = F_ptr[B_locations[3*nz_B + 1]]-1; k<F_ptr[B_locations[3*nz_B + 1]+1]-1; k++){
                                if(A_locations[3*nz_A] < F_row[k])
                                    break;
                                else if(A_locations[3*nz_A] == F_row[k]){
                                    C_i_block[block_nz] = A_locations[3*nz_A]; // == i_A
                                    C_j_block[block_nz] = B_locations[3*nz_B + 1]; // == j_B
                                    block_nz+=1;
                                }
                            }
                        }
                    }
                }
            }
        }
        #pragma omp critical
        {
            for(int k = nz_counter; k<(nz_counter+block_nz); k++){
                C_i_process[k] = C_i_block[k - nz_counter];
                C_j_process[k] = C_j_block[k - nz_counter];
            }
            nz_counter += block_nz;
        }
    }
    C_i_process = (unsigned int *)realloc(C_i_process, nz_counter*sizeof(unsigned int));
    C_j_process = (unsigned int *)realloc(C_j_process, nz_counter*sizeof(unsigned int));
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Request request;
    if(world_rank != 0){
        // tags: 0 -> C_i, 1 -> C_j, 2 -> nz_counter
        MPI_Isend(C_i_process, nz_counter, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &request);
        MPI_Isend(C_j_process, nz_counter, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD, &request);
        MPI_Isend(&nz_counter, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &request);
    }
    else {
        // add process' own NZs
        for(int i = *C_nz; i<(*C_nz+nz_counter); i++){
            C_i[i] = C_i_process[i - *C_nz];
            C_j[i] = C_j_process[i - *C_nz];
        }
        *C_nz += nz_counter;
        int processes_count;
        MPI_Comm_size(MPI_COMM_WORLD, &processes_count);

        // receive other processes' NZs
        unsigned int *buffer_i;
        unsigned int *buffer_j;
        for(int i = 1; i<processes_count; i++){

            MPI_Recv(&nz_counter, 1, MPI_INT, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            buffer_i = (unsigned int *)malloc(nz_counter*sizeof(unsigned int));
            MPI_Recv(buffer_i, nz_counter, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            buffer_j = (unsigned int *)malloc(nz_counter*sizeof(unsigned int));
            MPI_Recv(buffer_j, nz_counter, MPI_UNSIGNED, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for(int k = *C_nz; k<(*C_nz+nz_counter); k++){
                C_i[k] = buffer_i[k - *C_nz];
                C_j[k] = buffer_j[k - *C_nz];
            }
            *C_nz += nz_counter;
        }
        C_i = (int *)realloc(C_i, *C_nz*sizeof(int));
        C_j = (int *)realloc(C_j, *C_nz*sizeof(int));
        printf("C has %d non-zeros.\n", *C_nz);
    }
}


void BMM(int world_rank, int* chunk_offset, int* A_block_IDs, int* A_locations, int* A_nz_ptr, int* B_block_IDs, int* B_locations, 
            int* B_nz_ptr, int Blocks_per_row, int* A_nz_blocks_ptr, int* B_nz_blocks_ptr, int* C_i, int* C_j, int* C_nz){
    
	/* 
	Iterate every block of A and multiply it with every block of B. Depending on
	the "matched" bits coordinates activate the corresponding bit of output matrix C.
	The multiplication is done by iterating every non-zero of block_A and checking
	every non-zero of the corresponding block_B. Finally, the coordinates are checked.
	*/
	
    int total_blocks = (n/b)*(n/b);
    int blocks_per_row = n / b;
    int nz_counter = 0;

    unsigned int* C_i_process = (unsigned int *)malloc(b*b*(chunk_offset[world_rank+1]-chunk_offset[world_rank])*sizeof(unsigned int));
    if(C_i_process==NULL) exit(-1);
    unsigned int* C_j_process = (unsigned int *)malloc(b*b*(chunk_offset[world_rank+1]-chunk_offset[world_rank])*sizeof(unsigned int));
    if(C_j_process==NULL) exit(-1);

    #pragma omp parallel for
    for (int block = chunk_offset[world_rank]; block<chunk_offset[world_rank+1]; block++){
        int blocks_row = block / blocks_per_row;
        int blocks_col = block % blocks_per_row;
        int block_nz = 0;

        unsigned int* C_i_block = (unsigned int *)malloc(b*b*sizeof(unsigned int));
        if(C_i_block==NULL) exit(-1);
        unsigned int* C_j_block = (unsigned int *)malloc(b*b*sizeof(unsigned int));
        if(C_j_block==NULL) exit(-1);

        // for each NZ block of A in "blocks_row" row
        for(int i = A_nz_blocks_ptr[blocks_row]; i<A_nz_blocks_ptr[blocks_row+1]; i++){

            int current_A_block = A_block_IDs[i];

            // for each NZ block of B in "blocks_col" column
            for(int j = B_nz_blocks_ptr[blocks_col]; j<B_nz_blocks_ptr[blocks_col+1]; j++){

                int current_B_block = B_block_IDs[j];

                if ((current_A_block%blocks_per_row) != (current_B_block%blocks_per_row)){
                    if ((current_A_block%blocks_per_row) < (current_B_block%blocks_per_row))
                        break;
                    continue;
                }
                    
                // for each NZ bit in NZ block of A
                for(int nz_A = A_nz_ptr[i]; nz_A<A_nz_ptr[i+1]; nz_A++){
                    
                    int offset_A = A_locations[3*nz_A + 2];

                    // for each NZ bit in NZ block of B
                    for(int nz_B = B_nz_ptr[j]; nz_B<B_nz_ptr[j+1]; nz_B++){
                        
                        int offset_B = B_locations[3*nz_B + 2];
                        
                        if(offset_A == offset_B){
                            C_i_block[block_nz] = A_locations[3*nz_A]; // == i_A
                            C_j_block[block_nz] = B_locations[3*nz_B + 1]; // == j_B
                            block_nz+=1;
                        }
                    }
                }
            }
        }
        #pragma omp critical
        {
            for(int i = nz_counter; i<(nz_counter + block_nz); i++){
                C_i_process[i] = C_i_block[i - nz_counter];
                C_j_process[i] = C_j_block[i - nz_counter];
            }
            nz_counter += block_nz;
        }
    }
    C_i_process = (unsigned int *)realloc(C_i_process, nz_counter*sizeof(unsigned int));
    C_j_process = (unsigned int *)realloc(C_j_process, nz_counter*sizeof(unsigned int));
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Request request;
    if(world_rank != 0){
        // tags: 0 -> C_i, 1 -> C_j, 2 -> nz_counter
        MPI_Isend(C_i_process, nz_counter, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &request);
        MPI_Isend(C_j_process, nz_counter, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD, &request);
        MPI_Isend(&nz_counter, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &request);
    }
    else {
        // add process' own NZs
        for(int i = *C_nz; i<(*C_nz+nz_counter); i++){
            C_i[i] = C_i_process[i - *C_nz];
            C_j[i] = C_j_process[i - *C_nz];
        }
        *C_nz += nz_counter;
        int processes_count;
        MPI_Comm_size(MPI_COMM_WORLD, &processes_count);

        // receive other processes' NZs
        unsigned int *buffer_i;
        unsigned int *buffer_j;
        for(int i = 1; i<processes_count; i++){

            MPI_Recv(&nz_counter, 1, MPI_INT, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            buffer_i = (unsigned int *)malloc(nz_counter*sizeof(unsigned int));
            MPI_Recv(buffer_i, nz_counter, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            buffer_j = (unsigned int *)malloc(nz_counter*sizeof(unsigned int));
            MPI_Recv(buffer_j, nz_counter, MPI_UNSIGNED, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for(int k = *C_nz; k<(*C_nz+nz_counter); k++){
                C_i[k] = buffer_i[k - *C_nz];
                C_j[k] = buffer_j[k - *C_nz];
            }
            *C_nz += nz_counter;
        }
        C_i = (int *)realloc(C_i, *C_nz*sizeof(int));
        C_j = (int *)realloc(C_j, *C_nz*sizeof(int));
        printf("C has %d non-zeros.\n", *C_nz);
    }
}



void mergeSort(int *arr, int l, int r, int *mirror){
    if (l < r) {
        int m = l + (r - l) / 2;
        
        #pragma omp task
        mergeSort(arr, l, m, mirror);
        
        mergeSort(arr, m + 1, r,mirror);

        #pragma omp taskwait
        merge(arr, l, m, r,mirror);
    }
}

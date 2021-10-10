#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <pthread.h>
#include "auxiliary.h" 

#define n 500000
#define d 6 // non-zero elements per col
#define b 5000 // block size = bxb
#define FILTERED 1

void sparseBoolean(unsigned int* row, int nz_per_col){
    unsigned int* rand_rows; 
    for (int i=0; i<n; i++){
        rand_rows = (unsigned int *)calloc(nz_per_col, sizeof(unsigned int));
        if (rand_rows==NULL) exit(-1);

        for (int j=0; j<nz_per_col; j++){
            do{
                rand_rows[j] = (unsigned int)(rand() % n);
            } while(duplicateExists(rand_rows, nz_per_col, j));
        }
        // sort before adding, for easier search
        bubblesort(rand_rows);
        for (int j=0; j<nz_per_col; j++){
            row[i*nz_per_col + j] = rand_rows[j];
        }
    }
    free(rand_rows);
}


void bubblesort(unsigned int* array){
    int i, j;
    for (i = 0; i < d-1; i++){
        // Last i elements are already in place
        for (j = 0; j < d-i-1; j++){
            if (array[j] > array[j+1])
                swap(&array[j], &array[j+1]);
        }
    }
}


void swap(unsigned int *xp, unsigned int *yp){
    unsigned int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

// Searches for duplicate element in an array
bool duplicateExists(unsigned int* array, int length, int current){  
    for(int i = 0; i < current; i++) {    
        if(array[i] == array[current])    
            return true;      
    }
    return false; 
}


void merge(int *arr, int l, int m, int r,int *mirror){
    int k;
    int n1 = m - l + 1;
    int n2 = r - m;

    int *L,*R,*L_mirror,*R_mirror;
    L=(int *) malloc(n1 * sizeof(int));
    R =(int *) malloc(n2 * sizeof(int));
    L_mirror=(int *) malloc(n1 * sizeof(int));
    R_mirror =(int *) malloc(n2 * sizeof(int));


    for (int i = 0; i < n1; i++){
        L[i] = arr[l + i];
        L_mirror[i] = mirror[l+i];
    }

    for (int j = 0; j < n2; j++){
        R[j] = arr[m + 1 + j];
        R_mirror[j] = mirror[m+1+j];
    }

    int i = 0;
    int j = 0;
    k = l;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            mirror[k] = L_mirror[i];
            i++;
        }
        else {
            arr[k] = R[j];
            mirror[k] = R_mirror[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        arr[k] = L[i];
        mirror[k] = L_mirror[i];
        i++;
        k++;
    }
    while (j < n2) {
        arr[k] = R[j];
        mirror[k] = R_mirror[j];
        j++;
        k++;
    }
}



// fuction to calculate elapsed time between two timespecs (same as the one in the PDS deliverables)
double time_spent(struct timespec start,struct timespec end){
        struct timespec temp;
        if ((end.tv_nsec - start.tv_nsec) < 0){
            temp.tv_sec = end.tv_sec - start.tv_sec - 1;
            temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
        }
        else {
            temp.tv_sec = end.tv_sec - start.tv_sec;
            temp.tv_nsec = end.tv_nsec - start.tv_nsec;
        }
        return (double)temp.tv_sec +(double)((double)temp.tv_nsec/(double)1000000000);
}
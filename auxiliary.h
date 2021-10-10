#ifndef AUXILIARY_H_INCLUDED
#define AUXILIARY_H_INCLUDED


bool duplicateExists(unsigned int* array, int length, int current);
void sparseBoolean(unsigned int* row, int nz_per_col);
void swap(unsigned int *xp, unsigned int *yp);
void bubblesort(unsigned int* array);
double elapsed_time(struct timespec start, struct timespec end);

#endif // AUXILIARY_H_INCLUDED
CC=gcc
CFLAGS=-O3
MPICC=mpicc

default: all

BMM-serial:
	$(CC) $(CFLAGS) -o BMM-serial BMM-serial.c -lm 
	
BMM-openmp:
	$(CC) $(CFLAGS) -o BMM-openmp BMM-openmp.c -lm -fopenmp	

BMM-MPI-openmp:
	$(MPICC) $(CFLAGS) -o BMM-MPI-openmp BMM-MPI-openmp.c -lm -fopenmp
	

.PHONY: clean

all: BMM-serial BMM-openmp BMM-MPI-openmp

test:
	@printf "\n** Testing serial ΒΜΜ **\n\n"
	./BMM-serial
	@printf "\n** Testing ΒΜΜ with openmp **\n\n"
	./BMM-openmp
	@printf "\n** Testing ΒΜΜ with openmp and MPI **\n\n"
	mpiexec -np 8 ./BMM-MPI-openmp


clean:
	rm -f BMM-serial BMM-openmp BMM-MPI-openmp

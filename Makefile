##############################################################################
#	BY MINGMING LI, ARIZONA STATE UNIVERSITY, SINCE APRIL 2021 
##############################################################################
LIB= -lm
CC=mpicc
OBJFLAG=-c
OPTIM=-O2 -Wall -Wno-return-type

CFILES= Main.c\
	Parsing.c\
	Input.c\
	Mesh.c\
	Tree.c\
	Output.c\
	Construct_ien_id_lm.c\
	Boundary_condition.c\
	Array_work.c\
	Parallel.c\
	Shape_function.c\
	Initial_condition.c\
	Buoyancy.c\
	Viscosity.c\
	Build_matrix.c\
	Solver.c\
	Connection.c\
	Update_temperature.c\
	Matrix_transformation.c\
	Grain_size_evolution.c

HEADER = global_variables.h\
	 function_def.h
	
OBJFILES=$(CFILES:.c=.o)

default: MCC.mpi

MCC.mpi: $(OBJFILES) $(HEADER) Makefile
	$(CC) $(OPTIM) -Wall -g -Wcast-qual -Winline -Wpointer-arith -o MCC.mpi $(OBJFILES) $(LIB)

clean:
	rm -f *.o

clean-all:
	rm -f *.o MCC.mpi

Main.o: $(HEADER) Main.c
	$(CC) $(OPTIM) $(OBJFLAG)  Main.c
Parsing.o: $(HEADER) Parsing.c
	$(CC) $(OPTIM) $(OBJFLAG)  Parsing.c
Input.o: $(HEADER) Input.c
	$(CC) $(OPTIM) $(OBJFLAG)  Input.c
Mesh.o: $(HEADER) Mesh.c
	$(CC) $(OPTIM) $(OBJFLAG)  Mesh.c
Tree.o: $(HEADER) Tree.c
	$(CC) $(OPTIM) $(OBJFLAG)  Tree.c
Output.o: $(HEADER) Output.c
	$(CC) $(OPTIM) $(OBJFLAG)  Output.c
Construct_ien_id_lm.o: $(HEADER) Construct_ien_id_lm.c
	$(CC) $(OPTIM) $(OBJFLAG)  Construct_ien_id_lm.c
Boundary_condition.o: $(HEADER) Boundary_condition.c
	$(CC) $(OPTIM) $(OBJFLAG)  Boundary_condition.c
Array_work.o: $(HEADER) Array_work.c
	$(CC) $(OPTIM) $(OBJFLAG)  Array_work.c
Parallel.o: $(HEADER) Parallel.c
	$(CC) $(OPTIM) $(OBJFLAG)  Parallel.c
Shape_function.o: $(HEADER) Shape_function.c
	$(CC) $(OPTIM) $(OBJFLAG)  Shape_function.c
Initial_condition.o: $(HEADER) Initial_condition.c
	$(CC) $(OPTIM) $(OBJFLAG)  Initial_condition.c
Buoyancy.o: $(HEADER) Buoyancy.c
	$(CC) $(OPTIM) $(OBJFLAG)  Buoyancy.c
Viscosity.o: $(HEADER) Viscosity.c
	$(CC) $(OPTIM) $(OBJFLAG)  Viscosity.c
Build_matrix.o: $(HEADER) Build_matrix.c
	$(CC) $(OPTIM) $(OBJFLAG)  Build_matrix.c
Solver.o: $(HEADER) Solver.c
	$(CC) $(OPTIM) $(OBJFLAG)  Solver.c
Connection.o: $(HEADER) Connection.c
	$(CC) $(OPTIM) $(OBJFLAG)  Connection.c
Update_temperature.o: $(HEADER) Update_temperature.c
	$(CC) $(OPTIM) $(OBJFLAG)  Update_temperature.c
Matrix_transformation.o: $(HEADER) Matrix_transformation.c
	$(CC) $(OPTIM) $(OBJFLAG)  Matrix_transformation.c
Grain_size_evolution.o: $(HEADER) Grain_size_evolution.c
	$(CC) $(OPTIM) $(OBJFLAG)  Grain_size_evolution.c

FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: cluster
clean:
	rm -rf *.o cluster

cluster: cluster.o group.o  set.o  spmat.o vector_operations.o linked_list_s.o division.o mod_matrix.o
	gcc cluster.o group.o  set.o spmat.o vector_operations.o linked_list_s.o division.o mod_matrix.o -o cluster $(LIBS)

cluster.o: cluster.c division.h
	gcc $(FLAGS) -c cluster.c

group.o: group.c group.h spmat.h
	gcc $(FLAGS) -c group.c
	
set.o: set.c set.h linked_list_s.h
	gcc $(FLAGS) -c set.c

spmat.o: spmat.c spmat.h 
	gcc $(FLAGS) -c spmat.c

vector_operations.o: vector_operations.h vector_operations.c
	gcc $(FLAGS) -c vector_operations.c

linked_list_s.o: linked_list_s.h linked_list_s.c group.h
	gcc $(FLAGS) -c linked_list_s.c

division.o: division.h division.c  mod_matrix.h
	gcc $(FLAGS) -c division.c

mod_matrix.o: mod_matrix.h mod_matrix.c  set.h vector_operations.h
	gcc $(FLAGS) -c mod_matrix.c



 
 

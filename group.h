/*
 * group.h
 */

#ifndef GROUP_H_
#define GROUP_H_
#include "spmat.h"
#include <stdio.h>
#include <stdlib.h>

/* Struct group represents set of vertices- group g from the algorithms.
 * The struct contains the size of the group, an array g_nodes of the vertices,
 * an array kg of number of neighbours of each node, and Sub-sparse matrix Ag of A representing the group.
 *
 * allocate_group - Allocates new group
 * free_group - Frees the group
 * create_kg1_kg2 - Divides vector kg
 * create_first_group - Creates the initial group for the division.
 *
 * */

typedef struct _group{
	int size;
	int* g_nodes;
	double* kg;
	spmat* Ag;
}group;

/*Alloocates new group of size ng.
 * If error -function does nothing and program is terminated.
 * @param ng - size of the group
 * @param nnz - number of non zero values in Ag
 * @return pointer to the new group.
 * Error may occur because of malloc.
 */
group* allocate_group(int ng,int nnz);

/* Frees all resources used by g
 * @param g - group to be freed.
 * */
void free_group(group* g);


/* Creates kg1, kg2 from k according to vector s
 * @param k - Vector kg of original group
 * @param kg1 - Vector kg of first group
 * @param kg2 - Vector kg of the second group
 * @param ng - Size of the original group
 * @param s - Vector describes the division of the original group into 2 groups
 * s[i]=1 - Node i of g is in the first group
 * s[i]=-1 - Node i of g is in the second group
 * */
void create_kg1_kg2(double* k, double* kg1, double* kg2, int ng,double* s);


/*Creates and returns first group that contains all nodes from input file
 *Calculates M and puts into M_[0].
 *@param input - Input file
 *@param M_ - Address for the sum of degress of all vertices.
 *@return - Pointer to the initial group.
 *
 *Error may occur because of malloc.
 **/
group* create_first_group(FILE* input, double* M_);

#endif /* GROUP_H_ */

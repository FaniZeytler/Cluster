/*
 * spmat.h
 */

#ifndef SPMAT_H_
#define SPMAT_H_
#include <stdio.h>
#include <stdlib.h>

/*Struct for spmat as arrays */

 typedef struct arrays{
	int nnz;
	int* colind;
	int* rowptr;

}arrays;


/*Struct that  represents sparse matrix that is part of B_bar[g] for every group g.
 *Sparse matrix is represented with arrays.
 *
 * spmat_allocate_array - Allocates new sparse matrix.
 * free_spmat_array - Frees the sparse matrix.
 * mult_Ag_array  - Mults Ag with a vector.
 * sum_rowi_g_array - Calculate the sum of row i in the sparse matrix.
 * sum_Bg_col - Calculates sum of absolute for every column in B_bar[g].
 * update_sigma - Updates vector sigma that is used in calculation the score in modularity maximization.
 * calc_nnz_Ag_1_2 - Calculates nnz for Ag1 or Ag2
 * create_Ag1_Ag2 - Builds sub-sparse-matricies Ag1 and Ag2 according to the division of original group g.
 * createA - Builds the initial sparse matrix of the graph.
 */
typedef struct _spmat {
	/* Matrix size (n*n) */
	int		n;
	void	*private;
} spmat;


/* Allocates  new arrays sparse matrix .
 * @param n - Size of the sparse matrix
 * @param nnz - Number of non zero values in the sparse matrix
 * @return - Pointer to the new sparse matrix
 *
 * Error may occur because of malloc.
 *  */
spmat* spmat_allocate_array(int n, int nnz);

/*Frees all resources used by A
 */
void free_spmat_array(spmat* A);

/* Multiplies matrix Ag by vector bk, into result.
 * @param Ag - Sparse matrix
 * @param bk - Vector to mult with Ag
 * @param result - Vector that wiil contain the result of the mult. result is pre-allocated.
 *  */
void mult_Ag_array(spmat* Ag, const double *bk, double *result);

/*Returns sum of i-th row in Ag.
 * @param Ag - sparse matrix
 * @param i - index of row
 * @return sum of row i in Ag
 * */
double sum_rowi_g_array(spmat* Ag ,int i);

/*Calculates sum of absolute column i in B_bar[g] for every i into result.
 * @param Ag- Sparse matrix of group
 * @param ng - Size of the group
 * @param kg - Vector k of a group
 * @param f - Vector f that is part of B_bar[g]
 * @param M - Varible M of the original graph
 * @param result - Vector that contains the result for each column i
 * result is pre-allocated.
 *
 */
void sum_Bg_col(spmat* Ag, int ng, double* kg, double* f, double M,double* result);

/*Updates sigma according to change of index max_score_i in s
 * @param Ag- Sparse matrix of group
 * @param ng - Size of the group
 * @param kg - Vector k of a group
 * @param M - Varible M of the original graph
 * @param max_score_i - Index that changes in s in the modularity_maximization
 * @param ki - kg[max_score_i]
 * @param s_max - s[max_score_i]
 * @param sigma - The vector to be updated
 * The update of sigma - sigma[i]=sigma[i]+2*s[max_score_i]*B[g]_(max_score_i,i)
 * */
void update_sigma(spmat* Ag,int ng,double* kg,double M, int max_score_i,double ki,double s_max, double* sigma);

/*Calculates nnz (non-zero-values) in Ag1/2 (the new groups to create from g) according to x in s.
 * @param Ag - sparse matrix of original group g
 * @param s - Vector of division of g
 * s[i]=1 - Node i in g is in group1
 * s[i]=-1 - Node i in g is in group2
 * @param ng - Size of original group g
 * @param x - represents whether to calculate nnz1 or nnz2.
 * x=1 - nnz1
 * x=-1 - nnz2
 *  */
int calc_nnz_Ag_1_2(spmat* Ag, double* s,int ng,double x);
/*Calculates a row of Ag1/Ag2 according to x.
 * @param rowptr_g - rowptr of g1 or g2
 * @param colind_g - colind og g1 or g2
 * @param rowptr - rowptr of original Ag
 * @param colind - colind of original Ag
 * @param - represents whether to calculate for Ag1 or Ag2
 * @param s - Vector of division of g
 * @param ng - size of g
 * returns pointer to colind_g
 */
int* calc_row_Ag_1_2(int* rowptr_g,int* colind_g,int* rowptr, int* colind, int x,double* s,int ng);

/*Builds sub-sparse-matrices Ag1 and Ag2 from Ag according to s
 * @param s - Vector of division of g
 * @param Ag - spmat matrix of group g
 * @param ng - Size of group g
 * @param Ag1 - sub-sparse-matrix for group1
 * @param Ag2 - sub-sparse-matrix for group2
 * Ag1 and Ag2 are pre-allocated.
 * */
void create_Ag1_Ag2( double* s,spmat* Ag,int ng,spmat* Ag1, spmat* Ag2);


/*Creates sparse matrix A from input File and puts into vector k right values.
 * @param input - Input file
 * @param k - k[i] will be degree of node i
 * @param A - sparse matrix of the graph in the input file
 * A and K are pre-allocated.
 *
 *Error may occur because of malloc.
 */
void createA(FILE* input,double* k,spmat* A);



#endif /* SPMAT_H_ */

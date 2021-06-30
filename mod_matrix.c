/*
 * mod_matrix.c
 */
#include "mod_matrix.h"


/* Mults matrix B_bar[g] on vector b_k of size ng into result.
 * @param g - Group that currently in proccess
 * @param f - Vector f that realted to B_Bar[g]
 * @param bk - Vector to mult with B_bar[g]
 * @param result - Vector that will contain the result of the mult
 * @param k_kt_x - Help vector for the calculation: k_kt_x=k*k_transpose*bk
 * @param M - Sum of degrees of the original graph
 * Memory is pre-allocated to result*/
void mult_Bbar(group* g,double* f,double* bk,double* result,double* k_kt_x,double M){

	int ng=g->size;
	double* bk_=bk;
	double* result_=result;
	double kt_x=0;/*mult (k_transpose)*x*/
	spmat* Ag=g->Ag;
	double* kg=g->kg;
	int i=0;

	bk_=bk;
	mult_Ag_array(Ag,  bk, result);/*mult result= Ag* bk*/
	kt_x=mult_vectors( kg,bk,ng);/* kt_x= (k^t)*bk */
	mult_vector_scalar(kg,kt_x,k_kt_x,ng);/* k_kt_x= k(k^t)bk */
	Divide_M(k_kt_x,M,k_kt_x,ng);/* k_kt_x=k_kt_x/M*/
	sub_vectors(result,k_kt_x,result,ng);/* result = result- k_kt_x*/
	result_=result;

	for(i=0;i<ng;i++){
		result_[0]=result_[0]-f[0]*bk_[0];
		bk_++;
		f++;
		result_++;
	}
}

/* Creates vector f according to group g
 * @param g - Group that currently in proccess
 * @param f - f[i] will be the sum of the i-th row of B[g]
 * @param M - Sum of degrees of the original graph
 *
 * f is pre-allocated
 */
void createf(group* g,double* f ,double M){
	int ng=g->size;
	double* f_=f;
	double* kg=g->kg;
	double sum_k=0;
	double sum=0;
	spmat* Ag=g->Ag;
	int i;

	/* calculate sum of vector kg*/
	for(i=0;i<ng;i++){
		sum_k+=kg[0];
		kg++;
	}

	kg=g->kg;
	for(i=0;i<ng;i++){
		sum= sum_rowi_g_array( Ag ,i);/*Calculates the sum of row i of Ag*/
		f_[0]=sum-(kg[0]*sum_k/M);
		f_++;
		kg++;
	}
}


/* Returns norma1 of B_bar[g]
 * @param M - Sum of degrees of the original graph
 * @param g - Group that currently in proccess
 * @param f - f[i] is the sum of the i-th row of B[g]
 * @param help - Help vector for calculation
 * @return norma1 of B_bar[g]
 *
 * help is pre-allocated
 *
 *  */
double norma1_Bbar(double M,group* g,double* f,double* help){

	double max=0;
	int i=0;
	int ng=g->size;
	double* sum_col_=help;
	spmat* Ag=g->Ag;
	double* kg=g->kg;
	sum_Bg_col( Ag,ng,  kg,  f, M,help);
	max=help[0];
	sum_col_++;
	for(i=1;i<ng;i++){
		if(sum_col_[0]>max+EPSILON){
			max=sum_col_[0];
		}
		sum_col_++;
	}

	return max;
}


/* Substructs norma1 of B_bar from vector of f into f_shift
 * f_shift will be used in find eigen_value and find_eigen_vector
 * @param f  - Vector to be shifted
 * @param f_shift - Vector that will contain the result of the shift
 * @param norma1 - Norma1 that belongs to B_bar[g] that is in proccess
 * @param ng - size of f and f_shift and of group g taht is in proccess
 * f_shift is pre-allocated
 * */
void calc_fshift(double* f, double* f_shift, double norma1,int ng){
	int i=0;
	for(;i<ng;i++){
		f_shift[0]=f[0]-norma1;
		f_shift++;
		f++;
	}
}


/* Returns leading eigen value of B_bar[g].
 * @param M - Sum of degrees of the original graph
 * @param f - Shifted f of B_[g]. f[i] will be the sum of the i-th row of B[g] minus norma1
 * @param g - Group that currently in proccess
 * @param bk - Eigen vector of B_bar[g]
 * @param norma1 - Norma1 of B_bar[g]
 * @param help - Help vector for calculations
 * @param k_kt_x - Help vector for mult_B_Bar
 *
 * Error may occur because of division by zero.
 */
double find_eigenvalue(double M,double* f,group* g,double* bk,double norma1,double* help,double* k_kt_x){
	int ng=g->size;
	double result=0;
	double bk2=mult_vectors(bk,bk,ng); /* bk2=bk*bk */
	if(fabs(bk2)<EPSILON){
		printf("Division by zero in find_eigenvalue \n");
		exit(EXIT_FAILURE);
	}

	mult_Bbar(g,f,bk,help,k_kt_x,M); /* help=B_bar[g]*bk */
	result=mult_vectors(bk,help,ng)/bk2; /* result=(bk*help)/bk2 */
	result=result-norma1;

	return result;
}


/* Calculates eigen vector of B_bar[g] into bk.
  * @param M - Sum of degrees of the original graph
  * @param g - Group that currently in proccess
  * @param f - Shifted f of B_[g]. f[i] will be the sum of the i-th row of B[g] minus norma1.
  * @param ng  - size of g
  * @param k_kt_x - Help vector for mult_B_Bar
  * @param b0 - Intiail vector of the power iteration
  * @param bk - Vector that will contain the eigen_vector of B_bar[g]
 * memory has been allocated to bk
 *
 * Error may occur because of division by zero.
 */
void eigen_vector_B_bar(double M,group* g,double* f,int ng,double* k_kt_x,double* b0,double* bk){
	double* b_i;
	int i;
	unsigned long int iter=1;
	unsigned long int max_iter =PI_MAX*ng+200000;
	double norma;

	/* create random starting vector */
	b_i=b0;
	srand(time(NULL));
	for(i=0;i<ng;i++){
		b_i[0]=rand();
		b_i++;
	}

	mult_Bbar(g,f,b0,bk,k_kt_x,M);/* bk=B_bar*b0 */


	norma=calc_norma(bk,ng);/* claculate norma of B_bar*b0 */
	if(fabs(norma)<EPSILON){
		printf("Division by zero in eigen_vector_B_bar \n");
		exit(EXIT_FAILURE);
	}

	/* calculate bk=B_bar*b0/norma */
	b_i=bk;
	for(i=0;i<ng;i++){
		b_i[0]=b_i[0]/norma;
		b_i++;
	}

	while(check(b0,bk,ng)==0){
		iter++;
		if(iter>max_iter){
			printf("Infinite loop in power iteration \n");
			exit(EXIT_FAILURE);
		}
		swap(&bk,&b0);
		mult_Bbar(g,f,b0,bk,k_kt_x,M);/* bk=B_bar*b0 */

		/* claculate norma of B_bar*b0 */
		norma=calc_norma(bk,ng);
		if(fabs(norma)<0.00001){
			printf("Division by zero in calculate eigen vector \n");
			exit(EXIT_FAILURE);
		}

	    /* bk=B_bar*b0/norma */
		b_i=bk;
		for(i=0;i<ng;i++){
		   b_i[0]=b_i[0]/norma;
		   b_i++;
		}
	}


}


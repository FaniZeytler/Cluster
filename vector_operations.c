/*
 * vector_operations.c
 */
#include "vector_operations.h"
/*Turns vector into zero vector
 * @param a - Vector that will be zero vector.
 * @param size - Size of vector a
 * a is pre-allocated*/
void zero_vector(double* a, int size){
	int i;
	for(i=0;i<size;i++){
		a[0]=0;
		a++;
	}
}

/*Mults vector x of size n  with scalar t and puts into result.
 * result is pre-allocatsed*/
void mult_vector_scalar(double* x, double t, double* result,int n){
	int i=0;
	for(;i<n;i++){
		result[0]=x[0]*t;
		result++;
		x++;
	}
}

/*Divides vector x of size n in M into result
 * result is pre-allocatsed
 */
void Divide_M(double* x, double M, double* result,int n){
	int i=0;
	for(;i<n;i++){
		result[0]=x[0]/M;
		result++;
		x++;
	}
}

/*Subs x-y into result
 * result is pre-allocatsed
 */
void sub_vectors(double* x, double* y, double* result,int n){
	int i;
	for(i=0;i<n;i++){
		result[0]=x[0]-y[0];
		result++;
		x++;
		y++;
	}
}

/*Swaps pointers Ab_k and b_k*/
void swap(double** Ab_k, double** b_k){
	double* c=*b_k;
	*b_k=*Ab_k;
	*Ab_k=c;
}

/*Returns a*b
*/
double mult_vectors(double* a,double* b,int size){
	double result=0;
	int i=0;
	for(;i<size;i++){
		result+=a[0]*b[0];
		a++;
		b++;
	}
	return result;
}

/* Returns norma of vector bk of size ng*/
double calc_norma(double* bk,int ng){
	int i;
	double norma=0;
	for(i=0;i<ng;i++){
		norma+=bk[0]*bk[0];
		bk++;
	}
	norma=sqrt(norma);
	return norma;
}
/*Returns index with maximum value in array */
int find_max_index_array(long double* arr,int size){

	int i=0;
	long double temp;
	long double max=arr[0];
	int i_max=0;
	arr++;
	for(i=1;i<size;i++){
		temp=arr[0];
		if(arr[0]>max-EPSILON){
			i_max=i;
			max=temp;
		}
		arr++;
	}
	return i_max;
}

/*Checks if the vector b_k is close enough to Ab_k
 * used in find_eigenvalue*/
int check(double* b_k,double* Ab_k,int n){
	int i;
	for(i=0;i<n;i++){

		if(fabs(b_k[0]-Ab_k[0])>EPSILON){
			return 0;
		}
		b_k++;
		Ab_k++;
	}
	return 1;
}



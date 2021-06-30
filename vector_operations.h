/*
 * vector_operations.h
 */

#ifndef VECTOR_OPERATIONS_H_
#define VECTOR_OPERATIONS_H_
#include <math.h>
#define EPSILON 0.00001

/*Library of vector operations
 * zero_vector - Turns a vector into zero vector
 * mult_vector_scalar - Mults vector with scalar
 * Divide_M - Divide the vector with M
 * sub_vectors - Sub operation between two vectors
 * swap - Swap pointers of two vectors
 * mult_vectors - Mult operation between two vectors
 * calc_norma - Calculates norma of a vector
 * find_max_index_array - Finds the index of the maximum value in an array
 * check - Checks if two vectors are close to each other.
 *
 * */

/*Turns vector into zero vector
 * @param a - Vector that will be zero vector.
 * @param size - Size of vector a
 * a is pre-allocated*/
void zero_vector(double* a, int size);

/*Mults vector x of size n  with scalar t and puts into result.
 * result is pre-allocatsed*/
void mult_vector_scalar(double* x, double t, double* result,int n);

/*Divides vector x of size n in M into result
 * result is pre-allocatsed
 */
void Divide_M(double* x, double M, double* result,int n);

/*Subs x-y into result
 * result is pre-allocatsed
 */
void sub_vectors(double* x, double* y, double* result,int n);

/*Swaps pointers Ab_k and b_k*/
void swap(double** Ab_k, double** b_k);

/*Returns a*b
*/
double mult_vectors(double* a,double* b,int size);

/* Returns norma of vector bk of size ng*/
double calc_norma(double* bk,int ng);

/*Returns index with maximum value in array */
int find_max_index_array(long double* arr,int size);

/*Checks if the vector b_k is close enough to Ab_k
 * used in find_eigenvalue*/
int check(double* b_k,double* Ab_k,int n);

#endif /* VECTOR_OPERATIONS_H_ */

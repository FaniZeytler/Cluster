/*
 * mod_matrix.h
 */

#ifndef MOD_MATRIX_H_
#define MOD_MATRIX_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "vector_operations.h"
#include "set.h"
#define PI_MAX 50000

/* Library of functions that related to calculations with modularity matrix B_bar[g] or B[g].
 * The functions are used in the main algorithms of dividing a network.
 *
 * mult_Bbar - Mults matrix B_bar[g] of a group g with vector
 * createf - Calculates vector f of group g
 * norma1_Bbar - Calculates norma1 of matrix B_bar[g]
 * calc_fshift - Shifts the vector f with norma1 of B_bar[g]
 * find_eigenvalue - Finds the eigen value of B_bar[g]
 * eigen_vector_B_bar - Finds eigen vector of B_bar[g]
 */


/* Mults matrix B_bar[g] on vector b_k of size ng into result.
 * @param g - Group that currently in proccess
 * @param f - Vector f that realted to B_Bar[g]
 * @param bk - Vector to mult with B_bar[g]
 * @param result - Vector that will contain the result of the mult
 * @param k_kt_x - Help vector for the calculation: k_kt_x=k*k_transpose*bk
 * @param M - Sum of degrees of the original graph
 * Memory is pre-allocated to result*/
void mult_Bbar(group* g,double* f,double* bk,double* result,double* k_kt_x,double M);


/* Creates vector f according to group g
 * @param g - Group that currently in proccess
 * @param f - f[i] will be the sum of the i-th row of B[g]
 * @param M - Sum of degrees of the original graph
 *
 * f is pre-allocated
 */
void createf(group* g,double* f ,double M);

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
double norma1_Bbar(double M,group* g,double* f,double* help);

/* Substructs norma1 of B_bar from vector of f into f_shift
 * f_shift will be used in find eigen_value and find_eigen_vector
 * @param f  - Vector to be shifted
 * @param f_shift - Vector that will contain the result of the shift
 * @param norma1 - Norma1 that belongs to B_bar[g] that is in proccess
 * @param ng - size of f and f_shift and of group g taht is in proccess
 * f_shift is pre-allocated
 * */
void calc_fshift(double* f, double* f_shift, double norma1,int ng);

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
double find_eigenvalue(double M,double* f,group* g,double* bk,double norma1,double* help,double* k_kt_x);

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
void eigen_vector_B_bar(double M,group* g,double* f,int ng,double* k_kt_x,double* b0,double* bk);



#endif /* MOD_MATRIX_H_ */

/*
 * division.h
 */

#ifndef DIVISION_H_
#define DIVISION_H_
#include "mod_matrix.h"





/*Main algorithms for dividing a network and help functions of them.
 *
 * Algothims:
 * divide_network - Divides the network into groups.
 * Divide_2 - Divides one group into two groups or the trivial division of one group.
 * modularity_maximization - Maximizes the modularity of the division of a group
 *
 * Help functions:
 * find_max_score_i - Returns the index of a node that gives the maximium score if being moved
 * create_s_unit_vector - Turns s into unit vector
 * create_s -  Creates vector s according to eigen vector
 *  */




/* Returns the index of a node that gives the maximium score if being moved and puts
 * the maximum score into max_score.
 * Used in modularity maximization.
 * @param unmoved - unmoved[i]=0 if node i in g didnt move yet and unmoved[i]=1 otherwise
 * @param sigma - Help vector to calculate score for each k
 * @param g - Group that is currently being divided
 * @param s - Current division of g
 * s[i]=1 if node i of g  in group1
 * s[i]=-1 if node i of g in group2
 * @param M - Sum of degrees of the original graph
 * @param max_score - Address to put the max_score
 * @return index of a node in g that gives the max score
 */
int find_max_score_i(double* unmoved,double* sigma,group* g,double* s,double M, long double* max_score);


/* Maximizes the modulaity of given division s.
 * Changes s according to the new division. If no improvment has been found, the division s stays the same.
 * @param M - Sum of degrees of the original graph
 * @param g - Group that is currently being divided
 * @param s - Current division of g
 * @param k_kt_x - Help vector for mult_B_bar
 * @param help - Help vector for the calculations
 *
 * Error may occur because of malloc
 */
void modularity_maximization(double M, group* g, double* s,double* k_kt_x,double* help);




/*Turns s into unit vector
 * s is pre-allocated */
void create_s_unit_vector(double* s, int size);

/* Creates vector s according to eigen vector
 * s is pre-allocated
 * @param e_vector - Eigen vector
 * @param s - Vector that will represent the division of a group
 * @param ng - size of e_vector and s
 *  */
void create_s(double* e_vector,double* s,int ng);

/* Divdes group g into 2 groups g1 and g2 and puts into result. If the division is trival, g is added to result.
 * result is pre-allocated
 * @param g - Group that is currently being divided
 * @param result - Set that will conatian the division of g into groups
 * @param M - Sum of degrees of the original graph
 *
 * Error may occur because of malloc.
 */
void Divide_2(group* g,set* result,double M);

/* Returns the division of the network
 * Algorithm 3
 * @param g - Intial group that contains all nodes from the input graph
 * @param M - Sum of degrees of the original graph
 */
set* divide_network( group* g ,double M);



#endif /* DIVISION_H_ */

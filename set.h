/*
 * set.h
 */

#ifndef SET_H_
#define SET_H_
#include "linked_list_s.h"

/* Struct representing a set of groups contains size of the set in field n and the groups in linked list head.
 * The set used to P and O in algorithm 3.
 *
 * add_group - Adds new gorup to the set.
 * delete_first_group - Deletes the first group in the set.
 * allocate_set - Allocates new set.
 * free_set - Frees the set.
 * update_set - Inserts to given set the division into groups of given group.
 *
 */
typedef struct set{
	int n;
	linked_list_s* head;

}set;

/*Adds group x into the beggining of the linked list of s.
 * @param s - Set that the group x will be added to. s is pre-allocated.
 * @ param x - The group to add to s.
 *
 * Error may occur because of malloc.
 */
void add_group( set* s, group* x);

/*Returns first group of s, and deletes it from s
 **/
group* delete_first_group( set* s);

/*Allocates new set and returns the pointer to the new set
 * Error may occur because of malloc.*/
set* allocate_set();

/* Frees all resources used by a set.
 * @param s - set to be freed. */
void free_set(set* s);

/*Cerates division of g according to s and puts the new groups into result
 * @param s - Vector that represents the division of g.
 * @param g - The group to be divided.
 * @param result - The set to be updated with the division of g.result is pre-allocated.
 * If the division of g is one group - n of result will be 1.
 * If the division of g is two groups - n of result will be 2.
 * */
void update_set( double* s,group* g,set* result );

#endif /* SET_H_ */

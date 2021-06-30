/*
 * linked_list_s.h
 */

#ifndef LINKED_LIST_S_H_
#define LINKED_LIST_S_H_
#include "group.h"

/*Linked list is a node in set of groups containg one group and pointer to the next.
 *
 * free_linked_list_s - Frees all the items in head .
 * */

typedef struct linked_list_s
{
	group* g;
	struct linked_list_s *next;

}linked_list_s;

/* Frees all the items in head
 * @param head - the head of the linked_list*/
void free_linked_list_s(linked_list_s* head);

#endif /* LINKED_LIST_S_H_ */

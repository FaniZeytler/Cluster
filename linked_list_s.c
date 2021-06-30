/*
 * linked_list_s.c
 */

#include "linked_list_s.h"

/* Frees all the items in head
 * @param head - the head of the linked_list*/
void free_linked_list_s(linked_list_s* head){
	if (head!=NULL)
		{

		free_linked_list_s(head->next);
		free_group(head->g);
		free(head);
			}
}

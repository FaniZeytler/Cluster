/*
 * set.c
 */
#include "set.h"

/*Adds group x into the beggining of the linked list of s.
 * @param s - Set that the group x will be added to. s is pre-allocated.
 * @ param x - The group to add to s.
 *
 * Error may occur because of malloc.
 */
void add_group( set* s, group* x){
	linked_list_s* head;
	linked_list_s* new;
	if(s->n==0){
		s->head=(linked_list_s*)malloc(sizeof(linked_list_s));
		if(s->head==NULL){
			printf("Error in malloc in add_group \n" );
			exit(EXIT_FAILURE);
		}
		head=s->head;
		head->g=x;
		head->next=NULL;
	}
	else{
		new=(linked_list_s*)malloc(sizeof(linked_list_s));
		if(new==NULL){
			printf("Error in malloc in add_group \n" );
			exit(EXIT_FAILURE);
		}
		new->g=x;
		new->next=s->head;
		s->head=new;
	}
	s->n++;
}

/*Returns first group of s, and deletes it from s
 **/
group* delete_first_group( set* s){
	linked_list_s* head=s->head;
	group* g;
	if(s->n!=0){
		s->n--;
		g=head->g;
		s->head=head->next;
		free(head);
	}
	return g;
}

/*Allocates new set and returns the pointer to the new set
 * Error may occur because of malloc.*/
set* allocate_set(){
	set* s=(set*)malloc(sizeof(set));
	if(s==NULL){
		printf("Error in malloc in allocate_set \n" );
		exit(EXIT_FAILURE);
	}
	s->n=0;
	return s;
}

/* Frees all resources used by a set.
 * @param s - set to be freed. */
void free_set(set* s){
	free_linked_list_s(s->head);
	free(s);
}

/*Cerates division of g according to s and puts the new groups into result
 * @param s - Vector that represents the division of g.
 * @param g - The group to be divided.
 * @param result - The set to be updated with the division of g.result is pre-allocated.
 * If the division of g is one group - n of result will be 1.
 * If the division of g is two groups - n of result will be 2.
 * */
void update_set( double* s,group* g,set* result ){

	group* g1;
	group* g2;
	spmat* Ag;
	double* s_;
	int * g_nodes;
	int* g_nodes_1;
	int* g_nodes_2;
	int ng;
	int i;
	double x;
	int nnz1;
	int nnz2;
	int ng1=0;
	int ng2=0;
	ng=g->size;
	Ag=g->Ag;
	s_=s;
    /*calc ng1,ng2 */
	for(i=0;i<ng;i++){
		if(s_[0]==1)
		{
			ng1++;
		}
		else{
			ng2++;

		}
		s_++;
		}
	/*Size of one of the groups in the division is zero*/
	if(ng1==0 || ng2==0){
		add_group(result,g);
	}
	/*Two sizes of the group are not zero */
	else{
	x=1;
	nnz1=calc_nnz_Ag_1_2( Ag, s, ng, x);
	x=-1;
	nnz2=calc_nnz_Ag_1_2( Ag, s, ng, x);

    g1=allocate_group(ng1,nnz1);
    g2=allocate_group(ng2,nnz2);

    /*update g_nodes of g1 and g2 according to s and g*/
    g_nodes=g->g_nodes;
	g_nodes_1=g1->g_nodes;
	g_nodes_2=g2->g_nodes;
	s_=s;

	for(i=0;i<ng;i++){
		if(s_[0]==1){
			g_nodes_1[0]=g_nodes[0];
			g_nodes_1++;
		}
		else{
			g_nodes_2[0]=g_nodes[0];
			g_nodes_2++;
		}
		s_++;
		g_nodes++;
	}

	create_Ag1_Ag2(  s, Ag, ng,g1->Ag,g2->Ag);
	create_kg1_kg2(g->kg, g1->kg, g2->kg,  ng, s);
	add_group(result,g1);
	add_group(result,g2);
	free_group(g);
	}
}



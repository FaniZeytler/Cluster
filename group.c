/*
 * group.c
 */
#include "group.h"


/*Alloocates new group of size ng.
 * If error -function does nothing and program is terminated.
 * @param ng - size of the group
 * @param nnz - number of non zero values in Ag
 * @return pointer to the new group.
 * Error may occur because of malloc.
 */
group* allocate_group(int ng,int nnz){
	spmat* Ag=spmat_allocate_array(ng,nnz);
	group* g=(group*)malloc(sizeof(group));
	int* g_nodes=(int*)malloc(ng*sizeof(int));
	double* kg=(double*)malloc(ng*sizeof(double));
	if(g==NULL || g_nodes==NULL || kg==NULL){
		printf("Error in malloc in allocate_group \n" );
		exit(EXIT_FAILURE);
	}


	g->size=ng;
	g->g_nodes=g_nodes;
	g->kg=kg;
	g->Ag=Ag;

	return g;
}


/* Frees all resources used by g
 * @param g - group to be freed.
 * */
void free_group(group* g){
	free(g->g_nodes);
	free_spmat_array(g->Ag);
	free(g->kg);
	free(g);
}

/* Creates kg1, kg2 from k according to vector s
 * @param k - Vector kg of original group
 * @param kg1 - Vector kg of first group
 * @param kg2 - Vector kg of the second group
 * @param ng - Size of the original group
 * @param s - Vector describes the division of the original group into 2 groups
 * s[i]=1 - Node i of g is in the first group
 * s[i]=-1 - Node i of g is in the second group
 * */
void create_kg1_kg2(double* k, double* kg1, double* kg2, int ng,double* s){
	int i=0;
	for(i=0;i<ng;i++){
		if(s[0]==1){
			kg1[0]=k[0];
			kg1++;
		}
		else{
			kg2[0]=k[0];
			kg2++;
		}
		k++;
		s++;
	}

}



/*Creates and returns first group that contains all nodes from input file
 *Calculates M and puts into M_[0].
 *@param input - Input file
 *@param M_ - Address for the sum of degress of all vertices.
 *@return - Pointer to the initial group.
 *
 *Error may occur because of malloc.
 **/
group* create_first_group(FILE* input, double* M_){
	int n_node;/* number of neighbors of i*/
	int* line;
	int i;
	group* g;
	int n;
	int nnz=0; /* number of non-zero values in spmat */
	int* g_nodes;
	int temp;

	temp=fread(&n,sizeof(int),1,input);
	if(temp!=1){
		printf("Error in reading the Input File \n");
		exit(EXIT_FAILURE);
	}

	 /* calc nnz */
	line=(int*)malloc(sizeof(int)*n);
	if(line==NULL){
		printf("Error in malloc in create_first_group \n" );
		exit(EXIT_FAILURE);
	}
	for(i=0;i<n;i++){
		 temp=fread(&n_node,sizeof(int),1,input);
			if(temp!=1){
				printf("Error in reading the Input File \n");
				exit(EXIT_FAILURE);
			}
		 nnz+=n_node;
		 if(n_node!=0){
		 temp=fread(line,sizeof(int),n_node,input);
		 if(temp!=n_node){
			 printf("Error in reading the Input File \n");
			 exit(EXIT_FAILURE);
		 }
		 }
	}
	g=allocate_group(n,nnz);

	/*Update g_nodes of the first group*/
	g_nodes=g->g_nodes;
	for(i=0;i<n;i++){
		g_nodes[0]=i;
		g_nodes++;
	}
	rewind(input);
	temp=fread(&n,sizeof(int),1,input);
	if(temp!=1){
		printf("Error in reading the Input File \n");
		exit(EXIT_FAILURE);
	}
	createA( input, g->kg, g->Ag);/* Creates the sparse matrix A of the original graph and updates g->kg */
	M_[0]=nnz;

	free(line);
	return g;
}

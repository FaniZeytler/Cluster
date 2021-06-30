/*
 * cluster.c
 */
#include "division.h"


/*Creates vector representing group to be written into output file.
 * @param g_nodes - the nodes in the group
 * @param ng - size of the group
 * @param res - result vector to be written in the output file.
 * */
void create_vector_g(int* g_nodes, int ng,int* res){
	int* res_;
	int i;
	res[0]=ng;
	res_=res;
	res_++;

	for(i=0;i<ng;i++){
		res_[0]=g_nodes[0];
		res_++;
		g_nodes++;
	}
}
/*Writes the division of groups in the output file.
 * @param set_i - Head of linked list of the groups of the division
 * @param n - Size of the original graph
 * @param out - Output File
 *
 * Error mat occur because of malloc or writing into output file.
 */
void write_output(linked_list_s* set_i,int n,FILE* out){
	int temp;
	int * vec_g=(int*)malloc((n+1)*sizeof(int));
	 if(vec_g==NULL){
		printf("Error in malloc in main \n");
		exit(EXIT_FAILURE);
	}
	 while(set_i!=NULL){
		group* g_i=set_i->g;
		create_vector_g(g_i->g_nodes, g_i->size,vec_g);
		temp=fwrite(vec_g,sizeof(int),g_i->size+1,out);/* Writes the information about every group in the division*/
		 if(temp!=(g_i->size+1)){
			 printf("Error in writing into the Output File \n");
			 exit(EXIT_FAILURE);
		 }
		set_i=set_i->next;
	}
	 free(vec_g);
}

/* Main function for the clustering proccess.
 * @param argv[1]- input file.
 * @param argv[2] - output file.
 * @ return
 *  0 -success of the clustering.
 *  other then 0 - Error has occurred -unsuccessfull clustering.
 *
 */
int main(int argc, char* argv[])
{
	FILE* input;
	FILE* out;
	group* g;
	set* result;
	linked_list_s* set_i;


	int n;
	double M;
	int temp;

	if(argc!=3){
		printf("Error in Input \n");
		exit(EXIT_FAILURE);
	}

	input=fopen(argv[1],"rb");
	out=fopen(argv[2],"wb");
	if(input==NULL){
		printf("Input File is invalid \n");
		exit(EXIT_FAILURE);
	}
	if(out==NULL){
		printf("Output File is invalid \n");
		exit(EXIT_FAILURE);
	}
	 g=create_first_group(input,&M);
	 if(M==0){
		 printf("M=0->can't divide by zero \n");
		 exit(EXIT_FAILURE);
	 }

	 n=g->Ag->n;
	 result=divide_network( g,M);/*Divides the network*/
	 set_i=result->head;
	 temp=fwrite(&(result->n),sizeof(int),1,out);/*Writes number of groups*/
	 if(temp!=1){
		 printf("Error in writing the into Output File \n");
		 exit(EXIT_FAILURE);
	 }
	 write_output( set_i, n, out);
	free_set(result);
	fclose(input);
	fclose(out);
	return 0;
}

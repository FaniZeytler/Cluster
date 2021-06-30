/*
 * division.c
 */
#include "division.h"


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
int find_max_score_i(double* unmoved,double* sigma,group* g,double* s,double M, long double* max_score){
	int k=0;
	long double score;
	int ng=g->size;
	int count=0;
	double* kg=g->kg;
	int max_score_i=0;
	for(k=0;k<ng;k++){
		if(unmoved[0]==0){
			long double temp=(kg[0]*kg[0])/M;
			score=-4*s[0]*sigma[0]-4*temp;
            s[0]=-s[0];
			if(count==0){
				max_score[0]=score;
				max_score_i=k;
				count=1;
			}

			else{
				if(max_score[0]+EPSILON<score){
					max_score[0]=score;
					max_score_i=k;
				}
			}
			s[0]=-s[0];
		}
		s++;
		sigma++;
		unmoved++;
		kg++;

	}
	return max_score_i;


}

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
void modularity_maximization(double M, group* g, double* s,double* k_kt_x,double* help){

	long double max_score;
	int max_score_i;
	int max_improve_i;
	long double* improve_;
	double delta_Q2=0;
	int ng=g->size;
	double* unmoved=(double*)calloc(ng,sizeof(double));
	double* sigma=( double*)malloc(ng*sizeof( double));
	int* indices=(int*)malloc(ng*sizeof(int));
	int* indices_=indices;
	long double* improve=(long double*)malloc(ng*sizeof(long  double));
	double* kg=g->kg;
	spmat* Ag=g->Ag;
	double* zero=help;

	zero_vector(zero,ng);
	improve_=improve;

	if(unmoved==NULL  || sigma==NULL || indices==NULL || improve==NULL  ){
		printf("Error in malloc in modularity_maximization \n" );
		exit(EXIT_FAILURE);
	}

	do{
		int i;

		zero_vector(unmoved,ng);
		for(i=0;i<ng;i++){
			if(i==0){
				mult_Bbar( g, zero,s, sigma, k_kt_x,M);/*sigma=B[g]*s,used in calcution of score of each k*/
			}

            max_score_i= find_max_score_i( unmoved, sigma, g, s, M,  &max_score);/*Calculation of score to each k in unmoved*/
			s[max_score_i]=-s[max_score_i];

			/*Updating  sigma according to change of index max_score_i in s*/
			update_sigma(Ag, ng, kg,M, max_score_i,kg[max_score_i],s[max_score_i], sigma);

			indices_[0]=max_score_i;
			indices_++;
			if(i==0){
				improve_[0]=max_score;
			}
			else{
				improve_[0]=max_score+improve_[-1];
			}
			improve_++;
			unmoved[max_score_i]=1;
		}

		improve_=improve;
		improve_=improve;
		indices_=indices;
		max_improve_i=find_max_index_array(improve,ng);/*Finds the maximum state of improve*/
		indices_+=ng-1;

		/*Updates the division according to the maximum state*/
		for(i=ng-1;i>=max_improve_i+1;i--){
			s[indices_[0]]=-s[indices_[0]];
			indices_--;
		}
		indices_=indices;
		if(max_improve_i==ng-1){
			delta_Q2=0;
		}
		else{
			delta_Q2=improve[max_improve_i];
		}

	}while(delta_Q2>EPSILON);

	free(unmoved);
	free(indices);
	free(improve);
	free(sigma);
}


/*Turns s into unit vector
 * s is pre-allocated */
void create_s_unit_vector(double* s, int size){
	int i;
	for(i=0;i<size;i++){
		s[0]=1;
		s++;
	}
}


/* Creates vector s according to eigen vector
 * s is pre-allocated
 * @param e_vector - Eigen vector
 * @param s - Vector that will represent the division of a group
 * @param ng - size of e_vector and s
 *  */
void create_s(double* e_vector,double* s,int ng){
	int i;


	for(i=0;i<ng;i++){

		if(e_vector[0]>-EPSILON){
			s[0]=1;
		}
		else{
			s[0]=-1;
		}
		s++;
		e_vector++;

	}
}


/* Divdes group g into 2 groups g1 and g2 and puts into result. If the division is trival, g is added to result.
 * result is pre-allocated
 * @param g - Group that is currently being divided
 * @param result - Set that will conatian the division of g into groups
 * @param M - Sum of degrees of the original graph
 *
 * Error may occur because of malloc.
 */
void Divide_2(group* g,set* result,double M){

	int ng=g->size;
	double* f=(double*)malloc(ng*sizeof(double));
	double* k_kt_x=(double*)malloc(ng*sizeof(double));
	double* bk=(double*)malloc(ng*sizeof(double));
	double* b0=(double*)malloc(ng*sizeof(double));
	double* help=(double*)malloc(ng*sizeof(double));
	double* s=(double*)malloc(ng*sizeof(double));

	double* e_vector=bk;

	double e_value;
	double sBs;
	double norma1;
	if(f==NULL ||k_kt_x==NULL ||bk==NULL ||b0==NULL ||help==NULL || s==NULL){
		printf("Error in malloc in Divide Network \n" );
		exit(EXIT_FAILURE);
	}

	createf(g,f,M);
	norma1=norma1_Bbar( M,g, f,help);
	calc_fshift( f, f,  norma1, ng);/* To shift f with norma1 */
	eigen_vector_B_bar( M, g,f,ng,k_kt_x,b0,bk);
	e_value=find_eigenvalue( M,f,g,e_vector, norma1,help,k_kt_x);
	calc_fshift( f,f,(-1)*norma1, ng);/*To unshift f with norma1 */

	if(e_value<EPSILON){
		create_s_unit_vector(s,ng);
	}
	else{
		create_s( e_vector, s, ng);

		mult_Bbar(  g, f, s,help,k_kt_x,M);
		sBs=mult_vectors( s,help,ng);
		if(sBs<EPSILON){
			create_s_unit_vector(s,ng);

		}
	}
	      modularity_maximization(M, g,s,k_kt_x,help);
	      update_set(  s, g,result );/*Puts the division into one group or two of g into result*/

		free(f);
		free(help);
		free(k_kt_x);
		free(b0);
		free(bk);
		free(s);
}


/* Returns the division of the network
 * Algorithm 3
 * @param g - Intial group that contains all nodes from the input graph
 * @param M - Sum of degrees of the original graph
 */
set* divide_network( group* g ,double M){


	set* P=allocate_set();
	set* O=allocate_set();
	set* result=allocate_set();
	add_group(P,g);

	while(P->n!=0){
		group* g=delete_first_group(P);
		Divide_2(g,result,M);

		if(result->n==1){
			add_group(O,delete_first_group(result));
		}
		else{
			while(result->n!=0){
				group* temp= delete_first_group(result);
				if(temp->size==1){
					add_group(O,temp);
				}
				else{
					add_group(P,temp);
				}
			}
		}
	}
	free_set(P);
	free_set(result);
	return O;
}

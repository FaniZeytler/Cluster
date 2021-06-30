/*
 * spmat.c
 */


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "spmat.h"
#include "math.h"

/*Frees all resources used by A
 */
void free_spmat_array(spmat* A){
	arrays* a;
	int* rowptr;
	int* colind;
	a=(arrays*)(A->private);
	rowptr=(int*)(a->rowptr);
	colind=(int*)(a->colind);
	free(rowptr);
	free(colind);
	free(a);
	free(A);

}

/* Multiplies matrix Ag by vector bk, into result.
 * @param Ag - Sparse matrix
 * @param bk - Vector to mult with Ag
 * @param result - Vector that wiil contain the result of the mult. result is pre-allocated.
 *  */
void mult_Ag_array(spmat* Ag, const double *bk, double *result){
	arrays* a;
	int* rowptr;
	int* colind;
	int i;
	int j;
	double sum;
	int nnz;
	int ng=Ag->n;
	 a=(arrays*)(Ag->private);

	 rowptr=(int*)(a->rowptr);
	 colind=(int*)(a->colind);
	nnz=a->nnz;
	for(i=0;i<ng;i++){
		sum=0;

		for(j=rowptr[0];j<rowptr[1];j++){

			sum+=bk[colind[0]];
			if(j!=(nnz-1)){
			colind++;
			}
		}
		result[0]=sum;
		result++;
		rowptr++;
	}
}


/* Allocates  new arrays sparse matrix .
 * @param n - Size of the sparse matrix
 * @param nnz - Number of non zero values in the sparse matrix
 * @return - Pointer to the new sparse matrix
 *
 * Error may occur because of malloc.
 *  */
spmat* spmat_allocate_array(int n, int nnz){
	arrays* a;
	spmat *mat=(spmat*)malloc(sizeof(spmat));
	mat->private=(arrays*)malloc(sizeof(arrays));
	a=(arrays*)mat->private;
	a->colind=(int*)malloc(nnz*sizeof(int));
	a->rowptr=(int*)malloc((n+1)*sizeof(int));
	if(mat==NULL ||mat->private==NULL  || a->colind==NULL || a->rowptr==NULL){
		printf("error in malloc in spmat_allocate_array \n" );
		exit(EXIT_FAILURE);
	}
	a->rowptr[n]=nnz;
	a->nnz=nnz;
	mat->n=n;
	return mat;

}

/*Returns sum of i-th row in Ag.
 * @param Ag - sparse matrix
 * @param i - index of row
 * @return sum of row i in Ag
 * */
double sum_rowi_g_array(spmat* Ag,int i){
    arrays* a;
	int* rowptr;
	int* colind;
	int count;
	int end;

	a=(arrays*)(Ag->private);
	rowptr=(int*)(a->rowptr);
	colind=(int*)(a->colind);

	rowptr+=i;
	colind+=rowptr[0];
	count=rowptr[0];
	end=rowptr[1];

	return (end-count);
}


/*Calculates sum of absolute column i in B_bar[g] for every i into result.
 * @param Ag- Sparse matrix of group
 * @param ng - Size of the group
 * @param kg - Vector k of a group
 * @param f - Vector f that is part of B_bar[g]
 * @param M - Varible M of the original graph
 * @param result - Vector that contains the result for each column i
 * result is pre-allocated.
 *
 */
void sum_Bg_col(spmat* Ag, int ng, double* kg, double* f, double M,double* result){

    arrays* a;
	int* rowptr;
	int* colind;
	double sum_i;
	double* result_=result;
	double* kg_1=kg;
	double* kg_2=kg;
	double* f_=f;
	int i;
	int j;
	 a=(arrays*)(Ag->private);
	 rowptr=(int*)(a->rowptr);
	 colind=(int*)(a->colind);

	 for(i=0;i<ng;i++){
		 int count;
		 int end;
		 int ki=kg_1[0];

		 j=0;
		 count=rowptr[0];
		 end=rowptr[1];
		 kg_2=kg;
		 sum_i=0;

		 /*Go through row i of Ag and according to the cell in Ag calculate the absoulte of the cell B_bar[g]_ij
		  * and add it to sum of row i*/
		 while(j<ng){
			if(count!=end){
			 if(j==colind[0]){
				 if(i==j){
					 sum_i+=fabs(1-((ki*kg_2[0])/M)-f_[0]);
				 }
				 else{
				 sum_i+=fabs(1-((ki*kg_2[0])/M));
				 }
				 j++;
				 kg_2++;
				 colind++;
				 count++;
			 }
			 else{
				 if(colind[0]>j){
					 if(i==j){
						 sum_i+=fabs(0-(ki*kg_2[0])/M-f_[0]);
					 }
					 else{
					 sum_i+=fabs(0-(ki*kg_2[0])/M);
					 }

					 j++;
					 kg_2++;
				 }
				 else{
					 colind++;
					 count++;
				 }

			 }

		 }
			 else{
				 if(i==j){
					 sum_i+=fabs(0-(ki*kg_2[0])/M-f_[0]);
				 }
				 else{
				 sum_i+=fabs(0-(ki*kg_2[0])/M);
				 }
				 j++;
				 kg_2++;
			 }
		 }
		 result_[0]=sum_i;/* Puts sum of absoulute of row i of B_bar[g] in result[i]*/
		 f_++;
		 kg_1++;
		 result_++;
		 rowptr++;
	 }
}

/*Updates sigma according to change of index max_score_i in s
 * @param Ag- Sparse matrix of group
 * @param ng - Size of the group
 * @param kg - Vector k of a group
 * @param M - Varible M of the original graph
 * @param max_score_i - Index that changes in s in the modularity_maximization
 * @param ki - kg[max_score_i]
 * @param s_max - s[max_score_i]
 * @param sigma - The vector to be updated
 * The update of sigma - sigma[i]=sigma[i]+2*s[max_score_i]*B[g]_(max_score_i,i)
 * */
void update_sigma(spmat* Ag,int ng,double* kg,double M, int max_score_i,double ki,double s_max, double* sigma){
	    arrays* a;
		int* rowptr;
		int* colind;
		double* kg_2=kg;
		int j;
		int count;
	    int end;

		a=(arrays*)(Ag->private);
		rowptr=(int*)(a->rowptr);
		colind=(int*)(a->colind);
		rowptr+=max_score_i;
		colind+=rowptr[0];
		j=0;
		count=rowptr[0];
		end=rowptr[1];
		kg_2=kg;

		/*Adds to sigma the change needed according to change of index max_score_i in s(from mosularity_maximization).
		 * sigma[i]=sigma[i]+2*s[max_score_i]*B[g]_(max_score_i,i)
		 */
		while(j<ng){
		  if(count!=end){
			if(j==colind[0]){
				sigma[0]=sigma[0]+2*s_max*(1-((ki*kg_2[0])/M));
				j++;
				sigma++;
				kg_2++;
				colind++;
				count++;
				 }
			else{
				if(colind[0]>j){
					sigma[0]=sigma[0]+2*s_max*(0-((ki*kg_2[0])/M));

					j++;
					kg_2++;
					sigma++;
				}
				else{
					colind++;
					count++;
				}
			}
			}
		  else{
				sigma[0]=sigma[0]+2*s_max*(0-((ki*kg_2[0])/M));

				j++;
				sigma++;
				kg_2++;
			}
	 }
}


/*Calculates nnz (non-zero-values) in Ag1/2 (the new groups to create from g) according to x in s.
 * @param Ag - sparse matrix of original group g
 * @param s - Vector of division of g
 * s[i]=1 - Node i in g is in group1
 * s[i]=-1 - Node i in g is in group2
 * @param ng - Size of original group g
 * @param x - represents whether to calculate nnz1 or nnz2.
 * x=1 - nnz1
 * x=-1 - nnz2
 *  */
int calc_nnz_Ag_1_2(spmat* Ag, double* s,int ng,double x){
	arrays* a;
	int* rowptr;
	int* colind;
	int nnz=0;
	int i=0;
	double* s_1=s;
	double* s_2=s;

	 a=(arrays*)(Ag->private);
	 rowptr=(int*)(a->rowptr);
	 colind=(int*)(a->colind);
	 for(i=0;i<ng;i++){
		 if(s_1[0]==x){
			 int count=rowptr[0];
			 int end=rowptr[1];
			 int prev=0;
			 s_2=s;
			 while(count<end){
				 s_2+=colind[0]-prev;
				 prev=colind[0];
				 if(s_2[0]==x){
					 nnz++;
				 }
				 colind++;
				 count++;
			 }
		 }
		 else{
			 colind+=rowptr[1]-rowptr[0];
		 }
		 rowptr++;
		 s_1++;
	 }
	 return nnz;
}

/*Calculates a row of Ag1/Ag2 according to x.
 * @param rowptr_g - rowptr of g1 or g2
 * @param colind_g - colind og g1 or g2
 * @param rowptr - rowptr of original Ag
 * @param colind - colind of original Ag
 * @param - represents whether to calculate for Ag1 or Ag2
 * @param s - Vector of division of g
 * @param ng - size of g
 * returns pointer to colind_g
 */
int* calc_row_Ag_1_2(int* rowptr_g,int* colind_g,int* rowptr, int* colind, int x,double* s,int ng){
	int j=0;
	int k=0;
	int size_i=0;

	int count=rowptr[0];
	int end=rowptr[1];
	double* s_2=s;
	while(count<end && k<ng){
		if(colind[0]==k){
			if(s_2[0]==x){
				colind_g[0]=j;
				colind_g++;
				j++;
				size_i++;
			}
			colind++;
			count++;
			s_2++;
			k++;
		}
		else{
			if(colind[0]>k){
				if(s_2[0]==x){
					j++;
				}
				s_2++;
				k++;
			}
			else{
				colind++;
				count++;
			}
		}

	}
	rowptr_g[0]=rowptr_g[-1]+size_i;
	return colind_g;
}

/*Builds sub-sparse-matrices Ag1 and Ag2 from Ag according to s
 * @param s - Vector of division of g
 * @param Ag - spmat matrix of group g
 * @param ng - Size of group g
 * @param Ag1 - sub-sparse-matrix for group1
 * @param Ag2 - sub-sparse-matrix for group2
 * Ag1 and Ag2 are pre-allocated.
 * */
void create_Ag1_Ag2( double* s,spmat* Ag,int ng,spmat* Ag1, spmat* Ag2){
	double x;
    arrays* a_g1;
	int* rowptr_g1;
	int* colind_g1;
	arrays* a_g2;
	int* rowptr_g2;
	int* colind_g2;
	int i;
	double* s_1;
	arrays* a;
	int* rowptr;
	int* colind;

	 a=(arrays*)(Ag->private);
	 rowptr=(int*)(a->rowptr);
	 colind=(int*)(a->colind);

	 a_g1=(arrays*)(Ag1->private);
	 rowptr_g1=(int*)(a_g1->rowptr);
	 colind_g1=(int*)(a_g1->colind);

	 a_g2=(arrays*)(Ag2->private);
	 rowptr_g2=(int*)(a_g2->rowptr);
	 colind_g2=(int*)(a_g2->colind);

	s_1=s;
	rowptr_g1[0]=0;
	rowptr_g2[0]=0;
	rowptr_g1++;
	rowptr_g2++;
	for(i=0;i<ng;i++){
		if(s_1[0]==1){
			x=1;
			colind_g1=calc_row_Ag_1_2( rowptr_g1, colind_g1, rowptr, colind,  x, s, ng);/*Updates colind_g1 according to s and colind_g*/
			rowptr_g1++;
		}
		else{
			x=-1;
			colind_g2=calc_row_Ag_1_2( rowptr_g2, colind_g2, rowptr, colind,  x, s, ng);/*Updates colind_g2 according to s and colind_g*/
			rowptr_g2++;

		}
		colind+=rowptr[1]-rowptr[0];

		rowptr++;
		s_1++;

	}


}

/*Creates sparse matrix A from input File and puts into vector k right values.
 * @param input - Input file
 * @param k - k[i] will be degree of node i
 * @param A - sparse matrix of the graph in the input file
 * A and K are pre-allocated.
 *
 *Error may occur because of malloc.
 */
void createA(FILE* input,double* k,spmat* A){

	int j;
	int* rowptr;
	int* colind;
	int* line=(int*)malloc(sizeof(int)*A->n);/*neighbors of i*/
	int* line_;
	double* k_=k;
	int n_node;
	int temp;
	int i;
	int n=A->n;

	arrays* a=(arrays*)(A->private);
	if(line==NULL){
		   printf("Error in malloc in createB \n ");
		   exit(EXIT_FAILURE);
	}
	rowptr=(int*)(a->rowptr);
	rowptr[0]=0;
	rowptr++;
	colind=(int*)(a->colind);

/*Reads the neighbours of each node in the input and updates rowptr and colind of A*/
	for(i=0;i<n;i++){
			temp= fread(&n_node,sizeof(int),1,input);
			if(temp!=1){
				printf("Error in reading the Input File \n");
				 exit(EXIT_FAILURE);
			}
			k_[0]=n_node;
			k_++;
			if(n_node!=0){
			temp=fread(line,sizeof(int),n_node,input);
			 if(temp!=n_node){
				 printf("Error in reading the Input File \n");
				 exit(EXIT_FAILURE);
			 }
			line_=line;

			for(j=0;j<n_node;j++){
				colind[0]=line_[0];
				colind++;
				line_++;
			}
			}
			rowptr[0]=rowptr[-1]+n_node;
			rowptr++;

	}
	free(line);


}




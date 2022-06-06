/*
 *  smoother.c
 *
 *  Created by Xuehai Huang on 10/30/2008
 *  Modified by Chensong Zhang on 03/27/2009
 *  Copyright 2008 PSU. All rights reserved.
 *
 */
 
/*! \file smoother.c
 *  \brief Smoothers for multigrid method
 */
  
#include <stdio.h>
#include <math.h>
#include "header.h"
#include "matvec.h"

/**
 * \fn void gs(dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, int L)
 * \brief Gauss-Seidel method  as the smoother in solving Au=b with multigrid method
 * \param u initial guess and the new approximation to the solution obtained after L Gauss-Seidel steps
 * \param i_1 the index to begin with
 * \param i_n the index to end
 * \param s the step (s=1: forward, s=-1: backward)
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side
 * \param L number of iterations
 * \return void
 *
 * Gauss-Seidel smoother (Smoother_Type = 1)
 */
void gs(dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, int L)
{
	int i,j,k;
	double t,d;
	int *ia=A->IA,*ja=A->JA;
	double *aj=A->val,*bval=b->val,*uval=u->val;
	
	if (s > 0) {
		while (L--) {
			for (i=i_1;i<=i_n;i+=s) {
				t=bval[i];
				int begin_row=ia[i],end_row=ia[i+1]-1;
				for (k=begin_row;k<=end_row;k++) {
					j=ja[k];
					if (i!=j) t-=aj[k]*uval[j];
					else d=aj[k];
				} // end for k
				
				uval[i]=t/d;
			} // end for i
		} // end while		
	} // if s
	else {		
		while (L--) {
			for (i=i_1;i>=i_n;i+=s) {
				t=bval[i];
				int begin_row=ia[i],end_row=ia[i+1]-1;
				for (k=begin_row;k<=end_row;k++) {
					j=ja[k];
					if (i!=j) t-=aj[k]*uval[j];
					else d=aj[k];
				} // end for k
				
				uval[i]=t/d;
			} // end for i
		} // end while		
  } // end if
	return;

/*	int i,j,k,l=0;
	double t, d;
	
	for (l=0;l<L;l++)
	{
		for (i=i_1;i<=i_n;i+=s)
		{
			t=b->val[i];
			for (k=A->IA[i];k<A->IA[i+1];k++)
			{
				j=A->JA[k];
				if (i!=j)
					t-=A->val[k]*u->val[j];
				else
					d=A->val[k];
			}
			u->val[i]=t/d;
		}
	}*/
}

/**
 * \fn void jacobi(dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, int L)
 * \brief Jacobi method as the smoother in solving Au=b with multigrid method
 * \param u initial guess and the new approximation to the solution obtained after L Gauss-Seidel steps
 * \param i_1 the index to begin with
 * \param i_n the index to end
 * \param s the step
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side
 * \param L number of iterations
 * \return void
 *
 * Jacobi smoother (Smoother_Type = 2)
 */
void jacobi(dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, int L)
{
	int i,j,k,l=0;
	dvector t, d;
	create_dvector(u->row, &t);
	create_dvector(u->row, &d);

	for (l=0;l<L;l++)
	{
		for (i=i_1;i<=i_n;i+=s)
		{
			t.val[i]=b->val[i];
			for (k=A->IA[i];k<A->IA[i+1];k++)
			{
				j=A->JA[k];
				if(i!=j)
 					t.val[i]-=A->val[k]*u->val[j];
				else
					d.val[i]=A->val[k];
			}
		}
		for (i=i_1;i<=i_n;i+=s) u->val[i]=t.val[i]/d.val[i];
	}

	free_dvector(&t);
	free_dvector(&d);

}

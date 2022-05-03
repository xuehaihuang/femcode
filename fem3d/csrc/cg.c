/*
 *  cg.c
 *
 *  Created by Chensong Zhang on 03/28/2009.
 *  Copyright 2009 PSU. All rights reserved.
 *
 */

/*! \file cg.c
 *  \brief Preconditioned Conjugate Gradient Method
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "precond.h"
#include "matvec.h"

/**
 * \fn int pcg(dCSRmat *A, dvector *b, dvector *u, int MaxIt, double tol, precond *pre, int print_level)
 *	 \brief A preconditioned conjugate gradient (CG) method for solving Au=b 
 *	 \param *A	 pointer to the coefficient matrix
 *	 \param *b	 pointer to the dvector of right hand side
 *	 \param *u	 pointer to the dvector of DOFs
 *	 \param MaxIt integer, maximal number of iterations
 *	 \param tol double float, the tolerance for stopage
 *	 \param *pre pointer to the structure of precondition (precond) 
 * \param print_level how much information to print out
 *	 \return the number of iterations
 */
int pcg(dCSRmat *A, dvector *b, dvector *u, int MaxIt, double tol, precond *pre, int print_level)
{
	int iter=0,m=A->row;
	double alpha, beta, error, temp1, temp2, tempb;
	double *p, *z, *r, *t;

	p=(double *)calloc(m,sizeof(double));
	z=(double *)calloc(m,sizeof(double));
	r=(double *)calloc(m,sizeof(double));
	t=(double *)calloc(m,sizeof(double));

	// (b,b)
	tempb=dot_array(m,b->val,b->val);
	
	// r = b-A*u
	copy_array(m,b->val,r);
	sparse_mv(-1.0,A,u->val,r);
	
	// z = B*r
	if (pre == NULL)
		copy_array(m,r,z); /* No preconditioner, B=I */
	else
		pre->fct(A,r,z,pre->data); /* Preconditioning */
	
	// p = z
	copy_array(m,z,p);
	
	// temp1=(z_{k-1},r_{k-1})
	temp1=dot_array(m,z,r);
	
	while(iter<MaxIt)
	{
		iter++;
		
		// t=A*p
		init_array(m,t,0.0);
		sparse_mv(1.0,A,p,t);

		// comupte alpha_k=(z_{k-1},r_{k-1})/(A*p_{k-1},p_{k-1})
		temp2=dot_array(m,t,p);
		alpha=temp1/temp2;
		
		// compute u_k=u_{k-1} + alpha_k*p_{k-1}
		axpy_array(m,alpha,p,u->val);
		
		// compute r_k=r_{k-1} - alpha_k*A*p_{k-1}
		sparse_mv(-alpha,A,p,r);
		
		temp2=dot_array(m,r,r);
		if(temp2<1e-30) {
			iter=iter*(-1)-1; 
			break;
		}
		
		// relative residual = ||b-Au||_2/||b||_2=||r||_2/||b||_2
		error=sqrt(temp2/tempb);		
		if (print_level>1)
			printf("Iteration %3d: relative residual = %e\n",iter,error);
		if (error<tol) break;
		
		// z_k = B*r_k
		if (pre == NULL)
			copy_array(m,r,z);	 /* No preconditioner, B=I */
		else
			pre->fct(A,r,z,pre->data); /* preconditioning */
		
		// compute beta_k = (z_k, r_k)/(z_{k-1}, r_{k-1})
		temp2=dot_array(m,z,r);
		beta=temp2/temp1;
		temp1=temp2;
		
		// compute p_k = z_k + beta_k*p_{k-1}
		axpby_array(m,1.0,z,beta,p);
	}
	
	if(iter<0) {
		iter=iter*(-1)-1;
		free(p);
		free(r);
		free(z);
		free(t);
		return iter;
	}
	
	if (print_level>0) {
		if (iter>MaxIt)
			printf("Maximal iteration %d exceeded with relative residual %e.\n", MaxIt, error);
		else
			printf("Number of iterations = %d with relative residual %e.\n", iter, error);
	}
	
	free(p);
	free(r);
	free(z);
	free(t);
	
	return iter;
}


/**
 * \fn int den_pcg(ddenmat *A, dvector *b, dvector *u, int MaxIt, double tol, den_precond *pre, int print_level)
 *	 \brief A preconditioned conjugate gradient (CG) method for solving Au=b 
 *	 \param *A	 pointer to the coefficient matrix
 *	 \param *b	 pointer to the dvector of right hand side
 *	 \param *u	 pointer to the dvector of DOFs
 *	 \param MaxIt integer, maximal number of iterations
 *	 \param tol double float, the tolerance for stopage
 *	 \param *pre pointer to the structure of precondition (den_precond) 
 * \param print_level how much information to print out
 *	 \return the number of iterations
 */
int den_pcg(ddenmat *A, dvector *b, dvector *u, int MaxIt, double tol, den_precond *pre, int print_level)
{
	int iter=0,m=A->row;
	double alpha, beta, error, temp1, temp2, tempb;
	double *p, *z, *r, *t;
	
	p=(double *)calloc(m,sizeof(double));
	z=(double *)calloc(m,sizeof(double));
	r=(double *)calloc(m,sizeof(double));
	t=(double *)calloc(m,sizeof(double));
	
	// (b,b)
	tempb=dot_array(m,b->val,b->val);
	
	// r = b-A*u
	copy_array(m,b->val,r);
	denmat_mv(-1.0,A,u->val,r);
	
	temp2=dot_array(m,r,r);
	if(temp2<1e-30) {
		free(p);
		free(r);
		free(z);
		free(t);
		return 0;
	 }
	
	// z = B*r
	if (pre == NULL)
		copy_array(m,r,z); /* No preconditioner, B=I */
	else
		pre->fct(A,r,z,pre->data); /* Preconditioning */
	
	// p = z
	copy_array(m,z,p);
	
	// temp1=(z_{k-1},r_{k-1})
	temp1=dot_array(m,z,r);
	
	while(iter<MaxIt)
	{
		iter++;
		
		// t=A*p
		init_array(m,t,0.0);
		denmat_mv(1.0,A,p,t);
		
		// comupte alpha_k=(z_{k-1},r_{k-1})/(A*p_{k-1},p_{k-1})
		temp2=dot_array(m,t,p);
		alpha=temp1/temp2;
		
		// compute u_k=u_{k-1} + alpha_k*p_{k-1}
		axpy_array(m,alpha,p,u->val);
		
		// compute r_k=r_{k-1} - alpha_k*A*p_{k-1}
		denmat_mv(-alpha,A,p,r);
		
		temp2=dot_array(m,r,r);
		/*	if(temp2<1e-30) {
		 iter=iter*(-1)-1; break;
		 }*/
		
		// relative residual = ||b-Au||_2/||b||_2=||r||_2/||b||_2
		error=sqrt(temp2/tempb);		
		if (print_level>1)
			printf("Iteration %3d: relative residual = %e\n",iter,error);
		if (error<tol) break;
		
		// z_k = B*r_k
		if (pre == NULL)
			copy_array(m,r,z);	 /* No preconditioner, B=I */
		else
			pre->fct(A,r,z,pre->data); /* preconditioning */
		
		// compute beta_k = (z_k, r_k)/(z_{k-1}, r_{k-1})
		temp2=dot_array(m,z,r);
		beta=temp2/temp1;
		temp1=temp2;
		
		// compute p_k = z_k + beta_k*p_{k-1}
		axpby_array(m,1.0,z,beta,p);
	}
	
	if(iter<0) {
		iter=iter*(-1)-1;
		free(p);
		free(r);
		free(z);
		free(t);
		return iter;
	}
	
	if (print_level>0) {
		if (iter>MaxIt)
			printf("Maximal iteration %d exceeded with relative residual %e.\n", MaxIt, error);
		else
			printf("Number of iterations = %d with relative residual %e.\n", iter, error);
	}
	
	free(p);
	free(r);
	free(z);
	free(t);
	
	return iter;
}


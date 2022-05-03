/*
 *  multigrid.c
 *  classicAMG
 *
 *  Created by Xuehai Huang on 2/2/09.
 *  Modified by Chensong Zhang on 03/30/2009
 *  Copyright 2009 PSU. All rights reserved.
 *
 */

/*! \file multigrid.c
 *  \brief Abstract multigrid cycle
 */
 
#include <stdlib.h>

#include "precond.h"
#include "matvec.h"

/**
 * \fn void multigrid(dCSRmat *A, dvector *b, dvector *x, dCSRmat *R, dCSRmat *P,
 *                    int level, int levelNum, int smoother, int m1, int m2, int m0)
 * \brief Solve Ax=b by abstract multigrid cycle
 *
 * \param *A pointer to stiffness matrix of levelNum levels
 * \param *b pointer to the dvector of right hand side term
 * \param *x pointer to the dvector of dofs
 * \param *R pointer to resriction operator array of levelNum levels
 * \param *P pointer to interpolation operator array of levelNum levels
 * \param level current level
 * \param levelNum total level num of grid
 * \param smoother smoother type
 * \param m1 pre-smoothing times
 * \param m2 post-smoothing times
 * \param m0 times of correction on coarse grid
 */
void multigrid(dCSRmat *A, dvector *b, dvector *x, dCSRmat *R, dCSRmat *P, 
							 int level, int levelNum, int smoother, int m1, int m2, int m0)
{
	if(level<levelNum-1)
	{
		int i, j;
		dvector newb, e;
		int Rsize=R[level].row, Asize=A[level].row;
		double *r = (double*)calloc(Asize, sizeof(double));
		
		create_dvector(Rsize,&newb);
		create_dvector(Rsize,&e);

		/** pre smoothing */
		if (smoother == GS) {
			gs(x, 0, A[level].row-1, 1, &A[level], b, m1);
		}
		else if (smoother == JACOBI) {
			jacobi(x, 0, A[level].row-1, 1, &A[level], b, m1);
		}
		else if (smoother == SGS) {
			gs(x, 0, A[level].row-1, 1,  &A[level], b, m1);
			gs(x, A[level].row-1, 0, -1, &A[level], b, m2);
		}
		
		/** form residual r = b - A x */
		copy_array(Asize,b->val,r);
		sparse_mv(-1.0,&A[level],x->val,r);
		
		/** restriction */
		for(i=0;i<R[level].row;i++) {
			for(j=R[level].IA[i];j<R[level].IA[i+1];j++) {
				newb.val[i]+=R[level].val[j]*r[R[level].JA[j]];
			}
		}
				
		/** call MG recursively: m0 = 1 for V cycle, m0 = 2 for W cycle */
		for(i=0; i<m0; i++)
			multigrid(A, &newb, &e, R, P, level+1, levelNum, smoother, m1, m2, m0);
		
		/** prolongation */
		for(i=0;i<P[level].row;i++) {
			for(j=P[level].IA[i];j<P[level].IA[i+1];j++) {
				x->val[i]+=P[level].val[j]*e.val[P[level].JA[j]];
			}
		}
		
		/** post smoothing */
		if (smoother == GS) {
			gs(x, A[level].row-1, 0, -1, &A[level], b, m2);
		}
		else if (smoother == JACOBI) {
			jacobi(x, A[level].row-1, 0, -1, &A[level], b, m2);
		}
		else if (smoother == SGS) {
			gs(x, 0, A[level].row-1, 1,  &A[level], b, m1);
			gs(x, A[level].row-1, 0, -1, &A[level], b, m2);
		}
		
		free_dvector(&newb);
		free_dvector(&e);
		free(r);
	}
	else // coarsest level solver
	{
		int MaxNumIt = 100; 
		double CoarseTol = 1e-10;
		pcg(&A[level], b, x, MaxNumIt, CoarseTol, NULL, 0);
	}

}

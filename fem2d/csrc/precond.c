/*
 *  precond.c
 *  
 *
 *  Created by Chensong on 3/28/09.
 *  Copyright 2009 PSU. All rights reserved.
 *
 */

/*! \file precond.c
 *  \brief Preconditioners
 */
 
#include <stdio.h> 
#include <stdlib.h>

#include "precond.h"
#include "matvec.h"
 
/**
 * \fn void precond_null(dCSRmat *A, double *r, double *z, void *data)
 * \brief Do nothing preconditioner z=I*r
 * \param *A pointer to the stiffness matrix
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 */
void precond_null(dCSRmat *A, double *r, double *z, void *data)
{
	int m=A->row;	
	copy_array(m,r,z);
}

/**
 * \fn void precond_diag(dCSRmat *A, double *r, double *z, void *data)
 * \brief Diagonal preconditioner z=inv(D)*r
 * \param *A pointer to the stiffness matrix
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 */
void precond_diag(dCSRmat *A, double *r, double *z, void *data)
{
	int i, m=A->row;	
	dvector *diag = (dvector *)data;
	for (i=0;i<m;i++) {
		if (abs(diag->val[i])>1e-16) {
			z[i]=r[i]/diag->val[i];
		}
		else {
			z[i]=r[i];
		}
	}	
}

/**
 * \fn void precond_classicAMG(dCSRmat *A, double *r, double *z, void *data)
 * \brief get z from r by classic AMG
 * \param *A pointer to the stiffness matrix
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 */
void precond_classicAMG(dCSRmat *A, double *r, double *z, void *data)
{
	int i, m=A->row;
	dvector rr, zz;
	precond_data *amgdata=data;
	
	int smoother = amgdata->smoother;
	int pre = amgdata->presmooth_iter;
	int post = amgdata->postsmooth_iter;
	int maxit = amgdata->max_iter;
	int max_levels = amgdata->max_levels;
	
	zz.row=rr.row=m;
	rr.val=r; zz.val=z;
	init_dvector(&zz,0.0);
	
	for(i=0;i<maxit;i++)
		multigrid(amgdata->Aarray[0], &rr, &zz, amgdata->Aarray[1], 
						  amgdata->Aarray[2], 0, max_levels, smoother, pre, post, 1);	
}


/**
* \fn void precond_asP1ElasDG(dCSRmat *A, double *r, double *z, void *data)
* \brief get z from r by auxiliary space preconditioner of linear elasticity
* \param *A pointer to the stiffness matrix
* \param *r pointer to residual
* \param *z pointer to preconditioned residual
* \param *data pointer to precondition data
*/
void precond_asP1ElasDG(dCSRmat *A, double *r, double *z, void *data)
{
	int i, m = A->row;
	dvector rr, zz, z1, z2, r2;
	dvector z2s[2], r2s[2];
	precond_data *aspdata = data;

	int smoother = aspdata->smoother;
	int smooth_iter = aspdata->smooth_iter;
	int mg_smoother = aspdata->mg_smoother;
	int mg_smooth_iter = aspdata->mg_smooth_iter;
	int maxit = aspdata->max_iter;
	int levelNum = aspdata->max_levels;
	dCSRmat *R = aspdata->R;
	dCSRmat *P = aspdata->P;
	dCSRmat *As = aspdata->As;
	dCSRmat *Rs = aspdata->Rs;
	dCSRmat *Ps = aspdata->Ps;

	EDGE *edges = aspdata->edges;
	iCSRmat *edgesTran = aspdata->edgesTran;
	ivector *nodeCEdge = aspdata->nodeCEdge;
	// ivector *isInNode = aspdata->isInNode;
	// ivector *nondirichlet = aspdata->nondirichlet;
	// ivector *index = aspdata->index;
	int precond_type = aspdata->precond_type; // 1 additive; 2 multiplicative  

	zz.row = rr.row = m;
	rr.val = r; zz.val = z;

	init_array(m, z, 0);

	create_dvector(m, &z1);
	create_dvector(R->row, &z2);
	create_dvector(R->row, &r2);
	z2s[0].row = R->row / 2;
	z2s[1].row = R->row / 2;
	r2s[0].row = R->row / 2;
	r2s[1].row = R->row / 2;
	z2s[0].val = z2.val;
	z2s[1].val = z2.val + z2s[0].row;
	r2s[0].val = r2.val;
	r2s[1].val = r2.val + r2s[0].row;

	// additive
	if(precond_type ==1)
	{
		/** smoothing */
		if (smoother == JACOBI) {
			jacobi(&z1, 0, m - 1, 1, A, &rr, smooth_iter);
		}
		else if (smoother == SGS) {
			gs(&z1, 0, m - 1, 1, A, &rr, smooth_iter);
			gs(&z1, m - 1, 0, -1, A, &rr, smooth_iter);
		}

		sparse_mv0(1.0, R, r, r2.val);
		for (i = 0; i<maxit; i++)
			multigrid(As, &r2s[0], &z2s[0], Rs, Ps, 0, levelNum, mg_smoother, mg_smooth_iter, mg_smooth_iter, 1);
		for (i = 0; i<maxit; i++)
			multigrid(As, &r2s[1], &z2s[1], Rs, Ps, 0, levelNum, mg_smoother, mg_smooth_iter, mg_smooth_iter, 1);
		sparse_mv0(1.0, P, z2.val, zz.val);

		axpy_dvector(1.0, &z1, &zz);
	}
	else // multiplicative
	{
		/** smoothing */
		if (smoother == JACOBI) {
			jacobi(&zz, 0, m - 1, 1, A, &rr, smooth_iter);
		}
		else if (smoother == GS) {
			gs(&zz, 0, m - 1, 1, A, &rr, smooth_iter);
		}
		else if (smoother == SGS) {
			gs(&zz, 0, m - 1, 1, A, &rr, smooth_iter);
			gs(&zz, m - 1, 0, -1, A, &rr, smooth_iter);
		}

		/** form residual z1 = rr - A zz */
		copy_dvector(&rr, &z1);
		sparse_mv(-1.0, A, zz.val, z1.val);

		sparse_mv0(1.0, R, z1.val, r2.val);

		for (i = 0; i<maxit; i++)
			multigrid(As, &r2s[0], &z2s[0], Rs, Ps, 0, levelNum, mg_smoother, mg_smooth_iter, mg_smooth_iter, 1);
		for (i = 0; i<maxit; i++)
			multigrid(As, &r2s[1], &z2s[1], Rs, Ps, 0, levelNum, mg_smoother, mg_smooth_iter, mg_smooth_iter, 1);
		
		sparse_mv(1.0, P, z2.val, zz.val);

		/** smoothing */
		if (smoother == JACOBI) {
			jacobi(&zz, 0, m - 1, 1, A, &rr, smooth_iter);
		}
		else if (smoother == GS) {
			gs(&zz, m - 1, 0, -1, A, &rr, smooth_iter);
		}
		else if (smoother == SGS) {
			gs(&zz, 0, m - 1, 1, A, &rr, smooth_iter);
			gs(&zz, m - 1, 0, -1, A, &rr, smooth_iter);
		}
	}
	
	free_dvector(&z1);
	free_dvector(&r2);
	free_dvector(&z2);
}

/**
* \fn void precond_TriAsP1ElasDG(dCSRmat *A, double *r, double *z, void *data)
* \brief get z from r by block triangular preconditioner with auxiliary space method of linear elasticity
* discretized from mixed FEM for linear elasticity
*
*                 |M  B'| |sigma| = |f|
*                 |B -C | |u|     = |g|
*
* Use |D  B'| as the preconditioner in gmres and compute the inverse by
*     |B -C |
*
* the factorization
*
*                 |D  B'| |I Dinv*B'| = |D      0      |
*                 |B -C | |0   -I   |   |B B*Dinv*B'+C |
* More generally, we have
*                 |A11 A12| |I A11inv*A12| = |A11           0         |
*                 |A21 A22| |0     -I    |   |A21  A21*A11inv*A12-A22 |
* \param *A pointer to the stiffness matrix
* \param *r pointer to residual
* \param *z pointer to preconditioned residual
* \param *data pointer to precondition data
*/
void precond_TriAsP1ElasDG(dvector *r, dvector *z, void *data)
{
	int i, m = r[1].row;
	precond_data *aspdata = data;
	dvector *diag = aspdata->diag;
	dCSRmat *M = aspdata->precA[0];
	dCSRmat *Bt = aspdata->precA[1];
	dCSRmat *B = aspdata->precA[2];
	dCSRmat *Adg = aspdata->precA[3];
	double *scale = aspdata->precond_scale;

	// preconditioning r[0] by diagonal preconditioner
	for (i = 0; i < z[0].row; i++)
		z[0].val[i] = r[0].val[i] / diag->val[i];

	dvector tempVec;
	create_dvector(r[1].row, &tempVec);
	// tempVec = r[1] - B*z[0]
	copy_dvector(&r[1], &tempVec);
	sparse_mv(-1.0, B, z[0].val, tempVec.val);


	if (dot_dvector(&tempVec, &tempVec) < 1e-50)
		init_dvector(&z[1], 0);
	else
	{
		precond_asP1ElasDG(Adg, tempVec.val, z[1].val, data);
		/****  asP1ElasDG_PCG  ****
		precond *prec = (precond *)malloc(sizeof(precond));
		prec->data = data;
		prec->fct = precond_asP1ElasDG;

		// solver part
		int iter = pcg(Adg, &tempVec, &z[1], 5, 1e-8, prec, 0);
		//		printf("iter=%d\n",iter);
		/****  asP1ElasDG_PCG  ****/
	}

	free_dvector(&tempVec);

	int tri = 0;
	if(tri)
	{
//		axy_dvector(-1.0, &z[1], &z[1]);/////////////
	}
	else
	{
		// z[0] = z[0] + Dinv*Bt*z[1]
		create_dvector(r[0].row, &tempVec);
		sparse_mv(1.0, Bt, z[1].val, tempVec.val);
		for (i = 0; i < tempVec.row; i++)
			tempVec.val[i] = tempVec.val[i] / diag->val[i];
		axpy_dvector(1.0, &tempVec, &z[0]);
		free_dvector(&tempVec);

		// z[1] = -z[1]
		axy_dvector(-1.0, &z[1], &z[1]);
	}
	

	//	printf("%lf, %lf\n", scale[0],scale[1]);
	axy_dvector(scale[0], &z[0], &z[0]);
	axy_dvector(scale[1], &z[1], &z[1]);

}

/**
* \fn void precond_DiagAsP1ElasDG(dCSRmat *A, double *r, double *z, void *data)
* \brief get z from r by block diagonal preconditioner with auxiliary space method of linear elasticity
* \param *A pointer to the stiffness matrix
* \param *r pointer to residual
* \param *z pointer to preconditioned residual
* \param *data pointer to precondition data
*/
void precond_DiagAsP1ElasDG(dvector *r, dvector *z, void *data)
{
	int i, m = r[1].row;
	precond_data *aspdata = data;
	dvector *diag = aspdata->diag;
	dCSRmat *M = aspdata->precA[0];
	dCSRmat *Adg = aspdata->precA[3];
	double *scale = aspdata->precond_scale;

	// preconditioning r[0] by diagonal preconditioner
	for (i = 0; i < z[0].row; i++)
		z[0].val[i] = r[0].val[i] / diag->val[i];

//	gs(&z[0], 0, z[0].row - 1, 1, M, &r[0], 3);
//	gs(&z[0], z[0].row - 1, 0, -1, M, &r[0], 3);

	// preconditioning r[0] by ASP
//	copy_dvector(&r[1], &z[1]);
//	precond_asP1ElasDG(Adg, r[1].val, z[1].val, data);

//	for (i = 0; i < z[1].row; i++)
//		z[1].val[i] *= -1;
//	init_dvector(&z[1], 0);
	if(dot_dvector(&r[1], &r[1])<1e-50)
		init_dvector(&z[1], 0);
	else
	{
		precond_asP1ElasDG(Adg, r[1].val, z[1].val, data);

		/****  asP1ElasDG_PCG  ****
		precond *prec = (precond *)malloc(sizeof(precond));
		prec->data = data;
		prec->fct = precond_asP1ElasDG;

		// solver part
		int iter = pcg(Adg, &r[1], &z[1], 2, 1e-8, prec, 0);
//		printf("iter=%d\n",iter);
		/****  asP1ElasDG_PCG  ****/
	}
//	init_dvector(&z[1], 0);

//	printf("%lf, %lf\n", scale[0],scale[1]);
	axy_dvector(scale[0], &z[0], &z[0]);
	axy_dvector(scale[1], &z[1], &z[1]);

	/*printf("r[1].:\n");
	for (i = 0; i < r[1].row; i++)
		printf("%lf, ", r[1].val[i]);
	printf("\n");
	printf("z[1].:\n");
	for (i = 0; i < z[1].row; i++)
		printf("%lf, ", z[1].val[i]);
	printf("\n");
*/
}


/**
* \fn void precond_DiagAsP1ElasDGtemp(dCSRmat *A, double *r, double *z, void *data)
* \brief get z from r by block diagonal preconditioner with auxiliary space method of linear elasticity
* \param *A pointer to the stiffness matrix
* \param *r pointer to residual
* \param *z pointer to preconditioned residual
* \param *data pointer to precondition data
*/
void precond_DiagAsP1ElasDGtemp(dvector *r, dvector *z, void *data)
{
	int i, m = r[1].row;
	dvector rr, zz, z1, z2, r2;
	precond_data *aspdata = data;

	int smoother = aspdata->smoother;
	int smooth_iter = aspdata->smooth_iter;
	int mg_smoother = aspdata->mg_smoother;
	int mg_smooth_iter = aspdata->mg_smooth_iter;
	int maxit = aspdata->max_iter;
	int levelNum = aspdata->max_levels;
	dCSRmat *R = aspdata->R;
	dCSRmat *P = aspdata->P;
	dCSRmat *As = aspdata->As;
	dCSRmat *Adg = aspdata->precA[3];
	dvector *diag = aspdata->diag;

	// preconditioning r[0] by diagonal preconditioner
	for (i = 0; i < r[0].row; i++)
		z[0].val[i] = r[0].val[i] / diag->val[i];

	/*	printf("z[0].:\n");
	for (i = 0; i < z[0].row; i++)
	printf("%lf, ", z[0].val[i]);
	printf("\n");*/

	// preconditioning r[1] by ASP
	EDGE *edges = aspdata->edges;
	iCSRmat *edgesTran = aspdata->edgesTran;
	ivector *nodeCEdge = aspdata->nodeCEdge;
	ivector *isInNode = aspdata->isInNode;
	ivector *nondirichlet = aspdata->nondirichlet;
	ivector *index = aspdata->index;
	int precond_type = aspdata->precond_type; // 1 additive; 2 multiplicative  

	zz.row = rr.row = m;
	rr.val = r[1].val; zz.val = z[1].val;

	init_array(m, z[1].val, 0);

	//	if(dot_dvector(&rr, &rr)<1e-60)
	//		return;///////////////////////////////////////////////

	create_dvector(m, &z1);
	create_dvector(R->row, &z2);
	create_dvector(R->row, &r2);


	// additive
	if (precond_type == 1)
	{
		/** smoothing */
		if (smoother == JACOBI) {
			jacobi(&z1, 0, m - 1, 1, Adg, &rr, smooth_iter);
		}
		else if (smoother == SGS) {
			gs(&z1, 0, m - 1, 1, Adg, &rr, smooth_iter);
			gs(&z1, m - 1, 0, -1, Adg, &rr, smooth_iter);
		}

		sparse_mv0(1.0, R, rr.val, r2.val);
		for (i = 0; i<maxit; i++)
			mgvVectorP1_solve(As, &r2, &z2, edges, edgesTran, nodeCEdge, isInNode, nondirichlet, index, levelNum, mg_smoother, mg_smooth_iter);
		sparse_mv0(1.0, P, z2.val, zz.val);

		axpy_dvector(1.0, &z1, &zz);
	}
	else // multiplicative
	{
		/*		printf("rr1.:\n");
		for (i = 0; i < rr.row; i++)
		printf("%lf, ", rr.val[i]);
		printf("\n");
		printf("zz1.:\n");
		for (i = 0; i < zz.row; i++)
		printf("%lf, ", zz.val[i]);
		printf("\n");*/
		/** smoothing */
		if (smoother == JACOBI) {
			jacobi(&zz, 0, m - 1, 1, Adg, &rr, smooth_iter);
		}
		else if (smoother == GS) {
			gs(&zz, 0, m - 1, 1, Adg, &rr, smooth_iter);
		}
		else if (smoother == SGS) {
			gs(&zz, 0, m - 1, 1, Adg, &rr, smooth_iter);
			gs(&zz, m - 1, 0, -1, Adg, &rr, smooth_iter);
		}

		/*		printf("rra:\n");
		for (i = 0; i < 10; i++)
		printf("%lf, ", rr.val[i]);
		printf("\n");

		printf("zza:\n");
		for (i = 0; i < 10; i++)
		printf("%lf, ", zz.val[i]);
		printf("\n");*/

		/** form residual z1 = rr - A zz */
		copy_dvector(&rr, &z1);
		sparse_mv(-1.0, Adg, zz.val, z1.val);

		/*		printf("rr2.:\n");
		for (i = 0; i < rr.row; i++)
		printf("%lf, ", rr.val[i]);
		printf("\n");
		printf("zz2.:\n");
		for (i = 0; i < zz.row; i++)
		printf("%lf, ", zz.val[i]);
		printf("\n");
		printf("z1.:\n");
		for (i = 0; i < z1.row; i++)
		printf("%lf, ", z1.val[i]);
		printf("\n");*/

		/*		printf("z1a:\n");
		for (i = 0; i < 10; i++)
		printf("%lf, ", z1.val[i]);
		printf("\n");*/

		sparse_mv0(1.0, R, z1.val, r2.val);

		/*		printf("r2a:\n");
		for (i = 0; i < 10; i++)
		printf("%lf, ", r2.val[i]);
		printf("\n");
		printf("z2a:\n");
		for (i = 0; i < 10; i++)
		printf("%lf, ", z2.val[i]);
		printf("\n");*/
		for (i = 0; i<maxit; i++)
			mgvVectorP1_solve(As, &r2, &z2, edges, edgesTran, nodeCEdge, isInNode, nondirichlet, index, levelNum, mg_smoother, mg_smooth_iter);

		/*		printf("r2b:\n");
		for (i = 0; i < 10; i++)
		printf("%lf, ", r2.val[i]);
		printf("\n");
		printf("z2b:\n");
		for (i = 0; i < 10; i++)
		printf("%lf, ", z2.val[i]);
		printf("\n");*/


		sparse_mv(1.0, P, z2.val, zz.val);

		/** smoothing */
		if (smoother == JACOBI) {
			jacobi(&zz, 0, m - 1, 1, Adg, &rr, smooth_iter);
		}
		else if (smoother == GS) {
			gs(&zz, m - 1, 0, -1, Adg, &rr, smooth_iter);
		}
		else if (smoother == SGS) {
			gs(&zz, 0, m - 1, 1, Adg, &rr, smooth_iter);
			gs(&zz, m - 1, 0, -1, Adg, &rr, smooth_iter);
		}
	}

	free_dvector(&z1);
	free_dvector(&r2);
	free_dvector(&z2);
}

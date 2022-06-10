/*
 *  asElasDG.c
 *  auxiliary space preconditioner for linear elasticity discretized by dg
 *
 *  Created by Xuehai Huang on 02/29/2016.
 *  Copyright 2016 WZU. All rights reserved.
 *
 */

 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "precond.h"
#include "matvec.h"

 /**
 * \fn int mgvVectorP1AsElasDG(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
 * \brief Solve Ax=b by preconditioned conjugate gradient method (PCG),
 * with Auxiliary space method as precondition
 * \param *A	pointer to the dCSRmat matrix
 * \param *b	pointer to the dvector of right hand side
 * \param *x	pointer to the dvector of dofs
 * \param *param pointer to ASP parameters
 * \param print_level how much information to print out
 * \return the number of iterations
 */
int mgvVectorP1AsElasDG(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
{
	int levelNum = param->levelNum;
	int i;

	dCSRmat As[levelNum], tempA;
	dCSRmat P, PT;
	ELEMENT *elements = param->elements;
	idenmat *elementEdge = param->elementEdge;
	EDGE *edges = param->edges;
	dennode *nodes = param->nodes;
	iCSRmat *edgesTran = param->edgesTran;
	ivector *nodeCEdge = param->nodeCEdge;
	ELEMENT_DOF *elementDOFipdg = param->elementDOF;
	ELEMENT_DOF elementDOFas[levelNum];
	iCSRmat elementdofTranas[levelNum];
	ivector isInNode[levelNum], dirichlet[levelNum], nondirichlet[levelNum], index[levelNum];

	double lambda = param->lambda;
	double mu = param->mu;

	// setup preconditioner
	for (i = 0; i < levelNum; i++)
	{
		getElementDOF_Lagrange2d(&elementDOFas[i], &elements[i], &elementEdge[i], &edges[i], nodes[i].row, 1);
		getTransposeOfelementDoF(&elementDOFas[i], &elementdofTranas[i], 0);
		getBoundaryInfoVector2d(&edges[i], &nodes[i], elementDOFas[i].dof, elementDOFas[i].dop, &isInNode[i], &dirichlet[i], &nondirichlet[i], &index[i]);

		assembleStiffmatrixElasLagrange(&tempA, &elements[i], &elementEdge[i], &edges[i], &nodes[i], &elementDOFas[i], &elementdofTranas[i], mu);
		extractNondirichletMatrix11(&tempA, &As[i], &isInNode[i], &dirichlet[i], &nondirichlet[i], &index[i]);
		free_csr_matrix(&tempA);
	}
	interpVecP1toDG2d(&tempA, &elementDOFas[levelNum - 1], elementDOFipdg);
	extractNondirichletMatrix1c(&tempA, &P, &isInNode[levelNum - 1], &dirichlet[levelNum - 1], &index[levelNum - 1]);
	free_csr_matrix(&tempA);
	getTransposeOfSparse(&P, &PT);

	precond_data aspData;
	aspData.max_levels = levelNum;
	aspData.max_iter = param->mg_max_iter;
	aspData.tol = param->mg_tol;
	aspData.precond_type = param->precond_type;
	aspData.smoother = param->smoother;
	aspData.smooth_iter = param->smooth_iter;
	aspData.mg_smoother = param->mg_smoother;
	aspData.mg_smooth_iter = param->mg_smooth_iter;
	aspData.As = As;
	aspData.R = &PT;
	aspData.P = &P;

	aspData.edges = edges;
	aspData.edgesTran = edgesTran;
	aspData.nodeCEdge = nodeCEdge;
	aspData.isInNode = isInNode;
	aspData.nondirichlet = nondirichlet;
	aspData.index = index;

	// solver part
	int iter = 0;
	int MaxIt = param->max_iter;
	double tol = param->tol;
	double relres;

	clock_t solve_start, solve_end;
	solve_start = clock();
	while (++iter <= MaxIt) // MG solver here
	{
		relres = mgvVectorP1As_solve(A, b, x, &aspData);		
		printf("Iteration %3d: relative residual = %e\n", iter, relres);
		if (relres<tol) break;
	}
	solve_end = clock();

	double solveduration = (double)(solve_end - solve_start) / (double)(CLOCKS_PER_SEC);
	printf("Multigrid solve costs %f seconds.\n", solveduration);
	if (iter>MaxIt)
		printf("Maximal iteration %d exceeded with relative residual %e.\n", MaxIt, relres);
	else
		printf("Number of iterations = %d with relative residual %e.\n", iter, relres);

	for (i = 0; i<levelNum; i++)
	{
		free_csr_matrix(&As[i]);
		free_elementDOF(&elementDOFas[i]);
		free_icsr_matrix(&elementdofTranas[i]);
		free_ivector(&isInNode[i]);
		free_ivector(&dirichlet[i]);
		free_ivector(&nondirichlet[i]);
		free_ivector(&index[i]);
	}

	free_csr_matrix(&P);
	free_csr_matrix(&PT);

	return iter;
}

 /**
 * \fn double mgvVectorP1As_solve(dCSRmat *A, dvector *b, dvector *x, void *data)
 * \brief V-cycle geometric multigrid method for linear elasticity discretzed by conforming linear finite element method
 *        There is no Lame constants for this linear elasticity, i.e. lambda=0, mu=0.5
 *
 * \param *A pointer to stiffness matrix of levelNum levels
 * \param *b pointer to the dvector of right hand side term
 * \param *x pointer to the dvector of dofs
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
 * \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
 * \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
 * \param *isInNode pointer to boundary information of nodes: if the node is interior node, it will be 0; if the node is on the boundary, it will be -1
 * \param *nondirichlet pointer to the index of nondirichlet nodes
 * \param *index pointer to the transpose of dirichlet and nondirichlet nodes
 * \param levelNum total level num of grid
 * \param smoother smoother type
 * \param m smoothing times
 */
double mgvVectorP1As_solve(dCSRmat *A, dvector *b, dvector *x, void *data)
{
	int i, m = A->row;
	dvector r[2], e[2];
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

	EDGE *edges = aspdata->edges;
	iCSRmat *edgesTran = aspdata->edgesTran;
	ivector *nodeCEdge = aspdata->nodeCEdge;
	ivector *isInNode = aspdata->isInNode;
	ivector *nondirichlet = aspdata->nondirichlet;
	ivector *index = aspdata->index;


	create_dvector(m, &r[1]);
	create_dvector(m, &e[1]);
	create_dvector(R->row, &r[0]);
	create_dvector(R->row, &e[0]);


	/** form residual r = b - A x */
	copy_dvector(b, &r[1]);
	sparse_mv(-1.0, A, x->val, r[1].val);

	// Pre-smoothing
	if (smoother == GS) {
		gs(&e[1], 0, m - 1, 1, A, &r[1], smooth_iter);
	}
	else if (smoother == JACOBI) {
		jacobi(&e[1], 0, m - 1, 1, A, &r[1], smooth_iter);
	}
	else if (smoother == SGS) {
		gs(&e[1], 0, m - 1, 1, A, &r[1], smooth_iter);
		gs(&e[1], m - 1, 0, -1, A, &r[1], smooth_iter);
	}

	// Restrict the residual and keep it
	dvector temp;
	create_dvector(r[1].row, &temp);
	/** form residual temp = r - A e */
	copy_dvector(&r[1], &temp);
	sparse_mv(-1.0, A, e[1].val, temp.val);

	sparse_mv0(1.0, R, temp.val, r[0].val);
	free_dvector(&temp);

	for (i = 0; i<maxit; i++)
		mgvVectorP1_solve(As, &r[0], &e[0], edges, edgesTran, nodeCEdge, isInNode, nondirichlet, index, levelNum, mg_smoother, mg_smooth_iter);

	sparse_mv(1.0, P, e[0].val, e[1].val);

	// Post smoothing
	if (smoother == GS) {
		gs(&e[1], m - 1, 0, -1, A, &r[1], smooth_iter);
	}
	else if (smoother == JACOBI) {
		jacobi(&e[1], m - 1, 0, -1, A, &r[1], smooth_iter);
	}
	else if (smoother == SGS) {
		gs(&e[1], 0, m - 1, 1, A, &r[1], smooth_iter);
		gs(&e[1], m - 1, 0, -1, A, &r[1], smooth_iter);
	}

	// Update
	axpy_dvector(1.0, &e[1], x);

	// computer error
	copy_dvector(b, &r[1]);
	sparse_mv(-1.0, A, x->val, r[1].val);

	double aa = dot_dvector(&r[1], &r[1]);
	double bb = dot_dvector(b, b);
	double relres = sqrt(aa / bb);

	for (i = 0; i < 2; i++)
	{
		free_dvector(&r[i]);
		free_dvector(&e[i]);
	}

	return relres;
}


 /**
 * \fn void interpVecP1toDG2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFdg)
 * \brief the vector-version interpolation matrix from the 1st order Lagrange element to piecewise kth order polynomial in 2d
 * \param *P pointer to the vector-version interpolation matrix
 * \param *elementDOFp1 pointer to the relation between elements and degrees of freedom of the 1st order Lagrange element
 * \param *elementDOFdg pointer to the relation between elements and degrees of freedom of the piecewise kth order polynomial
 */
void interpVecP1toDG2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFdg)
{
	int i, j, ie, ii, k;
	int curnode[2];

	if(elementDOFp1->dop!=1)
	{
		P = NULL;
		return;
	}

	P->row = elementDOFdg->dof * 2;
	P->col = elementDOFp1->dof * 2;
	P->IA = (int*)calloc(P->row + 1, sizeof(int));
	P->JA = NULL;
	P->val = NULL;

	// step 1P: Find first the structure IA of the interpolation matrix P
	for (k = 0; k < elementDOFdg->row; k++)
	{
		for (i = 0; i < elementDOFdg->col; i++) //  for each node
		{
			curnode[0] = elementDOFdg->val[k][i];
			curnode[1] = curnode[0] + elementDOFdg->dof;
			P->IA[curnode[0] + 1] += elementDOFp1->col;
			P->IA[curnode[1] + 1] += elementDOFp1->col;
		}
	}

	for (i = 0; i<P->row; i++)
		P->IA[i + 1] += P->IA[i];

	P->nnz = P->IA[P->row];

	// step 2P: Find the structure JA of the interpolation matrix P
	P->JA = (int*)calloc(P->nnz, sizeof(int));
	for (k = 0; k < elementDOFdg->row; k++)
	{
		for (i = 0; i < elementDOFdg->col; i++) //  for each node
		{
			curnode[0] = elementDOFdg->val[k][i];
			curnode[1] = curnode[0] + elementDOFdg->dof;
			for (j = 0; j < elementDOFp1->col; j++)
			{
				P->JA[P->IA[curnode[0]] + j] = elementDOFp1->val[k][j];
				P->JA[P->IA[curnode[1]] + j] = elementDOFp1->val[k][j] + elementDOFp1->dof;
			}
		}
	}

	// step 3P: Loop element by element and compute the actual entries storing them in P
	P->val = (double*)calloc(P->nnz, sizeof(double));
	for (k = 0; k < elementDOFdg->row; k++)
	{
		if (elementDOFdg->dop == 0)
		{
			curnode[0] = elementDOFdg->val[k][0];
			curnode[1] = curnode[0] + elementDOFdg->dof;
			for (j = 0; j < 3; j++)
			{
				P->val[P->IA[curnode[0]] + j] = 1.0 / 3.0;
				P->val[P->IA[curnode[1]] + j] = 1.0 / 3.0;
			}
			continue;
		}

		for (i = 0; i < 3; i++) //  for each vertex
		{
			curnode[0] = elementDOFdg->val[k][i];
			curnode[1] = curnode[0] + elementDOFdg->dof;
			P->val[P->IA[curnode[0]] + i] = 1.0;
			P->val[P->IA[curnode[1]] + i] = 1.0;
		}

		for (ie = 0; ie < 3; ie++) //  for each dof in edge
		{
			for (ii = 0; ii < elementDOFdg->dop - 1; ii++)
			{
				curnode[0] = elementDOFdg->val[k][3 + ie*(elementDOFdg->dop - 1) + ii];
				curnode[1] = curnode[0] + elementDOFdg->dof;

				P->val[P->IA[curnode[0]] + ie] = 0;
				P->val[P->IA[curnode[0]] + (ie + 1) % 3] = ((double)elementDOFdg->dop - 1 - ii) / (double)elementDOFdg->dop;
				P->val[P->IA[curnode[0]] + (ie + 2) % 3] = (1.0 + ii) / (double)elementDOFdg->dop;

				for (j = 0; j < 3; j++)
					P->val[P->IA[curnode[1]] + j] = P->val[P->IA[curnode[0]] + j];
			}
		}

		if (elementDOFdg->dop > 2) //  for each dof in element
		{
			if (elementDOFdg->dop == 3)
			{
				curnode[0] = elementDOFdg->val[k][3*elementDOFdg->dop];
				curnode[1] = curnode[0] + elementDOFdg->dof;

				for (j = 0; j < 3; j++)
				{
					P->val[P->IA[curnode[0]] + j] = 1.0 / 3.0;
					P->val[P->IA[curnode[1]] + j] = 1.0 / 3.0;
				}
			}
			else if (elementDOFdg->dop == 4)
			{ 
				for (i = 0; i < 3; i++)
				{
					curnode[0] = elementDOFdg->val[k][3 * elementDOFdg->dop + i];
					curnode[1] = curnode[0] + elementDOFdg->dof;

					for (j = 0; j < 3; j++)
						P->val[P->IA[curnode[0]] + j] = 1.0 / 4.0;

					P->val[P->IA[curnode[0]] + i] = 2.0 / 4.0;

					for (j = 0; j < 3; j++)
						P->val[P->IA[curnode[1]] + j] = P->val[P->IA[curnode[0]] + j];
				}
			}
			else if (elementDOFdg->dop == 5)
			{
				for (i = 0; i < 3; i++)
				{
					curnode[0] = elementDOFdg->val[k][3 * elementDOFdg->dop + i];
					curnode[1] = curnode[0] + elementDOFdg->dof;

					for (j = 0; j < 3; j++)
						P->val[P->IA[curnode[0]] + j] = 1.0 / 5.0;

					P->val[P->IA[curnode[0]] + i] = 3.0 / 5.0;

					for (j = 0; j < 3; j++)
						P->val[P->IA[curnode[1]] + j] = P->val[P->IA[curnode[0]] + j];
				}

				for (ie = 0; ie < 3; ie++)
				{
					curnode[0] = elementDOFdg->val[k][3 * elementDOFdg->dop + 3 + ie];
					curnode[1] = curnode[0] + elementDOFdg->dof;

					for (j = 0; j < 3; j++)
						P->val[P->IA[curnode[0]] + j] = 2.0 / 5.0;

					P->val[P->IA[curnode[0]] + ie] = 1.0 / 5.0;

					for (j = 0; j < 3; j++)
						P->val[P->IA[curnode[1]] + j] = P->val[P->IA[curnode[0]] + j];
				}
			}
		} // elementDOFdg->dop > 2

	} // k
}
/*
 *  post_processing.c
 *  FEM
 *
 *  Created by Xuehai Huang on 03/13/12.
 *  Copyright 2012 WZU. All rights reserved.
 *
 */

/*! \file post_processing.c
 *  \brief Post-processing: Project the discontinuous Galerkin finite element solution to Raviart-Thomas element space
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "precond.h"
#include "matvec.h"


 /**
 * \fn void projL2PkVec2d(dvector *Qhu, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, void (*u)(double *, double *, double *), double *paras)
 * \brief L2 projection of the exact solution to the piecewise Lagrangian element space
 * \param *Qhu pointer to the L2 projection of u to the piecewise Lagrangian element space
 * \param *elements pointer to the structure of the triangulation
 * \param *nodes pointer to the nodes location of the triangulation
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
 */
void projL2PkVec2d(dvector *Qhu, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, void (*u)(double *, double *, double *), double *paras)
{
	int dop = elementDOF->dop;
	int dof = elementDOF->dof;
	int nt = elementDOF->row;
	ddenmat A;
	dvector b[2], w;
	create_dden_matrix(elementDOF->col, elementDOF->col, &A);
	create_dvector(A.col, &b[0]);
	create_dvector(A.col, &b[1]);
	create_dvector(A.col, &w);

	double phi[2], val[2];
	int i, j, k, i1;
	double x[2], **vertices, s;

	int num_qp, num_qp0;
	double lambdas[100][3], weight[100];
	double lambdas0[100][3], weight0[100];

	num_qp=getNumQuadPoints(dop * 2, 2); // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);

	num_qp0=49; // the number of numerical intergation points
	init_Gauss2d(num_qp0, lambdas0, weight0);

	for (k = 0; k < nt; k++)
	{
		// set parameters
		s = elements->vol[k];
		vertices = elements->vertices[k];
		// end set parameters

		init_dden_matrix(&A, 0.0);
		for (i = 0; i < A.row; i++)
		{
			for (j = 0; j < A.col; j++)
			{
				// A.val[i][j] = 0;
				for (i1 = 0; i1 < num_qp; i1++)
				{
					lagrange_basis(lambdas[i1], i, dop, &phi[0]);
					lagrange_basis(lambdas[i1], j, dop, &phi[1]);
					A.val[i][j] += s*weight[i1] * phi[0] * phi[1];
				}
			}
		}

		for (i = 0; i<A.row; i++)
		{
			b[0].val[i] = 0;
			b[1].val[i] = 0;
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lagrange_basis(lambdas0[i1], i, dop, &phi[0]);
				baryToCart2d(lambdas0[i1], x, vertices);
				u(x, val, paras);
				b[0].val[i] += s*weight0[i1] * val[0]*phi[0];
				b[1].val[i] += s*weight0[i1] * val[1]*phi[0];
			} // i1
		} // k1

		init_dvector(&w, 0.0);
		den_pcg(&A, &b[0], &w, 10000, 1e-14, NULL, 0);
		for (i = 0; i < A.row; i++)
			Qhu->val[elementDOF->val[k][i]] = w.val[i];
		init_dvector(&w, 0.0);
		den_pcg(&A, &b[1], &w, 10000, 1e-14, NULL, 0);
		for (i = 0; i < A.row; i++)
			Qhu->val[elementDOF->val[k][i] + dof] = w.val[i];
	}

	free_dden_matrix(&A);
	free_dvector(&b[0]);
	free_dvector(&b[1]);
	free_dvector(&w);
}

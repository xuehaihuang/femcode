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
* \fn void postprocess2newDisplacement(dvector *uhstar, dvector *sigmah, dvector *uh, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
* \brief L2 projection of the exact solution to the piecewise Lagrangian element space
* \param *uhstar pointer to postprocessed displacement
* \param *sigmah pointer to numerical solution
* \param *uh pointer to numerical solution
* \param *elements pointer to the structure of the triangulation
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void postprocess2newDisplacement(dvector *uhstar, dvector *sigmah, dvector *uh, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
{
	int dop1 = elementDOF[1].dop;
	int dop2 = elementDOF[2].dop;
	int ldof1 = elementDOF[1].col;
	int ldof2 = elementDOF[2].col;
	int nt = elements->row;

	ddenmat A;
	dvector b, x;
	create_dden_matrix((ldof1 + ldof2) * 2, (ldof1 + ldof2) * 2, &A);
	create_dvector(A.row, &b);
	create_dvector(A.row, &x);

	create_dvector(elementDOF[2].dof * 2, uhstar);

	double phi[2], phi1[3], phi2[3];
	int i, j, k, i1, j1, k1;
	double **gradLambda, s;
	double value[4];

	int num_qp, num_qp0;
	double lambdas[100][3], weight[100];
	double lambdas0[100][3], weight0[100];

	num_qp=getNumQuadPoints(dop1 + dop2, 2); // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);

	num_qp0=49; // the number of numerical intergation points
	init_Gauss2d(num_qp0, lambdas0, weight0);

	for (k = 0; k < nt; k++)
	{
		// set parameters
		s = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		// end set parameters

		// initial
		init_dden_matrix(&A, 0.0);
		init_dvector(&b, 0.0);

		// A
		for (i = 0; i < ldof2; i++)
		{
			// assemble top-left block of A
			for (j = 0; j < ldof2; j++)
			{
				for (i1 = 0; i1 < num_qp; i1++)
				{
					lagrange_basis1(lambdas[i1], gradLambda, i, dop2, phi1);
					lagrange_basis1(lambdas[i1], gradLambda, j, dop2, phi2);
					A.val[i][j] += s*weight[i1] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] / 2);
					A.val[i][j + ldof2] += s*weight[i1] * (phi1[1] * phi2[0] / 2);
					A.val[i + ldof2][j] += s*weight[i1] * (phi1[0] * phi2[1] / 2);
					A.val[i + ldof2][j + ldof2] += s*weight[i1] * (phi1[0] * phi2[0] / 2 + phi1[1] * phi2[1]);
				}
			}

			for (j = 0; j < ldof1; j++)
			{
				// assemble top-right block of A
				for (i1 = 0; i1 < num_qp; i1++)
				{
					lagrange_basis(lambdas[i1], i, dop2, &phi[0]);
					lagrange_basis(lambdas[i1], j, dop1, &phi[1]);
					A.val[i][j + ldof2 * 2] += s*weight[i1] * phi[0] * phi[1];
					A.val[i + ldof2][j + ldof1 + ldof2 * 2] += s*weight[i1] * phi[0] * phi[1];
				}
				// assemble bottom-left block of A
				A.val[j + ldof2 * 2][i] = A.val[i][j + ldof2 * 2];
				A.val[j + ldof1 + ldof2 * 2][i + ldof2] = A.val[i + ldof2][j + ldof1 + ldof2 * 2];
			}
		}

		// the first part of b
		for (i = 0; i<ldof2; i++)
		{
			for (i1 = 0; i1<num_qp0; i1++)
			{
				value[0] = 0;
				value[1] = 0;
				value[2] = 0;
				lagrange_basis1(lambdas0[i1], gradLambda, i, dop2, phi2);
				for (k1 = 0; k1<elementDOF[0].col; k1++)
				{
					huzhang_basis(lambdas0[i1], elements->nvector[k], elements->tvector[k], k1, elementDOF[0].dop, phi1);
					j1 = elementDOF[0].val[k][k1];
					value[0] += phi1[0] * sigmah->val[j1];
					value[1] += phi1[1] * sigmah->val[j1];
					value[2] += phi1[2] * sigmah->val[j1];
				}
				b.val[i] += s*weight0[i1] * (value[0] * phi2[0] + value[2] * phi2[1] - (value[0] + value[1]) * phi2[0] * lambda / (2 * lambda + 2 * mu)) / (2 * mu);
				b.val[i + ldof2] += s*weight0[i1] * (value[1] * phi2[1] + value[2] * phi2[0] - (value[0] + value[1]) * phi2[1] * lambda / (2 * lambda + 2 * mu)) / (2 * mu);
			} // i1
		} // i

		  // the second part of b
		for (i = 0; i<ldof1; i++)
		{
			for (i1 = 0; i1<num_qp; i1++)
			{
				value[0] = 0;
				value[1] = 0;
				lagrange_basis(lambdas[i1], i, dop1, &phi[1]);
				for (k1 = 0; k1<ldof1; k1++)
				{
					lagrange_basis(lambdas[i1], k1, dop1, &phi[0]);
					j1 = elementDOF[1].val[k][k1];
					value[0] += phi[0] * uh->val[j1];
					j1 += elementDOF[1].dof;
					value[1] += phi[0] * uh->val[j1];
				}
				b.val[i + ldof2 * 2] += s*weight[i1] * value[0] * phi[1];
				b.val[i + ldof1 + ldof2 * 2] += s*weight[i1] * value[1] * phi[1];
			} // i1
		} // i

		fgmres_den(&A, &b, &x, 80, 15000, 1e-15, NULL, 0);

		for (i = 0; i < ldof2; i++)
		{
			uhstar->val[elementDOF[2].val[k][i]] = x.val[i];
			uhstar->val[elementDOF[2].val[k][i] + elementDOF[2].dof] = x.val[i + ldof2];
		}
	}

	free_dden_matrix(&A);
	free_dvector(&b);
	free_dvector(&x);
}

/*
 *  error.c
 *  FEM
 *
 *  Created by Xuehai Huang on 06/10/09.
 *  Copyright 2009 PSU. All rights reserved.
 *
 */

/*! \file error.c
 *  \brief compute error between numerical solution and exact solution
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "header.h"
#include "matvec.h"

/**
 * \fn void geterrorslinearElasHuZhang2d(double *errors, dvector *sigmah, dvector *uh, dvector *Qhu, dvector *uhstar, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
 * \brief compute error between numerical solution and exact solution: L2 norm, Energy norm
 * \param *errors pointer to error between numerical solution and exact solution: L2 norm, H1 norm, Energy norm
 * \param *sigmah pointer to numerical solution
 * \param *uh pointer to numerical solution
 * \param *Qhu pointer to the L2 projection of u to the piecewise Lagrangian element space
 * \param *uhstar pointer to postprocessed displacement
* \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
 the fourth column stores -1 if the edge is on boundary
 * \param *nodes pointer to the nodes location of the triangulation
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param lambda Lame constant
 * \param mu Lame constant or Poisson ratio of plate
 */
void geterrorslinearElasHuZhang2d(double *errors, dvector *sigmah, dvector *uh, dvector *Qhu, dvector *uhstar, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
{
	int Nt = elements->row;

	int i, j, k, l;

	for (i = 0; i < 10; i++)
		errors[i] = 0;

	double phi0, phi1[3], phi2[3], val[4], paras[2];
	int k1, i1, j1, l1, l2;
	double x[2], s, **vertices, **gradLambda, **nv, **tv;
	double value[4];

	paras[0] = lambda;
	paras[1] = mu;

	int num_qp[5], num_qp1;
	double lambdas[5][100][3], weight[5][100];
	double lambdas1[100][2], weight1[100];

	num_qp[0]=49;
	num_qp[1]=49;
	num_qp[2]=49;
	num_qp[3]=getNumQuadPoints(elementDOF[1].dop * 2 - 2, 2);
	num_qp[4]=getNumQuadPoints(elementDOF[1].dop * 2, 2);
	for(i=0;i<5;i++)
		init_Gauss2d(num_qp[i], lambdas[i], weight[i]);

	num_qp1=getNumQuadPoints(elementDOF[1].dop * 2, 1);
	if (num_qp1>5) num_qp1 = 5;
	init_Gauss1d(num_qp1, lambdas1, weight1);

	for (k = 0; k<Nt; k++)
	{
		// set parameters
		s = elements->vol[k];
		vertices = elements->vertices[k];
		gradLambda = elements->gradLambda[k];
		nv = elements->nvector[k];
		tv = elements->tvector[k];
		// end set parameters

		// Energy norm and L2 norm of sigma-sigma_h
		for (i1 = 0; i1<num_qp[0]; i1++)
		{
			value[0] = 0;
			value[1] = 0;
			value[2] = 0;
			for (k1 = 0; k1<elementDOF[0].col; k1++)
			{
				huzhang_basis(lambdas[0][i1], nv, tv, k1, elementDOF[0].dop, phi1);
				j1 = elementDOF[0].val[k][k1];
				value[0] += phi1[0] * sigmah->val[j1];
				value[1] += phi1[1] * sigmah->val[j1];
				value[2] += phi1[2] * sigmah->val[j1];
			}
			baryToCart2d(lambdas[0][i1], x, vertices);
			linearElas2d_sigma(x, val, paras);
			value[0] -= val[0];
			value[1] -= val[1];
			value[2] -= val[2];
			value[3] = value[0] + value[1];
			value[0] = value[0] * value[0];
			value[1] = value[1] * value[1];
			value[2] = value[2] * value[2];
			value[3] = value[3] * value[3];
			if (lambda>-0.5)
				errors[0] += s*weight[0][i1] * (value[0] + value[1] + 2 * value[2] - value[3] * lambda / (2 * lambda + 2 * mu)) / (2 * mu);
			else
				errors[0] += s*weight[0][i1] * (value[0] + value[1] + 2 * value[2] - value[3] * 1.0 / 2.0) / (2 * mu);

			errors[9] += s*weight[0][i1] * (value[0] + value[1] + 2 * value[2]);
		}

		// L2 norm of divergence of sigma-sigma_h
		for (i1 = 0; i1<num_qp[1]; i1++)
		{
			value[0] = 0;
			value[1] = 0;
			for (k1 = 0; k1<elementDOF[0].col; k1++)
			{
				huzhang_basisDIV(lambdas[1][i1], gradLambda, nv, tv, k1, elementDOF[0].dop, phi1);
				j1 = elementDOF[0].val[k][k1];
				value[0] += phi1[0] * sigmah->val[j1];
				value[1] += phi1[1] * sigmah->val[j1];
			}
			baryToCart2d(lambdas[1][i1], x, vertices);
			linearElas2d_f(x, val, paras);
			value[0] += val[0];
			value[1] += val[1];
			errors[1] += s*weight[1][i1] * lpnormp_array(2, value, 2);
		}

		// L2 norm of u-u_h
		for (i1 = 0; i1<num_qp[2]; i1++)
		{
			value[0] = 0;
			value[1] = 0;
			for (k1 = 0; k1<elementDOF[1].col; k1++)
			{
				lagrange_basis(lambdas[2][i1], k1, elementDOF[1].dop, &phi0);
				j1 = elementDOF[1].val[k][k1];
				value[0] += phi0*uh->val[j1];
				j1 += elementDOF[1].dof;
				value[1] += phi0*uh->val[j1];
			}
			baryToCart2d(lambdas[2][i1], x, vertices);
			linearElas2d_u(x, val, paras);
			value[0] -= val[0];
			value[1] -= val[1];
			errors[3] += s*weight[2][i1] * lpnormp_array(2, value, 2);
		}

		// H1 semi-norm of Qhu-u_h
		for (i1 = 0; i1<num_qp[3]; i1++)
		{
			value[0] = 0;
			value[1] = 0;
			value[2] = 0;
			value[3] = 0;
			for (k1 = 0; k1<elementDOF[1].col; k1++)
			{
				lagrange_basis1(lambdas[3][i1], gradLambda, k1, elementDOF[1].dop, phi1);
				j1 = elementDOF[1].val[k][k1];
				value[0] += phi1[0] * (Qhu->val[j1] - uh->val[j1]);
				value[1] += phi1[1] * (Qhu->val[j1] - uh->val[j1]);
				j1 += elementDOF[1].dof;
				value[2] += phi1[0] * (Qhu->val[j1] - uh->val[j1]);
				value[3] += phi1[1] * (Qhu->val[j1] - uh->val[j1]);
			}
			errors[5] += s*weight[3][i1] * lpnormp_array(4, value, 2);
		}

		// L2 norm of Qhu-u_h
		for (i1 = 0; i1<num_qp[4]; i1++)
		{
			value[0] = 0;
			value[1] = 0;
			for (k1 = 0; k1<elementDOF[1].col; k1++)
			{
				lagrange_basis(lambdas[4][i1], k1, elementDOF[1].dop, &phi0);
				j1 = elementDOF[1].val[k][k1];
				value[0] += phi0*(Qhu->val[j1]- uh->val[j1]);
				j1 += elementDOF[1].dof;
				value[1] += phi0*(Qhu->val[j1] - uh->val[j1]);
			}
			errors[6] += s*weight[4][i1] * lpnormp_array(2, value, 2);
		}

		// H1 semi-norm of u-uhstar
		for (i1 = 0; i1<num_qp[2]; i1++)
		{
			value[0] = 0;
			value[1] = 0;
			value[2] = 0;
			value[3] = 0;
			for (k1 = 0; k1<elementDOF[2].col; k1++)
			{
				lagrange_basis1(lambdas[2][i1], gradLambda, k1, elementDOF[2].dop, phi1);
				j1 = elementDOF[2].val[k][k1];
				value[0] += phi1[0] * uhstar->val[j1];
				value[1] += phi1[1] * uhstar->val[j1];
				j1 += elementDOF[2].dof;
				value[2] += phi1[0] * uhstar->val[j1];
				value[3] += phi1[1] * uhstar->val[j1];
			}
			baryToCart2d(lambdas[2][i1], x, vertices);
			linearElas2d_gradu(x, val, paras);
			value[0] -= val[0];
			value[1] -= val[1];
			value[2] -= val[2];
			value[3] -= val[3];
			errors[7] += s*weight[2][i1] * lpnormp_array(4, value, 2);
		}

		// L2 norm of u-uhstar
		for (i1 = 0; i1<num_qp[2]; i1++)
		{
			value[0] = 0;
			value[1] = 0;
			for (k1 = 0; k1<elementDOF[2].col; k1++)
			{
				lagrange_basis(lambdas[2][i1], k1, elementDOF[2].dop, &phi0);
				j1 = elementDOF[2].val[k][k1];
				value[0] += phi0*uhstar->val[j1];
				j1 += elementDOF[2].dof;
				value[1] += phi0*uhstar->val[j1];
			}
			baryToCart2d(lambdas[2][i1], x, vertices);
			linearElas2d_u(x, val, paras);
			value[0] -= val[0];
			value[1] -= val[1];
			errors[8] += s*weight[2][i1] * lpnormp_array(2, value, 2);
		}
	}

	double elen, C11;
	int element[2], edgeNode[2];
	int count;
	int patchnodes[100];
	int edge;
	// error of lambdah for energy norm
	for (edge = 0; edge<edges->row; edge++)
	{
		edgeNode[0] = edges->val[edge][0];
		edgeNode[1] = edges->val[edge][1];
		element[0] = edges->val[edge][2];
		element[1] = edges->val[edge][3];
		elen = edges->length[edge];
		C11 = elen;

		count = 0;
		if (elementDOF[1].dop == 0)
		{
			patchnodes[count] = elementDOF[1].val[element[0]][0];
			count++;
			if (element[1] != -1)
			{
				patchnodes[count] = elementDOF[1].val[element[1]][0];
				count++;
			}
		} // if(elementDOF[1].dop==0)
		else
		{
			for (i = 0; i<3; i++)
			{
				if (elementEdge->val[element[0]][i] == edge)
					break;
			}
			patchnodes[count] = elementDOF[1].val[element[0]][(i + 1) % 3];
			count++;
			patchnodes[count] = elementDOF[1].val[element[0]][(i + 2) % 3];
			count++;
			for (j = 0; j<elementDOF[1].dop - 1; j++)
			{
				patchnodes[count] = elementDOF[1].val[element[0]][3 + i*(elementDOF[1].dop - 1) + j];
				count++;
			}

			if (element[1] != -1)
			{
				for (i = 0; i<3; i++)
				{
					if (elementEdge->val[element[1]][i] == edge)
						break;
				}
				patchnodes[count] = elementDOF[1].val[element[1]][(i + 1) % 3];
				count++;
				patchnodes[count] = elementDOF[1].val[element[1]][(i + 2) % 3];
				count++;
				for (j = 0; j<elementDOF[1].dop - 1; j++)
				{
					patchnodes[count] = elementDOF[1].val[element[1]][3 + i*(elementDOF[1].dop - 1) + j];
					count++;
				}
			}
		} // if(elementDOF[1].dop==0) else

		for (i1 = 0; i1<num_qp1; i1++)
		{
			for (i = 0; i<3; i++)
				phi2[i] = 0;

			for (k1 = 0; k1<count; k1++)
			{
				i = patchnodes[k1];
				jumpOperatorVector(lambdas1[i1][0], lambdas1[i1][1], edge, elements, elementEdge, edges, &elementDOF[1], i, phi1);

				for (j = 0; j<2; j++)
					phi2[j] += uh->val[i] * phi1[j];

				jumpOperatorVector(lambdas1[i1][0], lambdas1[i1][1], edge, elements, elementEdge, edges, &elementDOF[1], i + elementDOF[1].dof, phi1);
				for (j = 0; j<2; j++)
					phi2[j] += uh->val[i + elementDOF[1].dof] * phi1[j];
			}

			errors[2] += elen*weight1[i1] * C11* lpnormp_array(2, phi2, 2);
		}
	} // e

	errors[4] = errors[0] + errors[1] + errors[2];

	for (i = 0; i<10; i++)
		errors[i] = sqrt(errors[i]);
}


/**
* \fn void getposteriorierrorslinearElasHuZhang2d(double *errors, dvector *sigmah, dvector *uh, dvector *uhstar, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
* \brief compute error between numerical solution and exact solution: L2 norm, Energy norm
* \param *errors pointer to error between numerical solution and exact solution: L2 norm, H1 norm, Energy norm
* \param *sigmah pointer to numerical solution
* \param *uh pointer to numerical solution
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param nu Poisson ratio of plate
*/
void getposteriorierrorslinearElasHuZhang2d(double *errors, dvector *sigmah, dvector *uh, dvector *uhstar, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
{
	double paras[2];
	paras[0] = lambda;
	paras[1] = mu;
	dvector Qhf;
	create_dvector(uh->row, &Qhf);
	projL2PkVec2d(&Qhf, elements, nodes, &elementDOF[1], linearElas2d_f, paras);

	int Nt = elements->row;

	int i, j, k;

	for (i = 0; i<7; i++)
		errors[i] = 0;

	double phi0, phi1[3], phi2[3], val[4];
	int k1, i1, j1, l1, l2;
	double x[2], s, **vertices, **gradLambda, **nv, **tv;
	double value[4];

	int num_qp[3], num_qp1[2];
	double lambdas[3][100][3], weight[3][100];
	double lambdas1[2][100][2], weight1[2][100];

	num_qp[0]=getNumQuadPoints((elementDOF[0].dop - 2) * 2, 2);
	num_qp[1]=getNumQuadPoints(elementDOF[0].dop * 2, 2);
	num_qp[2]=49;
	for(i=0;i<3;i++)
		init_Gauss2d(num_qp[i], lambdas[i], weight[i]);

	num_qp1[0]=getNumQuadPoints(elementDOF[0].dop * 2, 1);
	if (num_qp1[0]>5) num_qp1[0] = 5;
	num_qp1[1]=getNumQuadPoints((elementDOF[0].dop - 1) * 2, 1);
	if (num_qp1[1]>5) num_qp1[1] = 5;
	init_Gauss1d(num_qp1[0], lambdas1[0], weight1[0]);
	init_Gauss1d(num_qp1[1], lambdas1[1], weight1[1]);

	for (k = 0; k<Nt; k++)
	{
		// set parameters
		s = elements->vol[k];
		vertices = elements->vertices[k];
		gradLambda = elements->gradLambda[k];
		nv = elements->nvector[k];
		tv = elements->tvector[k];
		// end set parameters

		// rotrot(Asigmah)
		for (i1 = 0; i1<num_qp[0]; i1++)
		{
			value[0] = 0;
			value[1] = 0;
			if (elementDOF[0].dop>1)
			{
				for (k1 = 0; k1<elementDOF[0].col; k1++)
				{
					j1 = elementDOF[0].val[k][k1];
					huzhang_basisROTROT(lambdas[0][i1], gradLambda, nv, tv, k1, elementDOF[0].dop, &phi0);
					value[0] += phi0 * sigmah->val[j1];
					huzhang_basisLaplaceTrace(lambdas[0][i1], gradLambda, nv, tv, k1, elementDOF[0].dop, &phi0);
					value[1] += phi0 * sigmah->val[j1];
				}
			}
			value[2] = (value[0] - value[1] * lambda / (2 * lambda + 2 * mu)) / (2 * mu);
			value[2] = value[2] * value[2];
			errors[0] += s*weight[0][i1] * value[2] * s*s;
		}

		// Asigmah - \varepsilon_h(uhstar)
		for (i1 = 0; i1<num_qp[1]; i1++)
		{
			value[0] = 0;
			value[1] = 0;
			value[2] = 0;
			// sigmah
			for (k1 = 0; k1<elementDOF[0].col; k1++)
			{
				huzhang_basis(lambdas[1][i1], nv, tv, k1, elementDOF[0].dop, phi1);
				j1 = elementDOF[0].val[k][k1];
				value[0] += phi1[0] * sigmah->val[j1];
				value[1] += phi1[1] * sigmah->val[j1];
				value[2] += phi1[2] * sigmah->val[j1];
			}
			// Asigmah
			value[3] = value[0] + value[1];
			value[0] = (value[0] - value[3] * lambda / (2 * lambda + 2 * mu)) / (2 * mu);
			value[1] = (value[1] - value[3] * lambda / (2 * lambda + 2 * mu)) / (2 * mu);
			value[2] = value[2] / (2 * mu);
			// Asigmah - \varepsilon_h(uhstar)
			for (k1 = 0; k1<elementDOF[2].col; k1++)
			{
				lagrange_basis1(lambdas[1][i1], gradLambda, k1, elementDOF[2].dop, phi1);
				j1 = elementDOF[2].val[k][k1];
				value[0] -= phi1[0] * uhstar->val[j1];
				value[2] -= phi1[1] * uhstar->val[j1] / 2;
				j1 += elementDOF[2].dof;
				value[1] -= phi1[1] * uhstar->val[j1];
				value[2] -= phi1[0] * uhstar->val[j1] / 2;
			}
			value[0] = value[0] * value[0];
			value[1] = value[1] * value[1];
			value[2] = value[2] * value[2];
			errors[1] += s*weight[1][i1] * (value[0] + value[1] + 2 * value[2]);
		}

		// oscillation: f-Qhf
		for (i1 = 0; i1<num_qp[2]; i1++)
		{
			value[0] = 0;
			value[1] = 0;
			for (k1 = 0; k1<elementDOF[1].col; k1++)
			{
				lagrange_basis(lambdas[2][i1], k1, elementDOF[1].dop, &phi0);
				j1 = elementDOF[1].val[k][k1];
				value[0] += phi0*Qhf.val[j1];
				j1 += elementDOF[1].dof;
				value[1] += phi0*Qhf.val[j1];
			}
			baryToCart2d(lambdas[2][i1], x, vertices);
			linearElas2d_f(x, val, paras);
			value[0] -= val[0];
			value[1] -= val[1];
			errors[2] += s*weight[2][i1] * lpnormp_array(2, value, 2)*s;
		}
	} // k

	int element[2];
	int count, istart;
	int *index;
	index = (int*)calloc(elementDOF[0].dof, sizeof(int));
	for (i = 0; i<elementDOF[0].dof; i++)
		index[i] = -1;

	int patchnodes[200];
	int edge;
	double elen;

	for (edge = 0; edge<edges->row; edge++)
	{
		element[0] = edges->val[edge][2];
		element[1] = edges->val[edge][3];
		elen = edges->length[edge];

		istart = -2;
		count = 0;
		for (i = 0; i<elementDOF[0].col; i++)
		{
			j = elementDOF[0].val[element[0]][i];
			patchnodes[count] = j;
			count++;
			index[j] = istart;
			istart = j;
		}

		if (element[1] != -1)
		{
			for (i = 0; i<elementDOF[0].col; i++)
			{
				j = elementDOF[0].val[element[1]][i];
				if (index[j] == -1)
				{
					patchnodes[count] = j;
					count++;
					index[j] = istart;
					istart = j;

				}
			}
		}

		for (j = 0; j<count; j++)
		{
			j1 = istart;
			istart = index[j1];
			index[j1] = -1;
		}

		// M_{tt}(Asigmah)
		for (i1 = 0; i1<num_qp1[0]; i1++)
		{
			value[0] = 0;

			for (k1 = 0; k1<count; k1++)
			{
				i = patchnodes[k1];
				jumpOperatorATensorTangent2(lambdas1[0][i1][0], lambdas1[0][i1][1], edge, elements, elementEdge, edges, &elementDOF[0], i, lambda, mu, &phi0);
				value[0] += sigmah->val[i] * phi0;
			}

			errors[3] += elen*weight1[0][i1] * value[0] * value[0] * elen;
		}

		// rot(Asigmah) t - \partial_t(M_{nt}(Asigmah))
		for (i1 = 0; i1<num_qp1[1]; i1++)
		{
			value[0] = 0;

			for (k1 = 0; k1<count; k1++)
			{
				i = patchnodes[k1];
				jumpOperatorRotATensorTangentPt(lambdas1[1][i1][0], lambdas1[1][i1][1], edge, elements, elementEdge, edges, &elementDOF[0], i, lambda, mu, &phi0);
				value[0] += sigmah->val[i] * phi0;
			}
			
			errors[4] += elen*weight1[1][i1] * value[0] * value[0] * elen * elen * elen;
		}
	} // edge
	free(index);

	free_dvector(&Qhf);

	errors[5] = errors[0] + errors[1] + errors[3] + errors[4];
	errors[6] = errors[5] + errors[2];

	for (i = 0; i<7; i++)
		errors[i] = sqrt(errors[i]);
}

/**
 * \fn void geterrorslinearElasHuangZhou2d(double *errors, dvector *sigmah, dvector *uh, dvector *Qhu, dvector *uhstar, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
 * \brief compute error between numerical solution and exact solution: L2 norm, Energy norm
 * \param *errors pointer to error between numerical solution and exact solution: L2 norm, H1 norm, Energy norm
 * \param *sigmah pointer to numerical solution
 * \param *uh pointer to numerical solution
 * \param *Qhu pointer to the L2 projection of u to the piecewise Lagrangian element space
 * \param *uhstar pointer to postprocessed displacement
* \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
 the fourth column stores -1 if the edge is on boundary
 * \param *nodes pointer to the nodes location of the triangulation
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param lambda Lame constant
 * \param mu Lame constant or Poisson ratio of plate
 */
void geterrorslinearElasHuangZhou2d(double *errors, dvector *sigmah, dvector *uh, dvector *Qhu, dvector *uhstar, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
{
	int Nt = elements->row;

	int i, j, k, l;

	for (i = 0; i < 10; i++)
		errors[i] = 0;

	double phi0, phi1[3], phi2[3], val[4], paras[2];
	int k1, i1, j1, l1, l2;
	double x[2], s, **vertices, **gradLambda, **nv, **tv, **tij;
	double value[4];

	paras[0] = lambda;
	paras[1] = mu;

	int num_qp[5], num_qp1;
	double lambdas[5][100][3], weight[5][100];
	double lambdas1[100][2], weight1[100];

	num_qp[0]=49;
	num_qp[1]=49;
	num_qp[2]=49;
	num_qp[3]=getNumQuadPoints(elementDOF[1].dop * 2 - 2, 2);
	num_qp[4]=getNumQuadPoints(elementDOF[1].dop * 2, 2);
	for(i=0;i<5;i++)
		init_Gauss2d(num_qp[i], lambdas[i], weight[i]);

	num_qp1=getNumQuadPoints(elementDOF[1].dop * 2, 1);
	if (num_qp1>5) num_qp1 = 5;
	init_Gauss1d(num_qp1, lambdas1, weight1);

	for (k = 0; k<Nt; k++)
	{
		// set parameters
		s = elements->vol[k];
		vertices = elements->vertices[k];
		gradLambda = elements->gradLambda[k];
		tij = elements->tij[k];
		nv = elements->nvector[k];
		tv = elements->tvector[k];
		// end set parameters

		// Energy norm and L2 norm of sigma-sigma_h
		for (i1 = 0; i1<num_qp[0]; i1++)
		{
			value[0] = 0;
			value[1] = 0;
			value[2] = 0;
			for (k1 = 0; k1<elementDOF[0].col; k1++)
			{
				// huzhang_basis(lambdas[0][i1], elements->nvector[k], elements->tvector[k], k1, elementDOF[0].dop, phi1);
				divS_huangzhou_basis(lambdas[0][i1], s, nv, tv, tij, k1, phi1);
				j1 = elementDOF[0].val[k][k1];
				value[0] += phi1[0] * sigmah->val[j1];
				value[1] += phi1[1] * sigmah->val[j1];
				value[2] += phi1[2] * sigmah->val[j1];
			}
			baryToCart2d(lambdas[0][i1], x, vertices);
			linearElas2d_sigma(x, val, paras);
			value[0] -= val[0];
			value[1] -= val[1];
			value[2] -= val[2];
			value[3] = value[0] + value[1];
			value[0] = value[0] * value[0];
			value[1] = value[1] * value[1];
			value[2] = value[2] * value[2];
			value[3] = value[3] * value[3];
			if (lambda>-0.5)
				errors[0] += s*weight[0][i1] * (value[0] + value[1] + 2 * value[2] - value[3] * lambda / (2 * lambda + 2 * mu)) / (2 * mu);
			else
				errors[0] += s*weight[0][i1] * (value[0] + value[1] + 2 * value[2] - value[3] * 1.0 / 2.0) / (2 * mu);

			errors[9] += s*weight[0][i1] * (value[0] + value[1] + 2 * value[2]);
		}

		// L2 norm of divergence of sigma-sigma_h
		for (i1 = 0; i1<num_qp[1]; i1++)
		{
			value[0] = 0;
			value[1] = 0;
			for (k1 = 0; k1<elementDOF[0].col; k1++)
			{
				divS_huangzhou_basisDIV(lambdas[1][i1], gradLambda, s, nv, tv, tij, k1, phi1);
				j1 = elementDOF[0].val[k][k1];
				value[0] += phi1[0] * sigmah->val[j1];
				value[1] += phi1[1] * sigmah->val[j1];
			}
			baryToCart2d(lambdas[1][i1], x, vertices);
			linearElas2d_f(x, val, paras);
			value[0] += val[0];
			value[1] += val[1];
			errors[1] += s*weight[1][i1] * lpnormp_array(2, value, 2);
		}

		// L2 norm of u-u_h
		for (i1 = 0; i1<num_qp[2]; i1++)
		{
			value[0] = 0;
			value[1] = 0;
			for (k1 = 0; k1<elementDOF[1].col; k1++)
			{
				lagrange_basis(lambdas[2][i1], k1, elementDOF[1].dop, &phi0);
				j1 = elementDOF[1].val[k][k1];
				value[0] += phi0*uh->val[j1];
				j1 += elementDOF[1].dof;
				value[1] += phi0*uh->val[j1];
			}
			baryToCart2d(lambdas[2][i1], x, vertices);
			linearElas2d_u(x, val, paras);
			value[0] -= val[0];
			value[1] -= val[1];
			errors[3] += s*weight[2][i1] * lpnormp_array(2, value, 2);
		}

		// L2 norm of symgrad(Qhu-u_h)
		for (i1 = 0; i1<num_qp[3]; i1++)
		{
			value[0] = 0;
			value[1] = 0;
			value[2] = 0;
			value[3] = 0;
			for (k1 = 0; k1<elementDOF[1].col; k1++)
			{
				lagrange_basis1(lambdas[3][i1], gradLambda, k1, elementDOF[1].dop, phi1);
				j1 = elementDOF[1].val[k][k1];
				value[0] += phi1[0] * (Qhu->val[j1] - uh->val[j1]);
				value[1] += phi1[1] * (Qhu->val[j1] - uh->val[j1]);
				j1 += elementDOF[1].dof;
				value[2] += phi1[0] * (Qhu->val[j1] - uh->val[j1]);
				value[3] += phi1[1] * (Qhu->val[j1] - uh->val[j1]);
			}
			value[1] = (value[1]+value[2])/2;
			value[2] = value[1];
			errors[5] += s*weight[3][i1] * lpnormp_array(4, value, 2);
		}

		// L2 norm of Qhu-u_h
		for (i1 = 0; i1<num_qp[4]; i1++)
		{
			value[0] = 0;
			value[1] = 0;
			for (k1 = 0; k1<elementDOF[1].col; k1++)
			{
				lagrange_basis(lambdas[4][i1], k1, elementDOF[1].dop, &phi0);
				j1 = elementDOF[1].val[k][k1];
				value[0] += phi0*(Qhu->val[j1]- uh->val[j1]);
				j1 += elementDOF[1].dof;
				value[1] += phi0*(Qhu->val[j1] - uh->val[j1]);
			}
			errors[6] += s*weight[4][i1] * lpnormp_array(2, value, 2);
		}

		// H1 semi-norm of u-uhstar
		// for (i1 = 0; i1<num_qp[2]; i1++)
		// {
		// 	value[0] = 0;
		// 	value[1] = 0;
		// 	value[2] = 0;
		// 	value[3] = 0;
		// 	for (k1 = 0; k1<elementDOF[2].col; k1++)
		// 	{
		// 		lagrange_basis1(lambdas[2][i1], gradLambda, k1, elementDOF[2].dop, phi1);
		// 		j1 = elementDOF[2].val[k][k1];
		// 		value[0] += phi1[0] * uhstar->val[j1];
		// 		value[1] += phi1[1] * uhstar->val[j1];
		// 		j1 += elementDOF[2].dof;
		// 		value[2] += phi1[0] * uhstar->val[j1];
		// 		value[3] += phi1[1] * uhstar->val[j1];
		// 	}
		// 	baryToCart2d(lambdas[2][i1], x, vertices);
		// 	linearElas2d_gradu(x, val, paras);
		// 	value[0] -= val[0];
		// 	value[1] -= val[1];
		// 	value[2] -= val[2];
		// 	value[3] -= val[3];
		// 	errors[7] += s*weight[2][i1] * lpnormp_array(4, value, 2);
		// }

		// L2 norm of u-uhstar
		// for (i1 = 0; i1<num_qp[2]; i1++)
		// {
		// 	value[0] = 0;
		// 	value[1] = 0;
		// 	for (k1 = 0; k1<elementDOF[2].col; k1++)
		// 	{
		// 		lagrange_basis(lambdas[2][i1], k1, elementDOF[2].dop, &phi0);
		// 		j1 = elementDOF[2].val[k][k1];
		// 		value[0] += phi0*uhstar->val[j1];
		// 		j1 += elementDOF[2].dof;
		// 		value[1] += phi0*uhstar->val[j1];
		// 	}
		// 	baryToCart2d(lambdas[2][i1], x, vertices);
		// 	linearElas2d_u(x, val, paras);
		// 	value[0] -= val[0];
		// 	value[1] -= val[1];
		// 	errors[8] += s*weight[2][i1] * lpnormp_array(2, value, 2);
		// }
	}

	double elen, C11;
	int element[2], edgeNode[2];
	int count;
	int patchnodes[100];
	int edge;
	// error of lambdah for energy norm
	for (edge = 0; edge<edges->row; edge++)
	{
		edgeNode[0] = edges->val[edge][0];
		edgeNode[1] = edges->val[edge][1];
		element[0] = edges->val[edge][2];
		element[1] = edges->val[edge][3];
		elen = edges->length[edge];
		C11 = 1./elen;

		count = 0;
		if (elementDOF[1].dop == 0)
		{
			patchnodes[count] = elementDOF[1].val[element[0]][0];
			count++;
			if (element[1] != -1)
			{
				patchnodes[count] = elementDOF[1].val[element[1]][0];
				count++;
			}
		} // if(elementDOF[1].dop==0)
		else
		{
			for (i = 0; i<3; i++)
			{
				if (elementEdge->val[element[0]][i] == edge)
					break;
			}
			patchnodes[count] = elementDOF[1].val[element[0]][(i + 1) % 3];
			count++;
			patchnodes[count] = elementDOF[1].val[element[0]][(i + 2) % 3];
			count++;
			for (j = 0; j<elementDOF[1].dop - 1; j++)
			{
				patchnodes[count] = elementDOF[1].val[element[0]][3 + i*(elementDOF[1].dop - 1) + j];
				count++;
			}

			if (element[1] != -1)
			{
				for (i = 0; i<3; i++)
				{
					if (elementEdge->val[element[1]][i] == edge)
						break;
				}
				patchnodes[count] = elementDOF[1].val[element[1]][(i + 1) % 3];
				count++;
				patchnodes[count] = elementDOF[1].val[element[1]][(i + 2) % 3];
				count++;
				for (j = 0; j<elementDOF[1].dop - 1; j++)
				{
					patchnodes[count] = elementDOF[1].val[element[1]][3 + i*(elementDOF[1].dop - 1) + j];
					count++;
				}
			}
		} // if(elementDOF[1].dop==0) else

		for (i1 = 0; i1<num_qp1; i1++)
		{
			for (i = 0; i<3; i++)
				phi2[i] = 0;

			for (k1 = 0; k1<count; k1++)
			{
				i = patchnodes[k1];
				jumpOperatorVector(lambdas1[i1][0], lambdas1[i1][1], edge, elements, elementEdge, edges, &elementDOF[1], i, phi1);

				for (j = 0; j<2; j++)
					phi2[j] += (Qhu->val[i]-uh->val[i]) * phi1[j];

				jumpOperatorVector(lambdas1[i1][0], lambdas1[i1][1], edge, elements, elementEdge, edges, &elementDOF[1], i + elementDOF[1].dof, phi1);
				for (j = 0; j<2; j++)
					phi2[j] += (Qhu->val[i + elementDOF[1].dof]-uh->val[i + elementDOF[1].dof]) * phi1[j];
			}

			errors[2] += elen*weight1[i1] * C11* lpnormp_array(2, phi2, 2);
		}
	} // e

	// H1 semi-norm of Qhu-u_h
	errors[4] = errors[5] + errors[2];

	for (i = 0; i<10; i++)
		errors[i] = sqrt(errors[i]);
}

/**
 * \fn void geterrorslinearElasHuangZhou2d(double *errors, dvector *sigmah, dvector *uh, dvector *Qhu, dvector *uhstar, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
 * \brief compute error between numerical solution and exact solution: L2 norm, Energy norm
 * \param *errors pointer to error between numerical solution and exact solution: L2 norm, H1 norm, Energy norm
 * \param *sigmah pointer to numerical solution
 * \param *uh pointer to numerical solution
 * \param *Qhu pointer to the L2 projection of u to the piecewise Lagrangian element space
 * \param *uhstar pointer to postprocessed displacement
* \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
 the fourth column stores -1 if the edge is on boundary
 * \param *nodes pointer to the nodes location of the triangulation
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param lambda Lame constant
 * \param mu Lame constant or Poisson ratio of plate
 */
void geterrorslinearElasReducedHuangZhou2d(double *errors, dvector *sigmah, dvector *uh, dvector *Qhu, dvector *uhstar, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
{
	int Nt = elements->row;

	int i, j, k, l;

	for (i = 0; i < 10; i++)
		errors[i] = 0;

	double phi0, phi1[4], phi2[4], val[4], paras[2];
	int k1, i1, j1, l1, l2;
	double x[2], s, **vertices, **gradLambda, **nv, **tv, **tij, *edgeslength;
	double **lbc;
	short *eorien;
	double value[4];

	paras[0] = lambda;
	paras[1] = mu;

	int num_qp[5], num_qp1;
	double lambdas[5][100][3], weight[5][100];
	double lambdas1[100][2], weight1[100];

	num_qp[0]=49;
	num_qp[1]=49;
	num_qp[2]=49;
	num_qp[3]=getNumQuadPoints(elementDOF[1].dop * 2 - 2, 2);
	num_qp[4]=getNumQuadPoints(elementDOF[1].dop * 2, 2);
	for(i=0;i<5;i++)
		init_Gauss2d(num_qp[i], lambdas[i], weight[i]);

	num_qp1=getNumQuadPoints(elementDOF[1].dop * 2, 1);
	if (num_qp1>5) num_qp1 = 5;
	init_Gauss1d(num_qp1, lambdas1, weight1);

	for (k = 0; k<Nt; k++){
		// set parameters
		s = elements->vol[k];
		vertices = elements->vertices[k];
		gradLambda = elements->gradLambda[k];
		tij = elements->tij[k];
		nv = elements->nvector[k];
		tv = elements->tvector[k];
		edgeslength = elements->edgeslength[k];
		eorien = elements->eorien[k];
		lbc = basisCoeffs->val[k];
		// end set parameters

		// Energy norm and L2 norm of sigma-sigma_h
		for (i1 = 0; i1<num_qp[0]; i1++){
			value[0] = 0;
			value[1] = 0;
			value[2] = 0;
			for (k1 = 0; k1<elementDOF[0].col; k1++){
				divS_reducedhuangzhou_basis(lambdas[0][i1], lbc, s, nv, tv, tij, k1, phi1);
				j1 = elementDOF[0].val[k][k1];
				value[0] += phi1[0] * sigmah->val[j1];
				value[1] += phi1[1] * sigmah->val[j1];
				value[2] += phi1[2] * sigmah->val[j1];
			}
			baryToCart2d(lambdas[0][i1], x, vertices);
			linearElas2d_sigma(x, val, paras);
			value[0] -= val[0];
			value[1] -= val[1];
			value[2] -= val[2];
			value[3] = value[0] + value[1];
			value[0] = value[0] * value[0];
			value[1] = value[1] * value[1];
			value[2] = value[2] * value[2];
			value[3] = value[3] * value[3];
			if (lambda>-0.5)
				errors[0] += s*weight[0][i1] * (value[0] + value[1] + 2 * value[2] - value[3] * lambda / (2 * lambda + 2 * mu)) / (2 * mu);
			else
				errors[0] += s*weight[0][i1] * (value[0] + value[1] + 2 * value[2] - value[3] * 1.0 / 2.0) / (2 * mu);

			errors[9] += s*weight[0][i1] * (value[0] + value[1] + 2 * value[2]);
		}

		// L2 norm of divergence of sigma-sigma_h
		for (i1 = 0; i1<num_qp[1]; i1++){
			value[0] = 0;
			value[1] = 0;
			for (k1 = 0; k1<elementDOF[0].col; k1++)
			{
				divS_reducedhuangzhou_basisDIV(lambdas[1][i1], gradLambda, lbc, s, nv, tv, tij, k1, phi1);
				j1 = elementDOF[0].val[k][k1];
				value[0] += phi1[0] * sigmah->val[j1];
				value[1] += phi1[1] * sigmah->val[j1];
			}
			baryToCart2d(lambdas[1][i1], x, vertices);
			linearElas2d_f(x, val, paras);
			value[0] += val[0];
			value[1] += val[1];
			errors[1] += s*weight[1][i1] * lpnormp_array(2, value, 2);
		}

		// L2 norm of u-u_h
		for (i1 = 0; i1<num_qp[2]; i1++){
			value[0] = 0;
			value[1] = 0;
			for (k1 = 0; k1<elementDOF[1].col; k1++){
				nedelec1st_basis(lambdas[2][i1], gradLambda, edgeslength, eorien, k1, 1, phi1);
				j1 = elementDOF[1].val[k][k1];
				value[0] += phi1[0] * uh->val[j1];
				value[1] += phi1[1] * uh->val[j1];
			}
			baryToCart2d(lambdas[2][i1], x, vertices);
			linearElas2d_u(x, val, paras);
			value[0] -= val[0];
			value[1] -= val[1];
			errors[3] += s*weight[2][i1] * lpnormp_array(2, value, 2);
		}

		// L2 norm of symgrad(Qhu-u_h)
		for (i1 = 0; i1<num_qp[3]; i1++){
			value[0] = 0;
			value[1] = 0;
			value[2] = 0;
			value[3] = 0;
			for (k1 = 0; k1<elementDOF[1].col; k1++){
				nedelec1st_basis1(lambdas[3][i1], gradLambda, edgeslength, eorien, k1, 1, phi1);
				j1 = elementDOF[1].val[k][k1];
				value[0] += phi1[0] * (Qhu->val[j1]- uh->val[j1]);
				value[1] += phi1[1] * (Qhu->val[j1]- uh->val[j1]);
				value[2] += phi1[2] * (Qhu->val[j1]- uh->val[j1]);
				value[3] += phi1[3] * (Qhu->val[j1]- uh->val[j1]);
			}
			value[1] = (value[1]+value[2])/2;
			value[2] = value[1];
			errors[5] += s*weight[3][i1] * lpnormp_array(4, value, 2);
		}

		// L2 norm of Qhu-u_h
		for (i1 = 0; i1<num_qp[4]; i1++){
			value[0] = 0;
			value[1] = 0;
			for (k1 = 0; k1<elementDOF[1].col; k1++){
				nedelec1st_basis(lambdas[4][i1], gradLambda, edgeslength, eorien, k1, 1, phi1);
				j1 = elementDOF[1].val[k][k1];
				value[0] += phi1[0] * (Qhu->val[j1]- uh->val[j1]);
				value[1] += phi1[1] * (Qhu->val[j1]- uh->val[j1]);
			}
			errors[6] += s*weight[4][i1] * lpnormp_array(2, value, 2);
		}

		// H1 semi-norm of u-uhstar
		// for (i1 = 0; i1<num_qp[2]; i1++)
		// {
		// 	value[0] = 0;
		// 	value[1] = 0;
		// 	value[2] = 0;
		// 	value[3] = 0;
		// 	for (k1 = 0; k1<elementDOF[2].col; k1++)
		// 	{
		// 		lagrange_basis1(lambdas[2][i1], gradLambda, k1, elementDOF[2].dop, phi1);
		// 		j1 = elementDOF[2].val[k][k1];
		// 		value[0] += phi1[0] * uhstar->val[j1];
		// 		value[1] += phi1[1] * uhstar->val[j1];
		// 		j1 += elementDOF[2].dof;
		// 		value[2] += phi1[0] * uhstar->val[j1];
		// 		value[3] += phi1[1] * uhstar->val[j1];
		// 	}
		// 	baryToCart2d(lambdas[2][i1], x, vertices);
		// 	linearElas2d_gradu(x, val, paras);
		// 	value[0] -= val[0];
		// 	value[1] -= val[1];
		// 	value[2] -= val[2];
		// 	value[3] -= val[3];
		// 	errors[7] += s*weight[2][i1] * lpnormp_array(4, value, 2);
		// }

		// L2 norm of u-uhstar
		// for (i1 = 0; i1<num_qp[2]; i1++)
		// {
		// 	value[0] = 0;
		// 	value[1] = 0;
		// 	for (k1 = 0; k1<elementDOF[2].col; k1++)
		// 	{
		// 		lagrange_basis(lambdas[2][i1], k1, elementDOF[2].dop, &phi0);
		// 		j1 = elementDOF[2].val[k][k1];
		// 		value[0] += phi0*uhstar->val[j1];
		// 		j1 += elementDOF[2].dof;
		// 		value[1] += phi0*uhstar->val[j1];
		// 	}
		// 	baryToCart2d(lambdas[2][i1], x, vertices);
		// 	linearElas2d_u(x, val, paras);
		// 	value[0] -= val[0];
		// 	value[1] -= val[1];
		// 	errors[8] += s*weight[2][i1] * lpnormp_array(2, value, 2);
		// }
	}

	double elen, C11;
	int element[2], edgeNode[2];
	int count;
	int patchnodes[100];
	int edge;
	// error of lambdah for energy norm
	for (edge = 0; edge<edges->row; edge++){
		edgeNode[0] = edges->val[edge][0];
		edgeNode[1] = edges->val[edge][1];
		element[0] = edges->val[edge][2];
		element[1] = edges->val[edge][3];
		elen = edges->length[edge];
		C11 = 1./elen;

		count = 0;
		for(i=0;i<3;i++){
			patchnodes[count] = elementDOF[1].val[element[0]][i];
			count++;
		}
		if (element[1] != -1){
			for(i=0;i<3;i++){
				patchnodes[count] = elementDOF[1].val[element[1]][i];
				count++;
			}
		}

		for (i1 = 0; i1<num_qp1; i1++){
			for (i = 0; i<3; i++)
				phi2[i] = 0;

			for (k1 = 0; k1<count; k1++){
				i = patchnodes[k1];
				jumpOperatorRM(lambdas1[i1][0], lambdas1[i1][1], edge, elements, elementEdge, edges, &elementDOF[1], i, phi1);

				for (j = 0; j<2; j++)
					phi2[j] += (Qhu->val[i]-uh->val[i]) * phi1[j];
			}

			errors[2] += elen*weight1[i1] * C11* lpnormp_array(2, phi2, 2);
		}
	} // e

	// H1 semi-norm of Qhu-u_h
	errors[4] = errors[5] + errors[2];

	for (i = 0; i<10; i++)
		errors[i] = sqrt(errors[i]);
}

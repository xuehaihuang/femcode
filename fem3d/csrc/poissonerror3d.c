/*
 *  poissonerror.c
 *
 *  Created by Xuehai Huang on Apr 30, 2022.
 *  Copyright 2022 SUFE. All rights reserved.
 *
 */

/*! \file poissonerror.c
 *  \brief Assembling for stiffness matrix and solve it
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "header.h"
#include "matvec.h"
#include "checkmat.h"

/**
 * \fn void geterrorsPoissonLagrange3d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
 * \brief compute error between numerical solution and exact solution: L2 norm, Energy norm
 * \param *errors pointer to error between numerical solution and exact solution: L2 norm, H1 norm, Energy norm
 * \param *sigmah pointer to numerical solution
 * \param *uh pointer to numerical solution
 * \param *elements pointer to the structure of the triangulation
 * \param *elementFace pointer to relation between tetrahedrons and faces: each row stores 4 faces index
 * \param *faces the first three columns store the three vertices corresponding to the face; 
 *				 the 4th and 5th columns store the elements which the face belongs to;
 *				 if the face is a boundary, the 5th column will stores -1;
 *				 the first column is in ascend order.
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges stores the two vertice corresponding to the edge
 *				 the first column is in ascend order.
 * \param *nodes pointer to the nodes location of the triangulation
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param lambda Lame constant
 * \param mu Lame constant or Poisson ratio of plate
 */
void geterrorsPoissonLagrange3d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{	
	int Nt=elements->row;
	init_array(3, errors, 0);
	
	int i,j,k;
	int face, edge;
			
	double phi0, phi1[3], phi2[3], val[3];
	int k1,i1,i2,j1;
	double x[3], vol, **gradLambda, **vertices;

	int num_qp;
	double lambdas[100][4], weight[100];

	dvector luh;
	create_dvector(elementDOF[0].col, &luh);
	
	num_qp=getNumQuadPoints_ShunnWilliams(9, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
		
	for(k=0;k<Nt;k++)
	{
		// set parameters
		vol=elements->vol[k];
        gradLambda = elements->gradLambda[k];
		vertices = elements->vertices[k];
		// end set parameters
		
		for (i = 0; i<elementDOF[0].col; i++)
		{
			j = elementDOF[0].val[k][i];
			luh.val[i] = uh->val[j];
		}

		// L2 norm of u-u_h
		for(i1=0;i1<num_qp;i1++)
		{
            axy_array(3, lambdas[i1][3], vertices[3], x);
			for(i=0;i<3;i++)
				axpy_array(3, lambdas[i1][i], vertices[i], x);
			val[0] = -poisson3d_u(x);

			for(k1=0;k1<elementDOF->col;k1++)
			{
				lagrange3d_basis(lambdas[i1], k1, elementDOF->dop, &phi0);
				val[0]+=phi0*luh.val[k1];
			}
			errors[0]+=vol*weight[i1]*val[0]*val[0];
		}

		// H1 semi-norm of u-u_h
		for(i1=0;i1<num_qp;i1++)
		{
            axy_array(3, lambdas[i1][3], vertices[3], x);
			for(i=0;i<3;i++)
				axpy_array(3, lambdas[i1][i], vertices[i], x);
            poisson3d_gradu(x, val);

			for(k1=0;k1<elementDOF->col;k1++)
			{
				lagrange3d_basis1(lambdas[i1], gradLambda, k1, elementDOF->dop, phi1);
				axpy_array(3, -luh.val[k1], phi1, val);
			}
			errors[1]+=vol*weight[i1]*dot_array(3, val, val);
		}
	}
	free_dvector(&luh);

	errors[2]=errors[0]+errors[1];
	
	for(i=0;i<3;i++)
		errors[i]=sqrt(errors[i]);
}
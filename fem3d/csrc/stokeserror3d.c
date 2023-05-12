/*
 *  stokeserror.c
 *
 *  Created by Xuehai Huang on May 8, 2023.
 *  Copyright 2023 SUFE. All rights reserved.
 *
 */

/*! \file stokeserror.c
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
 * \fn void geterrorsStokesNcP1P03d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
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
void geterrorsStokesNcP1P03d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{	
	int Nt=elements->row;
	init_array(4, errors, 0);
	
	int i,j,k;
	int face, edge;
			
	double phi0, phi1[3], phi2[3], val[9];
	int k1,i1,i2,j1;
	double x[3], vol, **gradLambda, **vertices;

	int num_qp;
	double lambdas[100][4], weight[100];

	dvector luh;
	create_dvector(elementDOF[0].col*3, &luh);
	
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
			luh.val[i] = uh[0].val[j];
			luh.val[i+elementDOF[0].col] = uh[0].val[j+elementDOF[0].dof];
			luh.val[i+elementDOF[0].col*2] = uh[0].val[j+elementDOF[0].dof*2];
		}

		// L2 norm of u-u_h
		for(i1=0;i1<num_qp;i1++)
		{
            axy_array(3, lambdas[i1][3], vertices[3], x);
			for(i=0;i<3;i++)
				axpy_array(3, lambdas[i1][i], vertices[i], x);
			stokes3d_u(x, val);
			for(i=0;i<3;i++) val[i]*=-1.0;

			for(k1=0;k1<elementDOF[0].col;k1++)
			{
				ncp13d_basis(lambdas[i1], k1, &phi0);
				val[0]+=phi0*luh.val[k1];
				val[1]+=phi0*luh.val[k1+elementDOF[0].col];
				val[2]+=phi0*luh.val[k1+elementDOF[0].col*2];
			}
			errors[0]+=vol*weight[i1]*dot_array(3, val, val);
		}

		// H1 semi-norm of u-u_h
		for(i1=0;i1<num_qp;i1++)
		{
            axy_array(3, lambdas[i1][3], vertices[3], x);
			for(i=0;i<3;i++)
				axpy_array(3, lambdas[i1][i], vertices[i], x);
			stokes3d_gradu(x, val);
			for(i=0;i<9;i++) val[i]*=-1.0;

			for(k1=0;k1<elementDOF[0].col;k1++)
			{
				ncp13d_basis1(gradLambda, k1, phi1);
				axpy_array(3, luh.val[k1], phi1, val);
				axpy_array(3, luh.val[k1+elementDOF[0].col], phi1, val+3);
				axpy_array(3, luh.val[k1+elementDOF[0].col*2], phi1, val+6);
			}
			errors[1]+=vol*weight[i1]*dot_array(9, val, val);
		}

		// L2 norm of p-p_h
		for(i1=0;i1<num_qp;i1++)
		{
            axy_array(3, lambdas[i1][3], vertices[3], x);
			for(i=0;i<3;i++)
				axpy_array(3, lambdas[i1][i], vertices[i], x);
			val[0]=-stokes3d_p(x);
			
			j = elementDOF[1].val[k][0];
			val[0] += uh[1].val[j];
			errors[3]+=vol*weight[i1]*val[0]*val[0];
		}
	}
	free_dvector(&luh);

	errors[2]=errors[0]+errors[1];
	
	for(i=0;i<4;i++)
		errors[i]=sqrt(errors[i]);
}
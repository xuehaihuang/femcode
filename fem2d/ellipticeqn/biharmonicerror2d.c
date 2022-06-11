/*
 *  biharmonicerror2d.c
 *
 *  Created by Xuehai Huang on Apr 30, 2022.
 *  Copyright 2022 SUFE. All rights reserved.
 *
 */

/*! \file biharmonicerror2d.c
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
 * \fn void geterrorsBiharmonicMorley2d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
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
void geterrorsBiharmonicMorley2d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{	
	int Nt=elements->row;
	init_array(4, errors, 0);
	
	int i,j,k;
	int face, edge;
			
	double phi0, phi1[3], phi2[3], val[3];
	int k1,i1,i2,j1;
	double x[2], s, **gradLambda, *nve[3], **vertices;

	int num_qp;
	double lambdas[100][3], weight[100];

	dvector luh;
	create_dvector(elementDOF[0].col, &luh);
	
	num_qp=49; 
	init_Gauss2d(num_qp, lambdas, weight); 
		
	for(k=0;k<Nt;k++)
	{
		// set parameters
		vertices = elements->vertices[k];
		s=elements->vol[k];
        gradLambda = elements->gradLambda[k];
		for(i=0;i<3;i++){
			j = elementEdge->val[k][i];
			nve[i] = edges->nvector[j];
		}
		// end set parameters
		
		for (i = 0; i<elementDOF[0].col; i++)
		{
			j = elementDOF[0].val[k][i];
			luh.val[i] = uh->val[j];
		}

		// L2 norm of u-u_h
		for(i1=0;i1<num_qp;i1++)
		{
            axpbyz_array(2, lambdas[i1][0], vertices[0], lambdas[i1][1], vertices[1], x);
			axpy_array(2, lambdas[i1][2], vertices[2], x);
			val[0] = -biharmonic2d_u(x, NULL);

			for(k1=0;k1<elementDOF->col;k1++)
			{
				morley_basis(lambdas[i1], gradLambda, nve, k1, &phi0);
				val[0]+=phi0*luh.val[k1];
			}
			errors[0]+=s*weight[i1]*val[0]*val[0];
		}

		// H1 semi-norm of u-u_h
		for(i1=0;i1<num_qp;i1++)
		{
			axpbyz_array(2, lambdas[i1][0], vertices[0], lambdas[i1][1], vertices[1], x);
			axpy_array(2, lambdas[i1][2], vertices[2], x);
            biharmonic2d_gradu(x, val, NULL);

			for(k1=0;k1<elementDOF->col;k1++)
			{
				morley_basis1(lambdas[i1], gradLambda, nve, k1, phi1);
				axpy_array(2, -luh.val[k1], phi1, val);
			}
			errors[1]+=s*weight[i1]*dot_array(2, val, val);
		}

		// H2 semi-norm of u-u_h
		for(i1=0;i1<num_qp;i1++)
		{
			axpbyz_array(2, lambdas[i1][0], vertices[0], lambdas[i1][1], vertices[1], x);
			axpy_array(2, lambdas[i1][2], vertices[2], x);
            biharmonic2d_hessu(x, val, NULL);

			for(k1=0;k1<elementDOF->col;k1++)
			{
				morley_basis2(gradLambda, nve, k1, phi2);
				axpy_array(3, -luh.val[k1], phi2, val);
			}
			errors[2]+=s*weight[i1]*(dot_array(2, val, val) + val[2]*val[2]);
		}
	}
	free_dvector(&luh);

	errors[3]=errors[0]+errors[1]+errors[2];
	
	for(i=0;i<4;i++)
		errors[i]=sqrt(errors[i]);
}
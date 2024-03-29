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
 * \fn void geterrorsPoissonLagrange2d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
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
void geterrorsPoissonLagrange2d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{	
	int Nt=elements->row;
	init_array(3, errors, 0);
	
	int i,j,k;
	int face, edge;
			
	double phi0, phi1[3], phi2[3], val[3];
	int k1,i1,i2,j1;
	double x[2], s, **gradLambda, **vertices;

	int num_qp;
	double lambdas[100][3], weight[100];

	dvector luh;
	create_dvector(elementDOF[0].col, &luh);
	
	num_qp=49; 
	init_Gauss2d(num_qp, lambdas, weight); 
		
	for(k=0;k<Nt;k++)
	{
		// set parameters
		s=elements->vol[k];
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
			baryToCart2d(lambdas[i1], x, vertices);
			// val[0] = -poisson2d_u(x, NULL);
			poisson2d_u(x, &val[0], NULL);
			val[0] *= -1;

			for(k1=0;k1<elementDOF->col;k1++)
			{
				lagrange_basis(lambdas[i1], k1, elementDOF->dop, &phi0);
				val[0]+=phi0*luh.val[k1];
			}
			errors[0]+=s*weight[i1]*val[0]*val[0];
		}

		// H1 semi-norm of u-u_h
		for(i1=0;i1<num_qp;i1++)
		{
			baryToCart2d(lambdas[i1], x, vertices);
            poisson2d_gradu(x, val, NULL);

			for(k1=0;k1<elementDOF->col;k1++)
			{
				lagrange_basis1(lambdas[i1], gradLambda, k1, elementDOF->dop, phi1);
				axpy_array(2, -luh.val[k1], phi1, val);
			}
			errors[1]+=s*weight[i1]*dot_array(2, val, val);
		}
	}
	free_dvector(&luh);

	errors[2]=errors[0]+errors[1];
	
	for(i=0;i<3;i++)
		errors[i]=sqrt(errors[i]);
}

/**
 * \fn void geterrorsPoissonCrouzeixRaviart2d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
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
void geterrorsPoissonCrouzeixRaviart2d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{	
	int Nt=elements->row;
	init_array(3, errors, 0);
	
	int i,j,k;
	int face, edge;
			
	double phi0, phi1[3], phi2[3], val[3];
	int k1,i1,i2,j1;
	double x[2], s, **gradLambda, **vertices;

	int num_qp;
	double lambdas[100][3], weight[100];

	dvector luh;
	create_dvector(elementDOF[0].col, &luh);
	
	num_qp=49; 
	init_Gauss2d(num_qp, lambdas, weight); 
		
	for(k=0;k<Nt;k++)
	{
		// set parameters
		s=elements->vol[k];
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
			baryToCart2d(lambdas[i1], x, vertices);
			// val[0] = -poisson2d_u(x, NULL);
			poisson2d_u(x, &val[0], NULL);
			val[0] *= -1;

			for(k1=0;k1<elementDOF->col;k1++)
			{
				cr_basis(2, lambdas[i1], k1, &phi0);
				val[0]+=phi0*luh.val[k1];
			}
			errors[0]+=s*weight[i1]*val[0]*val[0];
		}

		// H1 semi-norm of u-u_h
		for(i1=0;i1<num_qp;i1++)
		{
			baryToCart2d(lambdas[i1], x, vertices);
            poisson2d_gradu(x, val, NULL);

			for(k1=0;k1<elementDOF->col;k1++)
			{
				cr_basis1(2, gradLambda, k1, phi1);
				axpy_array(2, -luh.val[k1], phi1, val);
			}
			errors[1]+=s*weight[i1]*dot_array(2, val, val);
		}
	}
	free_dvector(&luh);

	errors[2]=errors[0]+errors[1];
	
	for(i=0;i<3;i++)
		errors[i]=sqrt(errors[i]);
}

/**
 * \fn void geterrorsPoissonRaviartThomas2d(double *errors, dvector *uh, dvector *Qhu, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
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
void geterrorsPoissonRaviartThomas2d(double *errors, dvector *uh, dvector *Qhu, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{	
	int Nt=elements->row;
	init_array(6, errors, 0);
	
	int i,j,k;
	int face, edge;
			
	double phi0, phi1[3], phi2[3], val[3];
	int k1,i1,i2,j1;
	double x[2], s, **gradLambda, **vertices, **tij, *height;
	short *eorien;
	
	int num_qp;
	double lambdas[100][3], weight[100];

	dvector lsigmah, luh, lQhu;
	create_dvector(elementDOF[0].col, &lsigmah);
	create_dvector(elementDOF[1].col, &luh);
	create_dvector(elementDOF[1].col, &lQhu);
	
	num_qp=49; 
	init_Gauss2d(num_qp, lambdas, weight); 
		
	for(k=0;k<Nt;k++)
	{
		// set parameters
		s=elements->vol[k];
        gradLambda = elements->gradLambda[k];
		vertices = elements->vertices[k];
		tij = elements->tij[k];
		height = elements->height[k];
		eorien = elements->eorien[k];
		// end set parameters

		for (i = 0; i<elementDOF[0].col; i++){
			j = elementDOF[0].val[k][i];
			lsigmah.val[i] = uh[0].val[j];
		}

		for (i = 0; i<elementDOF[1].col; i++){
			j = elementDOF[1].val[k][i];
			luh.val[i] = uh[1].val[j];
			lQhu.val[i] = Qhu->val[j];
		}

		// L2 norm of sigma - sigma_h
		for(i1=0;i1<num_qp;i1++){
			baryToCart2d(lambdas[i1], x, vertices);
			poisson2d_gradu(x, val, NULL);
			// ax_array(2, -1, val);

			for(k1=0;k1<elementDOF[0].col;k1++){
				rt_basis(lambdas[i1], height, tij, eorien, k1, elementDOF[0].dop, phi1);
				axpy_array(2, -lsigmah.val[k1], phi1, val);
			}
			errors[0]+=s*weight[i1]*dot_array(2, val, val);
		}

		// L2 norm of div(sigma - sigma_h)
		for(i1=0;i1<num_qp;i1++){
			baryToCart2d(lambdas[i1], x, vertices);
			poisson2d_f(x, &val[0], NULL);
			val[0] *= -1;

			for(k1=0;k1<elementDOF[0].col;k1++){
				rt_basisDIV(lambdas[i1], gradLambda, height, tij, eorien, k1, elementDOF[0].dop, &phi0);
				val[0] -= phi0*lsigmah.val[k1];
			}
			errors[1]+=s*weight[i1]*val[0]*val[0];
		}

		// L2 norm of u-u_h
		for(i1=0;i1<num_qp;i1++){
			baryToCart2d(lambdas[i1], x, vertices);
			poisson2d_u(x, &val[0], NULL);
			
			for(k1=0;k1<elementDOF[1].col;k1++){
				lagrange_basis(lambdas[i1], k1, elementDOF[1].dop, &phi0);
				val[0] -= phi0*luh.val[k1];
			}
			errors[2]+=s*weight[i1]*val[0]*val[0];
		}

		// H1 semi-norm of u-u_h
		for(i1=0;i1<num_qp;i1++){
			baryToCart2d(lambdas[i1], x, vertices);
            poisson2d_gradu(x, val, NULL);

			for(k1=0;k1<elementDOF[1].col;k1++){
				lagrange_basis1(lambdas[i1], gradLambda, k1, elementDOF[1].dop, phi1);
				axpy_array(2, -luh.val[k1], phi1, val);
			}
			errors[3]+=s*weight[i1]*dot_array(2, val, val);
		}

		// L2 norm of Qhu-u_h
		for(i1=0;i1<num_qp;i1++){
			val[0] = 0;
			for(k1=0;k1<elementDOF[1].col;k1++){
				lagrange_basis(lambdas[i1], k1, elementDOF[1].dop, &phi0);
				val[0] += phi0*(lQhu.val[k1]-luh.val[k1]);
			}
			errors[4]+=s*weight[i1]*val[0]*val[0];
		}

		// H1 semi-norm of Qhu-u_h
		for(i1=0;i1<num_qp;i1++){
			val[0] = 0; val[1] = 0;
			for(k1=0;k1<elementDOF[1].col;k1++){
				lagrange_basis1(lambdas[i1], gradLambda, k1, elementDOF[1].dop, phi1);
				axpy_array(2, lQhu.val[k1]-luh.val[k1], phi1, val);
			}
			errors[5]+=s*weight[i1]*dot_array(2, val, val);
		}

	}

	free_dvector(&lsigmah);
	free_dvector(&luh);
	free_dvector(&lQhu);
	
	for(i=0;i<6;i++)
		errors[i]=sqrt(errors[i]);

}

/**
 * \fn void geterrorsPoissonBrezziDouglasMarini2d(double *errors, dvector *uh, dvector *Qhu, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
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
void geterrorsPoissonBrezziDouglasMarini2d(double *errors, dvector *uh, dvector *Qhu, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{	
	int Nt=elements->row;
	init_array(6, errors, 0);
	
	int i,j,k;
	int face, edge;
			
	double phi0, phi1[3], phi2[3], val[3];
	int k1,i1,i2,j1;
	double x[2], s, **gradLambda, **vertices, **tij, *height, **nv, **tv;
	short *eorien;
	
	int num_qp;
	double lambdas[100][3], weight[100];

	dvector lsigmah, luh, lQhu;
	create_dvector(elementDOF[0].col, &lsigmah);
	create_dvector(elementDOF[1].col, &luh);
	create_dvector(elementDOF[1].col, &lQhu);
	
	num_qp=49; 
	init_Gauss2d(num_qp, lambdas, weight); 
		
	for(k=0;k<Nt;k++){
		// set parameters
		s=elements->vol[k];
        gradLambda = elements->gradLambda[k];
		vertices = elements->vertices[k];
		tij = elements->tij[k];
		height = elements->height[k];
		nv = elements->nvector[k];
		tv = elements->tvector[k];
		eorien = elements->eorien[k];
		// end set parameters

		for (i = 0; i<elementDOF[0].col; i++){
			j = elementDOF[0].val[k][i];
			lsigmah.val[i] = uh[0].val[j];
		}

		for (i = 0; i<elementDOF[1].col; i++){
			j = elementDOF[1].val[k][i];
			luh.val[i] = uh[1].val[j];
			lQhu.val[i] = Qhu->val[j];
		}

		// L2 norm of sigma - sigma_h
		for(i1=0;i1<num_qp;i1++){
			baryToCart2d(lambdas[i1], x, vertices);
			poisson2d_gradu(x, val, NULL);
			// ax_array(2, -1, val);

			for(k1=0;k1<elementDOF[0].col;k1++){
				bdm_basis(lambdas[i1], height, tij, nv, tv, eorien, k1, elementDOF[0].dop, phi1);
				axpy_array(2, -lsigmah.val[k1], phi1, val);
			}
			errors[0]+=s*weight[i1]*dot_array(2, val, val);
		}

		// L2 norm of div(sigma - sigma_h)
		for(i1=0;i1<num_qp;i1++){
			baryToCart2d(lambdas[i1], x, vertices);
			poisson2d_f(x, &val[0], NULL);
			val[0] *= -1;

			for(k1=0;k1<elementDOF[0].col;k1++){
				bdm_basisDIV(lambdas[i1], gradLambda, height, tij, nv, tv, eorien, k1, elementDOF[0].dop, &phi0);
				val[0] -= phi0*lsigmah.val[k1];
			}
			errors[1]+=s*weight[i1]*val[0]*val[0];
		}

		// L2 norm of u-u_h
		for(i1=0;i1<num_qp;i1++){
			baryToCart2d(lambdas[i1], x, vertices);
			poisson2d_u(x, &val[0], NULL);
			
			for(k1=0;k1<elementDOF[1].col;k1++){
				lagrange_basis(lambdas[i1], k1, elementDOF[1].dop, &phi0);
				val[0] -= phi0*luh.val[k1];
			}
			errors[2]+=s*weight[i1]*val[0]*val[0];
		}

		// H1 semi-norm of u-u_h
		for(i1=0;i1<num_qp;i1++){
			baryToCart2d(lambdas[i1], x, vertices);
            poisson2d_gradu(x, val, NULL);

			for(k1=0;k1<elementDOF[1].col;k1++){
				lagrange_basis1(lambdas[i1], gradLambda, k1, elementDOF[1].dop, phi1);
				axpy_array(2, -luh.val[k1], phi1, val);
			}
			errors[3]+=s*weight[i1]*dot_array(2, val, val);
		}

		// L2 norm of Qhu-u_h
		for(i1=0;i1<num_qp;i1++){
			val[0] = 0;
			for(k1=0;k1<elementDOF[1].col;k1++){
				lagrange_basis(lambdas[i1], k1, elementDOF[1].dop, &phi0);
				val[0] += phi0*(lQhu.val[k1]-luh.val[k1]);
			}
			errors[4]+=s*weight[i1]*val[0]*val[0];
		}

		// H1 semi-norm of Qhu-u_h
		for(i1=0;i1<num_qp;i1++){
			val[0] = 0; val[1] = 0;
			for(k1=0;k1<elementDOF[1].col;k1++){
				lagrange_basis1(lambdas[i1], gradLambda, k1, elementDOF[1].dop, phi1);
				axpy_array(2, lQhu.val[k1]-luh.val[k1], phi1, val);
			}
			errors[5]+=s*weight[i1]*dot_array(2, val, val);
		}

	}

	free_dvector(&lsigmah);
	free_dvector(&luh);
	free_dvector(&lQhu);
	
	for(i=0;i<6;i++)
		errors[i]=sqrt(errors[i]);

}

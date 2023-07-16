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
			baryToCart2d(lambdas[i1], x, vertices);
			biharmonic2d_u(x, &val[0], NULL);
			val[0] *= -1;

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
			baryToCart2d(lambdas[i1], x, vertices);
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
			baryToCart2d(lambdas[i1], x, vertices);
            biharmonic2d_hessu(x, val, NULL);

			for(k1=0;k1<elementDOF->col;k1++)
			{
				morley_basis2(gradLambda, nve, k1, phi2);
				axpy_array(3, -luh.val[k1], phi2, val);
			}
			errors[2]+=s*weight[i1]*(dot_array(3, val, val) + val[2]*val[2]);
		}
	}
	free_dvector(&luh);

	errors[3]=errors[0]+errors[1]+errors[2];
	
	for(i=0;i<4;i++)
		errors[i]=sqrt(errors[i]);
}

/**
 * \fn void geterrorsBiharmonicC0ipdg2d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
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
void geterrorsBiharmonicC0ipdg2d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{	
	int Nt=elements->row;
	init_array(6, errors, 0);
	
	int i,j,k;
	// int face, edge;
			
	double phi0, phi1[3], phi2[3], val[3];
	int k1,i1,i2,j1;
	double x[2], s, **gradLambda, *nve[3], **vertices;

	int dop = elementDOF->dop;

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
			baryToCart2d(lambdas[i1], x, vertices);
			biharmonic2d_u(x, &val[0], NULL);
			val[0] *= -1;

			for(k1=0;k1<elementDOF->col;k1++)
			{
				lagrange_basis(lambdas[i1], k1, dop, &phi0);
				val[0]+=phi0*luh.val[k1];
			}
			errors[0]+=s*weight[i1]*val[0]*val[0];
		}

		// H1 semi-norm of u-u_h
		for(i1=0;i1<num_qp;i1++)
		{
			baryToCart2d(lambdas[i1], x, vertices);
			biharmonic2d_gradu(x, val, NULL);

			for(k1=0;k1<elementDOF->col;k1++)
			{
				lagrange_basis1(lambdas[i1], gradLambda, k1, dop, phi1);
				axpy_array(2, -luh.val[k1], phi1, val);
			}
			errors[1]+=s*weight[i1]*dot_array(2, val, val);
		}

		// H2 semi-norm of u-u_h
		for(i1=0;i1<num_qp;i1++)
		{
			baryToCart2d(lambdas[i1], x, vertices);
            biharmonic2d_hessu(x, val, NULL);

			for(k1=0;k1<elementDOF->col;k1++)
			{
				lagrange_basis2(lambdas[i1], gradLambda, k1, dop, phi2);
				axpy_array(3, -luh.val[k1], phi2, val);
			}
			errors[2]+=s*weight[i1]*(dot_array(3, val, val) + val[2]*val[2]);
		}
	}
	free_dvector(&luh);

	int num_qp1;
	double lambdas1[100][2], weight1[100];
	num_qp1=getNumQuadPoints((dop - 1) * 2, 1);
	if (num_qp1>5) num_qp1 = 5;
	init_Gauss1d(num_qp1, lambdas1, weight1);

	int element[2];
	int count, istart;
	int *index;
	index = (int*)calloc(elementDOF->dof, sizeof(int));
	for (i = 0; i<elementDOF->dof; i++)
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
		for (i = 0; i<elementDOF->col; i++)
		{
			j = elementDOF->val[element[0]][i];
			patchnodes[count] = j;
			count++;
			index[j] = istart;
			istart = j;
		}

		if (element[1] != -1)
		{
			for (i = 0; i<elementDOF->col; i++)
			{
				j = elementDOF->val[element[1]][i];
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

		// jump of normal derivative
		for (i1 = 0; i1<num_qp1; i1++)
		{
			val[0] = 0;

			for (k1 = 0; k1<count; k1++)
			{
				i = patchnodes[k1];
				jumpNormalDerivativeLagrange2d(lambdas1[i1], edge, elements, elementEdge, edges, elementDOF, i, &phi0);
				val[0] += uh->val[i] * phi0;
			}

			errors[4] += weight1[i1] * val[0] * val[0];
		}
	} // edge
	free(index);

	errors[3] = errors[0]+errors[1]+errors[2];
	errors[5] = errors[3]+errors[4];
	
	for(i=0;i<6;i++)
		errors[i]=sqrt(errors[i]);
}

/**
 * \fn void geterrorsBiharmonicC0ipdgExtrap2d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran)
 * \brief compute error between numerical solution and exact solution: L inf norm, 1 inf norm, 2 inf norm
 * \param *errors pointer to error between numerical solution and exact solution: L inf norm, 1 inf norm, 2 inf norm
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
void geterrorsBiharmonicC0ipdgExtrap2d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran)
{	
	// int Nt=elements->row;
	init_array(9, errors, 0);
	
	int i,j,k,l;
	// int face, edge;
			
	double phi0, phi1[3], phi2[3], val[3];
	int k1; //,i1,i2,j1;
	double *x, s, **gradLambda, *nve[3], **vertices;

	int dop = elementDOF->dop;

	int num_qp=3;
	double lambdas[100][3], weight[100];
	lambdas[0][0] = 1; lambdas[0][1] = 0; lambdas[0][2] = 0;
	lambdas[1][0] = 0; lambdas[1][1] = 1; lambdas[1][2] = 0;
	lambdas[2][0] = 0; lambdas[2][1] = 0; lambdas[2][2] = 1;

	ddenmat elementuh[2];
	ddenmat3 elementGraduh[2], elementHessuh[2];
	ddenmat Gradwh[2], Hesswh[2];
	dvector luh, wh[2];

	for(l=0;l<2;l++){
		create_dden_matrix(elements[l].row, 3, elementuh+l);
		create_dden_matrix3(elements[l].row, 3, 2, elementGraduh+l);
		create_dden_matrix3(elements[l].row, 3, 3, elementHessuh+l);
		
		create_dvector(nodes[l].row, wh+l);
		create_dden_matrix(nodes[l].row, 2, Gradwh+l);
		create_dden_matrix(nodes[l].row, 3, Hesswh+l);
	}

	for(l=0;l<2;l++){	
		create_dvector(elementDOF[l].col, &luh);
		for(k=0;k<elements[l].row;k++){
			// set parameters
			vertices = elements[l].vertices[k];
			s=elements[l].vol[k];
        	gradLambda = elements[l].gradLambda[k];
			for(i=0;i<3;i++){
				j = elementEdge[l].val[k][i];
				nve[i] = edges[l].nvector[j];
			}
			// end set parameters
		
			for (i = 0; i<elementDOF[l].col; i++)
			{
				j = elementDOF[l].val[k][i];
				luh.val[i] = uh[l].val[j];
			}

			// u_h, grad u_h, hess u_h
			for(i=0;i<3;i++)
			{	
				for(k1=0;k1<elementDOF[l].col;k1++)
				{
					lagrange_basis(lambdas[i], k1, dop, &phi0);
					elementuh[l].val[k][i] += phi0*luh.val[k1];

					lagrange_basis1(lambdas[i], gradLambda, k1, dop, phi1);
					axpy_array(2, luh.val[k1], phi1, elementGraduh[l].val[k][i]);

					lagrange_basis2(lambdas[i], gradLambda, k1, dop, phi2);
					axpy_array(3, luh.val[k1], phi2, elementHessuh[l].val[k][i]);
				}
			}
		}
		free_dvector(&luh);
	}
	
	int element, idx, nn;
	for(l=0;l<2;l++){
		for(i=0;i<nodes[l].row;i++){
			nn = elementdofTran[l].IA[i + 1] - elementdofTran[l].IA[i];
			for (j = elementdofTran[l].IA[i]; j<elementdofTran[l].IA[i + 1]; j++){
				element = elementdofTran[l].JA[j];
				for(idx=0;idx<3;idx++){
					if(elementDOF[l].val[element][idx]==i)
						break;
				}
				wh[l].val[i] += elementuh[l].val[element][idx] / nn;
				axpy_array(2, 1./nn, elementGraduh[l].val[element][idx], Gradwh[l].val[i]);
				axpy_array(3, 1./nn, elementHessuh[l].val[element][idx], Hesswh[l].val[i]);
			}
		}
	}

	for(l=0;l<2;l++){
		free_dden_matrix(elementuh+l);
		free_dden_matrix3(elementGraduh+l);
		free_dden_matrix3(elementHessuh+l);
	}

	double c[3][2];
	for(i=0;i<3;i++){
		c[i][0]=0; c[i][1]=1;
	}

	if(dop == 2){
		c[0][0] = -1.0/3.0; c[0][1] = 4.0/3.0;
		c[1][0] = -1.0/3.0; c[1][1] = 4.0/3.0;
		c[2][0] = -1.0; c[2][1] = 2.0;
	}
	else if(dop == 3){
		c[0][0] = -1.0/15.0; c[0][1] = 16.0/15.0;
		c[1][0] = -1.0/7.0; c[1][1] = 8.0/7.0;
		c[2][0] = -1.0/3.0; c[2][1] = 4.0/3.0;
	}
	
	for(i=0;i<nodes[0].row;i++){
		x = nodes[0].val[i];
		if(fabs(x[0]-0.5)<0.25+1e-10 && fabs(x[1]-0.5)<0.25+1e-10){
			phi0 = c[0][0]*wh[0].val[i] + c[0][1]*wh[1].val[i];
			biharmonic2d_u(x, &val[0], NULL);
			phi0 -= val[0];
			phi0 = fabs(phi0);
			if(phi0>errors[0]) errors[0]=phi0;

			axpbyz_array(2, c[1][0], Gradwh[0].val[i], c[1][1], Gradwh[1].val[i], phi1);
			biharmonic2d_gradu(x, val, NULL);
			axpy_array(2, -1.0, val, phi1);
			phi0 = sqrt(phi1[0]*phi1[0]+phi1[1]*phi1[1]);
			if(phi0>errors[1]) errors[1]=phi0;

			axpbyz_array(3, c[2][0], Hesswh[0].val[i], c[2][1], Hesswh[1].val[i], phi2);
			biharmonic2d_hessu(x, val, NULL);
			axpy_array(3, -1.0, val, phi2);
			phi0 = sqrt(phi2[0]*phi2[0]+phi2[1]*phi2[1]+2*phi2[2]*phi2[2]);
			if(phi0>errors[2]) errors[2]=phi0;
		}
	}

	for(l=0;l<2;l++){
		for(i=0;i<nodes[l].row;i++){
			x = nodes[l].val[i];
			if(fabs(x[0]-0.5)<0.25+1e-10 && fabs(x[1]-0.5)<0.25+1e-10){
				phi0 = wh[l].val[i];
				biharmonic2d_u(x, &val[0], NULL);
				phi0 -= val[0];
				phi0 = fabs(phi0);
				if(phi0>errors[3+3*l]) errors[3+3*l]=phi0;

				copy_array(2, Gradwh[l].val[i], phi1);
				biharmonic2d_gradu(x, val, NULL);
				axpy_array(2, -1.0, val, phi1);
				phi0 = sqrt(phi1[0]*phi1[0]+phi1[1]*phi1[1]);
				if(phi0>errors[3+3*l+1]) errors[3+3*l+1]=phi0;

				copy_array(3, Hesswh[l].val[i], phi2);
				biharmonic2d_hessu(x, val, NULL);
				axpy_array(3, -1.0, val, phi2);
				phi0 = sqrt(phi2[0]*phi2[0]+phi2[1]*phi2[1]+2*phi2[2]*phi2[2]);
				if(phi0>errors[3+3*l+2]) errors[3+3*l+2]=phi0;
			}
		}
	}
	
	for(l=0;l<2;l++){
		free_dvector(wh+l);
		free_dden_matrix(Gradwh+l);
		free_dden_matrix(Hesswh+l);
	}
}

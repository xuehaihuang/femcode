/*
 *  post_processing.c
 *  FEM
 *
 *  Created by Xuehai Huang on 06/08/23.
 *  Copyright 2023 SUFE. All rights reserved.
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
 * \fn void postprocessQuadcurlDistribMixedFEM3d(dvector *uhstar, dvector *sigmah, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, ELEMENT_DOF *elementDOFuhstar, int curlfemtype)
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
 * \param curlfemtype type of Nedelec element
 */
void postprocessQuadcurlDistribMixedFEM3d(dvector *uhstar, dvector *sigmah, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, ELEMENT_DOF *elementDOFuhstar, int curlfemtype)
{	
	int Nt=elements->row;
	int i, j, k;
	
	int face, edge;

	double phi0[9], phi1[9], phi2[9], val[9], val1;
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], **grd_lambda, **nv, *nvf[4], *tve, **vertices;
	double vol, s, elen;
	int ei[2];
	int count;
	short *forien, *eorien;
	int **fpermi, **eperm;
	double *lambdaConst;
			
	int num_qp, num_qp2, num_qp1;
	double lambdas[100][4], lambdas2[100][4], lambdas1[100][2], weight[100], weight2[100], weight1[100];
	double lambdasT[4]; 
	int dop = elementDOFuhstar->dop;
	
	dvector lsigmah, luh;
	create_dvector(8*elementDOF[0].col, &lsigmah);
	create_dvector(elementDOF[1].col, &luh);
	
	num_qp=getNumQuadPoints_ShunnWilliams(2*(dop-2), 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	num_qp2=getNumQuadPoints_ShunnWilliams(2*dop-1, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp2, lambdas2, weight2); // Shunn-Williams intergation initial
	num_qp1=getNumQuadPoints(dop, 1); // the number of numerical intergation points
	init_Gauss1d(num_qp1, lambdas1, weight1); // Gauss intergation initial
		
	ddenmat lA, lb; // local A
	int row = elementDOFuhstar->col+(dop+1)*(dop+2)*(dop+3)/6+2;
	create_dden_matrix(row, row, &lA);
	create_dden_matrix(lA.row, 1, &lb);
	for(k=0;k<Nt;k++){
		// set parameters
		vol = elements->vol[k];
		grd_lambda = elements->gradLambda[k];
		vertices = elements->vertices[k];
		nv = elements->nvector[k];
		forien = elements->forien[k];
		eorien = elements->eorien[k];
		eperm = elements->eperm[k];
		fpermi = elements->fpermi[k];
		lambdaConst = elements->lambdaConst[k];
		// end set parameters
		
		init_dden_matrix(&lA, 0.0);
		init_dden_matrix(&lb, 0.0);

		for (i = 0; i<8*elementDOF[0].col; i++){
			j = elementDOF[0].col *8*k +i;
			lsigmah.val[i] = sigmah->val[j];
		}

		for (i = 0; i<elementDOF[1].col; i++){
			j = elementDOF[1].val[k][i];
			luh.val[i] = uh->val[j];
		}

		for (k1 = 0; k1<elementDOFuhstar->col; k1++){
			for (k2 = 0; k2<elementDOFuhstar->col; k2++){
				val1 = 0;
				for (i1 = 0; i1<num_qp; i1++){
					nedelec1st3d_basisGradCurl(lambdas[i1], grd_lambda, eorien, fpermi, k1, dop, phi1);
					nedelec1st3d_basisGradCurl(lambdas[i1], grd_lambda, eorien, fpermi, k2, dop, phi2);
					val1 += vol*weight[i1] * dot_array(9, phi1, phi2);
				}
				lA.val[k1][k2] += val1;
			} // k2
			for (i1 = 0; i1<num_qp; i1++){
				nedelec1st3d_basisGradCurl(lambdas[i1], grd_lambda, eorien, fpermi, k1, dop, phi1);
				init_array(9, val, 0);
				for (k2 = 0; k2<elementDOF[0].col; k2++){
					lagrange3d_basis(lambdas[i1], k2, elementDOF[0].dop, phi0);
					val[1] += lsigmah.val[k2] * phi0[0];
					val[2] += lsigmah.val[elementDOF[0].col+k2] * phi0[0];
					val[3] += lsigmah.val[2*elementDOF[0].col+k2] * phi0[0];
					val[5] += lsigmah.val[3*elementDOF[0].col+k2] * phi0[0];
					val[6] += lsigmah.val[4*elementDOF[0].col+k2] * phi0[0];
					val[7] += lsigmah.val[5*elementDOF[0].col+k2] * phi0[0];
					val[0] += lsigmah.val[6*elementDOF[0].col+k2] * phi0[0]/sqrt(2);
					val[4] -= lsigmah.val[6*elementDOF[0].col+k2] * phi0[0]/sqrt(2);
					val[0] += lsigmah.val[7*elementDOF[0].col+k2] * phi0[0]/sqrt(6);
					val[4] += lsigmah.val[7*elementDOF[0].col+k2] * phi0[0]/sqrt(6);
					val[8] -= lsigmah.val[7*elementDOF[0].col+k2] * phi0[0]*2/sqrt(6);
				} // k2
				lb.val[k1][0] += vol*weight[i1] * dot_array(9, phi1, val);
			} // i1
		} // k1

		for (k1 = 4; k1<(dop+1)*(dop+2)*(dop+3)/6; k1++){
			for (k2 = 0; k2<elementDOFuhstar->col; k2++){
				val1 = 0;
				for (i1 = 0; i1<num_qp2; i1++){
					bernstein3d_basis1(lambdas2[i1], grd_lambda, k1, dop, phi1);
					nedelec1st3d_basis(lambdas2[i1], grd_lambda, eorien, fpermi, k2, dop, phi2);
					val1 += vol*weight2[i1] * dot_array(3, phi1, phi2);
				}
				lA.val[k1-4+elementDOFuhstar->col][k2] += val1;
				lA.val[k2][k1-4+elementDOFuhstar->col] += val1;
			} // k2
			for (i1 = 0; i1<num_qp2; i1++){
				bernstein3d_basis1(lambdas2[i1], grd_lambda, k1, dop, phi1);
				init_array(3, val, 0);
				for (k2 = 0; k2<elementDOF[1].col; k2++){
					if(curlfemtype == 1)
						nedelec1st3d_basis(lambdas2[i1], grd_lambda, eorien, fpermi, k2, elementDOF[1].dop, phi2);
					else
						nedelec2nd3d_basis(lambdas2[i1], grd_lambda, eperm, k2, elementDOF[1].dop, phi2);
					axpy_array(3, luh.val[k2], phi2, val);
				} // k2
				lb.val[k1-4+elementDOFuhstar->col][0] += vol*weight2[i1] * dot_array(3, phi1, val);
			} // i1
		} // k1

		for (k1 = 0; k1<6; k1++){
			edge = elementEdge->val[k][k1];
			tve = edges->tvector[edge];
			elen = edges->length[edge];
			edge2vv3d(k1, ei);
			for(i=0;i<4;i++) lambdasT[i]=0;
			for (k2 = 0; k2<elementDOFuhstar->col; k2++){
				val1 = 0;
				for (i1 = 0; i1<num_qp1; i1++){
					lambdasT[ei[0]] = lambdas1[i1][0];
					lambdasT[ei[1]] = lambdas1[i1][1];
					nedelec1st3d_basis(lambdasT, grd_lambda, eorien, fpermi, k2, dop, phi2);
					val1 += elen*weight1[i1] * dot_array(3, tve, phi2);
				}
				lA.val[row-6+k1][k2] += val1;
				lA.val[k2][row-6+k1] += val1;
			} // k2
			for (i1 = 0; i1<num_qp1; i1++){
				lambdasT[ei[0]] = lambdas1[i1][0];
				lambdasT[ei[1]] = lambdas1[i1][1];
				init_array(3, val, 0);
				for (k2 = 0; k2<elementDOF[1].col; k2++){
					if(curlfemtype == 1)
						nedelec1st3d_basis(lambdasT, grd_lambda, eorien, fpermi, k2, elementDOF[1].dop, phi2);
					else
						nedelec2nd3d_basis(lambdasT, grd_lambda, eperm, k2, elementDOF[1].dop, phi2);
					axpy_array(3, luh.val[k2], phi2, val);
				} // k2
				lb.val[row-6+k1][0] += elen*weight1[i1] * dot_array(3, tve, val);
			} // i1
		} // k1

		AxBrref(&lA, &lb);
		for (i = 0; i<elementDOFuhstar->col; i++){
			uhstar->val[elementDOFuhstar->val[k][i]] = lb.val[i][0];
		}

	}
// exit(1);////////////////////////
	free_dden_matrix(&lA);
	free_dden_matrix(&lb);
	free_dvector(&lsigmah);
	free_dvector(&luh);
}

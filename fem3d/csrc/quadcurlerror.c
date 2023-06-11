/*
 *  quadcurlerror.c
 *  FEM
 *
 *  Created by Xuehai Huang on May 3, 2022.
 *  Copyright 2022 SUFE. All rights reserved.
 *
 */

/*! \file quadcurlerror.c
 *  \brief compute error between numerical solution and exact solution
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <time.h>
#include "header.h"
#include "matvec.h"

/**
 * \fn void geterrorsQuadcurlHuang3d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
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
void geterrorsQuadcurlHuang3d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{	
	int Nt=elements->row;
	int i, j, k;

	init_array(4, errors, 0);
	
	int face, edge;

	double phi1[9], val[9];
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], *xK, **grd_lambda, **nv, *nvf[4], *etv[6], **vertices;
	double vol, s[4], h[6];
	int rowstart[3], row31[3];
	int count;
	short *forien, *eorien;
	int **fpermi;
	double *lambdaConst;
			
	int num_qp;
	double lambdas[100][4], lambdas2[100][3], weight[100];

	dvector luh;
	create_dvector(elementDOF[0].col, &luh);
	
	num_qp=getNumQuadPoints_ShunnWilliams(9, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
		
	for(k=0;k<Nt;k++)
	{
		// set parameters
		vol = elements->vol[k];
		xK = elements->barycenter[k];
		grd_lambda = elements->gradLambda[k];
		vertices = elements->vertices[k];
		nv = elements->nvector[k];
		forien = elements->forien[k];
		eorien = elements->eorien[k];
		fpermi = elements->fpermi[k];
		lambdaConst = elements->lambdaConst[k];
		for (i = 0; i < elementFace->col; i++)
		{
			face = elementFace->val[k][i];
			nvf[i] = faces->nvector[face];
			s[i] = faces->area[face];
		}
		for (i = 0; i<6; i++)
		{
			edge = elementEdge->val[k][i];
			etv[i] = edges->tvector[edge];
			h[i] = edges->length[edge];
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
			axy_array(3, lambdas[i1][3], vertices[3], x);
			for(i=0;i<3;i++)
				axpy_array(3, lambdas[i1][i], vertices[i], x);
			quadcurl3d_u(x, val);

			for(k1=0;k1<elementDOF[0].col;k1++)
			{
				huangQuadcurl3d_basis(x, xK, lambdas[i1], grd_lambda, vertices, nvf, eorien, fpermi, k1, phi1);
				axpy_array(3, -luh.val[k1], phi1, val);
			}
			errors[0] += vol*weight[i1] * dot_array(3, val, val);
		}

		// L2 norm of curl(u-u_h)
		for(i1=0;i1<num_qp;i1++)
		{
			axy_array(3, lambdas[i1][3], vertices[3], x);
			for(i=0;i<3;i++)
				axpy_array(3, lambdas[i1][i], vertices[i], x);
			quadcurl3d_curlu(x, val);

			for(k1=0;k1<elementDOF[0].col;k1++)
			{
				huangQuadcurl3d_basisCurl(lambdas[i1], grd_lambda, nvf, eorien, fpermi, k1, phi1);
				axpy_array(3, -luh.val[k1], phi1, val);
			}
			errors[1] += vol*weight[i1] * dot_array(3, val, val);
		}

		// H1 semi-norm of curl(u-u_h)
		for(i1=0;i1<num_qp;i1++)
		{
			axy_array(3, lambdas[i1][3], vertices[3], x);
			for(i=0;i<3;i++)
				axpy_array(3, lambdas[i1][i], vertices[i], x);
			quadcurl3d_gradcurlu(x, val);

			for(k1=0;k1<elementDOF[0].col;k1++)
			{
				huangQuadcurl3d_basisGradCurl(grd_lambda, nvf, eorien, fpermi, k1, phi1);
				axpy_array(9, -luh.val[k1], phi1, val);
			}
			errors[2] += vol*weight[i1] * dot_array(9, val, val);
		}
	}
	free_dvector(&luh);

	errors[3] = errors[0] + errors[1] + errors[2];

	for (i = 0; i < 4; i++)
		errors[i] = sqrt(errors[i]);
}

/**
 * \fn void geterrorsQuadcurlHuangZhang3d(double *errors, dvector *uh, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
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
void geterrorsQuadcurlHuangZhang3d(double *errors, dvector *uh, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{	
	int Nt=elements->row;
	int i, j, k;

	init_array(4, errors, 0);
	
	int face, edge;

	double phi1[9], val[9];
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], *xK, **grd_lambda, **nv, *nvf[4], *etv[6], **vertices;
	double vol, s[4], h[6];
	int rowstart[3], row31[3];
	int count;
	short *forien, *eorien;
	int **fpermi;
	double *lambdaConst;
			
	int num_qp;
	double lambdas[100][4], lambdas2[100][3], weight[100];

	dvector luh, lvh;
	create_dvector(elementDOF[0].col, &luh);
	create_dvector(elementDOF[0].col, &lvh);
	
	num_qp=getNumQuadPoints_ShunnWilliams(9, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
		
	for(k=0;k<Nt;k++)
	{
		// set parameters
		vol = elements->vol[k];
		xK = elements->barycenter[k];
		grd_lambda = elements->gradLambda[k];
		vertices = elements->vertices[k];
		nv = elements->nvector[k];
		forien = elements->forien[k];
		eorien = elements->eorien[k];
		fpermi = elements->fpermi[k];
		lambdaConst = elements->lambdaConst[k];
		for (i = 0; i < elementFace->col; i++)
		{
			face = elementFace->val[k][i];
			nvf[i] = faces->nvector[face];
			s[i] = faces->area[face];
		}
		for (i = 0; i<6; i++)
		{
			edge = elementEdge->val[k][i];
			etv[i] = edges->tvector[edge];
			h[i] = edges->length[edge];
		}
		// end set parameters

		for (i = 0; i<elementDOF[0].col; i++)
		{
			j = elementDOF[0].val[k][i];
			lvh.val[i] = uh->val[j];
		}

		ddenmat lbC;
		lbC.row=basisCoeffs->row; lbC.col=basisCoeffs->col; lbC.val=basisCoeffs->val[k];
		Atxy_ddenmat(1.0, &lbC, lvh.val, luh.val);

		// L2 norm of u-u_h
		for(i1=0;i1<num_qp;i1++)
		{
			axy_array(3, lambdas[i1][3], vertices[3], x);
			for(i=0;i<3;i++)
				axpy_array(3, lambdas[i1][i], vertices[i], x);
			quadcurl3d_u(x, val);

			for(k1=0;k1<elementDOF[0].col;k1++)
			{
				huangzhang03d_basis(lambdas[i1], grd_lambda, k1, phi1);
				axpy_array(3, -luh.val[k1], phi1, val);
			}
			errors[0] += vol*weight[i1] * dot_array(3, val, val);
		}

		// L2 norm of curl(u-u_h)
		for(i1=0;i1<num_qp;i1++)
		{
			axy_array(3, lambdas[i1][3], vertices[3], x);
			for(i=0;i<3;i++)
				axpy_array(3, lambdas[i1][i], vertices[i], x);
			quadcurl3d_curlu(x, val);

			for(k1=0;k1<elementDOF[0].col;k1++)
			{
				huangzhang03d_basisCurl(lambdas[i1], grd_lambda, k1, phi1);
				axpy_array(3, -luh.val[k1], phi1, val);
			}
			errors[1] += vol*weight[i1] * dot_array(3, val, val);
		}

		// H1 semi-norm of curl(u-u_h)
		for(i1=0;i1<num_qp;i1++)
		{
			axy_array(3, lambdas[i1][3], vertices[3], x);
			for(i=0;i<3;i++)
				axpy_array(3, lambdas[i1][i], vertices[i], x);
			quadcurl3d_gradcurlu(x, val);

			for(k1=0;k1<elementDOF[0].col;k1++)
			{
				huangzhang03d_basisGradCurl(lambdas[i1], grd_lambda, k1, phi1);
				axpy_array(9, -luh.val[k1], phi1, val);
			}
			errors[2] += vol*weight[i1] * dot_array(9, val, val);
		}
	}
	free_dvector(&luh);
	free_dvector(&lvh);

	errors[3] = errors[0] + errors[1] + errors[2];

	for (i = 0; i < 4; i++)
		errors[i] = sqrt(errors[i]);
}

/**
 * \fn void geterrorsQuadcurlDistribMixedFEM3d(double *errors, dvector *sigmah, dvector *uh, dvector *lambdah, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, int curlfemtype)
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
void geterrorsQuadcurlDistribMixedFEM3d(double *errors, dvector *sigmah, dvector *uh, dvector *lambdah, dvector *uhstar, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, ELEMENT_DOF *elementDOFstar, int curlfemtype)
{	
	int Nt=elements->row;
	int i, j, k;

	init_array(7, errors, 0);
	
	int face, edge;

	double phi0[9], phi1[9], phi2[9], val[9];
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], **grd_lambda, **nv, *nvf[4], *etv[6], **vertices;
	double vol, s[4], h[6];
	int count;
	short *forien, *eorien;
	int **fpermi, **eperm;
	double *lambdaConst;
			
	int num_qp;
	double lambdas[100][4], lambdas2[100][3], weight[100];

	dvector lsigmah, luh, luhstar;
	create_dvector(8*elementDOF[0].col, &lsigmah);
	create_dvector(elementDOF[1].col, &luh);
	create_dvector(elementDOFstar->col, &luhstar);
	
	num_qp=getNumQuadPoints_ShunnWilliams(9, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
		
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
		// for (i = 0; i < elementFace->col; i++){
		// 	face = elementFace->val[k][i];
		// 	nvf[i] = faces->nvector[face];
		// 	s[i] = faces->area[face];
		// }
		// for (i = 0; i<6; i++){
		// 	edge = elementEdge->val[k][i];
		// 	etv[i] = edges->tvector[edge];
		// 	h[i] = edges->length[edge];
		// }
		// end set parameters

		for (i = 0; i<8*elementDOF[0].col; i++){
			j = elementDOF[0].col *8*k +i;
			lsigmah.val[i] = sigmah->val[j];
		}

		for (i = 0; i<elementDOF[1].col; i++){
			j = elementDOF[1].val[k][i];
			luh.val[i] = uh->val[j];
		}

		for (i = 0; i<elementDOFstar->col; i++){
			j = elementDOFstar->val[k][i];
			luhstar.val[i] = uhstar->val[j];
		}

		// L2 norm of sigma-sigma_h
		for(i1=0;i1<num_qp;i1++){
			axy_array(3, lambdas[i1][3], vertices[3], x);
			for(i=0;i<3;i++)
				axpy_array(3, lambdas[i1][i], vertices[i], x);
			quadcurl3d_gradcurlu(x, val);

			for(k1=0;k1<elementDOF[0].col;k1++){
				lagrange3d_basis(lambdas[i1], k1, elementDOF[0].dop, phi0);

				val[1] -= lsigmah.val[k1] * phi0[0];
				val[2] -= lsigmah.val[elementDOF[0].col+k1] * phi0[0];
				val[3] -= lsigmah.val[2*elementDOF[0].col+k1] * phi0[0];
				val[5] -= lsigmah.val[3*elementDOF[0].col+k1] * phi0[0];
				val[6] -= lsigmah.val[4*elementDOF[0].col+k1] * phi0[0];
				val[7] -= lsigmah.val[5*elementDOF[0].col+k1] * phi0[0];
				val[0] -= lsigmah.val[6*elementDOF[0].col+k1] * phi0[0]/sqrt(2);
				val[4] += lsigmah.val[6*elementDOF[0].col+k1] * phi0[0]/sqrt(2);
				val[0] -= lsigmah.val[7*elementDOF[0].col+k1] * phi0[0]/sqrt(6);
				val[4] -= lsigmah.val[7*elementDOF[0].col+k1] * phi0[0]/sqrt(6);
				val[8] += lsigmah.val[7*elementDOF[0].col+k1] * phi0[0]*2/sqrt(6);
			}
			errors[0] += vol*weight[i1] * dot_array(9, val, val);
		}

		// L2 norm of u-u_h
		for(i1=0;i1<num_qp;i1++){
			axy_array(3, lambdas[i1][3], vertices[3], x);
			for(i=0;i<3;i++)
				axpy_array(3, lambdas[i1][i], vertices[i], x);
			quadcurl3d_u(x, val);

			for(k1=0;k1<elementDOF[1].col;k1++){
				if(curlfemtype == 1)
					nedelec1st3d_basis(lambdas[i1], grd_lambda, eorien, fpermi, k1, elementDOF[1].dop, phi1);
				else
					nedelec2nd3d_basis(lambdas[i1], grd_lambda, eperm, k1, elementDOF[1].dop, phi1);
				axpy_array(3, -luh.val[k1], phi1, val);
			}
			errors[1] += vol*weight[i1] * dot_array(3, val, val);
		}

		// L2 norm of curl(u-u_h)
		for(i1=0;i1<num_qp;i1++){
			axy_array(3, lambdas[i1][3], vertices[3], x);
			for(i=0;i<3;i++)
				axpy_array(3, lambdas[i1][i], vertices[i], x);
			quadcurl3d_curlu(x, val);

			for(k1=0;k1<elementDOF[1].col;k1++){
				if(curlfemtype == 1)
					nedelec1st3d_basisCurl(lambdas[i1], grd_lambda, eorien, fpermi, k1, elementDOF[1].dop, phi1);
				else
					nedelec2nd3d_basisCurl(lambdas[i1], grd_lambda, eperm, k1, elementDOF[1].dop, phi1);
				axpy_array(3, -luh.val[k1], phi1, val);
			}
			errors[2] += vol*weight[i1] * dot_array(3, val, val);
		}

		// H1 seminorm of curl(u-u_h)
		for(i1=0;i1<num_qp;i1++){
			axy_array(3, lambdas[i1][3], vertices[3], x);
			for(i=0;i<3;i++)
				axpy_array(3, lambdas[i1][i], vertices[i], x);
			quadcurl3d_gradcurlu(x, val);

			for(k1=0;k1<elementDOF[1].col;k1++){
				if(curlfemtype == 1)
					nedelec1st3d_basisGradCurl(lambdas[i1], grd_lambda, eorien, fpermi, k1, elementDOF[1].dop, phi1);
				else
					nedelec2nd3d_basisGradCurl(lambdas[i1], grd_lambda, eperm, k1, elementDOF[1].dop, phi1);
				axpy_array(9, -luh.val[k1], phi1, val);
			}
			errors[3] += vol*weight[i1] * dot_array(9, val, val);
		}

		// L2 norm of u-u_h*
		for(i1=0;i1<num_qp;i1++){
			axy_array(3, lambdas[i1][3], vertices[3], x);
			for(i=0;i<3;i++)
				axpy_array(3, lambdas[i1][i], vertices[i], x);
			quadcurl3d_u(x, val);

			for(k1=0;k1<elementDOFstar->col;k1++){
				nedelec1st3d_basis(lambdas[i1], grd_lambda, eorien, fpermi, k1, elementDOFstar->dop, phi1);
				axpy_array(3, -luhstar.val[k1], phi1, val);
			}
			errors[4] += vol*weight[i1] * dot_array(3, val, val);
		}

		// L2 norm of curl(u-u_h*)
		for(i1=0;i1<num_qp;i1++){
			axy_array(3, lambdas[i1][3], vertices[3], x);
			for(i=0;i<3;i++)
				axpy_array(3, lambdas[i1][i], vertices[i], x);
			quadcurl3d_curlu(x, val);

			for(k1=0;k1<elementDOFstar->col;k1++){
				nedelec1st3d_basisCurl(lambdas[i1], grd_lambda, eorien, fpermi, k1, elementDOFstar->dop, phi1);
				axpy_array(3, -luhstar.val[k1], phi1, val);
			}
			errors[5] += vol*weight[i1] * dot_array(3, val, val);
		}

		// H1 seminorm of curl(u-u_h*)
		for(i1=0;i1<num_qp;i1++){
			axy_array(3, lambdas[i1][3], vertices[3], x);
			for(i=0;i<3;i++)
				axpy_array(3, lambdas[i1][i], vertices[i], x);
			quadcurl3d_gradcurlu(x, val);

			for(k1=0;k1<elementDOFstar->col;k1++){
				nedelec1st3d_basisGradCurl(lambdas[i1], grd_lambda, eorien, fpermi, k1, elementDOFstar->dop, phi1);
				axpy_array(9, -luhstar.val[k1], phi1, val);
			}
			errors[6] += vol*weight[i1] * dot_array(9, val, val);
		}
	}
	free_dvector(&lsigmah);
	free_dvector(&luh);
	free_dvector(&luhstar);

	for (i = 0; i < 7; i++)
		errors[i] = sqrt(errors[i]);
}

/*
 *  basiscoeff.c
 *
 *  Created by Xuehai Huang on 5/10/2022.
 *  Copyright 2022 SUFE. All rights reserved.
 *
 */

/*! \file basiscoeff.c
 *  \brief coefficients of Basis functions
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"
#include "matvec.h"


/**
 * \fn void getHuangZhangBasisCoeffs(ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges)
 * \brief generate the coefficients of basis functions
 * \param *basisCoeffs pointer to the coefficients of basis functions
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param *nodes pointer to the nodes location of the triangulation
 * \return void
 */
void getHuangZhangBasisCoeffs(ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges)
{
	int i, j, k;
	int face;
	ddenmat A, B;
	double **C;
	short *eorien;
	double v[4][3], **vertices, **bcFace, **gradLambda, gradCross[4][4][3], nxv[4][3], nxvT[4][3], nxvB[4][3];
	double *nv, *tv[2], *t, s;
	double c[4], cc[4][4], val[3];

	create_dden_matrix(32, 32, &A);
	create_dden_matrix(32, 32, &B);
	for(k=0;k<elements->row;k++)
	{
		eorien = elements->eorien[k];
		vertices = elements->vertices[k];
		bcFace = elements->bcFace[k];
		gradLambda = elements->gradLambda[k];
		for(i=0;i<3;i++)
		{
			for(j=i+1;j<4;j++)
			{
				cross_array(gradLambda[i], gradLambda[j], gradCross[i][j]);
				axy_array(3, -1.0, gradCross[i][j], gradCross[j][i]);
			}
		}

		init_dden_matrix(&A, 0.0);
		// edge e01
		A.val[0][0] = 1.0/3.0 * eorien[0]; A.val[0][1] = 1.0/6.0 * eorien[0];
		A.val[1][0] = 1.0/6.0 * eorien[0]; A.val[1][1] = 1.0/3.0 * eorien[0];
		// edge e12
		A.val[2][2] = 1.0/3.0 * eorien[1]; A.val[2][3] = 1.0/6.0 * eorien[1];
		A.val[3][2] = 1.0/6.0 * eorien[1]; A.val[3][3] = 1.0/3.0 * eorien[1];
		// edge e23
		A.val[4][4] = 1.0/3.0 * eorien[2]; A.val[4][5] = 1.0/6.0 * eorien[2];
		A.val[5][4] = 1.0/6.0 * eorien[2]; A.val[5][5] = 1.0/3.0 * eorien[2];
		// edge e30
		A.val[6][6] = 1.0/3.0 * eorien[3]; A.val[6][7] = 1.0/6.0 * eorien[3];
		A.val[7][6] = 1.0/6.0 * eorien[3]; A.val[7][7] = 1.0/3.0 * eorien[3];
		// edge e02
		A.val[8][8] = 1.0/3.0 * eorien[4]; A.val[8][9] = 1.0/6.0 * eorien[4];
		A.val[9][8] = 1.0/6.0 * eorien[4]; A.val[9][9] = 1.0/3.0 * eorien[4];
		// edge e13
		A.val[10][10] = 1.0/3.0 * eorien[5]; A.val[10][11] = 1.0/6.0 * eorien[5];
		A.val[11][10] = 1.0/6.0 * eorien[5]; A.val[11][11] = 1.0/3.0 * eorien[5];
		
		// face F0
		face = elementFace->val[k][0];
		nv = faces->nvector[face];
		tv[0] = faces->t1vector[face];
		tv[1] = faces->t2vector[face];
		s = faces->area[face];
		axpyz_array(3, -1.0, bcFace[0], vertices[1], v[1]);
		axpyz_array(3, -1.0, bcFace[0], vertices[2], v[2]);
		axpyz_array(3, -1.0, bcFace[0], vertices[3], v[3]);
		cross_array(nv, v[1], nxv[1]);
		cross_array(nv, v[2], nxv[2]);
		cross_array(nv, v[3], nxv[3]);
		axpyz_array(3, 1.0, nxv[1], nxv[2], nxv[0]);
		axpy_array(3, 1.0, nxv[3], nxv[0]);
		axpbyz_array(3, 1.0/4.0, nxv[0], 1.0/4.0, nxv[1], nxvT[1]);
		axpbyz_array(3, 1.0/4.0, nxv[0], 1.0/4.0, nxv[2], nxvT[2]);
		axpbyz_array(3, 1.0/4.0, nxv[0], 1.0/4.0, nxv[3], nxvT[3]);
		axpyz_array(3, 2.0, nxv[0], nxv[1], nxvB[1]);
		axpyz_array(3, 2.0, nxv[0], nxv[2], nxvB[2]);
		axpyz_array(3, 2.0, nxv[0], nxv[3], nxvB[3]);
		// N12, N13
		for(i=0;i<2;i++)
		{
			j = 12 + i;
			t = tv[i];
			c[1] = dot_array(3, t, gradLambda[1]);
			c[2] = dot_array(3, t, gradLambda[2]);
			c[3] = dot_array(3, t, gradLambda[3]);
			A.val[j][2] = c[2]/6 - c[1]/12; 
			A.val[j][3] = c[2]/12 - c[1]/6;
			A.val[j][4] = c[3]/6 - c[2]/12; 
			A.val[j][5] = c[3]/12 - c[2]/6;
			A.val[j][10] = c[3]/6 - c[1]/12; 
			A.val[j][11] = c[3]/12 - c[1]/6;
			A.val[j][12] = c[2]/12 - c[1]/12; 
			A.val[j][13] = c[3]/12 - c[1]/12;
		}
		// N14, N15
		for(i=0;i<2;i++)
		{
			j = 14 + i;
			t = tv[i];
			cc[0][1] = dot_array(3, t, gradCross[0][1]); cc[1][0] = -cc[0][1]; 
			cc[0][2] = dot_array(3, t, gradCross[0][2]); cc[2][0] = -cc[0][2];
			cc[0][3] = dot_array(3, t, gradCross[0][3]); cc[3][0] = -cc[0][3];
			cc[1][2] = dot_array(3, t, gradCross[1][2]); cc[2][1] = -cc[1][2];
			cc[1][3] = dot_array(3, t, gradCross[1][3]); cc[3][1] = -cc[1][3];
			cc[2][3] = dot_array(3, t, gradCross[2][3]); cc[3][2] = -cc[2][3];
			A.val[j][1] = cc[0][1];
			A.val[j][2] = A.val[j][3] = cc[1][2];
			A.val[j][4] = A.val[j][5] = cc[2][3];
			A.val[j][6] = cc[3][0];
			A.val[j][9] = cc[0][2];
			A.val[j][10] = A.val[j][11] = cc[1][3];
			A.val[j][12] = (cc[0][1] - cc[0][2]) / 3.0; 
			A.val[j][13] = (cc[0][1] - cc[0][3]) / 3.0; 
			A.val[j][14] = -A.val[j][13];
			A.val[j][15] = -A.val[j][12];
			A.val[j][16] = A.val[j][12];
			A.val[j][17] = (cc[0][3] - cc[0][2]) / 3.0; 
			A.val[j][18] = -A.val[j][17];
			A.val[j][19] = A.val[j][13];
			cross_array(t, gradLambda[0], val);
			A.val[j][21] = A.val[j][22] = A.val[j][23] = val[0]/180.0;
			A.val[j][25] = A.val[j][26] = A.val[j][27] = val[1]/180.0;
			A.val[j][29] = A.val[j][30] = A.val[j][31] = val[2]/180.0;
		}
		// N16
		A.val[16][1] = dot_array(3, gradCross[0][1], nxvT[1]);
		A.val[16][2] = dot_array(3, gradCross[1][2], nxvT[1]);
		A.val[16][3] = dot_array(3, gradCross[1][2], nxvT[2]);
		A.val[16][4] = dot_array(3, gradCross[2][3], nxvT[2]);
		A.val[16][5] = dot_array(3, gradCross[2][3], nxvT[3]);
		A.val[16][6] = dot_array(3, gradCross[3][0], nxvT[3]);
		A.val[16][9] = dot_array(3, gradCross[0][2], nxvT[2]);
		A.val[16][10] = dot_array(3, gradCross[1][3], nxvT[1]);
		A.val[16][11] = dot_array(3, gradCross[1][3], nxvT[3]);
		A.val[16][12] = dot_array(3, gradCross[1][2], nxvT[3]) * 2.0 / 3.0 - dot_array(3, gradCross[2][3], nxvT[1]) / 3.0 - dot_array(3, gradCross[3][1], nxvT[2]) / 3.0;
		A.val[16][13] = dot_array(3, gradCross[1][3], nxvT[2]) * 2.0 / 3.0 - dot_array(3, gradCross[3][2], nxvT[1]) / 3.0 - dot_array(3, gradCross[2][1], nxvT[3]) / 3.0;
		A.val[16][14] = dot_array(3, gradCross[0][2], nxvT[3]) / 3.0 + dot_array(3, gradCross[0][3], nxvT[2]) * 2.0 / 3.0;
		A.val[16][15] = dot_array(3, gradCross[0][3], nxvT[2]) / 3.0 + dot_array(3, gradCross[0][2], nxvT[3]) * 2.0 / 3.0;
		A.val[16][16] = dot_array(3, gradCross[0][3], nxvT[1]) / 3.0 + dot_array(3, gradCross[0][1], nxvT[3]) * 2.0 / 3.0;
		A.val[16][17] = dot_array(3, gradCross[0][1], nxvT[3]) / 3.0 + dot_array(3, gradCross[0][3], nxvT[1]) * 2.0 / 3.0;
		A.val[16][18] = dot_array(3, gradCross[0][1], nxvT[2]) / 3.0 + dot_array(3, gradCross[0][2], nxvT[1]) * 2.0 / 3.0;
		A.val[16][19] = dot_array(3, gradCross[0][2], nxvT[1]) / 3.0 + dot_array(3, gradCross[0][1], nxvT[2]) * 2.0 / 3.0;
		cross_array(nxvB[1], gradLambda[0], val);
		A.val[16][21] = val[0]/1260.0;
		A.val[16][25] = val[1]/1260.0;
		A.val[16][29] = val[2]/1260.0;
		cross_array(nxvB[2], gradLambda[0], val);
		A.val[16][22] = val[0]/1260.0;
		A.val[16][26] = val[1]/1260.0;
		A.val[16][30] = val[2]/1260.0;
		cross_array(nxvB[3], gradLambda[0], val);
		A.val[16][23] = val[0]/1260.0;
		A.val[16][27] = val[1]/1260.0;
		A.val[16][31] = val[2]/1260.0;
		// scale
		for(i=12;i<17;i++)
			ax_array(32, sqrt(s), A.val[i]);
		
		// face F1
		face = elementFace->val[k][1];
		nv = faces->nvector[face];
		tv[0] = faces->t1vector[face];
		tv[1] = faces->t2vector[face];
		s = faces->area[face];
		axpyz_array(3, -1.0, bcFace[1], vertices[0], v[0]);
		axpyz_array(3, -1.0, bcFace[1], vertices[2], v[2]);
		axpyz_array(3, -1.0, bcFace[1], vertices[3], v[3]);
		cross_array(nv, v[0], nxv[0]);
		cross_array(nv, v[2], nxv[2]);
		cross_array(nv, v[3], nxv[3]);
		axpyz_array(3, 1.0, nxv[0], nxv[2], nxv[1]);
		axpy_array(3, 1.0, nxv[3], nxv[1]);
		axpbyz_array(3, 1.0/4.0, nxv[1], 1.0/4.0, nxv[0], nxvT[0]);
		axpbyz_array(3, 1.0/4.0, nxv[1], 1.0/4.0, nxv[2], nxvT[2]);
		axpbyz_array(3, 1.0/4.0, nxv[1], 1.0/4.0, nxv[3], nxvT[3]);
		axpyz_array(3, 2.0, nxv[1], nxv[0], nxvB[0]);
		axpyz_array(3, 2.0, nxv[1], nxv[2], nxvB[2]);
		axpyz_array(3, 2.0, nxv[1], nxv[3], nxvB[3]);
		// N17, N18
		for(i=0;i<2;i++)
		{
			j = 17 + i;
			t = tv[i];
			c[0] = dot_array(3, t, gradLambda[0]);
			c[2] = dot_array(3, t, gradLambda[2]);
			c[3] = dot_array(3, t, gradLambda[3]);
			A.val[j][4] = c[3]/6 - c[2]/12; 
			A.val[j][5] = c[3]/12 - c[2]/6;
			A.val[j][6] = c[0]/6 - c[3]/12; 
			A.val[j][7] = c[0]/12 - c[3]/6;
			A.val[j][8] = c[2]/6 - c[0]/12; 
			A.val[j][9] = c[2]/12 - c[0]/6;
			A.val[j][14] = c[3]/12 - c[0]/12; 
			A.val[j][15] = c[2]/12 - c[0]/12;
		}
		// N19, N20
		for(i=0;i<2;i++)
		{
			j = 19 + i;
			t = tv[i];
			cc[0][1] = dot_array(3, t, gradCross[0][1]); cc[1][0] = -cc[0][1]; 
			cc[0][2] = dot_array(3, t, gradCross[0][2]); cc[2][0] = -cc[0][2];
			cc[0][3] = dot_array(3, t, gradCross[0][3]); cc[3][0] = -cc[0][3];
			cc[1][2] = dot_array(3, t, gradCross[1][2]); cc[2][1] = -cc[1][2];
			cc[1][3] = dot_array(3, t, gradCross[1][3]); cc[3][1] = -cc[1][3];
			cc[2][3] = dot_array(3, t, gradCross[2][3]); cc[3][2] = -cc[2][3];
			A.val[j][0] = cc[0][1];
			A.val[j][3] = cc[1][2];
			A.val[j][4] = A.val[j][5] = cc[2][3];
			A.val[j][6] = A.val[j][7] = cc[3][0];
			A.val[j][8] = A.val[j][9] = cc[0][2];
			A.val[j][11] = cc[1][3];
			A.val[j][12] = (cc[1][2] - cc[1][0]) / 3.0; 
			A.val[j][13] = (cc[1][3] - cc[1][0]) / 3.0; 
			A.val[j][14] = -A.val[j][13];
			A.val[j][15] = -A.val[j][12];
			A.val[j][16] = A.val[j][12];
			A.val[j][17] = A.val[j][13]; 
			A.val[j][18] = A.val[j][12];
			A.val[j][19] = A.val[j][13];
			cross_array(t, gradLambda[1], val);
			A.val[j][20] = A.val[j][22] = A.val[j][23] = val[0]/180.0;
			A.val[j][24] = A.val[j][26] = A.val[j][27] = val[1]/180.0;
			A.val[j][28] = A.val[j][30] = A.val[j][31] = val[2]/180.0;
		}
		// N21
		A.val[21][0] = dot_array(3, gradCross[0][1], nxvT[0]);
		A.val[21][3] = dot_array(3, gradCross[1][2], nxvT[2]);
		A.val[21][4] = dot_array(3, gradCross[2][3], nxvT[2]);
		A.val[21][5] = dot_array(3, gradCross[2][3], nxvT[3]);
		A.val[21][6] = dot_array(3, gradCross[3][0], nxvT[3]);
		A.val[21][7] = dot_array(3, gradCross[3][0], nxvT[0]);
		A.val[21][8] = dot_array(3, gradCross[0][2], nxvT[0]);
		A.val[21][9] = dot_array(3, gradCross[0][2], nxvT[2]);
		A.val[21][11] = dot_array(3, gradCross[1][3], nxvT[3]);
		A.val[21][12] = dot_array(3, gradCross[1][3], nxvT[2]) / 3.0 + dot_array(3, gradCross[1][2], nxvT[3]) * 2.0 / 3.0;
		A.val[21][13] = dot_array(3, gradCross[1][2], nxvT[3]) / 3.0 + dot_array(3, gradCross[1][3], nxvT[2]) * 2.0 / 3.0;
		A.val[21][14] = dot_array(3, gradCross[0][3], nxvT[2]) * 2.0 / 3.0 - dot_array(3, gradCross[3][2], nxvT[0]) / 3.0 - dot_array(3, gradCross[2][0], nxvT[3]) / 3.0;
		A.val[21][15] = dot_array(3, gradCross[0][2], nxvT[3]) * 2.0 / 3.0 - dot_array(3, gradCross[2][3], nxvT[0]) / 3.0 - dot_array(3, gradCross[3][0], nxvT[2]) / 3.0;
		A.val[21][16] = dot_array(3, gradCross[3][1], nxvT[0]) / 3.0 + dot_array(3, gradCross[0][1], nxvT[3]) * 2.0 / 3.0;
		A.val[21][17] = dot_array(3, gradCross[1][3], nxvT[0]) / 3.0 - dot_array(3, gradCross[1][0], nxvT[3]) / 3.0;
		A.val[21][18] = dot_array(3, gradCross[1][2], nxvT[0]) / 3.0 - dot_array(3, gradCross[1][0], nxvT[2]) / 3.0;
		A.val[21][19] = dot_array(3, gradCross[2][1], nxvT[0]) / 3.0 + dot_array(3, gradCross[0][1], nxvT[2]) * 2.0 / 3.0;
		cross_array(nxvB[0], gradLambda[1], val);
		A.val[21][20] = val[0]/1260.0;
		A.val[21][24] = val[1]/1260.0;
		A.val[21][28] = val[2]/1260.0;
		cross_array(nxvB[2], gradLambda[1], val);
		A.val[21][22] = val[0]/1260.0;
		A.val[21][26] = val[1]/1260.0;
		A.val[21][30] = val[2]/1260.0;
		cross_array(nxvB[3], gradLambda[1], val);
		A.val[21][23] = val[0]/1260.0;
		A.val[21][27] = val[1]/1260.0;
		A.val[21][31] = val[2]/1260.0;
		// scale
		for(i=17;i<22;i++)
			ax_array(32, sqrt(s), A.val[i]);
		
		// face F2
		face = elementFace->val[k][2];
		nv = faces->nvector[face];
		tv[0] = faces->t1vector[face];
		tv[1] = faces->t2vector[face];
		s = faces->area[face];
		axpyz_array(3, -1.0, bcFace[2], vertices[0], v[0]);
		axpyz_array(3, -1.0, bcFace[2], vertices[1], v[1]);
		axpyz_array(3, -1.0, bcFace[2], vertices[3], v[3]);
		cross_array(nv, v[0], nxv[0]);
		cross_array(nv, v[1], nxv[1]);
		cross_array(nv, v[3], nxv[3]);
		axpyz_array(3, 1.0, nxv[0], nxv[1], nxv[2]);
		axpy_array(3, 1.0, nxv[3], nxv[2]);
		axpbyz_array(3, 1.0/4.0, nxv[2], 1.0/4.0, nxv[0], nxvT[0]);
		axpbyz_array(3, 1.0/4.0, nxv[2], 1.0/4.0, nxv[1], nxvT[1]);
		axpbyz_array(3, 1.0/4.0, nxv[2], 1.0/4.0, nxv[3], nxvT[3]);
		axpyz_array(3, 2.0, nxv[2], nxv[0], nxvB[0]);
		axpyz_array(3, 2.0, nxv[2], nxv[1], nxvB[1]);
		axpyz_array(3, 2.0, nxv[2], nxv[3], nxvB[3]);
		// N22, N23
		for(i=0;i<2;i++)
		{
			j = 22 + i;
			t = tv[i];
			c[0] = dot_array(3, t, gradLambda[0]);
			c[1] = dot_array(3, t, gradLambda[1]);
			c[3] = dot_array(3, t, gradLambda[3]);
			A.val[j][0] = c[1]/6 - c[0]/12; 
			A.val[j][1] = c[1]/12 - c[0]/6;
			A.val[j][6] = c[0]/6 - c[3]/12; 
			A.val[j][7] = c[0]/12 - c[3]/6;
			A.val[j][10] = c[3]/6 - c[1]/12; 
			A.val[j][11] = c[3]/12 - c[1]/6;
			A.val[j][16] = c[1]/12 - c[0]/12; 
			A.val[j][17] = c[3]/12 - c[0]/12;
		}
		// N24, N25
		for(i=0;i<2;i++)
		{
			j = 24 + i;
			t = tv[i];
			cc[0][1] = dot_array(3, t, gradCross[0][1]); cc[1][0] = -cc[0][1]; 
			cc[0][2] = dot_array(3, t, gradCross[0][2]); cc[2][0] = -cc[0][2];
			cc[0][3] = dot_array(3, t, gradCross[0][3]); cc[3][0] = -cc[0][3];
			cc[1][2] = dot_array(3, t, gradCross[1][2]); cc[2][1] = -cc[1][2];
			cc[1][3] = dot_array(3, t, gradCross[1][3]); cc[3][1] = -cc[1][3];
			cc[2][3] = dot_array(3, t, gradCross[2][3]); cc[3][2] = -cc[2][3];
			A.val[j][0] = A.val[j][1] = cc[0][1];
			A.val[j][2] = cc[1][2];
			A.val[j][5] = cc[2][3];
			A.val[j][6] = A.val[j][7] = cc[3][0];
			A.val[j][8] = cc[0][2];
			A.val[j][10] = A.val[j][11] = cc[1][3];
			A.val[j][12] = (cc[2][0] - cc[2][1]) / 3.0; 
			A.val[j][13] = (cc[2][3] - cc[2][1]) / 3.0; 
			A.val[j][14] = (cc[2][3] - cc[2][0]) / 3.0;
			A.val[j][15] = -A.val[j][12];
			A.val[j][16] = A.val[j][12];
			A.val[j][17] = -A.val[j][14]; 
			A.val[j][18] = A.val[j][14];
			A.val[j][19] = -A.val[j][12];
			cross_array(t, gradLambda[2], val);
			A.val[j][20] = A.val[j][21] = A.val[j][23] = val[0]/180.0;
			A.val[j][24] = A.val[j][25] = A.val[j][27] = val[1]/180.0;
			A.val[j][28] = A.val[j][29] = A.val[j][31] = val[2]/180.0;
		}
		// N26
		A.val[26][0] = dot_array(3, gradCross[0][1], nxvT[0]);
		A.val[26][1] = dot_array(3, gradCross[0][1], nxvT[1]);
		A.val[26][2] = dot_array(3, gradCross[1][2], nxvT[1]);
		A.val[26][5] = dot_array(3, gradCross[2][3], nxvT[3]);
		A.val[26][6] = dot_array(3, gradCross[3][0], nxvT[3]);
		A.val[26][7] = dot_array(3, gradCross[3][0], nxvT[0]);
		A.val[26][8] = dot_array(3, gradCross[0][2], nxvT[0]);
		A.val[26][10] = dot_array(3, gradCross[1][3], nxvT[1]);
		A.val[26][11] = dot_array(3, gradCross[1][3], nxvT[3]);
		A.val[26][12] = dot_array(3, gradCross[3][2], nxvT[1]) / 3.0 + dot_array(3, gradCross[1][2], nxvT[3]) * 2.0 / 3.0;
		A.val[26][13] = dot_array(3, gradCross[2][3], nxvT[1]) / 3.0 - dot_array(3, gradCross[2][1], nxvT[3]) / 3.0;
		A.val[26][14] = dot_array(3, gradCross[2][3], nxvT[0]) / 3.0 - dot_array(3, gradCross[2][0], nxvT[3]) / 3.0;
		A.val[26][15] = dot_array(3, gradCross[3][2], nxvT[0]) / 3.0 + dot_array(3, gradCross[0][2], nxvT[3]) * 2.0 / 3.0;
		A.val[26][16] = dot_array(3, gradCross[0][1], nxvT[3]) * 2.0 / 3.0 - dot_array(3, gradCross[1][3], nxvT[0]) / 3.0 - dot_array(3, gradCross[3][0], nxvT[1]) / 3.0;
		A.val[26][17] = dot_array(3, gradCross[0][3], nxvT[1]) * 2.0 / 3.0 - dot_array(3, gradCross[3][1], nxvT[0]) / 3.0 - dot_array(3, gradCross[1][0], nxvT[3]) / 3.0;
		A.val[26][18] = dot_array(3, gradCross[1][2], nxvT[0]) / 3.0 + dot_array(3, gradCross[0][2], nxvT[1]) * 2.0 / 3.0;
		A.val[26][19] = dot_array(3, gradCross[2][1], nxvT[0]) / 3.0 - dot_array(3, gradCross[2][0], nxvT[1]) / 3.0;
		cross_array(nxvB[0], gradLambda[2], val);
		A.val[26][20] = val[0]/1260.0;
		A.val[26][24] = val[1]/1260.0;
		A.val[26][28] = val[2]/1260.0;
		cross_array(nxvB[1], gradLambda[2], val);
		A.val[26][21] = val[0]/1260.0;
		A.val[26][25] = val[1]/1260.0;
		A.val[26][29] = val[2]/1260.0;
		cross_array(nxvB[3], gradLambda[2], val);
		A.val[26][23] = val[0]/1260.0;
		A.val[26][27] = val[1]/1260.0;
		A.val[26][31] = val[2]/1260.0;
		// scale
		for(i=22;i<27;i++)
			ax_array(32, sqrt(s), A.val[i]);

		// face F3
		face = elementFace->val[k][3];
		nv = faces->nvector[face];
		tv[0] = faces->t1vector[face];
		tv[1] = faces->t2vector[face];
		s = faces->area[face];
		axpyz_array(3, -1.0, bcFace[3], vertices[0], v[0]);
		axpyz_array(3, -1.0, bcFace[3], vertices[1], v[1]);
		axpyz_array(3, -1.0, bcFace[3], vertices[2], v[2]);
		cross_array(nv, v[0], nxv[0]);
		cross_array(nv, v[1], nxv[1]);
		cross_array(nv, v[2], nxv[2]);
		axpyz_array(3, 1.0, nxv[0], nxv[1], nxv[3]);
		axpy_array(3, 1.0, nxv[2], nxv[3]);
		axpbyz_array(3, 1.0/4.0, nxv[3], 1.0/4.0, nxv[0], nxvT[0]);
		axpbyz_array(3, 1.0/4.0, nxv[3], 1.0/4.0, nxv[1], nxvT[1]);
		axpbyz_array(3, 1.0/4.0, nxv[3], 1.0/4.0, nxv[2], nxvT[2]);
		axpyz_array(3, 2.0, nxv[3], nxv[0], nxvB[0]);
		axpyz_array(3, 2.0, nxv[3], nxv[1], nxvB[1]);
		axpyz_array(3, 2.0, nxv[3], nxv[2], nxvB[2]);
		// N27, N28
		for(i=0;i<2;i++)
		{
			j = 27 + i;
			t = tv[i];
			c[0] = dot_array(3, t, gradLambda[0]);
			c[1] = dot_array(3, t, gradLambda[1]);
			c[2] = dot_array(3, t, gradLambda[2]);
			A.val[j][0] = c[1]/6 - c[0]/12; 
			A.val[j][1] = c[1]/12 - c[0]/6;
			A.val[j][2] = c[2]/6 - c[1]/12; 
			A.val[j][3] = c[2]/12 - c[1]/6;
			A.val[j][8] = c[2]/6 - c[0]/12; 
			A.val[j][9] = c[2]/12 - c[0]/6;
			A.val[j][18] = c[2]/12 - c[0]/12; 
			A.val[j][19] = c[1]/12 - c[0]/12;
		}
		// N29, N30
		for(i=0;i<2;i++)
		{
			j = 29 + i;
			t = tv[i];
			cc[0][1] = dot_array(3, t, gradCross[0][1]); cc[1][0] = -cc[0][1]; 
			cc[0][2] = dot_array(3, t, gradCross[0][2]); cc[2][0] = -cc[0][2];
			cc[0][3] = dot_array(3, t, gradCross[0][3]); cc[3][0] = -cc[0][3];
			cc[1][2] = dot_array(3, t, gradCross[1][2]); cc[2][1] = -cc[1][2];
			cc[1][3] = dot_array(3, t, gradCross[1][3]); cc[3][1] = -cc[1][3];
			cc[2][3] = dot_array(3, t, gradCross[2][3]); cc[3][2] = -cc[2][3];
			A.val[j][0] = A.val[j][1] = cc[0][1];
			A.val[j][2] = A.val[j][3] = cc[1][2];
			A.val[j][4] = cc[2][3];
			A.val[j][7] = cc[3][0];
			A.val[j][8] = A.val[j][9] = cc[0][2];
			A.val[j][10] = cc[1][3];
			A.val[j][12] = (cc[3][2] - cc[3][1]) / 3.0; 
			A.val[j][13] = (cc[3][0] - cc[3][1]) / 3.0; 
			A.val[j][14] = -A.val[j][13];
			A.val[j][15] = (cc[3][2] - cc[3][0]) / 3.0;
			A.val[j][16] = -A.val[j][13];
			A.val[j][17] = A.val[j][15]; 
			A.val[j][18] = -A.val[j][15];
			A.val[j][19] = A.val[j][13];
			cross_array(t, gradLambda[3], val);
			A.val[j][20] = A.val[j][21] = A.val[j][22] = val[0]/180.0;
			A.val[j][24] = A.val[j][25] = A.val[j][26] = val[1]/180.0;
			A.val[j][28] = A.val[j][29] = A.val[j][30] = val[2]/180.0;
		}
		// N31
		A.val[31][0] = dot_array(3, gradCross[0][1], nxvT[0]);
		A.val[31][1] = dot_array(3, gradCross[0][1], nxvT[1]);
		A.val[31][2] = dot_array(3, gradCross[1][2], nxvT[1]);
		A.val[31][3] = dot_array(3, gradCross[1][2], nxvT[2]);
		A.val[31][4] = dot_array(3, gradCross[2][3], nxvT[2]);
		A.val[31][7] = dot_array(3, gradCross[3][0], nxvT[0]);
		A.val[31][8] = dot_array(3, gradCross[0][2], nxvT[0]);
		A.val[31][9] = dot_array(3, gradCross[0][2], nxvT[2]);
		A.val[31][10] = dot_array(3, gradCross[1][3], nxvT[1]);
		A.val[31][12] = dot_array(3, gradCross[3][2], nxvT[1]) / 3.0 - dot_array(3, gradCross[3][1], nxvT[2]) / 3.0;
		A.val[31][13] = dot_array(3, gradCross[2][3], nxvT[1]) / 3.0 + dot_array(3, gradCross[1][3], nxvT[2]) * 2.0 / 3.0;
		A.val[31][14] = dot_array(3, gradCross[2][3], nxvT[0]) / 3.0 + dot_array(3, gradCross[0][3], nxvT[2]) * 2.0 / 3.0;
		A.val[31][15] = dot_array(3, gradCross[3][2], nxvT[0]) / 3.0 - dot_array(3, gradCross[3][0], nxvT[2]) / 3.0;
		A.val[31][16] = dot_array(3, gradCross[3][1], nxvT[0]) / 3.0 - dot_array(3, gradCross[3][0], nxvT[1]) / 3.0;
		A.val[31][17] = dot_array(3, gradCross[1][3], nxvT[0]) / 3.0 + dot_array(3, gradCross[0][3], nxvT[1]) * 2.0 / 3.0;
		A.val[31][18] = dot_array(3, gradCross[0][2], nxvT[1]) * 2.0 / 3.0 - dot_array(3, gradCross[2][1], nxvT[0]) / 3.0 - dot_array(3, gradCross[1][0], nxvT[2]) / 3.0;
		A.val[31][19] = dot_array(3, gradCross[0][1], nxvT[2]) * 2.0 / 3.0 - dot_array(3, gradCross[1][2], nxvT[0]) / 3.0 - dot_array(3, gradCross[2][0], nxvT[1]) / 3.0;
		cross_array(nxvB[0], gradLambda[3], val);
		A.val[31][20] = val[0]/1260.0;
		A.val[31][24] = val[1]/1260.0;
		A.val[31][28] = val[2]/1260.0;
		cross_array(nxvB[1], gradLambda[3], val);
		A.val[31][21] = val[0]/1260.0;
		A.val[31][25] = val[1]/1260.0;
		A.val[31][29] = val[2]/1260.0;
		cross_array(nxvB[2], gradLambda[3], val);
		A.val[31][22] = val[0]/1260.0;
		A.val[31][26] = val[1]/1260.0;
		A.val[31][30] = val[2]/1260.0;
		// scale
		for(i=27;i<32;i++)
			ax_array(32, sqrt(s), A.val[i]);
		
		init_dden_matrix(&B, 0.0);
		for(i=0;i<B.row;i++)
			B.val[i][i]=1.0;

/*********************************************************************************************
		FILE *outputFile;
		if(k==0)
		{
			outputFile=fopen("output/basisCoeffA.dat", "w");
			for(i=0;i<A.row;i++)
			{
				for(j=0;j<A.col-1;j++)
					fprintf(outputFile,"%.16lf\t",A.val[i][j]);
				fprintf(outputFile,"%.16lf\n",A.val[i][j]);
			}
			fclose(outputFile);			
		}
		********************************************************************************************/

		AxBrref(&A, &B);
		C=basisCoeffs->val[k];
		for(i=0;i<basisCoeffs->row;i++)
		{
			for(j=0;j<basisCoeffs->col;j++)
				C[i][j]=B.val[j][i];
		}
		/*********************************************************************************************
		// FILE *outputFile;
		if(k==0)
		{
			outputFile=fopen("output/basisCoeffs.dat", "w");
			for(i=0;i<basisCoeffs->row;i++)
			{
				for(j=0;j<basisCoeffs->col-1;j++)
					fprintf(outputFile,"%.16lf\t",basisCoeffs->val[k][i][j]);
				fprintf(outputFile,"%.16lf\n",basisCoeffs->val[k][i][j]);
			}
			fclose(outputFile);
		}
		********************************************************************************************/
	} // k
	
	free_dden_matrix(&A);
	free_dden_matrix(&B);
}
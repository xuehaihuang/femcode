/*
 *  assemble.c
 *  DGFEM
 *
 *  Created by Xuehai Huang on 10/30/08.
 *  Modified by Xuehai Huang on 1/11/2015.
 *  Copyright 2015 WZU. All rights reserved.
 *
 */

/*! \file assemble.c
 *  \brief Assembling for stiffness matrix
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "header.h"
#include "matvec.h"

/**
 * \fn int getmesh(int domain_num, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, int levelNum)
 * \brief generate mesh information
 * \param domain_num number of domain
 * \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices, the last 3 columns store the indexes of edge's midpoints
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                       the fourth column stores -1 if the edge is on boundary
 * \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
                                       the third column stores the Dirichlet status of the points(0: nondirichlet, -1: dirichlet)
 * \param levelNum total level number of grid
 * \return 1 if succeed 0 if fail
 */
int getmesh(int domain_num, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, int levelNum)
{
	iCSRmat elementsTran[levelNum]; // the transpose of elements. 
	
	int IsExist = getCoarseInfo(domain_num, &nodes[0], &elements[0], &faces[0], &edges[0], &elementsTran[0], &elementEdge[0]);
	if(IsExist==0)
	{
		printf("Constructing coarse grid fails!\n");
		return 0;
	}

	int i,j,k,l,m;
	
	for(i=0;i<levelNum-1;i++)
		uniformrefine(&nodes[i], &elements[i], &faces[i], &edges[i], &elementsTran[i], &nodes[i+1], &elements[i+1], &faces[i+1], &edges[i + 1], &elementsTran[i+1], &elementEdge[i + 1]);

	for (l = 0; l<levelNum; l++)
	{
		free(elementsTran[l].IA);
		free(elementsTran[l].JA);
	}

	int L=levelNum-1;

	// generate elementFace			
	int point[3], element1, element2;
	for (m = 0; m < levelNum; m++)
	{
		create_iden_matrix(elements[m].row, 4, elementFace + m);
		for (i = 0; i < faces[m].row; i++)
		{
			point[0] = faces[m].val[i][0];
			point[1] = faces[m].val[i][1];
			point[2] = faces[m].val[i][2];
			element1 = faces[m].val[i][3];
			element2 = faces[m].val[i][4];

			for (j = 0; j < 4; j++)
			{
				if (elements[m].val[element1][j] == point[0])
					break;
			}
			for (k = 0; k < 4; k++)
			{
				if (elements[m].val[element1][k] == point[1])
					break;
			}
			for (l = 0; l < 4; l++)
			{
				if (elements[m].val[element1][l] == point[2])
					break;
			}
			elementFace[m].val[element1][6 - j - k - l] = i;

			if (element2 > -1)
			{
				for (j = 0; j < 4; j++)
				{
					if (elements[m].val[element2][j] == point[0])
						break;
				}
				for (k = 0; k < 4; k++)
				{
					if (elements[m].val[element2][k] == point[1])
						break;
				}
				for (l = 0; l < 4; l++)
				{
					if (elements[m].val[element2][l] == point[2])
						break;
				}
				elementFace[m].val[element2][6 - j - k - l] = i;
			}
		}
	}
	// end generate elementFace
		
	return 1;
}

/**
* \fn void getElementDOFdg(ELEMENT_DOF *elementDOF, int nt, int ldof, int dop)
 * \brief get the degrees of freedom of piecewise element
 * \param nt number of elements
 * \param ldof number of local degrees of freedom
 * \param dop degree of polynomial
*/
void getElementDOFdg(ELEMENT_DOF *elementDOF, int nt, int ldof, int dop)
{
	int i,k;
	int count;

	create_elementDOF(dop, nt * ldof, nt, ldof, elementDOF);

	for(k=0;k<nt;k++){
		count=k*elementDOF->col;
		for(i=0;i<elementDOF->col;i++)
			elementDOF->val[k][i]=count+i;
	}
}

/**
 * \fn void getElementDOF1d_Continue(ELEMENT_DOF *edgeDOF, ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
 * \brief get the degrees of freedom of Lagrange element on edges
 * \param *edgeDOF pointer to relation between edges and DOFs
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param nvertices number of vertices
 * \param dop degree of polynomial
 */
void getElementDOF1d_Continue(ELEMENT_DOF *edgeDOF, ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
{
	int i,j,k;
	int nt=elements->row;
	int ne=edges->row;
	int nn=nvertices;
	if(dop<2)
		dop=1;

	create_elementDOF(dop, nn + ne*(dop-1), ne, dop+1, edgeDOF);

	int node, edge;
	for(j=0;j<ne;j++)
	{
		for(i=0;i<2;i++)
		{
			node=edges->val[j][i];
			edgeDOF->val[j][i]=node;
		}
		
		for(i=0;i<dop-1;i++)
			edgeDOF->val[j][2+i] = nn + j*(dop-1) + i;
	}

	create_elementDOF(dop, nn + ne*(dop-1), nt, dop*3, elementDOF);

	int orient;
	for(k=0;k<nt;k++)
	{
		for(i=0;i<3;i++)
		{
			node=elements->val[k][i];
			elementDOF->val[k][i]=node;
		}
		
		for(j=0;j<3;j++)
		{
			edge=elementEdge->val[k][j];
			if(elements->val[k][(j+1)%3]== edges->val[edge][0])
				orient=1;
			else
				orient=0;
			
			if(orient==1)
			{
				for(i=0;i<dop-1;i++)
					elementDOF->val[k][3+(dop-1)*j+i] = nn + edge*(dop-1) + i;
			}
			else
			{
				for(i=0;i<dop-1;i++)
					elementDOF->val[k][3+(dop-1)*j+dop-2-i] = nn + edge*(dop-1) + i;
			}
		}
	}
}

/**
 * \fn void getElementDOF_Lagrange2d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, int nedges, int nvertices, int dop)
 * \brief get the degrees of freedom of Lagrange element
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertices, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param nvertices number of vertices
 * \param dop degree of polynomial
 */
void getElementDOF_Lagrange2d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
{
	int i,j,k;
	int nt=elements->row;
	int ne=edges->row;
	int nn=nvertices;
	if(dop<2)
		dop=1;

	create_elementDOF(dop, nn + ne*(dop-1) + nt*(dop-1)*(dop-2)/2, nt, (dop+1)*(dop+2)/2, elementDOF);

	int node, edge;
	int orient;
	for(k=0;k<nt;k++)
	{
		for(i=0;i<3;i++)
		{
			node=elements->val[k][i];
			elementDOF->val[k][i]=node;
		}
		
		for(j=0;j<3;j++)
		{
			edge=elementEdge->val[k][j];
			if(elements->val[k][(j+1)%3]== edges->val[edge][0])
				orient=1;
			else
				orient=0;
			
			if(orient==1)
			{
				for(i=0;i<dop-1;i++)
					elementDOF->val[k][3+(dop-1)*j+i] = nn + edge*(dop-1) + i;
			}
			else
			{
				for(i=0;i<dop-1;i++)
					elementDOF->val[k][3+(dop-1)*j+dop-2-i] = nn + edge*(dop-1) + i;
			}
		}

		for(i=0;i<(dop-1)*(dop-2)/2;i++)
		{
			elementDOF->val[k][3*dop+i] = nn + ne*(dop-1) + k*(dop-1)*(dop-2)/2 + i;
		}
	}
}

/**
 * \fn void getElementDOF_Lagrange3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
 * \brief get the degrees of freedom of Lagrange element in three dimensions
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *elements pointer to triangulation: the first 4 columns store the indexes of vertices
 * \param *elementFace pointer to relation between tetrahedrons and faces: each row stores 4 faces index
 * \param *faces pointer to faces: the first three columns store the three vertices, the fourth and fifth columns store the affiliated elements
                                   the fifth column stores -1 if the face is on boundary
 * \param *elementEdge pointer to relation between tetrahedrons and edges: each row stores 6 edges index
 * \param *edges pointer to edges: store the two vertice
 * \param nvertices number of vertices
 * \param dop degree of polynomial
 */
// void getElementDOF_Lagrange3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
void getElementDOF_Lagrange3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, idenmat *elementEdge, int nfaces, int nedges, int nvertices, int dop)
{
	int i,j,k,l,ii[3];
	int nt=elements->row;
	int nf=nfaces;
	int ne=nedges;
	int nn=nvertices;
	if(dop<2)
		dop=1;

	create_elementDOF(dop, nn + ne*(dop-1) + nf*(dop-1)*(dop-2)/2 + nt*(dop-1)*(dop-2)*(dop-3)/6, nt, (dop+1)*(dop+2)*(dop+3)/6, elementDOF);
	
	int node, edge, face, v[3];
	int *perm;
	// int localFaces[4][3];
	for(k=0;k<nt;k++){
		// for(i=0;i<4;i++)
		// {
		// 	face2vertices3d(i, v);
		// 	for(j=0;j<3;j++)
		// 		localFaces[i][j]=elements->val[k][v[j]];
		// }
		
		copy_arrayint(4, elements->val[k], elementDOF->val[k]);
		// for(i=0;i<4;i++)
		// {
		// 	node=elements->val[k][i];
		// 	elementDOF->val[k][i]=node;
		// }
		
		for(j=0;j<6;j++){
			// edge2vv3d(j, v);
			// edge=elementEdge->val[k][j];
			// if(elements->val[k][v[0]] == edges->val[edge][0])
			// {
			// 	perm[0]=0;
			// 	perm[1]=1;
			// }
			// else
			// {
			// 	perm[0]=1;
			// 	perm[1]=0;
			// }
			edge = elementEdge->val[k][j];
			perm = elements->eperm[k][j];

			for(ii[0]=0;ii[0]<dop-1;ii[0]++){
				ii[1]=dop-2-ii[0];
				elementDOF->val[k][4+(dop-1)*j+ ii[0]] = nn + edge*(dop-1) + ii[perm[0]];
			}
		}

		for(j=0;j<4;j++){
			face = elementFace->val[k][j];
			perm = elements->fperm[k][j];
			// getPermutation(localFaces[j], faces->val[face], perm, 3);

			if(dop==3)
				elementDOF->val[k][4+(dop-1)*6+(dop-1)*(dop-2)/2*j] = nn + ne*(dop-1) + face*(dop-1)*(dop-2)/2;

			// if(dop==4)
			// {
			// 	for(i=0;i<3;i++)
			// 		elementDOF->val[k][4+(dop-1)*6+(dop-1)*(dop-2)/2*j+ i] = nn + ne*(dop-1) + face*(dop-1)*(dop-2)/2 + perm[i];
			// }

			if (dop > 3){
				int curleft, curright, rotl, roti;
				curleft = 4 + (dop - 1) * 6 + (dop - 1) * (dop - 2) / 2 * j;
				curright = nn + ne * (dop - 1) + face * (dop - 1) * (dop - 2) / 2;
				for (l = 0; l <= dop - 3; l++)
				{
					for (i = 0; i <= l; i++)
					{
						v[perm[0]] = dop - 3 - l;
						v[perm[1]] = l - i;
						v[perm[2]] = i;
						rotl = dop - 3 - v[0];
						roti = v[2];
						elementDOF->val[k][curleft + (l + 1) * l / 2 + i] = curright + (rotl + 1) * rotl / 2 + roti;
						// elementDOF->val[k][curleft + (v[0] + 1) * v[0] / 2 + v[1]] = curright + (l + 1) * l / 2 + i;
					}
				}
			}
		}

		for (i = 0; i < (dop - 1) * (dop - 2) * (dop - 3) / 6; i++){
			elementDOF->val[k][2 * (dop * dop + 1) + i] = nn + ne * (dop - 1) + nf * (dop - 1) * (dop - 2) / 2 + k * (dop - 1) * (dop - 2) * (dop - 3) / 6 + i;
		}
	}
}

/**
* \fn void getElementDOF_NoncfmP13d(ELEMENT_DOF *elementDOF, idenmat *elementFace, int nf)
* \brief get the degrees of freedom of nonconforming P1 element in three dimensions
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param ne number of edges
*/
void getElementDOF_NoncfmP13d(ELEMENT_DOF *elementDOF, idenmat *elementFace, int nf)
{
	int i, k;
	int nt = elementFace->row;

	create_elementDOF(1, nf, nt, 4, elementDOF);

	for (k = 0; k<nt; k++){
		for (i = 0; i<4; i++)
			elementDOF->val[k][i] = elementFace->val[k][i];
	}
}

/**
* \fn void getElementDOF_Morley3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges)
* \brief get the degrees of freedom of Morley element in three dimensions
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elements pointer to triangulation: the first 4 columns store the indexes of vertices
* \param *elementFace pointer to relation between tetrahedrons and faces: each row stores 4 faces index
* \param *faces pointer to faces: the first three columns store the three vertices, the fourth and fifth columns store the affiliated elements
                                  the fifth column stores -1 if the face is on boundary
* \param *elementEdge pointer to relation between tetrahedrons and edges: each row stores 6 edges index
* \param *edges pointer to edges: store the two vertice
*/
void getElementDOF_Morley3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges)
{
	int i, j, k;
	int nt = elements->row;
	int nf = faces->row;
	int ne = edges->row;

	create_elementDOF(2, ne + nf, nt, 10, elementDOF);

	for (k = 0; k<nt; k++)
	{
		for (j = 0; j<6; j++)
			elementDOF->val[k][j] = elementEdge->val[k][j];

		for (j = 0; j<4; j++)
			elementDOF->val[k][6 + j] = ne + elementFace->val[k][j];
	}
}

/**
 * \fn void getElementDOF_HuZhang(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
 * \brief get the degrees of freedom of Hu-Zhang element
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param nvertices number of vertices
 * \param dop degree of polynomial
 */
void getElementDOF_HuZhang(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
{
	int i,j,k;
	int nt=elements->row;
	int ne=edges->row;
	int nn=nvertices;
	if(dop<2)
		dop=1;

	create_elementDOF(dop, nn*3 + ne*(dop-1)*2 + nt*((dop-1)*3 + (dop-1)*(dop-2)/2*3), nt, (dop+1)*(dop+2)/2*3, elementDOF);

	int node, edge;
	int orient;
	for(k=0;k<nt;k++)
	{
		for(i=0;i<3;i++)
		{
			node=elements->val[k][i];
			elementDOF->val[k][i]=node;
			elementDOF->val[k][i+3]=node+nn;
			elementDOF->val[k][i+6]=node+nn*2;
		}
		
		for(j=0;j<3;j++)
		{
			edge=elementEdge->val[k][j];
			if(elements->val[k][(j+1)%3]== edges->val[edge][0])
				orient=1;
			else
				orient=0;
			
			if(orient==1)
			{
				for(i=0;i<dop-1;i++)
				{
					elementDOF->val[k][9+(dop-1)*j+i] = nn*3 + edge*(dop-1) + i;
					elementDOF->val[k][9+(dop-1)*3+(dop-1)*j+i] = nn*3 + ne*(dop-1) + edge*(dop-1) + i;
				}
			}
			else
			{
				for(i=0;i<dop-1;i++)
				{
					elementDOF->val[k][9+(dop-1)*j+dop-2-i] = nn*3 + edge*(dop-1) + i;
					elementDOF->val[k][9+(dop-1)*3+(dop-1)*j+dop-2-i] = nn*3 + ne*(dop-1) + edge*(dop-1) + i;
				}
			}

			for(i=0;i<dop-1;i++)
				elementDOF->val[k][9+(dop-1)*6+(dop-1)*j+i] = nn*3 + ne*(dop-1)*2 + k*((dop-1)*3 + (dop-1)*(dop-2)/2*3) + (dop-1)*j + i;
	//			elementDOF->val[k][9+(dop-1)*6+(dop-1)*j+i] = nn*3 + ne*(dop-1)*2 + k*(dop-1)*3 + (dop-1)*j + i;
		} // j

		for(i=0;i<(dop-1)*(dop-2)/2*3;i++)
			elementDOF->val[k][9*dop+i] = nn*3 + ne*(dop-1)*2 + k*((dop-1)*3 + (dop-1)*(dop-2)/2*3) + (dop-1)*3 + i;
	
/*		for(i=0;i<(dop-1)*(dop-2)/2;i++)
		{
			elementDOF->val[k][9*dop+i] = nn*3 + ne*(dop-1)*2 + nt*(dop-1)*3 + k*(dop-1)*(dop-2)/2 + i;
			elementDOF->val[k][9*dop+(dop-1)*(dop-2)/2+i] = nn*3 + ne*(dop-1)*2 + nt*((dop-1)*3 + (dop-1)*(dop-2)/2) + k*(dop-1)*(dop-2)/2 + i;
			elementDOF->val[k][9*dop+(dop-1)*(dop-2)+i] = nn*3 + ne*(dop-1)*2 + nt*((dop-1)*3 + (dop-1)*(dop-2)) + k*(dop-1)*(dop-2)/2 + i;
		}*/
	} // k
}

/**
 * \fn void getElementDOF_HuZhang3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
 * \brief get the degrees of freedom of Hu-Zhang element in three dimensions
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *elements pointer to triangulation: the first 4 columns store the indexes of vertices
 * \param *elementFace pointer to relation between tetrahedrons and faces: each row stores 4 faces index
 * \param *faces pointer to faces: the first three columns store the three vertices, the fourth and fifth columns store the affiliated elements
                                   the fifth column stores -1 if the face is on boundary
 * \param *elementEdge pointer to relation between tetrahedrons and edges: each row stores 6 edges index
 * \param *edges pointer to edges: store the two vertice
 * \param nvertices number of vertices
 * \param dop degree of polynomial
 */
void getElementDOF_HuZhang3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
{
	int i,j,k,ii[3];
	int nt=elements->row;
	int nf=faces->row;
	int ne=edges->row;
	int nn=nvertices;
	if(dop<2)
		dop=1;

	// (dop-1)*6 + (dop-1)*(dop-2)/2*3*4 + (dop-1)*(dop-2)*(dop-3) = (dop*dop-1)*dop
	create_elementDOF(dop, nn*6 + ne*(dop-1)*5 + nf*(dop-1)*(dop-2)/2*3 + nt*(dop*dop-1)*dop, nt, (dop+1)*(dop+2)*(dop+3), elementDOF);

	int node, edge, face, v[3];
	int perm[3];
	int localFaces[4][3];
	for(k=0;k<nt;k++)
	{
		for(i=0;i<4;i++)
		{
			face2vertices3d(i, v);
			for(j=0;j<3;j++)
				localFaces[i][j]=elements->val[k][v[j]];
		}

		for(i=0;i<4;i++)
		{
			node=elements->val[k][i];
			for(j=0;j<6;j++)
				elementDOF->val[k][i+4*j]=node+nn*j;
		}

		for(j=0;j<6;j++)
		{
			edge2vv3d(j, v);
			edge=elementEdge->val[k][j];
			if(elements->val[k][v[0]] == edges->val[edge][0])
			{
				perm[0]=0;
				perm[1]=1;
			}
			else
			{
				perm[0]=1;
				perm[1]=0;
			}

			for(ii[0]=0;ii[0]<dop-1;ii[0]++)
			{
				ii[1]=dop-2-ii[0];
				for(i=0;i<5;i++)
					elementDOF->val[k][24+(dop-1)*j+ii[0] + 6*(dop-1)*i] = nn*6+edge*(dop-1)+ii[perm[0]] + ne*(dop-1)*i;
			}
		}

		for(j=0;j<4;j++)
		{
			face=elementFace->val[k][j];
			getPermutation(localFaces[j], faces->val[face], perm, 3);

			if(dop==3)
			{
				for(i=0;i<3;i++)
					elementDOF->val[k][24+(dop-1)*30+(dop-1)*(dop-2)/2*j + 2*(dop-1)*(dop-2)*i] = nn*6+ne*(dop-1)*5+face*(dop-1)*(dop-2)/2 + nf*(dop-1)*(dop-2)/2*i;
			}

			if(dop==4)
			{
				for(ii[0]=0;ii[0]<3;ii[0]++)
				{
					for(i=0;i<3;i++)
						elementDOF->val[k][24+(dop-1)*30+(dop-1)*(dop-2)/2*j+ii[0] + 2*(dop-1)*(dop-2)*i] = nn*6+ne*(dop-1)*5+face*(dop-1)*(dop-2)/2+perm[ii[0]] + nf*(dop-1)*(dop-2)/2*i;
				}
			}
		}

		for(i=0;i<(dop*dop-1)*dop;i++)
			elementDOF->val[k][24+(dop-1)*30+(dop-1)*(dop-2)*6+i] = nn*6 + ne*(dop-1)*5 + nf*(dop-1)*(dop-2)/2*3 + k*(dop*dop-1)*dop + i;

	} // k
}

/**
* \fn void getElementDOF_Nedelec1st3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, int dop)
* \brief get the degrees of freedom of the first kind Nedelec element in three dimensions
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elements pointer to triangulation: the first 4 columns store the indexes of vertices
* \param *elementFace pointer to relation between tetrahedrons and faces: each row stores 4 faces index
* \param *faces pointer to faces: the first three columns store the three vertices, the fourth and fifth columns store the affiliated elements
								  the fifth column stores -1 if the face is on boundary
* \param *elementEdge pointer to relation between tetrahedrons and edges: each row stores 6 edges index
* \param *edges pointer to edges: store the two vertice
* \param nvertices number of vertices
* \param dop degree of polynomial
*/
void getElementDOF_Nedelec1st3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, int dop)
{
	int i, j, k, l, ii[3];
	int curleft, curright, rotl, roti, step;
	int nt = elements->row;
	int nf = faces->row;
	int ne = edges->row;
//	int nn = nvertices;
	if (dop<2)
		dop = 1;

	create_elementDOF(dop, ne*dop + nf*(dop - 1)*dop + nt*(dop - 2)*(dop - 1)*dop / 2, nt, dop*(dop + 2)*(dop + 3) / 2, elementDOF);

	int node, edge, face, v[3];
	int *perm;
	// int localFaces[4][3];
	for (k = 0; k<nt; k++)
	{
		for (j = 0; j<6; j++)
		{
			edge = elementEdge->val[k][j];
			perm = elements->eperm[k][j];

			for (ii[0] = 0; ii[0]<dop; ii[0]++)
			{
				ii[1] = dop - 1 - ii[0];
				elementDOF->val[k][dop*j + ii[0]] = edge*dop + ii[perm[0]];
			}
		}

		if (dop < 2)
			continue;

		for (j = 0; j<4; j++)
		{
			face = elementFace->val[k][j];
			perm = elements->fperm[k][j];

			if (dop == 2)
			{
				for (i = 0; i<(dop - 1)*dop; i++)
				{
					elementDOF->val[k][dop * 6 + (dop - 1)*dop * j + i] = ne*dop + face*(dop - 1)*dop + i; //perm[i];
				}
			}

			if (dop == 3) // dop > 2 
			{
				curleft = dop * 6 + (dop - 1)*dop * j;
				curright = dop * ne + (dop - 1)*dop * face;
				step = (dop - 1)*dop/2;
				for (l = 0; l <= dop - 2; l++)
				{
					for (i = 0; i <= l; i++)
					{
						v[perm[0]] = dop - 2 - l;
						v[perm[1]] = l - i;
						v[perm[2]] = i;
						rotl = dop - 2 - v[0];
						roti = v[2];
						elementDOF->val[k][curleft + (l + 1) * l / 2 + i] = curright + (rotl + 1) * rotl / 2 + roti;
						elementDOF->val[k][curleft + step + (l + 1) * l / 2 + i] = curright + step + (rotl + 1) * rotl / 2 + roti;
					}
				}
			} // dop = 3
		} // j

		curleft = dop * 6 + (dop - 1)*dop * 4;
		curright = ne*dop + nf*(dop - 1)*dop + k*(dop - 2)*(dop - 1)*dop / 2;
		step = (dop - 2)*(dop - 1)*dop / 6;
		for (i = 0; i < (dop - 2)*(dop - 1)*dop / 2; i++)
		{
			elementDOF->val[k][curleft + i] = curright + i;
		}
	} // k
}

/**
* \fn void getElementDOF_Nedelec2nd3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, int dop)
* \brief get the degrees of freedom of the first kind Nedelec element in three dimensions
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elements pointer to triangulation: the first 4 columns store the indexes of vertices
* \param *elementFace pointer to relation between tetrahedrons and faces: each row stores 4 faces index
* \param *faces pointer to faces: the first three columns store the three vertices, the fourth and fifth columns store the affiliated elements
								  the fifth column stores -1 if the face is on boundary
* \param *elementEdge pointer to relation between tetrahedrons and edges: each row stores 6 edges index
* \param *edges pointer to edges: store the two vertice
* \param nvertices number of vertices
* \param dop degree of polynomial
*/
void getElementDOF_Nedelec2nd3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, int dop)
{
	int i, j, k, ii[3];
	int nt = elements->row;
	int nf = faces->row;
	int ne = edges->row;
	//	int nn = nvertices;
	if (dop<2)
		dop = 1;

	create_elementDOF(dop, ne*(dop + 1) + nf*(dop - 1)*(dop + 1) + nt*(dop - 2)*(dop - 1)*(dop + 1) / 2, nt, (dop + 1)*(dop + 2)*(dop + 3) / 2, elementDOF);

	int node, edge, face, v[3];
	int *perm, permp[2];
	// int localFaces[4][3];
	for (k = 0; k<nt; k++)
	{
		// for (i = 0; i<4; i++)
		// {
		// 	face2vertices3d(i, v);
		// 	for (j = 0; j<3; j++)
		// 		localFaces[i][j] = elements->val[k][v[j]];
		// }

		for (j = 0; j<6; j++)
		{
/*			edge = elementEdge->val[k][j];
			perm = elements->eperm[k][j];*/
			edge2vv3d(j, v);
			edge = elementEdge->val[k][j];
			if (elements->val[k][v[0]] == edges->val[edge][0])
			{
				permp[0] = 0;
				permp[1] = 1;
			}
			else
			{
				permp[0] = 1;
				permp[1] = 0;
			}
			perm = permp;

			for (ii[0] = 0; ii[0]<dop+1; ii[0]++)
			{
				ii[1] = dop - ii[0];
				elementDOF->val[k][(dop + 1)*j + ii[0]] = edge*(dop + 1) + ii[perm[0]];
			}
		}

		if (dop < 2)
			continue;

		for (j = 0; j<4; j++)
		{
			face = elementFace->val[k][j];
			perm = elements->fperm[k][j];

			if (dop == 2)
			{
				for (i = 0; i<(dop - 1)*(dop + 1); i++)
				{
					elementDOF->val[k][(dop + 1) * 6 + (dop - 1)*(dop + 1) * j + i] = ne*(dop + 1) + face*(dop - 1)*(dop + 1) + perm[i];
				}
			}
		}

		for (i = 0; i < (dop - 2)*(dop - 1)*(dop + 1) / 2; i++)
			elementDOF->val[k][(dop + 1) * 6 + (dop - 1)*(dop + 1) * 4 + i] = ne*(dop + 1) + nf*(dop - 1)*(dop + 1) + k*(dop - 2)*(dop - 1)*(dop + 1) / 2 + i;
	} // k
}

/**
* \fn void getElementDOF_CHHcurlHermite3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
* \brief get the degrees of freedom of Hermite-type Christiansen-Hu-Hu H(curl) element in three dimensions
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elements pointer to triangulation: the first 4 columns store the indexes of vertices
* \param *elementFace pointer to relation between tetrahedrons and faces: each row stores 4 faces index
* \param *faces pointer to faces: the first three columns store the three vertices, the fourth and fifth columns store the affiliated elements
								  the fifth column stores -1 if the face is on boundary
* \param *elementEdge pointer to relation between tetrahedrons and edges: each row stores 6 edges index
* \param *edges pointer to edges: store the two vertice
* \param nvertices number of vertices
* \param dop degree of polynomial
*/
void getElementDOF_CHHcurlHermite3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
{
	int i, j, k, ii[3];
	int nt = elements->row;
	int nf = faces->row;
	int ne = edges->row;
	int nn = nvertices;
	if (dop<2)
		dop = 1;

	create_elementDOF(dop, nn * 3 + ne*(dop - 1) + nf*(dop - 1)*(dop + 1) + nt*(dop - 2)*(dop - 1)*(dop + 1) / 2, nt, (dop + 1)*(dop + 2)*(dop + 3) / 2, elementDOF);

	int node, edge, face, v[3];
	int perm[3];
	int localFaces[4][3];
	for (k = 0; k<nt; k++)
	{
		for (i = 0; i<4; i++)
		{
			face2vertices3d(i, v);
			for (j = 0; j<3; j++)
				localFaces[i][j] = elements->val[k][v[j]];
		}

		for (i = 0; i<4; i++)
		{
			node = elements->val[k][i];
			for (j = 0; j<3; j++)
				elementDOF->val[k][i + 4 * j] = node + nn*j;
		}

		for (j = 0; j<6; j++)
		{
			edge2vv3d(j, v);
			edge = elementEdge->val[k][j];
			if (elements->val[k][v[0]] == edges->val[edge][0])
			{
				perm[0] = 0;
				perm[1] = 1;
			}
			else
			{
				perm[0] = 1;
				perm[1] = 0;
			}

			for (ii[0] = 0; ii[0]<dop - 1; ii[0]++)
			{
				ii[1] = dop - 2 - ii[0];
				elementDOF->val[k][12 + (dop - 1)*j + ii[0]] = nn * 3 + edge*(dop - 1) + ii[perm[0]];
			}
		}

		for (j = 0; j<4; j++)
		{
			face = elementFace->val[k][j];
			getPermutation(localFaces[j], faces->val[face], perm, 3);

			if (dop == 2)
			{
				for (ii[0] = 0; ii[0]<(dop - 1)*(dop + 1); ii[0]++)
				{
					elementDOF->val[k][12 + (dop - 1) * 6 + (dop - 1)*(dop + 1) * j + ii[0]] = nn * 3 + ne*(dop - 1) + face*(dop - 1)*(dop + 1) + perm[ii[0]];
				}
			}
		}

		for (i = 0; i < (dop - 2)*(dop - 1)*(dop + 1) / 2; i++)
			elementDOF->val[k][12 + (dop - 1) * 6 + (dop - 1)*(dop + 1) * 4 + i] = nn * 3 + ne*(dop - 1) + nf*(dop - 1)*(dop + 1) + k*(dop - 2)*(dop - 1)*(dop + 1) / 2 + i;
	} // k
}

/**
* \fn void getElementDOF_HuangGradcurl3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges)
* \brief get the degrees of freedom of the Huang element in three dimensions for grad-curl
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elements pointer to triangulation: the first 4 columns store the indexes of vertices
* \param *elementFace pointer to relation between tetrahedrons and faces: each row stores 4 faces index
* \param *faces pointer to faces: the first three columns store the three vertices, the fourth and fifth columns store the affiliated elements
								  the fifth column stores -1 if the face is on boundary
* \param *elementEdge pointer to relation between tetrahedrons and edges: each row stores 6 edges index
* \param *edges pointer to edges: store the two vertice
* \param nvertices number of vertices
* \param dop degree of polynomial
*/
void getElementDOF_HuangGradcurl3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges)
{
	int i, j, k;
	int nt = elements->row;
	int nf = faces->row;
	int ne = edges->row;
//	int nn = nvertices;

	create_elementDOF(2, ne + nf*2, nt, 14, elementDOF);

	int edge, face;
	for (k = 0; k<nt; k++)
	{
		for (j = 0; j<6; j++)
		{
			edge = elementEdge->val[k][j];
			elementDOF->val[k][j] = edge;
		}

		for (j = 0; j<4; j++)
		{
			face = elementFace->val[k][j];
			for (i = 0; i<2; i++)
				elementDOF->val[k][6 + 2 * j + i] = ne + face*2 + i;
		}
	} // k
}

/**
* \fn void getElementDOF_HuangZhang3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges)
* \brief get the degrees of freedom of the Huang-Zhang element in three dimensions for grad-curl
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elements pointer to triangulation: the first 4 columns store the indexes of vertices
* \param *elementFace pointer to relation between tetrahedrons and faces: each row stores 4 faces index
* \param *faces pointer to faces: the first three columns store the three vertices, the fourth and fifth columns store the affiliated elements
								  the fifth column stores -1 if the face is on boundary
* \param *elementEdge pointer to relation between tetrahedrons and edges: each row stores 6 edges index
* \param *edges pointer to edges: store the two vertice
* \param nvertices number of vertices
* \param dop degree of polynomial
*/
void getElementDOF_HuangZhang3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges)
{
	int i, j, k, ii[3];
	int nt = elements->row;
	int nf = faces->row;
	int ne = edges->row;
//	int nn = nvertices;
	int *perm;

	create_elementDOF(5, ne*2 + nf*5, nt, 32, elementDOF);

	int edge, face;
	for (k = 0; k<nt; k++)
	{
		for (j = 0; j<6; j++)
		{
			edge = elementEdge->val[k][j];
			perm = elements->eperm[k][j];

			for(ii[0]=0;ii[0]<2;ii[0]++)
			{
				ii[1]=1-ii[0];
				elementDOF->val[k][2*j+ ii[0]] = edge*2 + ii[perm[0]];
			}
		}

		for (j = 0; j<4; j++)
		{
			face = elementFace->val[k][j];
			for (i = 0; i<5; i++)
				elementDOF->val[k][12 + 5 * j + i] = ne*2 + face*5 + i;
		}
	} // k
}

/**
 * \fn void assembleMassmatrixTracelessDGPk3d(dBDmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF)
 * \brief assemble mass matrix (u, v)
 * \param *A pointer to mass matrix
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
 * \return void
 */
void assembleMassmatrixTracelessDGPk3d(dBDmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF)
{	
	int i,j,k,l;

	// int nvertices=nodes->row;
	// int nedges=edges->row;
	// int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[3], phi2[3], val1;
	int k1,k2,i1,j1,l1,l2,i2;
	double x[3];
	double vol;
	int count;
	// double *lambdaConst;

	int num_qp;
	double lambdas[100][4], weight[100];	

	/************************************************** stiffness matrix A *****************************************************************/
	int nb = elements->row * 8;
	create_dbd_matrix(elementDOF[0].col * nb, elementDOF[0].col * nb, nb, A);
	for(i=0;i<A->nb;i++)
		create_dden_matrix(elementDOF[0].col, elementDOF[0].col, A->blk+i);

	ddenmat *lA, *lB; // local A

	num_qp = getNumQuadPoints_ShunnWilliams(2*elementDOF->dop, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vol = elements->vol[k];
		// gradLambda = elements->gradLambda[k];
		// end set parameters

		// init_dden_matrix(&lA, 0.0);
		lA = A->blk + k*8;
		for (k1 = 0; k1<elementDOF->col; k1++){
			for (k2 = 0; k2<elementDOF->col; k2++){
				val1 = 0;
				for (i1 = 0; i1<num_qp; i1++){
					lagrange3d_basis(lambdas[i1], k1, elementDOF->dop, phi1);
					lagrange3d_basis(lambdas[i1], k2, elementDOF->dop, phi2);
					val1 += vol * weight[i1] * phi1[0] * phi2[0];
				}
				lA->val[k1][k2] += val1;
			} // k2
		} // k1

		for(i=1;i<8;i++){
			lB = lA+i;
			for (k1 = 0; k1<elementDOF->col; k1++){
				for (k2 = 0; k2<elementDOF->col; k2++)
				lB->val[k1][k2] = lA->val[k1][k2];
			}
		} // i
	} // k
}

/**
 * \fn void assembleMassmatrixLagrange3d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF)
 * \brief assemble mass matrix (u, v)
 * \param *A pointer to mass matrix
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
 * \return void
 */
void assembleMassmatrixLagrange3d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF)
{	
	int i,j,k,l;

	// int nvertices=nodes->row;
	// int nedges=edges->row;
	// int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[3], phi2[3], val1;
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3];
	double vol;
	int count;
	// double *lambdaConst;

	int num_qp;
	double lambdas[100][4], weight[100];	

	/************************************************** stiffness matrix A *****************************************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF[0].col*elementDOF[0].col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(elementDOF[0].col, elementDOF[0].col, &lA);

	num_qp = getNumQuadPoints_ShunnWilliams(2*elementDOF->dop, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vol = elements->vol[k];
		// gradLambda = elements->gradLambda[k];
		// end set parameters

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF->col; k1++){
			for (k2 = 0; k2<elementDOF->col; k2++){
				val1 = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					lagrange3d_basis(lambdas[i1], k1, elementDOF->dop, phi1);
					lagrange3d_basis(lambdas[i1], k2, elementDOF->dop, phi2);
					val1 += vol * weight[i1] * phi1[0] * phi2[0];
				}
				lA.val[k1][k2] += val1;
			} // k2
		} // k1

		l = elementDOF[0].col*elementDOF[0].col * k;
		for (i = 0; i<elementDOF[0].col; i++)
		{
			for (j = 0; j<elementDOF[0].col; j++)
			{
				ia[l] = elementDOF[0].val[k][i];
				ja[l] = elementDOF[0].val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // i
		} // j
	} // k
	free_dden_matrix(&lA);	

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	dIJtoCSReps(A, ia, ja, va, N, 0, 0, 0);
	
}

/**
* \fn void assembleMassmatrixP03d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF)
* \brief assemble mass matrix
* \param *M pointer to mass matrix
* \param *elements pointer to the structure of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void assembleMassmatrixP03d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF)
{
	int i, j, k;

	M->row = elementDOF->dof;
	M->col = elementDOF->dof;
	M->IA = (int*)calloc(M->row + 1, sizeof(int));
	M->JA = NULL;
	M->val = NULL;
	for (i = 0; i<M->row; i++)
		M->IA[i + 1] = 1;
	for (i = 0; i<M->row; i++)
		M->IA[i + 1] += M->IA[i];
	M->nnz = M->IA[M->row];
	M->JA = (int*)calloc(M->nnz, sizeof(int));
	for (i = 0; i<M->nnz; i++)
		M->JA[i] = i;
	M->val = (double*)calloc(M->nnz, sizeof(double));
	for (k = 0; k < elements->row; k++)
		M->val[k] = elements->vol[k];
}

/**
 * \fn void assembleBiGradLagrange3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
 * \brief assemble stiffness matrix (grad u, grad v)
 * \param *A pointer to stiffness matrix
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
 * \return void
 */
void assembleBiGradLagrange3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{	
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[3], phi2[3], val1;
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], **gradLambda, **nv, *nvf[4], *etv[6], **vertices;
	double vol, s[4], h[6];
	int rowstart[3], row31[3];
	int count;
	short *forien, *eorien;
	int **fpermi, **eperm;
	double *lambdaConst;

	int num_qp;
	double lambdas[100][4], lambdas2[100][3], weight[100];
	
	/************************************************** stiffness matrix A *****************************************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF[0].col*elementDOF[0].col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(elementDOF[0].col, elementDOF[0].col, &lA);
	
	num_qp = getNumQuadPoints_ShunnWilliams(2*(elementDOF[0].dop - 1), 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vol = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		// end set parameters

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				val1 = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					lagrange3d_basis1(lambdas[i1], gradLambda, k1, elementDOF[0].dop, phi1);
					lagrange3d_basis1(lambdas[i1], gradLambda, k2, elementDOF[0].dop, phi2);
					val1 += vol * weight[i1] * dot_array(3, phi1, phi2);
				}
				lA.val[k1][k2] += val1;
			} // k2
		} // k1
		
		l = elementDOF[0].col*elementDOF[0].col * k;
		for (i = 0; i<elementDOF[0].col; i++)
		{
			for (j = 0; j<elementDOF[0].col; j++)
			{
				ia[l] = elementDOF[0].val[k][i];
				ja[l] = elementDOF[0].val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // j
		} // i
	} // k	
	free_dden_matrix(&lA);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	dIJtoCSReps(A, ia, ja, va, N, 0, 0, 0);
}

/**
 * \fn void assembleRHSLagrange3d(dvector *b, ELEMENT *elements, ELEMENT_DOF *elementDOF, double (*f)(double *))
 * \brief assemble stiffness matrix (f, v)
 * \param *b pointer to right hand side
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
 * \return void
 */
void assembleRHSLagrange3d(dvector *b, ELEMENT *elements, ELEMENT_DOF *elementDOF, double (*f)(double *))
{
	int i,j,k,k1,i1;
	
	double phi;
	double x[3], **gradLambda, **vertices;
	double vol;

	int num_qp;
	double lambdas[100][4], weight[100];
			
	dvector lb;
	create_dvector(elementDOF->col, &lb);
	/************************************************** right hand side b *****************************************************************/
	create_dvector(elementDOF->dof, b);
	num_qp = getNumQuadPoints_ShunnWilliams(9, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vol = elements->vol[k];
		vertices = elements->vertices[k];
		// end set parameters

		init_dvector(&lb, 0.0);
		for (i = 0; i<elementDOF->col; i++)
		{
			for (i1 = 0; i1<num_qp; i1++)
			{
				lagrange3d_basis(lambdas[i1], i, elementDOF->dop, &phi);
				axy_array(3, lambdas[i1][3], vertices[3], x);
				for(j=0;j<3;j++)
					axpy_array(3, lambdas[i1][j], vertices[j], x);
					lb.val[i] += vol*weight[i1] * f(x)*phi;
				// b->val[i] += vol*weight[i1] * poisson3d_f(x)*phi;
			} // i1
		} // k1

		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			i = elementDOF->val[k][k1];
			b->val[i] += lb.val[k1];
		} // k1
	} // k
	free_dvector(&lb);
}

/**
 * \fn void assembleRHSdgPoly3d(dvector *b, ELEMENT *elements, ELEMENT_DOF *elementDOF, double (*f)(double *))
 * \brief assemble stiffness matrix (f, v)
 * \param *b pointer to right hand side
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
 * \return void
 */
void assembleRHSdgPoly3d(dvector *b, ELEMENT *elements, ELEMENT_DOF *elementDOF, double (*f)(double *))
{
	int i,j,k,k1,i1;
	
	double phi;
	double x[3], **gradLambda, **vertices;
	double vol;

	int num_qp;
	double lambdas[100][4], weight[100];
			
	dvector lb;
	create_dvector(elementDOF->col, &lb);
	/************************************************** right hand side b *****************************************************************/
	create_dvector(elementDOF->dof, b);
	num_qp = getNumQuadPoints_ShunnWilliams(9, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vol = elements->vol[k];
		vertices = elements->vertices[k];
		// end set parameters

		init_dvector(&lb, 0.0);
		for (i = 0; i<elementDOF->col; i++)
		{
			for (i1 = 0; i1<num_qp; i1++)
			{
				lagrange3d_basis(lambdas[i1], i, elementDOF->dop, &phi);
				axy_array(3, lambdas[i1][3], vertices[3], x);
				for(j=0;j<3;j++)
					axpy_array(3, lambdas[i1][j], vertices[j], x);
					lb.val[i] += vol*weight[i1] * f(x)*phi;
				// b->val[i] += vol*weight[i1] * poisson3d_f(x)*phi;
			} // i1
		} // k1

		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			i = elementDOF->val[k][k1];
			b->val[i] += lb.val[k1];
		} // k1
	} // k
	free_dvector(&lb);
}

/**
 * \fn void assembleBiGradVectorNcP13d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *BT pointer to stiffness matrix
 * \param *C pointer to stiffness matrix
 * \param *b pointer to right hand side
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
 * \return void
 */
void assembleBiGradVectorNcP13d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Cx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[3], phi2[3], val1;
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], **gradLambda, **nv, *nvf[4], *etv[6], **vertices;
	double vol, s[4], h[6];
	int rowstart[3], row31[3];
	int count;
	short *forien, *eorien;
	int **fpermi;
	double *lambdaConst;

	int num_qp;
	double lambdas[100][4], lambdas2[100][3], weight[100];
		
	/************************************************** stiffness matrix A *****************************************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF[0].col*elementDOF[0].col*elements->row;
	ia = (int*)malloc(3*N * sizeof(int));
	ja = (int*)malloc(3*N * sizeof(int));
	va = (double*)malloc(3*N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(elementDOF[0].col, elementDOF[0].col, &lA);
	
	num_qp = getNumQuadPoints_ShunnWilliams(2*(elementDOF[0].dop - 1), 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vol = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		// end set parameters

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				val1 = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					ncp13d_basis1(gradLambda, k1, phi1);
					ncp13d_basis1(gradLambda, k2, phi2);
					val1 += vol * weight[i1] * dot_array(3, phi1, phi2);
				}
				lA.val[k1][k2] += val1;
			} // k2
		} // k1
		
		l = elementDOF[0].col*elementDOF[0].col * k;
		for (i = 0; i<elementDOF[0].col; i++)
		{
			for (j = 0; j<elementDOF[0].col; j++)
			{
				ia[l] = elementDOF[0].val[k][i];
				ja[l] = elementDOF[0].val[k][j];
				va[l] = lA.val[i][j];
				ia[l+N] = ia[l] + elementDOF[0].dof;
				ja[l+N] = ja[l] + elementDOF[0].dof;
				va[l+N] = va[l];
				ia[l+N*2] = ia[l] + elementDOF[0].dof*2;
				ja[l+N*2] = ja[l] + elementDOF[0].dof*2;
				va[l+N*2] = va[l];
				l++;
			} // j
		} // i
	} // k	
	free_dden_matrix(&lA);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	dIJtoCSReps(A, ia, ja, va, N*3, 0, 0, 0);
}

/**
 * \fn void assembleDivNcP1P03d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *BT pointer to stiffness matrix
 * \param *C pointer to stiffness matrix
 * \param *b pointer to right hand side
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
 * \return void
 */
void assembleDivNcP1P03d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[3], phi2[3], val[3];
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], **gradLambda, **nv, *nvf[4], *etv[6], **vertices;
	double vol, s[4], h[6];
	int rowstart[3], row31[3];
	int count;
	short *forien, *eorien;
	int **fpermi;
	double *lambdaConst;

	int num_qp;
	double lambdas[100][4], lambdas2[100][3], weight[100];
	
	/************************************************** stiffness matrix A *****************************************************************/
	int *ia, *ja;
	double *va;
	int N = 3*elementDOF[0].col*elementDOF[1].col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA[3]; // local A
	create_dden_matrix(elementDOF[1].col, elementDOF[0].col, lA);
	create_dden_matrix(elementDOF[1].col, elementDOF[0].col, lA+1);
	create_dden_matrix(elementDOF[1].col, elementDOF[0].col, lA+2);

	num_qp = getNumQuadPoints_ShunnWilliams(elementDOF[0].dop+elementDOF[1].dop-1, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vol = elements->vol[k];
		gradLambda = elements->gradLambda[k];
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

		init_dden_matrix(&lA[0], 0.0); init_dden_matrix(&lA[1], 0.0); init_dden_matrix(&lA[2], 0.0);
		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				val[0] = 0; val[1] = 0; val[2] = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					// lagrange3d_basis1(lambdas[i1], gradLambda, k1, elementDOF[1].dop, phi1);
					ncp13d_basis1(gradLambda, k2, phi2);
					val[0] += vol*weight[i1] * phi2[0];
					val[1] += vol*weight[i1] * phi2[1];
					val[2] += vol*weight[i1] * phi2[2];
				}
				lA[0].val[k1][k2] += val[0];
				lA[1].val[k1][k2] += val[1];
				lA[2].val[k1][k2] += val[2];
			} // k2
		} // k1

		l = 3*elementDOF[0].col*elementDOF[1].col * k;
		for (i = 0; i<elementDOF[1].col; i++)
		{
			for (j = 0; j<elementDOF[0].col; j++)
			{
				ia[l] = elementDOF[1].val[k][i];
				ja[l] = elementDOF[0].val[k][j];
				va[l] = lA[0].val[i][j];
				l++;
				ia[l] = elementDOF[1].val[k][i];
				ja[l] = elementDOF[0].val[k][j]+elementDOF[0].dof;
				va[l] = lA[1].val[i][j];
				l++;
				ia[l] = elementDOF[1].val[k][i];
				ja[l] = elementDOF[0].val[k][j]+elementDOF[0].dof*2;
				va[l] = lA[2].val[i][j];
				l++;
			} // j
		} // i
	} // k
	free_dden_matrix(&lA[0]);
	free_dden_matrix(&lA[1]);
	free_dden_matrix(&lA[2]);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	dIJtoCSReps(A, ia, ja, va, N, 0, 0, 0);
}

/**
 * \fn void assembleRHSVectorNcP13d(dvector *b, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, void (*f)(double *, double *))
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *BT pointer to stiffness matrix
 * \param *C pointer to stiffness matrix
 * \param *b pointer to right hand side
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
 * \return void
 */
void assembleRHSVectorNcP13d(dvector *b, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, void (*f)(double *, double *))
{
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[3], phi2[3];
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], **grd_lambda, **nv, *nvf[4], *etv[6], **vertices;
	double vol, s[4], h[6];
	int rowstart[3], row31[3];
	int count;
	short *forien, *eorien;
	int **fpermi;
	double *lambdaConst;

	int num_qp;
	double lambdas[100][4], weight[100];
	
	dvector lb;
	create_dvector(3*elementDOF[0].col, &lb);
	/************************************************** right hand side b *****************************************************************/
	create_dvector(3*elementDOF[0].dof, b);
	num_qp = getNumQuadPoints_ShunnWilliams(9, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vol = elements->vol[k];
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

		init_dvector(&lb, 0.0);
		for (i = 0; i<elementDOF[0].col; i++)
		{
			for (i1 = 0; i1<num_qp; i1++)
			{
				ncp13d_basis(lambdas[i1], i, &phi);
				axy_array(3, lambdas[i1][3], vertices[3], x);
				for(j=0;j<3;j++)
					axpy_array(3, lambdas[i1][j], vertices[j], x);
				// maxwell3d_f(x, phi2);
				f(x, phi2);
				lb.val[i] += vol*weight[i1] * phi2[0] * phi;
				lb.val[i+elementDOF[0].col] += vol*weight[i1] * phi2[1] * phi;
				lb.val[i+2*elementDOF[0].col] += vol*weight[i1] * phi2[2] * phi;
			} // i1
		} // i

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			b->val[i] += lb.val[k1];
			b->val[i + elementDOF[0].dof] += lb.val[k1+elementDOF[0].col];
			b->val[i + elementDOF[0].dof*2] += lb.val[k1+elementDOF[0].col*2];
		} // k1
	} // k
	free_dvector(&lb);
}

/**
 * \fn void assembleBiCurlNedelec1st3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *BT pointer to stiffness matrix
 * \param *C pointer to stiffness matrix
 * \param *b pointer to right hand side
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
 * \return void
 */
void assembleBiCurlNedelec1st3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Cx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[3], phi2[3], val1;
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], **grd_lambda, **nv, *nvf[4], *etv[6], **vertices;
	double vol, s[4], h[6];
	int rowstart[3], row31[3];
	int count;
	short *forien, *eorien;
	int **fpermi;
	double *lambdaConst;

	int num_qp;
	double lambdas[100][4], lambdas2[100][3], weight[100];
		
	/************************************************** stiffness matrix A *****************************************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF[0].col*elementDOF[0].col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(elementDOF[0].col, elementDOF[0].col, &lA);

	num_qp = getNumQuadPoints_ShunnWilliams(2*(elementDOF[0].dop-1), 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vol = elements->vol[k];
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
		}
		// end set parameters

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				val1 = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					nedelec1st3d_basisCurl(lambdas[i1], grd_lambda, eorien, fpermi, k1, elementDOF[0].dop, phi1);
					nedelec1st3d_basisCurl(lambdas[i1], grd_lambda, eorien, fpermi, k2, elementDOF[0].dop, phi2);
					val1 += vol*weight[i1] * dot_array(3, phi1, phi2);
				}
				lA.val[k1][k2] += val1;
			} // k2
		} // k1

		l = elementDOF[0].col*elementDOF[0].col * k;
		for (i = 0; i<elementDOF[0].col; i++)
		{
			for (j = 0; j<elementDOF[0].col; j++)
			{
				ia[l] = elementDOF[0].val[k][i];
				ja[l] = elementDOF[0].val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // j
		} // i
	} // k
	free_dden_matrix(&lA);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	dIJtoCSReps(A, ia, ja, va, N, 0, 0, 0);
}

/**
 * \fn void assembleNedelec1stGradLagrange3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1)
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *BT pointer to stiffness matrix
 * \param *C pointer to stiffness matrix
 * \param *b pointer to right hand side
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
 * \return void
 */
void assembleNedelec1stGradLagrange3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Cx2   = 0
	where x1: u_h, x2: lambda_h
	**/

	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[3], phi2[3], val1;
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], **grd_lambda, **nv, *nvf[4], *etv[6], **vertices;
	double vol, s[4], h[6];
	int rowstart[3], row31[3];
	int count;
	short *forien, *eorien;
	int **fpermi;
	double *lambdaConst;

	int num_qp;
	double lambdas[100][4], weight[100];
		
	/************************************************** stiffness matrix A *****************************************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF0->col*elementDOF1->col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(elementDOF1->col, elementDOF0->col, &lA);

	num_qp = getNumQuadPoints_ShunnWilliams(elementDOF0->dop+elementDOF1->dop-1, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vol = elements->vol[k];
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

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF1->col; k1++)
		{
			for (k2 = 0; k2<elementDOF0->col; k2++)
			{
				val1 = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					lagrange3d_basis1(lambdas[i1], grd_lambda, k1, elementDOF1->dop, phi1);
					nedelec1st3d_basis(lambdas[i1], grd_lambda, eorien, fpermi, k2, elementDOF0->dop, phi2);
					val1 += vol*weight[i1] * dot_array(3, phi1, phi2);
				}
				lA.val[k1][k2] += val1;
			} // k2
		} // k1

		l = elementDOF0->col*elementDOF1->col * k;
		for (i = 0; i<elementDOF1->col; i++)
		{
			for (j = 0; j<elementDOF0->col; j++)
			{
				ia[l] = elementDOF1->val[k][i];
				ja[l] = elementDOF0->val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // j
		} // i
	} // k
	free_dden_matrix(&lA);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	dIJtoCSReps(A, ia, ja, va, N, 0, 0, 0);
}

/**
 * \fn void assembleRHSNedelec1st3d(dvector *b, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, void (*f)(double *, double *))
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *BT pointer to stiffness matrix
 * \param *C pointer to stiffness matrix
 * \param *b pointer to right hand side
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
 * \return void
 */
void assembleRHSNedelec1st3d(dvector *b, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, void (*f)(double *, double *))
{
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[3], phi2[3];
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], **grd_lambda, **nv, *nvf[4], *etv[6], **vertices;
	double vol, s[4], h[6];
	int rowstart[3], row31[3];
	int count;
	short *forien, *eorien;
	int **fpermi;
	double *lambdaConst;

	int num_qp;
	double lambdas[100][4], weight[100];
	
	int *index;
	int istart;
	
	dvector lb;
	create_dvector(elementDOF[0].col, &lb);
	/************************************************** right hand side b *****************************************************************/
	create_dvector(elementDOF[0].dof, b);
	num_qp = getNumQuadPoints_ShunnWilliams(9, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vol = elements->vol[k];
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

		init_dvector(&lb, 0.0);
		for (i = 0; i<elementDOF[0].col; i++)
		{
			for (i1 = 0; i1<num_qp; i1++)
			{
				nedelec1st3d_basis(lambdas[i1], grd_lambda, eorien, fpermi, i, elementDOF[0].dop, phi1);
				axy_array(3, lambdas[i1][3], vertices[3], x);
				for(j=0;j<3;j++)
					axpy_array(3, lambdas[i1][j], vertices[j], x);
				// maxwell3d_f(x, phi2);
				f(x, phi2);
				lb.val[i] += vol*weight[i1] * dot_array(3, phi1, phi2);
			} // i1
		} // i

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			b->val[i] += lb.val[k1];
		} // k1
	} // k
	free_dvector(&lb);
}

/**
 * \fn void assembleBiCurlNedelec2nd3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *BT pointer to stiffness matrix
 * \param *C pointer to stiffness matrix
 * \param *b pointer to right hand side
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
 * \return void
 */
void assembleBiCurlNedelec2nd3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[3], phi2[3], val1;
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], **grd_lambda, **nv, *nvf[4], *etv[6], **vertices;
	double vol, s[4], h[6];
	int rowstart[3], row31[3];
	int count;
	short *forien, *eorien;
	int **fpermi, **eperm;
	double *lambdaConst;

	int num_qp;
	double lambdas[100][4], weight[100];
	
	/************************************************** stiffness matrix A *****************************************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF[0].col*elementDOF[0].col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(elementDOF[0].col, elementDOF[0].col, &lA);

	num_qp = getNumQuadPoints_ShunnWilliams(2*(elementDOF[0].dop-1), 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
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
		for (i = 0; i < elementFace->col; i++)
		{
			face = elementFace->val[k][i];
			nvf[i] = faces->nvector[face];
		}
		// end set parameters

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				val1 = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					nedelec2nd3d_basisCurl(lambdas[i1], grd_lambda, eperm, k1, elementDOF[0].dop, phi1);
					nedelec2nd3d_basisCurl(lambdas[i1], grd_lambda, eperm, k2, elementDOF[0].dop, phi2);
					val1 += vol*weight[i1] * dot_array(3, phi1, phi2);
				}
				lA.val[k1][k2] += val1;
			} // k2
		} // k1

		l = elementDOF[0].col*elementDOF[0].col * k;
		for (i = 0; i<elementDOF[0].col; i++)
		{
			for (j = 0; j<elementDOF[0].col; j++)
			{
				ia[l] = elementDOF[0].val[k][i];
				ja[l] = elementDOF[0].val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // j
		} // i
	} // k	
	free_dden_matrix(&lA);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	dIJtoCSReps(A, ia, ja, va, N, 0, 0, 0);
}

/**
 * \fn void assembleNedelec2ndGradLagrange3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1)
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *BT pointer to stiffness matrix
 * \param *C pointer to stiffness matrix
 * \param *b pointer to right hand side
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
 * \return void
 */
void assembleNedelec2ndGradLagrange3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1)
{	
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[3], phi2[3], val1;
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], **grd_lambda, **nv, *nvf[4], *etv[6], **vertices;
	double vol, s[4], h[6];
	int rowstart[3], row31[3];
	int count;
	short *forien, *eorien;
	int **fpermi, **eperm;
	double *lambdaConst;

	int num_qp;
	double lambdas[100][4], weight[100];
		
	/************************************************** stiffness matrix A *****************************************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF0->col*elementDOF1->col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));

	ddenmat lA; // local A
	create_dden_matrix(elementDOF1->col, elementDOF0->col, &lA);

	num_qp = getNumQuadPoints_ShunnWilliams(elementDOF0->dop+elementDOF1->dop-1, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
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

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF1->col; k1++)
		{
			for (k2 = 0; k2<elementDOF0->col; k2++)
			{
				val1 = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					lagrange3d_basis1(lambdas[i1], grd_lambda, k1, elementDOF1->dop, phi1);
					nedelec2nd3d_basis(lambdas[i1], grd_lambda, eperm, k2, elementDOF0->dop, phi2);
					val1 += vol*weight[i1] * dot_array(3, phi1, phi2);
				}
				lA.val[k1][k2] += val1;
			} // k2
		} // k1

		l = elementDOF0->col*elementDOF1->col * k;
		for (i = 0; i<elementDOF1->col; i++)
		{
			for (j = 0; j<elementDOF0->col; j++)
			{
				ia[l] = elementDOF1->val[k][i];
				ja[l] = elementDOF0->val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // j
		} // i
	} // k
	free_dden_matrix(&lA);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	dIJtoCSReps(A, ia, ja, va, N, 0, 0, 0);
}

/**
 * \fn void assembleRHSNedelec2nd3d(dvector *b, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, void (*f)(double *, double *))
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *BT pointer to stiffness matrix
 * \param *C pointer to stiffness matrix
 * \param *b pointer to right hand side
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
 * \return void
 */
void assembleRHSNedelec2nd3d(dvector *b, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, void (*f)(double *, double *))
{
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[3], phi2[3];
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], **grd_lambda, **nv, *nvf[4], *etv[6], **vertices;
	double vol, s[4], h[6];
	int rowstart[3], row31[3];
	int count;
	short *forien, *eorien;
	int **fpermi, **eperm;
	double *lambdaConst;

	int num_qp;
	double lambdas[100][4], weight[100];
	
	int *index;
	int istart;
	
	dvector lb;
	create_dvector(elementDOF[0].col, &lb);
	 /************************************************** right hand side b *****************************************************************/
	create_dvector(elementDOF[0].dof, b);
	num_qp = getNumQuadPoints_ShunnWilliams(9, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
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

		init_dvector(&lb, 0.0);
		for (i = 0; i<elementDOF[0].col; i++)
		{
			for (i1 = 0; i1<num_qp; i1++)
			{
				nedelec2nd3d_basis(lambdas[i1], grd_lambda, eperm, i, elementDOF[0].dop, phi1);
				axy_array(3, lambdas[i1][3], vertices[3], x);
				for(j=0;j<3;j++)
					axpy_array(3, lambdas[i1][j], vertices[j], x);
				// maxwell3d_f(x, phi2);
				f(x, phi2);
				lb.val[i] += vol*weight[i1] * dot_array(3, phi1, phi2);
			} // i1
		} // k1 

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			b->val[i] += lb.val[k1];
		} // k1 
	} // k
	free_dvector(&lb);
}

/**
 * \fn void assembleRHSNedelec3d(dvector *b, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, int nrs, int curlfemtype, void (*f)(double *, double *))
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *BT pointer to stiffness matrix
 * \param *C pointer to stiffness matrix
 * \param *b pointer to right hand side
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
 * \param nrs extend vector b with additional nrs zeros, nrs = 0 means no extension
 * \param curlfemtype type of Nedelec element
 * \return void
 */
void assembleRHSNedelec3d(dvector *b, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, int nrs, int curlfemtype, void (*f)(double *, double *))
{
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[3], phi2[3];
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], **grd_lambda, **nv, *nvf[4], *etv[6], **vertices;
	double vol, s[4], h[6];
	int rowstart[3], row31[3];
	int count;
	short *forien, *eorien;
	int **fpermi, **eperm;
	double *lambdaConst;

	int num_qp;
	double lambdas[100][4], weight[100];
	
	int *index;
	int istart;
	
	dvector lb;
	create_dvector(elementDOF[0].col, &lb);
	 /************************************************** right hand side b *****************************************************************/
	create_dvector(elementDOF[0].dof+nrs, b);
	num_qp = getNumQuadPoints_ShunnWilliams(9, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
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

		init_dvector(&lb, 0.0);
		for (i = 0; i<elementDOF[0].col; i++)
		{
			for (i1 = 0; i1<num_qp; i1++)
			{
				if(curlfemtype == 1)
					nedelec1st3d_basis(lambdas[i1], grd_lambda, eorien, fpermi, i, elementDOF[0].dop, phi1);
				else
					nedelec2nd3d_basis(lambdas[i1], grd_lambda, eperm, i, elementDOF[0].dop, phi1);
				axy_array(3, lambdas[i1][3], vertices[3], x);
				for(j=0;j<3;j++)
					axpy_array(3, lambdas[i1][j], vertices[j], x);
				// maxwell3d_f(x, phi2);
				f(x, phi2);
				lb.val[i] += vol*weight[i1] * dot_array(3, phi1, phi2);
			} // i1
		} // k1 

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			b->val[i] += lb.val[k1];
		} // k1 
	} // k
	free_dvector(&lb);
}

/**
 * \fn void assembleBiCurlHuangZhang3d(dCSRmat *A, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *basisCoeffs pointer to coefficients of basis
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
 * \return void
 */
void assembleBiCurlHuangZhang3d(dCSRmat *A, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Cx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[3], phi2[3], val1;
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
	
	/************************************************** stiffness matrix A *****************************************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF[0].col*elementDOF[0].col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA, lB; // local A: lA = C lB C^T
	create_dden_matrix(elementDOF[0].col, elementDOF[0].col, &lA);
	create_dden_matrix(elementDOF[0].col, elementDOF[0].col, &lB);
	// step 3A: Loop element by element and compute the actual entries storing them in A
	num_qp = getNumQuadPoints_ShunnWilliams(2*(elementDOF[0].dop-1), 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vol = elements->vol[k];
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
		}
		// end set parameters

		init_dden_matrix(&lB, 0.0);
		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				val1 = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					huangzhang03d_basisCurl(lambdas[i1], grd_lambda, k1, phi1);
					huangzhang03d_basisCurl(lambdas[i1], grd_lambda, k2, phi2);
					val1 += vol*weight[i1] * dot_array(3, phi1, phi2);
				}
				lB.val[k1][k2] += val1;
			} // k2
		} // k1

		ddenmat lbC;
		lbC.row=basisCoeffs->row; lbC.col=basisCoeffs->col; lbC.val=basisCoeffs->val[k];
		ABAt_ddenmat(1.0, &lbC, &lB, &lA);

		l = elementDOF[0].col*elementDOF[0].col * k;
		for (i = 0; i<elementDOF[0].col; i++)
		{
			for (j = 0; j<elementDOF[0].col; j++)
			{
				ia[l] = elementDOF[0].val[k][i];
				ja[l] = elementDOF[0].val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // j
		} // i
	} // k
	free_dden_matrix(&lA);
	free_dden_matrix(&lB);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	dIJtoCSReps(A, ia, ja, va, N, 0, 0, 0);
}

/**
 * \fn void assembleBiGradcurlHuangZhang3d(dCSRmat *A, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *BT pointer to stiffness matrix
 * \param *C pointer to stiffness matrix
 * \param *b pointer to right hand side
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
 * \return void
 */
void assembleBiGradcurlHuangZhang3d(dCSRmat *A, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[9], phi2[9], val1;
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
	
	/************************************************** stiffness matrix A *****************************************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF[0].col*elementDOF[0].col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA, lB; // local A: lA = C lB C^T
	create_dden_matrix(elementDOF[0].col, elementDOF[0].col, &lA);
	create_dden_matrix(elementDOF[0].col, elementDOF[0].col, &lB);

	num_qp = getNumQuadPoints_ShunnWilliams(2*(elementDOF[0].dop-2), 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vol = elements->vol[k];
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
		}
		// end set parameters

		init_dden_matrix(&lB, 0.0);
		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{	
				val1 = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					huangzhang03d_basisGradCurl(lambdas[i1], grd_lambda, k1, phi1);
					huangzhang03d_basisGradCurl(lambdas[i1], grd_lambda, k2, phi2);
					val1 += vol*weight[i1] * dot_array(9, phi1, phi2);
				}
				lB.val[k1][k2] += val1;
			} // k2
		} // k1

		ddenmat lbC;
		lbC.row=basisCoeffs->row; lbC.col=basisCoeffs->col; lbC.val=basisCoeffs->val[k];
		ABAt_ddenmat(1.0, &lbC, &lB, &lA);

		l = elementDOF[0].col*elementDOF[0].col * k;
		for (i = 0; i<elementDOF[0].col; i++)
		{
			for (j = 0; j<elementDOF[0].col; j++)
			{
				ia[l] = elementDOF[0].val[k][i];
				ja[l] = elementDOF[0].val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // j
		} // i
	} // k	
	free_dden_matrix(&lA);
	free_dden_matrix(&lB);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	dIJtoCSReps(A, ia, ja, va, N, 0, 0, 0);
}

/**
 * \fn void assembleBiGradcurlperturbHuangZhang3d(dCSRmat *A, double paraeps, short nitsche, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *BT pointer to stiffness matrix
 * \param *C pointer to stiffness matrix
 * \param *b pointer to right hand side
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
 * \return void
 */
void assembleBiGradcurlperturbHuangZhang3d(dCSRmat *A, double paraeps, short nitsche, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi0[9], phi1[9], phi2[9], val1;
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], *xK, **grd_lambda, **nv, *nvf[4], *etv[6], **vertices;
	double vol, s, h[6], C11;
	int fi[3], rowstart[3], row31[3];
	int count;
	short *forien, *eorien;
	int **fpermi;
	double *lambdaConst;

	int num_qp, num_qp1, num_qp2;
	double lambdasT[4]; 
	double lambdas[100][4], weight[100];
	double lambdas1[100][4], weight1[100];
	double lambdas2[100][3], weight2[100];
	
	/************************************************** stiffness matrix A *****************************************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF[0].col*elementDOF[0].col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA, lB; // local A: lA = C lB C^T
	create_dden_matrix(elementDOF[0].col, elementDOF[0].col, &lA);
	create_dden_matrix(elementDOF[0].col, elementDOF[0].col, &lB);
	
	num_qp = getNumQuadPoints_ShunnWilliams(2*(elementDOF[0].dop-2), 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	num_qp1 = getNumQuadPoints_ShunnWilliams(2*(elementDOF[0].dop-1), 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp1, lambdas1, weight1); // Shunn-Williams intergation initial
	num_qp2=getNumQuadPoints_ShunnWilliams(2*(elementDOF[0].dop-1), 2); // the number of numerical intergation points
	init_ShunnWilliams2d(num_qp2, lambdas2, weight2); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vol = elements->vol[k];
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
		}
		// end set parameters

		init_dden_matrix(&lB, 0.0);
		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{	
				val1 =0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					huangzhang03d_basisGradCurl(lambdas[i1], grd_lambda, k1, phi1);
					huangzhang03d_basisGradCurl(lambdas[i1], grd_lambda, k2, phi2);
					val1 += vol*weight[i1] * dot_array(9, phi1, phi2);
				}
				if(nitsche > 0){
					for (i = 0; i < elementFace->col; i++){
						face = elementFace->val[k][i];
						s = faces->area[face];
						C11 = 5.0 / sqrt(s);
						if (faces->bdFlag[face] == 1 || faces->bdFlag[face] == 2 || faces->bdFlag[face] == 3 || faces->bdFlag[face] == 4){ // Dirichlet boundary
							face2vertices3d(i, fi);
							for (i1 = 0; i1<num_qp2; i1++){
								lambdasT[i]=0;
								for(j=0;j<3;j++) lambdasT[fi[j]] = lambdas2[i1][j];

								huangzhang03d_basisGradCurl(lambdasT, grd_lambda, k1, phi1);
								huangzhang03d_basisCurl(lambdasT, grd_lambda, k2, phi2);
								for(j=0;j<3;j++) phi0[j] = dot_array(3, nv[i], phi1+j);
								val1 -= s*weight2[i1] * dot_array(3, phi0, phi2);

								huangzhang03d_basisCurl(lambdasT, grd_lambda, k1, phi1);
								huangzhang03d_basisGradCurl(lambdasT, grd_lambda, k2, phi2);
								for(j=0;j<3;j++) phi0[j] = dot_array(3, nv[i], phi2+j);
								val1 -= s*weight2[i1] * dot_array(3, phi1, phi0);

								huangzhang03d_basisCurl(lambdasT, grd_lambda, k1, phi1);
								huangzhang03d_basisCurl(lambdasT, grd_lambda, k2, phi2);
								val1 += s*weight2[i1] * C11 * dot_array(3, phi1, phi2);
							}
						}
					}
				}

				lB.val[k1][k2] += val1*paraeps*paraeps;
				val1 = 0;
				for (i1 = 0; i1<num_qp1; i1++)
				{
					huangzhang03d_basisCurl(lambdas1[i1], grd_lambda, k1, phi1);
					huangzhang03d_basisCurl(lambdas1[i1], grd_lambda, k2, phi2);
					val1 += vol*weight1[i1] * dot_array(3, phi1, phi2);
				}
				lB.val[k1][k2] += val1;
			} // k2
		} // k1

		ddenmat lbC;
		lbC.row=basisCoeffs->row; lbC.col=basisCoeffs->col; lbC.val=basisCoeffs->val[k];
		ABAt_ddenmat(1.0, &lbC, &lB, &lA);

		l = elementDOF[0].col*elementDOF[0].col * k;
		for (i = 0; i<elementDOF[0].col; i++)
		{
			for (j = 0; j<elementDOF[0].col; j++)
			{
				ia[l] = elementDOF[0].val[k][i];
				ja[l] = elementDOF[0].val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // i
		} // j
	} // k	
	free_dden_matrix(&lA);
	free_dden_matrix(&lB);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	dIJtoCSReps(A, ia, ja, va, N, 0, 0, 0);
}

/**
 * \fn void assembleHuangZhangGradLagrange3d(dCSRmat *A, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *basisCoeffs pointer to coefficients of basis
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
 * \return void
 */
void assembleHuangZhangGradLagrange3d(dCSRmat *A, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{	
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[9], phi2[9], val1;
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], *xK, **grd_lambda, **nv, *nvf[4], *etv[6], **vertices;
	double vol, s[4], h[6];
	int rowstart[3], row31[3];
	int count;
	short *forien, *eorien;
	int **fpermi;
	double *lambdaConst;

	int num_qp;
	double lambdas[100][4], weight[100];
		
	/************************************************** stiffness matrix A *****************************************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF[0].col*elementDOF[1].col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA, lB; // local A: lA = lB C^T
	create_dden_matrix(elementDOF[1].col, elementDOF[0].col, &lA);
	create_dden_matrix(elementDOF[1].col, elementDOF[0].col, &lB);

	num_qp = getNumQuadPoints_ShunnWilliams(elementDOF[0].dop+elementDOF[1].dop-1, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
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

		init_dden_matrix(&lB, 0.0);
		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				val1 = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					lagrange3d_basis1(lambdas[i1], grd_lambda, k1, elementDOF[1].dop, phi1);
					huangzhang03d_basis(lambdas[i1], grd_lambda, k2, phi2);
					val1 += vol*weight[i1] * dot_array(3, phi1, phi2);
				}
				lB.val[k1][k2] += val1;
			} // k2
		} // k1

		ddenmat lbC;
		lbC.row=basisCoeffs->row; lbC.col=basisCoeffs->col; lbC.val=basisCoeffs->val[k];
		ABt_ddenmat(1.0, &lB, &lbC, &lA);

		l = elementDOF[0].col*elementDOF[1].col * k;
		for (i = 0; i<elementDOF[1].col; i++)
		{
			for (j = 0; j<elementDOF[0].col; j++)
			{
				ia[l] = elementDOF[1].val[k][i];
				ja[l] = elementDOF[0].val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // i
		} // j
	} // k
	free_dden_matrix(&lA);
	free_dden_matrix(&lB);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	dIJtoCSReps(A, ia, ja, va, N, 0, 0, 0);
}

/**
 * \fn void assembleRHSHuangZhang3d(dvector *b, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, void (*f)(double *, double *))
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *BT pointer to stiffness matrix
 * \param *C pointer to stiffness matrix
 * \param *b pointer to right hand side
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
 * \return void
 */
void assembleRHSHuangZhang3d(dvector *b, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, void (*f)(double *, double *))
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Cx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[9], phi2[9];
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
	
	int *index;
	int istart;

	dvector la, lb;
	create_dvector(elementDOF[0].col, &la);
	create_dvector(elementDOF[0].col, &lb);
	/************************************************** right hand side b *****************************************************************/
	create_dvector(elementDOF[0].dof, b);
	num_qp = getNumQuadPoints_ShunnWilliams(9, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
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

		init_dvector(&la, 0.0);
		for (i = 0; i<elementDOF[0].col; i++){	
			for (i1 = 0; i1<num_qp; i1++){
				axy_array(3, lambdas[i1][3], vertices[3], x);
				for(j=0;j<3;j++)
					axpy_array(3, lambdas[i1][j], vertices[j], x);
				huangzhang03d_basis(lambdas[i1], grd_lambda, i, phi1);
				f(x, phi2);
				la.val[i] += vol*weight[i1] * dot_array(3, phi1, phi2);
			} // i1
		} // i

		ddenmat lbC;
		lbC.row=basisCoeffs->row; lbC.col=basisCoeffs->col; lbC.val=basisCoeffs->val[k];
		Axy_ddenmat(1.0, &lbC, la.val, lb.val);

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			b->val[i] += lb.val[k1];
		} // k1 
	} // k
	free_dvector(&la);
	free_dvector(&lb);
}

/**
 * \fn void assembleBiCurlHuang3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
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
 * \return void
 */
void assembleBiCurlHuang3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Cx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[3], phi2[3], val1;
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
	
	/************************************************** stiffness matrix A *****************************************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF[0].col*elementDOF[0].col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(elementDOF[0].col, elementDOF[0].col, &lA);

	num_qp = getNumQuadPoints_ShunnWilliams(2*(elementDOF[0].dop-1), 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vol = elements->vol[k];
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
		}
		// end set parameters

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{	
				val1 = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					huangQuadcurl3d_basisCurl(lambdas[i1], grd_lambda, nvf, eorien, fpermi, k1, phi1);
					huangQuadcurl3d_basisCurl(lambdas[i1], grd_lambda, nvf, eorien, fpermi, k2, phi2);
					val1 += vol*weight[i1] * dot_array(3, phi1, phi2);
				}
				lA.val[k1][k2] += val1;
			} // k2
		} // k1

		l = elementDOF[0].col*elementDOF[0].col * k;
		for (i = 0; i<elementDOF[0].col; i++)
		{
			for (j = 0; j<elementDOF[0].col; j++)
			{
				ia[l] = elementDOF[0].val[k][i];
				ja[l] = elementDOF[0].val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // j
		} // i
	} // k
	free_dden_matrix(&lA);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	dIJtoCSReps(A, ia, ja, va, N, 0, 0, 0);
}

/**
 * \fn void assembleBiGradcurlHuang3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *BT pointer to stiffness matrix
 * \param *C pointer to stiffness matrix
 * \param *b pointer to right hand side
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
 * \return void
 */
void assembleBiGradcurlHuang3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[9], phi2[9], val1;
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
		
	
	/************************************************** stiffness matrix A *****************************************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF[0].col*elementDOF[0].col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(elementDOF[0].col, elementDOF[0].col, &lA);

	num_qp = getNumQuadPoints_ShunnWilliams(2*(elementDOF[0].dop-2), 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vol = elements->vol[k];
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
		}
		// end set parameters

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				val1 = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					huangQuadcurl3d_basisGradCurl(grd_lambda, nvf, eorien, fpermi, k1, phi1);
					huangQuadcurl3d_basisGradCurl(grd_lambda, nvf, eorien, fpermi, k2, phi2);
					val1 += vol*weight[i1] * dot_array(9, phi1, phi2);
				}
				lA.val[k1][k2] += val1;
			} // k2
		} // k1

		l = elementDOF[0].col*elementDOF[0].col * k;
		for (i = 0; i<elementDOF[0].col; i++)
		{
			for (j = 0; j<elementDOF[0].col; j++)
			{
				ia[l] = elementDOF[0].val[k][i];
				ja[l] = elementDOF[0].val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // j
		} // i
	} // k	
	free_dden_matrix(&lA);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	dIJtoCSReps(A, ia, ja, va, N, 0, 0, 0);
}

/**
 * \fn void assembleHuangGradLagrange3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *BT pointer to stiffness matrix
 * \param *C pointer to stiffness matrix
 * \param *b pointer to right hand side
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
 * \return void
 */
void assembleHuangGradLagrange3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{	
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[9], phi2[9], val1;
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], *xK, **grd_lambda, **nv, *nvf[4], *etv[6], **vertices;
	double vol, s[4], h[6];
	int rowstart[3], row31[3];
	int count;
	short *forien, *eorien;
	int **fpermi;
	double *lambdaConst;

	int num_qp;
	double lambdas[100][4], weight[100];
		
	/************************************************** stiffness matrix A *****************************************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF[0].col*elementDOF[1].col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(elementDOF[1].col, elementDOF[0].col, &lA);

	num_qp = getNumQuadPoints_ShunnWilliams(elementDOF[0].dop+elementDOF[1].dop-1, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
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

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				val1 = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					lagrange3d_basis1(lambdas[i1], grd_lambda, k1, elementDOF[1].dop, phi1);
					axy_array(3, lambdas[i1][3], vertices[3], x);
					for(i2=0;i2<3;i2++)
						axpy_array(3, lambdas[i1][i2], vertices[i2], x);
					huangQuadcurl3d_basis(x, xK, lambdas[i1], grd_lambda, vertices, nvf, eorien, fpermi, k2, phi2);
					val1 += vol*weight[i1] * dot_array(3, phi1, phi2);
				}
				lA.val[k1][k2] += val1;
			} // k2
		} // k1

		l = elementDOF[0].col*elementDOF[1].col * k;
		for (i = 0; i<elementDOF[1].col; i++)
		{
			for (j = 0; j<elementDOF[0].col; j++)
			{
				ia[l] = elementDOF[1].val[k][i];
				ja[l] = elementDOF[0].val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // j
		} // i
	} // k
	free_dden_matrix(&lA);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	dIJtoCSReps(A, ia, ja, va, N, 0, 0, 0);
}

/**
 * \fn void assembleRHSHuang3d(dvector *b, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, void (*f)(double *, double *))
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *BT pointer to stiffness matrix
 * \param *C pointer to stiffness matrix
 * \param *b pointer to right hand side
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
 * \return void
 */
void assembleRHSHuang3d(dvector *b, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, void (*f)(double *, double *))
{
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[9], phi2[9];
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
		
	dvector lb;
	create_dvector(elementDOF[0].col, &lb);
	/************************************************** right hand side b *****************************************************************/
	create_dvector(elementDOF[0].dof, b);
	num_qp = getNumQuadPoints_ShunnWilliams(9, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
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

		init_dvector(&lb, 0.0);
		for (i = 0; i<elementDOF[0].col; i++)
		{
			for (i1 = 0; i1<num_qp; i1++)
			{
				axy_array(3, lambdas[i1][3], vertices[3], x);
				for(j=0;j<3;j++)
					axpy_array(3, lambdas[i1][j], vertices[j], x);
				huangQuadcurl3d_basis(x, xK, lambdas[i1], grd_lambda, vertices, nvf, eorien, fpermi, i, phi1);
				// quadcurl3d_f(x, phi2);
				f(x, phi2);
				lb.val[i] += vol*weight[i1] * dot_array(3, phi1, phi2);
			} // i1
		} // k1

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			b->val[i] += lb.val[k1];
		} // k1 
	} // k
	free_dvector(&lb);
}

/**
 * \fn void assembleDistribCurldivTracelessDGNedelecHybrid3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, int curlfemtype)
 * \brief assemble stiffness matrix
 *  The basis of Tracelss tensor: e1 e2^T, e1 e3^T, e2 e1^T, e2 e3^T, e3 e1^T, e3 e2^T, (e1 e1^T - e2 e2^T)/sqrt(2), (e1 e1^T + e2 e2^T - 2 e3 e3^T)/sqrt(6) 
 * \param *A pointer to stiffness matrix
 * \param *BT pointer to stiffness matrix
 * \param *C pointer to stiffness matrix
 * \param *b pointer to right hand side
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
 * \return void
 */
void assembleDistribCurldivTracelessDGNedelecHybrid3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, int curlfemtype)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Cx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	
	int i,j,k,l,li;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[3], phi2[3], val[16], val1;
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], **grd_lambda, **nv, *nvf, *t1vf, *t2vf, *etv[6], **vertices;
	double vol, s, h[6];
	int fi[3];
	int count;
	short *forien, *eorien;
	int **fpermi, **eperm;
	double *lambdaConst;

	int num_qp;
	double lambdas[100][4], weight[100];
	int num_qp2;
	double lambdasT[4]; 
	double lambdas2[100][3], weight2[100];
		
	/************************************************** stiffness matrix A *****************************************************************/
	int *ia, *ja;
	double *va;
	int N1 = 8*elementDOF[0].col*elementDOF[1].col*elements->row;
	int N2 = 64*elementDOF[0].col*elementDOF[2].col*elements->row;
	int N = N1+N2;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA[16]; // local A
	for(i=0;i<8;i++)
		create_dden_matrix(elementDOF[1].col, elementDOF[0].col, lA+i);

	num_qp = getNumQuadPoints_ShunnWilliams(elementDOF[0].dop+elementDOF[1].dop-2, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	num_qp2 = getNumQuadPoints_ShunnWilliams(elementDOF[0].dop+elementDOF[1].dop-1, 2); // the number of numerical intergation points
	init_ShunnWilliams2d(num_qp2, lambdas2, weight2); // Shunn-Williams intergation initial

	for (k = 0; k<elements->row; k++){
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

		for(i=0;i<8;i++)
			init_dden_matrix(lA+i, 0.0);
		for (k1 = 0; k1<elementDOF[1].col; k1++){
			for (k2 = 0; k2<elementDOF[0].col; k2++){
				init_array(8, val, 0);
				for (i1 = 0; i1<num_qp; i1++){
					if(curlfemtype == 1)
						nedelec1st3d_basisCurl(lambdas[i1], grd_lambda, eorien, fpermi, k1, elementDOF[1].dop, phi1);
					else
						nedelec2nd3d_basisCurl(lambdas[i1], grd_lambda, eperm, k1, elementDOF[1].dop, phi1);
					lagrange3d_basis1(lambdas[i1], grd_lambda, k2, elementDOF[0].dop, phi2);
					val[0] += vol*weight[i1] * phi1[0]*phi2[1];
					val[1] += vol*weight[i1] * phi1[0]*phi2[2];
					val[2] += vol*weight[i1] * phi1[1]*phi2[0];
					val[3] += vol*weight[i1] * phi1[1]*phi2[2];
					val[4] += vol*weight[i1] * phi1[2]*phi2[0];
					val[5] += vol*weight[i1] * phi1[2]*phi2[1];
					val[6] += vol*weight[i1] * (phi1[0]*phi2[0] - phi1[1]*phi2[1]) / sqrt(2);
					val[7] += vol*weight[i1] * (phi1[0]*phi2[0] + phi1[1]*phi2[1] - 2*phi1[2]*phi2[2]) / sqrt(6);
				}

				for (i = 0; i < elementFace->col; i++){
					face = elementFace->val[k][i];
					// t1vf = faces->t1vector[face];
					// t2vf = faces->t2vector[face];
					s = faces->area[face];
					for(j=0;j<3;j++)
						fi[j]=find_iarray(4, elements->val[k], faces->val[face][j]);
					init_array(4, lambdasT, 0);
					for (i1 = 0; i1<num_qp2; i1++){
						for(j=0;j<3;j++) lambdasT[fi[j]] = lambdas2[i1][j];
						
						if(curlfemtype == 1)
							nedelec1st3d_basisCurl(lambdasT, grd_lambda, eorien, fpermi, k1, elementDOF[1].dop, phi1);
						else
							nedelec2nd3d_basisCurl(lambdasT, grd_lambda, eperm, k1, elementDOF[1].dop, phi1);
						lagrange3d_basis(lambdasT, k2, elementDOF[0].dop, phi2);
					
						val1 = - s*weight2[i1] * phi2[0] * dot_array(3, nv[i], phi1);
						val[0] += val1 * nv[i][0]*nv[i][1];
						val[1] += val1 * nv[i][0]*nv[i][2];
						val[2] += val1 * nv[i][1]*nv[i][0];
						val[3] += val1 * nv[i][1]*nv[i][2];
						val[4] += val1 * nv[i][2]*nv[i][0];
						val[5] += val1 * nv[i][2]*nv[i][1];
						val[6] += val1 * (nv[i][0]*nv[i][0] - nv[i][1]*nv[i][1]) / sqrt(2);
						val[7] += val1 * (nv[i][0]*nv[i][0] + nv[i][1]*nv[i][1] - 2*nv[i][2]*nv[i][2]) / sqrt(6);
					}
				}

				for(i=0;i<8;i++)
				 lA[i].val[k1][k2] += val[i];
			} // k2
		} // k1

		l = 8*elementDOF[0].col*elementDOF[1].col * k;
		for (k1 = 0; k1<elementDOF[1].col; k1++){
			for (k2 = 0; k2<elementDOF[0].col; k2++){
				for (i = 0; i<8; i++){
					ia[l] = elementDOF[1].val[k][k1];
					ja[l] = elementDOF[0].col*(8*k+i) + k2;
					va[l] = lA[i].val[k1][k2];
					l++;
				} // i
			} // k2
		} // k1
	} // k

	for(i=0;i<8;i++)
		free_dden_matrix(lA+i);

/************** Lagrange multiplier part **************/
	for(i=0;i<16;i++)
		create_dden_matrix(elementDOF[2].col, elementDOF[0].col, lA+i);

	for (k = 0; k<elements->row; k++){
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

		for (i = 0; i < elementFace->col; i++){
			face = elementFace->val[k][i];
			t1vf = faces->t1vector[face];
			t2vf = faces->t2vector[face];
			s = faces->area[face];
			for(j=0;j<3;j++)
				fi[j]=find_iarray(4, elements->val[k], faces->val[face][j]);
					
			for(j=0;j<16;j++)
				init_dden_matrix(lA+j, 0.0);
			init_array(4, lambdasT, 0);
			for (k1 = 0; k1<elementDOF[2].col; k1++){
				for (k2 = 0; k2<elementDOF[0].col; k2++){
					init_array(16, val, 0);
					for (i1 = 0; i1<num_qp2; i1++){
						// lambdasT[i]=0;
						for(j=0;j<3;j++) lambdasT[fi[j]] = lambdas2[i1][j];
						
						lagrange_basis(lambdas2[i1], k1, elementDOF[2].dop, phi1);
						lagrange3d_basis(lambdasT, k2, elementDOF[0].dop, phi2);
						val1 = - s*weight2[i1] * phi2[0] * phi1[0];
						val[0] += val1 * t1vf[0]*nv[i][1];
						val[1] += val1 * t1vf[0]*nv[i][2];
						val[2] += val1 * t1vf[1]*nv[i][0];
						val[3] += val1 * t1vf[1]*nv[i][2];
						val[4] += val1 * t1vf[2]*nv[i][0];
						val[5] += val1 * t1vf[2]*nv[i][1];
						val[6] += val1 * (t1vf[0]*nv[i][0] - t1vf[1]*nv[i][1]) / sqrt(2);
						val[7] += val1 * (t1vf[0]*nv[i][0] + t1vf[1]*nv[i][1] - 2*t1vf[2]*nv[i][2]) / sqrt(6);
						val[8] += val1 * t2vf[0]*nv[i][1];
						val[9] += val1 * t2vf[0]*nv[i][2];
						val[10] += val1 * t2vf[1]*nv[i][0];
						val[11] += val1 * t2vf[1]*nv[i][2];
						val[12] += val1 * t2vf[2]*nv[i][0];
						val[13] += val1 * t2vf[2]*nv[i][1];
						val[14] += val1 * (t2vf[0]*nv[i][0] - t2vf[1]*nv[i][1]) / sqrt(2);
						val[15] += val1 * (t2vf[0]*nv[i][0] + t2vf[1]*nv[i][1] - 2*t2vf[2]*nv[i][2]) / sqrt(6);
					}
					for(j=0;j<16;j++)
						lA[j].val[k1][k2] += val[j];
				} // k2
			} // k1
		
			l = N1+16*elementDOF[0].col*elementDOF[2].col * (4*k + i);
			for (k1 = 0; k1<elementDOF[2].col; k1++){
				for (k2 = 0; k2<elementDOF[0].col; k2++){
					for (j = 0; j<16; j++){
						ia[l] = elementDOF[1].dof + elementDOF[2].col*(2*face + j/8) + k1;
						ja[l] = elementDOF[0].col*(8*k + j%8) + k2;
						va[l] = lA[j].val[k1][k2];
						l++;
					} // j
				} // k2
			} // k1
		} // i, face
	} // k
	
	for(i=0;i<16;i++)
		free_dden_matrix(lA+i);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	dIJtoCSReps(A, ia, ja, va, N, 0, 0, 0);
}

/**
* \fn void jumpOperatorMorley3d(double *prt_lambdas, int face, ELEMENT *elements, idenmat *elementFace, FACE *faces, ELEMENT_DOF *elementDOF, int node, double *jump)
* \brief the jump of uh for Morley element in 3d
* \param prt_lambdas pointer to the area coordiante
* \param face the index of current face
* \param *elements pointer to the structure of the triangulation
* \param *elementFace pointer to relation between tetrahedrons and faces: each row stores 4 faces index
* \param *faces pointer to faces: the first three columns store the three vertices, the fourth and fifth columns store the affiliated elements
the fifth column stores -1 if the face is on boundary
* \param *elementDOF pointer to relation between elements and DOFs
* \param node index of current node variable
* \param *jump pointer to the result of jump operator
* \return void
*/
void jumpOperatorMorley3d(double *prt_lambdas, int face, ELEMENT *elements, idenmat *elementFace, FACE *faces, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i, j, k, l, m;
	int nodeindex, ei, fi;
	int element, v[3];
	double lambdas[4];
	double vol, **grd_lambda, **nv, *nvf[4];
	double phi;

	*jump = 0;

	if (node<0 && node >= elementDOF->dof)
		return;

	v[0] = faces->val[face][0];
	v[1] = faces->val[face][1];
	v[2] = faces->val[face][2];

	for (m = 0; m < 2; m++)
	{
		element = faces->val[face][3 + m];
		if (element == -1)
			continue;

		for (nodeindex = 0; nodeindex<elementDOF->col; nodeindex++)
		{
			if (elementDOF->val[element][nodeindex] == node)
				break;
		}
		if (nodeindex == elementDOF->col)
			continue;

		for (fi = 0; fi<4; fi++)
		{
			if (elementFace->val[element][fi] == face)
				break;
		}

		for (i = 0; i<4; i++)
		{
			if (elements->val[element][i] == v[0])
				break;
		}
		for (j = 0; j<4; j++)
		{
			if (elements->val[element][j] == v[1])
				break;
		}
		for (k = 0; k<4; k++)
		{
			if (elements->val[element][k] == v[2])
				break;
		}
		l = 6 - i - j - k;
		lambdas[i] = prt_lambdas[0];
		lambdas[j] = prt_lambdas[1];
		lambdas[k] = prt_lambdas[2];
		lambdas[l] = 0;

		vol = elements->vol[element];
		grd_lambda = elements->gradLambda[element];
		nv = elements->nvector[element];
		for (i = 0; i < elementFace->col; i++)
			nvf[i] = faces->nvector[elementFace->val[element][i]];

		morley3d_basis(lambdas, vol, grd_lambda , nv, nvf, nodeindex, &phi);

		*jump += phi*dot_array(3, nv[fi], faces->nvector[face]);
	}
}

/**
* \fn void averageOperatorNormalDerivativeMorley3d(double *prt_lambdas, int face, ELEMENT *elements, idenmat *elementFace, FACE *faces, ELEMENT_DOF *elementDOF, int node, double *jump)
* \brief the average of partial_{n_F} uh for Morley element in 3d
* \param prt_lambdas pointer to the area coordiante
* \param face the index of current face
* \param *elements pointer to the structure of the triangulation
* \param *elementFace pointer to relation between tetrahedrons and faces: each row stores 4 faces index
* \param *faces pointer to faces: the first three columns store the three vertices, the fourth and fifth columns store the affiliated elements
the fifth column stores -1 if the face is on boundary
* \param *elementDOF pointer to relation between elements and DOFs
* \param node index of current node variable
* \param *jump pointer to the result of jump operator
* \return void
*/
void averageOperatorNormalDerivativeMorley3d(double *prt_lambdas, int face, ELEMENT *elements, idenmat *elementFace, FACE *faces, ELEMENT_DOF *elementDOF, int node, double *ave)
{
	int i, j, k, l, m;
	int nodeindex, ei, fi;
	int element, v[3];
	double lambdas[4];
	double vol, **grd_lambda, **nv, *nvf[4];
	double phi[3];

	*ave = 0;

	if (node<0 && node >= elementDOF->dof)
		return;

	v[0] = faces->val[face][0];
	v[1] = faces->val[face][1];
	v[2] = faces->val[face][2];

	for (m = 0; m < 2; m++)
	{
		element = faces->val[face][3 + m];
		if (element == -1)
			continue;

		for (nodeindex = 0; nodeindex<elementDOF->col; nodeindex++)
		{
			if (elementDOF->val[element][nodeindex] == node)
				break;
		}
		if (nodeindex == elementDOF->col)
			continue;

		for (fi = 0; fi<4; fi++)
		{
			if (elementFace->val[element][fi] == face)
				break;
		}

		for (i = 0; i<4; i++)
		{
			if (elements->val[element][i] == v[0])
				break;
		}
		for (j = 0; j<4; j++)
		{
			if (elements->val[element][j] == v[1])
				break;
		}
		for (k = 0; k<4; k++)
		{
			if (elements->val[element][k] == v[2])
				break;
		}
		l = 6 - i - j - k;
		lambdas[i] = prt_lambdas[0];
		lambdas[j] = prt_lambdas[1];
		lambdas[k] = prt_lambdas[2];
		lambdas[l] = 0;

		vol = elements->vol[element];
		grd_lambda = elements->gradLambda[element];
		nv = elements->nvector[element];
		for (i = 0; i < elementFace->col; i++)
			nvf[i] = faces->nvector[elementFace->val[element][i]];

		morley3d_basis1(lambdas, vol, grd_lambda, nv, nvf, nodeindex, phi);

		*ave += dot_array(3, phi, faces->nvector[face]);
	} // m

	if (faces->val[face][4] != -1)
		*ave *= 0.5;
}

/**
 * \fn void jumpOperatorVectorTensor3d(double *prt_lambdas, int face, ELEMENT *elements, idenmat *elementFace, FACE *faces, ELEMENT_DOF *elementDOF, int node, double *jump)
 * \brief the jump operator for vector otimes normal derivative in 3d
 * \param prt_lambdas pointer to the area coordiante
 * \param face the index of current face
 * \param *elements pointer to the structure of the triangulation
 * \param *elementFace pointer to relation between tetrahedrons and faces: each row stores 4 faces index
 * \param *faces pointer to faces: the first three columns store the three vertices, the fourth and fifth columns store the affiliated elements
                                   the fifth column stores -1 if the face is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *jump pointer to the result of jump operator
 * \return void
 */
void jumpOperatorVectorTensor3d(double *prt_lambdas, int face, ELEMENT *elements, idenmat *elementFace, FACE *faces, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i,j,k,l;
	int nodeindex, ei, fi, ni3;
	int element, v[3];
	double lambdas[4];
	double *nv;
	double phi, vec[3];

	for(i=0;i<6;i++)
		jump[i]=0;

	if(node<0 && node>=elementDOF->dof*3)
		return;

	element=node/(elementDOF->col*3);
	ni3=node%(elementDOF->col*3);

	for(fi=0;fi<4;fi++)
	{
		if(elementFace->val[element][fi]==face)
			break;
	}

	if(fi==4)
		return;

	nv=elements->nvector[element][fi];
	vec[0]=0;vec[1]=0;vec[2]=0;
	
	if(elementDOF->dop==0)
	{
		if(ni3<elementDOF->col)
		{
			vec[0]=1;
		}
		else if(ni3<elementDOF->col*2)
		{
			vec[1]=1;
		}
		else
		{
			vec[2]=1;
		}
		tensorproduct3d_array(vec, nv, jump);
		return;
	}

	nodeindex=node%elementDOF->col;

	if(nodeindex<4 && nodeindex>=0)
	{
		if(nodeindex==fi)
			return;
	}
	else if(nodeindex<6*elementDOF->dop-2)
	{
		ei=(nodeindex-4)/(elementDOF->dop-1);
		edge2vv3d(ei, v);
		if(v[0]==fi || v[1]==fi)
			return;
	}
	else if(nodeindex<2*(elementDOF->dop*elementDOF->dop+1))
	{
		if((nodeindex-(6*elementDOF->dop-2))/((elementDOF->dop-1)*(elementDOF->dop-2)/2) != fi)
			return;
	}
	else
		return;

	v[0]=faces->val[face][0];
	v[1]=faces->val[face][1];
	v[2]=faces->val[face][2];
	for(i=0;i<4;i++)
	{
		if(elements->val[element][i]==v[0])
			break;
	}
	for(j=0;j<4;j++)
	{
		if(elements->val[element][j]==v[1])
			break;
	}
	for(k=0;k<4;k++)
	{
		if(elements->val[element][k]==v[2])
			break;
	}
	l=6-i-j-k;
	lambdas[i]=prt_lambdas[0];
	lambdas[j]=prt_lambdas[1];
	lambdas[k]=prt_lambdas[2];
	lambdas[l]=0;

	lagrange3d_basis(lambdas, nodeindex, elementDOF->dop, &phi);

	if(ni3<elementDOF->col)
	{
		vec[0]=phi;
	}
	else if(ni3<elementDOF->col*2)
	{
		vec[1]=phi;
	}
	else
	{
		vec[2]=phi;
	}
	tensorproduct3d_array(vec, nv, jump);
}

/**
* \fn void jumpOperatorVector3d(double *prt_lambdas, int face, ELEMENT *elements, idenmat *elementFace, FACE *faces, ELEMENT_DOF *elementDOF, int node, double *jump)
* \brief the jump operator for vector otimes normal derivative in 3d
* \param prt_lambdas pointer to the area coordiante
* \param face the index of current face
* \param *elements pointer to the structure of the triangulation
* \param *elementFace pointer to relation between tetrahedrons and faces: each row stores 4 faces index
* \param *faces pointer to faces: the first three columns store the three vertices, the fourth and fifth columns store the affiliated elements
the fifth column stores -1 if the face is on boundary
* \param *elementDOF pointer to relation between elements and DOFs
* \param node index of current node variable
* \param *jump pointer to the result of jump operator
* \return void
*/
void jumpOperatorVector3d(double *prt_lambdas, int face, ELEMENT *elements, idenmat *elementFace, FACE *faces, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i, j, k, l;
	int nodeindex, ei, fi;
	int element, v[3];
	double lambdas[4];
	double *nv, *nve;
	double phi;

	for (i = 0; i<3; i++)
		jump[i] = 0;

	if (node<0 && node >= elementDOF->dof * 3)
		return;

	element = (node % elementDOF->dof) / elementDOF->col;
	int li = node / elementDOF->dof;

	for (fi = 0; fi<4; fi++)
	{
		if (elementFace->val[element][fi] == face)
			break;
	}

	if (fi == 4)
		return;

	nv = elements->nvector[element][fi];
	nve = faces->nvector[face];
	double sgn = nve[0] * nv[0] + nve[1] * nv[1] + nve[2] * nv[2];
	
	if (elementDOF->dop == 0)
	{
		jump[li] = sgn;
		return;
	}

	nodeindex = node%elementDOF->col;

	if (nodeindex<4 && nodeindex >= 0)
	{
		if (nodeindex == fi)
			return;
	}
	else if (nodeindex<6 * elementDOF->dop - 2)
	{
		ei = (nodeindex - 4) / (elementDOF->dop - 1);
		edge2vv3d(ei, v);
		if (v[0] == fi || v[1] == fi)
			return;
	}
	else if (nodeindex<2 * (elementDOF->dop*elementDOF->dop + 1))
	{
		if ((nodeindex - (6 * elementDOF->dop - 2)) / ((elementDOF->dop - 1)*(elementDOF->dop - 2) / 2) != fi)
			return;
	}
	else
		return;

	v[0] = faces->val[face][0];
	v[1] = faces->val[face][1];
	v[2] = faces->val[face][2];
	for (i = 0; i<4; i++)
	{
		if (elements->val[element][i] == v[0])
			break;
	}
	for (j = 0; j<4; j++)
	{
		if (elements->val[element][j] == v[1])
			break;
	}
	for (k = 0; k<4; k++)
	{
		if (elements->val[element][k] == v[2])
			break;
	}
	l = 6 - i - j - k;
	lambdas[i] = prt_lambdas[0];
	lambdas[j] = prt_lambdas[1];
	lambdas[k] = prt_lambdas[2];
	lambdas[l] = 0;

	lagrange3d_basis(lambdas, nodeindex, elementDOF->dop, &phi);

	jump[li] = phi*sgn;
}

/**
* \fn void jumpOperatorVector3dOld(double *prt_lambdas, int face, ELEMENT *elements, idenmat *elementFace, FACE *faces, ELEMENT_DOF *elementDOF, int node, double *jump)
* \brief the jump operator for vector otimes normal derivative in 3d
* \param prt_lambdas pointer to the area coordiante
* \param face the index of current face
* \param *elements pointer to the structure of the triangulation
* \param *elementFace pointer to relation between tetrahedrons and faces: each row stores 4 faces index
* \param *faces pointer to faces: the first three columns store the three vertices, the fourth and fifth columns store the affiliated elements
the fifth column stores -1 if the face is on boundary
* \param *elementDOF pointer to relation between elements and DOFs
* \param node index of current node variable
* \param *jump pointer to the result of jump operator
* \return void
*/
void jumpOperatorVector3dOld(double *prt_lambdas, int face, ELEMENT *elements, idenmat *elementFace, FACE *faces, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i, j, k, l;
	int nodeindex, ei, fi, ni3;
	int element, v[3];
	double lambdas[4];
	double *nv, *nve;
	double phi;

	for (i = 0; i<3; i++)
		jump[i] = 0;

	if (node<0 && node >= elementDOF->dof * 3)
		return;

	element = node / (elementDOF->col * 3);
	ni3 = node % (elementDOF->col * 3);

	for (fi = 0; fi<4; fi++)
	{
		if (elementFace->val[element][fi] == face)
			break;
	}

	if (fi == 4)
		return;

	nv = elements->nvector[element][fi];
	nve = faces->nvector[face];
	double sgn = nve[0] * nv[0] + nve[1] * nv[1] + nve[2] * nv[2];
	
	if (elementDOF->dop == 0)
	{
		if (ni3<elementDOF->col)
		{
			jump[0] = sgn;
		}
		else if (ni3<elementDOF->col * 2)
		{
			jump[1] = sgn;
		}
		else
		{
			jump[2] = sgn;
		}
		return;
	}

	nodeindex = node%elementDOF->col;

	if (nodeindex<4 && nodeindex >= 0)
	{
		if (nodeindex == fi)
			return;
	}
	else if (nodeindex<6 * elementDOF->dop - 2)
	{
		ei = (nodeindex - 4) / (elementDOF->dop - 1);
		edge2vv3d(ei, v);
		if (v[0] == fi || v[1] == fi)
			return;
	}
	else if (nodeindex<2 * (elementDOF->dop*elementDOF->dop + 1))
	{
		if ((nodeindex - (6 * elementDOF->dop - 2)) / ((elementDOF->dop - 1)*(elementDOF->dop - 2) / 2) != fi)
			return;
	}
	else
		return;

	v[0] = faces->val[face][0];
	v[1] = faces->val[face][1];
	v[2] = faces->val[face][2];
	for (i = 0; i<4; i++)
	{
		if (elements->val[element][i] == v[0])
			break;
	}
	for (j = 0; j<4; j++)
	{
		if (elements->val[element][j] == v[1])
			break;
	}
	for (k = 0; k<4; k++)
	{
		if (elements->val[element][k] == v[2])
			break;
	}
	l = 6 - i - j - k;
	lambdas[i] = prt_lambdas[0];
	lambdas[j] = prt_lambdas[1];
	lambdas[k] = prt_lambdas[2];
	lambdas[l] = 0;

	lagrange3d_basis(lambdas, nodeindex, elementDOF->dop, &phi);

	if (ni3<elementDOF->col)
	{
		jump[0] = phi*sgn;
	}
	else if (ni3<elementDOF->col * 2)
	{
		jump[1] = phi*sgn;
	}
	else
	{
		jump[2] = phi*sgn;
	}
}

/**
 * \fn void getElementFaceEdgeGeoInfo(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes)
 * \brief compute the geometric information of elements, faces and edges
 * \param *elements stores 4 nodes corresponding to the element
 * \param *faces the first three columns store the three vertices corresponding to the face; 
 *				 the 4th and 5th columns store the elements which the face belongs to;
 *				 if the face is a boundary, the 5th column will stores -1;
 *				 the first column is in ascend order.
 * \param *edges stores the two vertice corresponding to the edge
 *				 the first column is in ascend order.
 * \param *nodes the first column stores the x coordinate of points, the second column stores the y coordinate of points, 
 *               the third column stores the z coordinate of points
 * \return void
 */
void getElementFaceEdgeGeoInfo(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes)
{
	int i,j,k;
	int vertex, v[3], face, lface[3], edge;
	double tet[4][3], a[3], b[3], len, vol;
	double t01[3], t12[3], t23[3], t30[3];
	int *perm, *permi;

	for(i=0;i<elements->row;i++)
	{
		for(j=0;j<4;j++)
		{
			vertex=elements->val[i][j];
			copy_array(3, nodes->val[vertex], tet[j]);
			copy_array(3, nodes->val[vertex], elements->vertices[i][j]); 
		}
		elements->vol[i]=volume(elements->vertices[i]);
		vol = elements->vol[i];

		for(j=0;j<3;j++)
			elements->barycenter[i][j] = (tet[0][j]+tet[1][j]+tet[2][j]+tet[3][j]) / 4.0;

		for(j=0;j<4;j++)
			axpbyz_array(3, 4.0/3.0, elements->barycenter[i], -1.0/3.0, tet[j], elements->bcFace[i][j]);

		if(elements->vol[i]<0)
		{
			printf("The volume of the %d-th tetrahedron is %e.\n", i, elements->vol[i]);
			exit(1);
		}

		axpyz_array(3, -1.0, tet[0], tet[1], t01);
		axpyz_array(3, -1.0, tet[1], tet[2], t12);
		axpyz_array(3, -1.0, tet[2], tet[3], t23);
		axpyz_array(3, -1.0, tet[3], tet[0], t30);

		cross_array(t23, t12, elements->gradLambda[i][0]);
		ax_array(3, 1.0/(6*vol), elements->gradLambda[i][0]);

		cross_array(t23, t30, elements->gradLambda[i][1]);
		ax_array(3, 1.0/(6*vol), elements->gradLambda[i][1]);
		
		cross_array(t01, t30, elements->gradLambda[i][2]);
		ax_array(3, 1.0/(6*vol), elements->gradLambda[i][2]);
		
		cross_array(t01, t12, elements->gradLambda[i][3]);
		ax_array(3, 1.0/(6*vol), elements->gradLambda[i][3]);

		for(j=0;j<4;j++)
		{
			elements->lambdaConst[i][j] = 1 - dot_array(3, tet[j], elements->gradLambda[i][j]);

			face2vertices3d(j, v);
			axpyz_array(3, -1.0, tet[v[0]], tet[v[1]], a);
			axpyz_array(3, -1.0, tet[v[0]], tet[v[2]], b);
			cross_array(a, b, elements->nvector[i][j]);
			len=sqrt(dot_array(3, elements->nvector[i][j], elements->nvector[i][j]));
			ax_array(3, 1.0/len, elements->nvector[i][j]);			

			face = elementFace->val[i][j];
			for (k = 0; k < 3; k++)
				lface[k] = elements->val[i][v[k]];
			getPermutation(lface, faces->val[face], elements->fperm[i][j], 3);
			perm = elements->fperm[i][j]; permi = elements->fpermi[i][j];
			for(k=0;k<3;k++) permi[perm[k]]=k;
		} // j

		for (j = 0; j < 6; j++)
		{
			edge2vv3d(j, v);
			edge = elementEdge->val[i][j];
			if (elements->val[i][v[0]] == edges->val[edge][0])
			{
				elements->eperm[i][j][0] = 0;
				elements->eperm[i][j][1] = 1;
				elements->eorien[i][j] = 1;
			}
			else
			{
				elements->eperm[i][j][0] = 1;
				elements->eperm[i][j][1] = 0;
				elements->eorien[i][j] = -1;
			}
		}
	} // elements i

	for(i=0;i<faces->row;i++)
	{
		for(j=0;j<3;j++)
		{
			vertex=faces->val[i][j];
			copy_array(3, nodes->val[vertex], tet[j]);
			// for(k=0;k<3;k++)
			// 	tet[j][k]=nodes->val[vertex][k];
		}

		axpbyz_array(3, 1.0/3.0, tet[0], 1.0/3.0, tet[1], faces->barycenter[i]);
		axpy_array(3, 1.0/3.0, tet[2], faces->barycenter[i]);
		
		axpyz_array(3, -1.0, tet[0], tet[1], faces->t01[i]);
		axpyz_array(3, -1.0, tet[0], tet[2], faces->t02[i]);
		axpyz_array(3, -1.0, tet[1], tet[2], faces->t12[i]);
		
		cross_array(faces->t01[i], faces->t02[i], faces->nvector[i]);
		len=sqrt(dot_array(3, faces->nvector[i], faces->nvector[i]));
		ax_array(3, 1.0/len, faces->nvector[i]);

		orthocomplement_array(faces->nvector[i], faces->t1vector[i], faces->t2vector[i]);

		faces->area[i] = len/2.0;		
	} // faces i
	
	for(i=0;i<edges->row;i++)
	{
		for(j=0;j<2;j++)
		{
			vertex=edges->val[i][j];
			copy_array(3, nodes->val[vertex], tet[j]);
			// for(k=0;k<3;k++)
			// 	tet[j][k]=nodes->val[vertex][k];
		}

		axpyz_array(3, -1.0, tet[0], tet[1], edges->tvector[i]);
		edges->length[i]=sqrt(dot_array(3, edges->tvector[i], edges->tvector[i]));
		ax_array(3, 1.0/edges->length[i], edges->tvector[i]);

		orthocomplement_array(edges->tvector[i], edges->n1vector[i], edges->n2vector[i]);
	} // edges i

	for(i=0;i<elements->row;i++)
	{
		for(j=0;j<4;j++)
		{
			face = elementFace->val[i][j];

			if(dot_array(3, elements->nvector[i][j], faces->nvector[face]) > 0)
				elements->forien[i][j] = 1;
			else
				elements->forien[i][j] = -1;
		}
	}
}

/**
* \fn void getFreenodesInfo(ELEMENT_DOF *elementDOF)
* \brief get freenodes information without boundary conditions
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfo(ELEMENT_DOF *elementDOF)
{
	int i, j, k, estride, fstride, nnf;

	int dof = elementDOF->dof;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(dof, nfFlag);
	create_ivector(dof, index);

	nnf = 0; // number of non-free nodes

	create_ivector(nnf, nfreenodes);
	create_ivector(dof - nnf, freenodes);

	for (i = 0; i<dof; i++){
		freenodes->val[i] = i;
		index->val[i] = i;
	}
}

/**
* \fn void getFreenodesInfoFaces(FACE *faces, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of piecewise polynomial on faces
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoFaces(FACE *faces, ELEMENT_DOF *elementDOF)
{
	int i, j, k, fstride, nnf;

	int nf = faces->row;
	int dof = elementDOF->dof;
	int dop = elementDOF->dop;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(dof, nfFlag);
	create_ivector(dof, index);

	nnf = 0; // number of non-free nodes
	fstride = (dop+1)*(dop+2)/2;
	for (k = 0; k<nf; k++)
	{
		if (faces->bdFlag[k] == 1 || faces->bdFlag[k] == 2 || faces->bdFlag[k] == 3 || faces->bdFlag[k] == 4) // Dirichlet boundary
		{
			for (i = 0; i < fstride; i++)
				nfFlag->val[k*fstride + i] = 1;
			nnf += fstride;
		}
	}


	create_ivector(nnf, nfreenodes);
	create_ivector(dof - nnf, freenodes);

	j = 0; k = 0;
	for (i = 0; i<dof; i++)
	{
		if (nfFlag->val[i] == 1) //  non-free node
		{
			nfreenodes->val[k] = i;
			index->val[i] = k;
			k++;
		}
		else // free variable
		{
			freenodes->val[j] = i;
			index->val[i] = j;
			j++;
		}
	}
}

/**
* \fn void getFreenodesInfoFacesVector(FACE *faces, ELEMENT_DOF *elementDOF, int n)
* \brief get freenodes information of piecewise polynomial on faces
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param *elementDOF pointer to relation between elements and DOFs
* \param n the dimension of vector
* \return void
*/
void getFreenodesInfoFacesVector(FACE *faces, ELEMENT_DOF *elementDOF, int n)
{
	int i, j, k, fstride, nnf;

	int nf = faces->row;
	int dof = elementDOF->dof*n;
	int dop = elementDOF->dop;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(dof, nfFlag);
	create_ivector(dof, index);

	nnf = 0; // number of non-free nodes
	fstride = (dop+1)*(dop+2)/2*n;
	for (k = 0; k<nf; k++){
		if (faces->bdFlag[k] == 1 || faces->bdFlag[k] == 2 || faces->bdFlag[k] == 3 || faces->bdFlag[k] == 4) // Dirichlet boundary
		{
			for (i = 0; i < fstride; i++)
				nfFlag->val[k*fstride + i] = 1;
			nnf += fstride;
		}
	}


	create_ivector(nnf, nfreenodes);
	create_ivector(dof - nnf, freenodes);

	j = 0; k = 0;
	for (i = 0; i<dof; i++){
		if (nfFlag->val[i] == 1) //  non-free node
		{
			nfreenodes->val[k] = i;
			index->val[i] = k;
			k++;
		}
		else // free variable
		{
			freenodes->val[j] = i;
			index->val[i] = j;
			j++;
		}
	}
}

/**
* \fn void getFreenodesInfoLagrange3d(FACE *faces, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of Lagrange element in three dimensions
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoLagrange3d(FACE *faces, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i, j, k, estride, fstride, nnf;

	int nn = nodes->row;
	int nf = faces->row;
	int ne = edges->row;
	int dof = elementDOF->dof;
	int dop = elementDOF->dop;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(dof, nfFlag);
	create_ivector(dof, index);

	nnf = 0; // number of non-free nodes
	estride = dop - 1;
	fstride = (dop-1)*(dop-2)/2;
	for (i = 0; i < nn; i++)
	{
		if (nodes->bdFlag[i] == 1 || nodes->bdFlag[i] == 2 || nodes->bdFlag[i] == 3 || nodes->bdFlag[i] == 4)
		{
			nfFlag->val[i] = 1;
			nnf++;
		}
	}

	for (j = 0; j<ne; j++)
	{
		if (edges->bdFlag[j] == 1 || edges->bdFlag[j] == 2 || edges->bdFlag[j] == 3 || edges->bdFlag[j] == 4) // Dirichlet boundary
		{
			for (i = 0; i < estride; i++)
			{
				nfFlag->val[nn + j*estride + i] = 1;
			}
			nnf += estride;
		}
	}

	for (k = 0; k<nf; k++)
	{
		if (faces->bdFlag[k] == 1 || faces->bdFlag[k] == 2 || faces->bdFlag[k] == 3 || faces->bdFlag[k] == 4) // Dirichlet boundary
		{
			for (i = 0; i < fstride; i++)
				nfFlag->val[nn + ne*estride + k*fstride + i] = 1;
			nnf += fstride;
		}
	}


	create_ivector(nnf, nfreenodes);
	create_ivector(dof - nnf, freenodes);

	j = 0; k = 0;
	for (i = 0; i<dof; i++)
	{
		if (nfFlag->val[i] == 1) //  non-free node
		{
			nfreenodes->val[k] = i;
			index->val[i] = k;
			k++;
		}
		else // free variable
		{
			freenodes->val[j] = i;
			index->val[i] = j;
			j++;
		}
	}
}

/**
* \fn void getFreenodesInfoMorley3d(FACE *faces, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of Morley element in three dimensions
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
*                                 the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoMorley3d(FACE *faces, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i, j, k, nnf;

	int nn = nodes->row;
	int nf = faces->row;
	int ne = edges->row;
	int dof = elementDOF->dof;
	//	int dop = elementDOF->dop;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(dof, nfFlag);
	create_ivector(dof, index);

	nnf = 0; // number of non-free nodes
			 	
/*	for (j = 0; j<ne; j++)
	{
		if (edges->bdFlag[j] == 1 || edges->bdFlag[j] == 2 || edges->bdFlag[j] == 3 || edges->bdFlag[j] == 4) // Dirichlet boundary
		{
			nfFlag->val[j] = 1;
			nnf += 1;
		}
	}*/
	for (j = 0; j<nf; j++)
	{
		if (faces->bdFlag[j] == 1 || faces->bdFlag[j] == 2 || faces->bdFlag[j] == 3 || faces->bdFlag[j] == 4) // Dirichlet boundary
		{
			nfFlag->val[ne + j] = 1;
			nnf += 1;
		}
	}

	create_ivector(nnf, nfreenodes);
	create_ivector(dof - nnf, freenodes);

	j = 0; k = 0;
	for (i = 0; i<dof; i++)
	{
		if (nfFlag->val[i] == 1) //  non-free node
		{
			nfreenodes->val[k] = i;
			index->val[i] = k;
			k++;
		}
		else // free variable
		{
			freenodes->val[j] = i;
			index->val[i] = j;
			j++;
		}
	}
}

/**
* \fn void getFreenodesInfoNoncfmP1Vector3d(FACE *faces, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of nonconforming P1 element in 3d-vector version
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								  the fourth column stores -1 if the edge is on boundary
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoNoncfmP1Vector3d(FACE *faces, ELEMENT_DOF *elementDOF)
{
	int i, j, k, nnf;

	int nf = faces->row;
	int dof = elementDOF->dof * 3;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(dof, nfFlag);
	create_ivector(dof, index);

	nnf = 0; // number of non-free nodes
	for (k = 0; k<nf; k++)
	{
		if (faces->bdFlag[k] == 1 || faces->bdFlag[k] == 2 || faces->bdFlag[k] == 3 || faces->bdFlag[k] == 4) // Dirichlet boundary
		{
			nfFlag->val[k] = 1;
			nfFlag->val[k+elementDOF->dof] = 1;
			nfFlag->val[k+elementDOF->dof*2] = 1;	
			nnf += 3;
		}
	}


	create_ivector(nnf, nfreenodes);
	create_ivector(dof - nnf, freenodes);

	j = 0; k = 0;
	for (i = 0; i<dof; i++)
	{
		if (nfFlag->val[i] == 1) //  non-free node
		{
			nfreenodes->val[k] = i;
			index->val[i] = k;
			k++;
		}
		else // free variable
		{
			freenodes->val[j] = i;
			index->val[i] = j;
			j++;
		}
	}
}

/**
* \fn void getFreenodesInfoNedelec1st3d(FACE *faces, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of the first kind Nedelec element in three dimensions
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoNedelec1st3d(FACE *faces, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i, j, k, nnf;

	int nn = nodes->row;
	int nf = faces->row;
	int ne = edges->row;
	int dof = elementDOF->dof;
	int dop = elementDOF->dop;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(dof, nfFlag);
	create_ivector(dof, index);

	nnf = 0; // number of non-free nodes
	for (j = 0; j<ne; j++)
	{
		if (edges->bdFlag[j] == 1 || edges->bdFlag[j] == 2 || edges->bdFlag[j] == 3 || edges->bdFlag[j] == 4) // Dirichlet boundary
		{
			for (i = 0; i < dop; i++)
			{
				nfFlag->val[j*dop + i] = 1;
			}
			nnf += dop;
		}
	}

	for (j = 0; j<nf; j++)
	{
		if (faces->bdFlag[j] == 1 || faces->bdFlag[j] == 2 || faces->bdFlag[j] == 3 || faces->bdFlag[j] == 4) // Dirichlet boundary
		{
			for (i = 0; i < (dop - 1)*dop; i++)
				nfFlag->val[ne*dop + j*(dop - 1)*dop + i] = 1;
			nnf += (dop - 1)*dop;
		}
	}


	create_ivector(nnf, nfreenodes);
	create_ivector(dof - nnf, freenodes);

	j = 0; k = 0;
	for (i = 0; i<dof; i++)
	{
		if (nfFlag->val[i] == 1) //  non-free node
		{
			nfreenodes->val[k] = i;
			index->val[i] = k;
			k++;
		}
		else // free variable
		{
			freenodes->val[j] = i;
			index->val[i] = j;
			j++;
		}
	}
}

/**
* \fn void getFreenodesInfoNedelec2nd3d(FACE *faces, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of the first kind Nedelec element in three dimensions
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoNedelec2nd3d(FACE *faces, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i, j, k, nnf;

	int nn = nodes->row;
	int nf = faces->row;
	int ne = edges->row;
	int dof = elementDOF->dof;
	int dop = elementDOF->dop;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(dof, nfFlag);
	create_ivector(dof, index);

	nnf = 0; // number of non-free nodes
	for (j = 0; j<ne; j++)
	{
		if (edges->bdFlag[j] == 1 || edges->bdFlag[j] == 2 || edges->bdFlag[j] == 3 || edges->bdFlag[j] == 4) // Dirichlet boundary
		{
			for (i = 0; i < (dop+1); i++)
			{
				nfFlag->val[j*(dop + 1) + i] = 1;
			}
			nnf += (dop + 1);
		}
	}

	for (j = 0; j<nf; j++)
	{
		if (faces->bdFlag[j] == 1 || faces->bdFlag[j] == 2 || faces->bdFlag[j] == 3 || faces->bdFlag[j] == 4) // Dirichlet boundary
		{
			for (i = 0; i < (dop - 1)*(dop + 1); i++)
				nfFlag->val[ne*(dop + 1) + j*(dop - 1)*(dop + 1) + i] = 1;
			nnf += (dop - 1)*(dop + 1);
		}
	}


	create_ivector(nnf, nfreenodes);
	create_ivector(dof - nnf, freenodes);

	j = 0; k = 0;
	for (i = 0; i<dof; i++)
	{
		if (nfFlag->val[i] == 1) //  non-free node
		{
			nfreenodes->val[k] = i;
			index->val[i] = k;
			k++;
		}
		else // free variable
		{
			freenodes->val[j] = i;
			index->val[i] = j;
			j++;
		}
	}
}

/**
* \fn void getFreenodesInfoCHHcurlHermite3d(FACE *faces, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of Hermite-type Christiansen-Hu-Hu H(curl) element in three dimensions
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoCHHcurlHermite3d(FACE *faces, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i, j, k, nnf;

	int nn = nodes->row;
	int nf = faces->row;
	int ne = edges->row;
	int dof = elementDOF->dof;
	int dop = elementDOF->dop;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(dof, nfFlag);
	create_ivector(dof, index);

	nnf = 0; // number of non-free nodes
	for (i = 0; i < nn; i++)
	{
		if (nodes->bdFlag[i] == 1 || nodes->bdFlag[i] == 2 || nodes->bdFlag[i] == 3 || nodes->bdFlag[i] == 4)
		{
			nfFlag->val[i] = 1;
			nfFlag->val[i + nn] = 1;
			nfFlag->val[i + nn * 2] = 1;
			nnf += 3;
		}
	}

	for (j = 0; j<ne; j++)
	{
		if (edges->bdFlag[j] == 1 || edges->bdFlag[j] == 2 || edges->bdFlag[j] == 3 || edges->bdFlag[j] == 4) // Dirichlet boundary
		{
			for (i = 0; i < dop - 1; i++)
			{
				nfFlag->val[nn * 3 + j*(dop - 1) + i] = 1;
			}
			nnf += (dop - 1);
		}
	}

	for (j = 0; j<nf; j++)
	{
		if (faces->bdFlag[j] == 1 || faces->bdFlag[j] == 2 || faces->bdFlag[j] == 3 || faces->bdFlag[j] == 4) // Dirichlet boundary
		{
			for (i = 0; i < (dop - 1)*(dop + 1); i++)
				nfFlag->val[nn * 3 + ne*(dop - 1) + j*(dop - 1)*(dop + 1) + i] = 1;
			nnf += (dop - 1)*(dop + 1);
		}
	}


	create_ivector(nnf, nfreenodes);
	create_ivector(dof - nnf, freenodes);

	j = 0; k = 0;
	for (i = 0; i<dof; i++)
	{
		if (nfFlag->val[i] == 1) //  non-free node
		{
			nfreenodes->val[k] = i;
			index->val[i] = k;
			k++;
		}
		else // free variable
		{
			freenodes->val[j] = i;
			index->val[i] = j;
			j++;
		}
	}
}

/**
* \fn void getFreenodesInfoHuangGradcurl3d(FACE *faces, EDGE *edges, dennode *nodes, int type, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of the Huang element in three dimensions for grad-curl
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param type the type of boundary condition: 0: free, 1: clamped, 2: simply supported
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoHuangGradcurl3d(FACE *faces, EDGE *edges, dennode *nodes, int type, ELEMENT_DOF *elementDOF)
{
	int i, j, k, nnf;

	int nn = nodes->row;
	int nf = faces->row;
	int ne = edges->row;
	int dof = elementDOF->dof;
	int dop = elementDOF->dop;

	if(type<0)
		type = 0;
	if(type>2)
		type = 2;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(dof, nfFlag);
	create_ivector(dof, index);

	nnf = 0; // number of non-free nodes
	if(type>0) // simply supported
	{
		for (j = 0; j<ne; j++)
		{
			if (edges->bdFlag[j] == 1 || edges->bdFlag[j] == 2 || edges->bdFlag[j] == 3 || edges->bdFlag[j] == 4) // Dirichlet boundary
			{
				nfFlag->val[j] = 1;
				nnf++;
			}
		}
	}
	
	if(type==1) // clamped
	{
		for (j = 0; j<nf; j++)
		{
			if (faces->bdFlag[j] == 1 || faces->bdFlag[j] == 2 || faces->bdFlag[j] == 3 || faces->bdFlag[j] == 4) // Dirichlet boundary
			{
				for (i = 0; i < 2; i++)
					nfFlag->val[ne + j*2 + i] = 1;
				nnf += 2;
			}
		}
	}


	create_ivector(nnf, nfreenodes);
	create_ivector(dof - nnf, freenodes);

	j = 0; k = 0;
	for (i = 0; i<dof; i++)
	{
		if (nfFlag->val[i] == 1) //  non-free node
		{
			nfreenodes->val[k] = i;
			index->val[i] = k;
			k++;
		}
		else // free variable
		{
			freenodes->val[j] = i;
			index->val[i] = j;
			j++;
		}
	}
}

/**
* \fn void getFreenodesInfoHuangZhang3d(FACE *faces, EDGE *edges, dennode *nodes, int type, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of the Huang-Zhang element in three dimensions for grad-curl
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param type the type of boundary condition: 0: free, 1: clamped, 2: simply supported
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoHuangZhang3d(FACE *faces, EDGE *edges, dennode *nodes, int type, ELEMENT_DOF *elementDOF)
{
	int i, j, k, nnf;

	int nn = nodes->row;
	int nf = faces->row;
	int ne = edges->row;
	int dof = elementDOF->dof;
	int dop = elementDOF->dop;

	if(type<0)
		type = 0;
	if(type>2)
		type = 2;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(dof, nfFlag);
	create_ivector(dof, index);

	nnf = 0; // number of non-free nodes
	if(type>0) // simply supported
	{
		for (j = 0; j<ne; j++)
		{
			if (edges->bdFlag[j] == 1 || edges->bdFlag[j] == 2 || edges->bdFlag[j] == 3 || edges->bdFlag[j] == 4) // Dirichlet boundary
			{
				nfFlag->val[j*2] = 1;
				nfFlag->val[j*2+1] = 1;
				nnf += 2;
			}
		}

		for (j = 0; j<nf; j++)
		{
			if (faces->bdFlag[j] == 1 || faces->bdFlag[j] == 2 || faces->bdFlag[j] == 3 || faces->bdFlag[j] == 4) // Dirichlet boundary
			{
				for (i = 0; i < 2; i++)
					nfFlag->val[ne*2 + j*5 + i] = 1;
				nnf += 2;
			}
		}
	}
	
	if(type==1) // clamped
	{
		for (j = 0; j<nf; j++)
		{
			if (faces->bdFlag[j] == 1 || faces->bdFlag[j] == 2 || faces->bdFlag[j] == 3 || faces->bdFlag[j] == 4) // Dirichlet boundary
			{
				for (i = 2; i < 5; i++)
					nfFlag->val[ne*2 + j*5 + i] = 1;
				nnf += 3;
			}
		}
	}


	create_ivector(nnf, nfreenodes);
	create_ivector(dof - nnf, freenodes);

	j = 0; k = 0;
	for (i = 0; i<dof; i++)
	{
		if (nfFlag->val[i] == 1) //  non-free node
		{
			nfreenodes->val[k] = i;
			index->val[i] = k;
			k++;
		}
		else // free variable
		{
			freenodes->val[j] = i;
			index->val[i] = j;
			j++;
		}
	}
}

/**
* \fn void combineFreenodesInfo(ivector *nfFlag1, ivector *nfFlag2, ELEMENT_DOF *elementDOF)
* \brief combines nfFlag1 and nfFlag2 into nfFlag
* \param *nfFlage1 pointer to the first freenode information
* \param *nfFlage2 pointer to the second freenode information
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void combineFreenodesInfo(ivector *nfFlag1, ivector *nfFlag2, ELEMENT_DOF *elementDOF)
{
	int i, j, k;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	int dof1 = nfFlag1->row;
	int dof2 = nfFlag2->row;
	int dof = dof1 + dof2;
	elementDOF->dof = dof;
	elementDOF->val = NULL;
	
	create_ivector(dof, nfFlag);
	create_ivector(dof, index);
	for(i=0;i<dof1;i++)
		nfFlag->val[i] = nfFlag1->val[i];
	for(i=0;i<dof2;i++)
		nfFlag->val[dof1+i] = nfFlag2->val[i];

	int nnf = 0;
	for(i=0;i<dof;i++){
		if(nfFlag->val[i] == 1) nnf++;
	}

	create_ivector(nnf, nfreenodes);
	create_ivector(dof - nnf, freenodes);

	j = 0; k = 0;
	for (i = 0; i<dof; i++)
	{
		if (nfFlag->val[i] == 1){ //  non-free node
			nfreenodes->val[k] = i;
			index->val[i] = k;
			k++;
		}
		else{ // free variable
			freenodes->val[j] = i;
			index->val[i] = j;
			j++;
		}
	}
}

void getFaceEdgeNTtensor(double **fnv, double **ft1v, double **ft2v, double **etv, double **en1v, double **en2v, double (*fnn)[6], double (*fnt1)[6], double (*fnt2)[6], double (*ft1t1)[6], double (*ft2t2)[6], double (*ft1t2)[6], double (*ett)[6], double (*etn1)[6], double (*etn2)[6], double (*en1n1)[6], double (*en2n2)[6], double (*en1n2)[6])
{
	int i;
	// double subscript of tensor: 00 11 22 12 20 01
	for(i=0;i<4;i++)
	{
		tensorproduct3d_array(fnv[i], fnv[i], fnn[i]);
		tensorproduct3d_array(fnv[i], ft1v[i], fnt1[i]);
		tensorproduct3d_array(fnv[i], ft2v[i], fnt2[i]);
		tensorproduct3d_array(ft1v[i], ft1v[i], ft1t1[i]);
		tensorproduct3d_array(ft2v[i], ft2v[i], ft2t2[i]);
		tensorproduct3d_array(ft1v[i], ft2v[i], ft1t2[i]);
	}

	for(i=0;i<6;i++)
	{
		tensorproduct3d_array(etv[i], etv[i], ett[i]);
		tensorproduct3d_array(etv[i], en1v[i], etn1[i]);
		tensorproduct3d_array(etv[i], en2v[i], etn2[i]);
		tensorproduct3d_array(en1v[i], en1v[i], en1n1[i]);
		tensorproduct3d_array(en2v[i], en2v[i], en2n2[i]);
		tensorproduct3d_array(en1v[i], en2v[i], en1n2[i]);
	}
}

/**
 * \fn void getBoundaryInfoNormalTrace(EDGE *edges, ELEMENT_DOF *elementDOF, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index)
 * \brief get information of triangulation on the boundary
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *isInNode pointer to boundary information of nodes: if the node is interior node, it will be 0; if the node is on the boundary, it will be -1
 * \param *dirichlet pointer to the index of dirichlet nodes
 * \param *nondirichlet pointer to the index of nondirichlet nodes
 * \param *index pointer to the transpose of dirichlet and nondirichlet nodes
 * \return void
 */
void getBoundaryInfoNormalTrace(EDGE *edges, ELEMENT_DOF *elementDOF, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index)
{
	int i,j,k,dr;

	int dof=elementDOF->dof;
	create_ivector(dof, isInNode);
	create_ivector(dof, index);
	
	for(j=0;j<edges->row;j++)
	{
		if(edges->val[j][3]==-1)
		{
			for(i=0;i<elementDOF->col;i++)
				isInNode->val[elementDOF->val[j][i]]=-1;
		}
	}
	
	dr=0;
	for(i=0;i<isInNode->row;i++)
	{
		if(isInNode->val[i]==-1)
			dr++;
	}
	create_ivector(dr, dirichlet);
	create_ivector(isInNode->row-dr, nondirichlet);

	j=0;k=0;
	for(i=0;i<isInNode->row;i++)
	{
		if(isInNode->val[i]==-1) //  Dirichlet boundary node
		{
			dirichlet->val[k]=i;
			index->val[i]=k;
			k++;
		}
		else // free variable
		{
			nondirichlet->val[j]=i;
			index->val[i]=j;
			j++;
		}
	}
}


/**
 * \fn void getBoundaryInfoEdge(EDGE *edges, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index)
 * \brief get information of triangulation on the boundary
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param *isInNode pointer to boundary information of nodes: if the node is interior node, it will be 0; if the node is on the boundary, it will be -1
 * \param *dirichlet pointer to the index of dirichlet nodes
 * \param *nondirichlet pointer to the index of nondirichlet nodes
 * \param *index pointer to the transpose of dirichlet and nondirichlet nodes
 * \return void
 */
void getBoundaryInfoEdge(EDGE *edges, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index)
{
	int i,j,k,dr;

	int ne=edges->row;
	create_ivector(ne, isInNode);
	create_ivector(ne, index);
	
	for(j=0;j<ne;j++)
	{
		if(edges->val[j][3]==-1)
			isInNode->val[j]=-1;
	}
	
	dr=0;
	for(i=0;i<isInNode->row;i++)
	{
		if(isInNode->val[i]==-1)
			dr++;
	}
	create_ivector(dr, dirichlet);
	create_ivector(isInNode->row-dr, nondirichlet);

	j=0;k=0;
	for(i=0;i<isInNode->row;i++)
	{
		if(isInNode->val[i]==-1) //  Dirichlet boundary node
		{
			dirichlet->val[k]=i;
			index->val[i]=k;
			k++;
		}
		else // free variable
		{
			nondirichlet->val[j]=i;
			index->val[i]=j;
			j++;
		}
	}
}
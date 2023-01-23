/*
 *  assemble.c
 *  DGFEM
 *
 *  Created by Xuehai Huang on 10/30/08.
 *  Modified by Xuehai Huang on 3/21/12.
 *  Copyright 2012 WZU. All rights reserved.
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
 * \fn int getmesh(int domain_num, ELEMENT *ptr_elements, idenmat *ptr_elementEdge, EDGE *ptr_edges, dennode *ptr_nodes, iCSRmat *edgesTran, ivector *nodeCEdge, int levelNum)
 * \brief generate mesh information
 * \param domain_num number of domain
 * \param *ptr_elements pointer to triangulation: the first 3 columns store the indexes of vertices, the last 3 columns store the indexes of edge's midpoints
 * \param *ptr_elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *ptr_edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                       the fourth column stores -1 if the edge is on boundary
 * \param *ptr_nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
                                       the third column stores the Dirichlet status of the points(0: nondirichlet, -1: dirichlet)
 * \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
 * \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
 * \param levelNum total level number of grid
 * \return 1 if succeed 0 if fail
 */
int getmesh(int domain_num, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, int levelNum)
{
	// ddenmat nodes; // the first column stores the x coordinate of points, the second column stores the y coordinate of points
	// idenmat elements[levelNum]; // triangulation: store 3 points corresponding to the element in each row
	// idenmat edges[levelNum]; // the first two columns store the two vertice, the third column stores the affiliated element
	iCSRmat elementsTran[levelNum]; // the transpose of elements. 
	// ivector isInNode; // if the node is interior node, it will be 0; if the node is on the boundary, it will be -1

	int IsExist=getCoarseInfo(domain_num, &nodes[0], &elements[0], &edges[0], &elementsTran[0], &edgesTran[0], nodeCEdge);
	if(IsExist==0)
	{
		printf("Constructing coarse grid fails!\n");
		return 0;
	}
	int i,j,k,l;
	int nvertices;

	for(l=0;l<levelNum-1;l++)
	{
		uniformrefine(&nodes[l], &elements[l], &edges[l], &elementsTran[l], &nodes[l+1], &elements[l+1], &edges[l+1], &elementsTran[l+1], &edgesTran[l+1], nodeCEdge);
	}

	for(l=0;l<levelNum;l++)
	{	
		free(elementsTran[l].IA);
		free(elementsTran[l].JA);
	}
		
	// generate elementEdge	
	int point1, point2, element1, element2;
	
	for (l = 0; l < levelNum; l++)
	{
		create_iden_matrix(elements[l].row, 3, elementEdge + l);

		for (i = 0; i<edges[l].row; i++)
		{
			point1 = edges[l].val[i][0];
			point2 = edges[l].val[i][1];
			element1 = edges[l].val[i][2];
			element2 = edges[l].val[i][3];

			for (j = 0; j<3; j++)
			{
				if (elements[l].val[element1][j] == point1)
					break;
			}
			for (k = 0; k<3; k++)
			{
				if (elements[l].val[element1][k] == point2)
					break;
			}
			elementEdge[l].val[element1][3 - j - k] = i;

			if (element2>-1)
			{
				for (j = 0; j<3; j++)
				{
					if (elements[l].val[element2][j] == point1)
						break;
				}
				for (k = 0; k<3; k++)
				{
					if (elements[l].val[element2][k] == point2)
						break;
				}
				elementEdge[l].val[element2][3 - j - k] = i;
			}
		}
	}
	// end generate ptr_elementEdge
	
	return 1;
}


/**
 * \fn void getElementDOF1D(ELEMENT_DOF *elementDOF, int ne, int dop)
 * \brief get the degrees of freedom of piecewise Lagrange element on edges
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param ne number of edges
 * \param dop degree of polynomial
 */
void getElementDOF1D(ELEMENT_DOF *elementDOF, int ne, int dop)
{
	int i,j;
	int count;
	if(dop<0)
	{
		elementDOF=NULL;
		return;
	}

	create_elementDOF(dop, ne*(dop+1), ne, dop+1, elementDOF);

	for(j=0;j<ne;j++)
	{
		count=j*elementDOF->col;
		for(i=0;i<elementDOF->col;i++)
			elementDOF->val[j][i]=count+i;
	}

}

/**
 * \fn void getElementDOF1D_Continue(ELEMENT_DOF *edgeDOF, ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
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
void getElementDOF1D_Continue(ELEMENT_DOF *edgeDOF, ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
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
 * \fn void getElementDOF(ELEMENT_DOF *elementDOF, int nt, int dop)
 * \brief get the degrees of freedom of piecewise Lagrange element
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param nt number of elements
 * \param dop degree of polynomial
 */
void getElementDOF(ELEMENT_DOF *elementDOF, int nt, int dop)
{
	int i,k;
	int count;
	if(dop<1)
		dop=0;

	create_elementDOF(dop, nt*(dop+1)*(dop+2)/2, nt, (dop+1)*(dop+2)/2, elementDOF);

	for(k=0;k<nt;k++)
	{
		count=k*elementDOF->col;
		for(i=0;i<elementDOF->col;i++)
			elementDOF->val[k][i]=count+i;
	}
}

/**
 * \fn void getElementDOF_Lagrange2d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, int nedges, int nvertices, int dop)
 * \brief get the degrees of freedom of Lagrange element
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param nvertices number of vertices
 * \param dop degree of polynomial
 */
void getElementDOF_Lagrange2d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
{
	int i,j,k,ii[3];
	int nt=elements->row;
	int ne=edges->row;
	int nn=nvertices;
	if(dop<2)
		dop=1;

	create_elementDOF(dop, nn + ne*(dop-1) + nt*(dop-1)*(dop-2)/2, nt, (dop+1)*(dop+2)/2, elementDOF);

	int node, edge;
	int *perm;
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
			perm = elements->eperm[k][j];

			for(ii[0]=0;ii[0]<dop-1;ii[0]++)
			{
				ii[1]=dop-2-ii[0];
				elementDOF->val[k][3+(dop-1)*j+ ii[0]] = nn + edge*(dop-1) + ii[perm[0]];
			}
		}

		for(i=0;i<(dop-1)*(dop-2)/2;i++)
		{
			elementDOF->val[k][3*dop+i] = nn + ne*(dop-1) + k*(dop-1)*(dop-2)/2 + i;
		}
	}
}

/**
 * \fn void getElementDOF_CrouzeixRaviart2d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges)
 * \brief get the degrees of freedom of Crouzeix-Raviart element
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements the fourth column stores -1 if the edge is on boundary
 */
void getElementDOF_CrouzeixRaviart2d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges)
{
	int i,j,k;
	int nt=elements->row;
	int ne=edges->row;
	// int nn=nvertices;
	
	create_elementDOF(1, ne, nt, 3, elementDOF);

	int edge;
	for(k=0;k<nt;k++){	
		for(j=0;j<3;j++){
			edge=elementEdge->val[k][j];
			elementDOF->val[k][j] = edge;
		}
	}
}

/**
* \fn void getElementDOF_MINI2d(ELEMENT_DOF *elementDOF, ELEMENT *elements, int nvertices)
* \brief get the degrees of freedom of MINI element for Stokes equation in two dimensions
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
* \param nvertices number of vertices
*/
void getElementDOF_MINI2d(ELEMENT_DOF *elementDOF, ELEMENT *elements, int nvertices)
{
	int i, j, k;
	int nt = elements->row;
	int nn = nvertices;

	create_elementDOF(3, nn + nt, nt, 4, elementDOF);

	int node, edge;
	for (k = 0; k<nt; k++){
		for (i = 0; i<3; i++){
			node = elements->val[k][i];
			elementDOF->val[k][i] = node;
		}
		elementDOF->val[k][3] = nn + k;
	}
}

/**
 * \fn void getElementDOF_Morley2d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices)
 * \brief get the degrees of freedom of Morley element
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param nvertices number of vertices
 */
void getElementDOF_Morley2d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices)
{
	int i,j,k,ii[3];
	int nt=elements->row;
	int ne=edges->row;
	int nn=nvertices;
	
	create_elementDOF(2, nn + ne, nt, 6, elementDOF);

	int node, edge;
	int *perm;
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
			elementDOF->val[k][3+j] = nn + edge;
		}
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
 * \fn void getElementDOF_HuangZhou(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices)
 * \brief get the degrees of freedom of Huang-Zhou element
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param nvertices number of vertices
 * \param dop degree of polynomial
 */
void getElementDOF_HuangZhou(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices)
{
	int i,j,k;
	int nt=elements->row;
	int ne=edges->row;
	int nn=nvertices;

	create_elementDOF(3, nn*3 + ne*3 + nt*3, nt, 21, elementDOF);

	int node, edge;
	int orient;
	for(k=0;k<nt;k++){
		for(i=0;i<3;i++){
			node=elements->val[k][i];
			elementDOF->val[k][i]=node;
			elementDOF->val[k][i+3]=node+nn;
			elementDOF->val[k][i+6]=node+nn*2;
		}
		
		for(j=0;j<3;j++){
			edge=elementEdge->val[k][j];
			if(elements->val[k][(j+1)%3]== edges->val[edge][0])
				orient=1;
			else
				orient=0;
			
			if(orient==1){
				for(i=0;i<2;i++)
					elementDOF->val[k][9+2*j+i] = nn*3 + edge*2 + i;
			}
			else{
				for(i=0;i<2;i++)
					elementDOF->val[k][9+2*j+1-i] = nn*3 + edge*2 + i;
			}
			
			elementDOF->val[k][15+j] = nn*3 + ne*2 + edge;
		} // j

		for(i=0;i<3;i++)
			elementDOF->val[k][18+i] = nn*3 + ne*3 + k*3 + i;
	} // k
}

/**
* \fn void getFreenodesInfoLagrange2d(EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of Lagrange element in two dimensions
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoLagrange2d(EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i, j, k, estride, fstride, nnf;

	int nn = nodes->row;
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
* \fn void getFreenodesInfoCrouzeixRaviart2d(EDGE *edges, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of Crouzeix-Raviart element in two dimensions
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoCrouzeixRaviart2d(EDGE *edges, ELEMENT_DOF *elementDOF)
{
	int i, j, k, estride, fstride, nnf;

	// int nn = nodes->row;
	int ne = edges->row;
	int dof = elementDOF->dof;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(dof, nfFlag);
	create_ivector(dof, index);

	nnf = 0; // number of non-free nodes
	for (j = 0; j<ne; j++){
		if (edges->bdFlag[j] == 1 || edges->bdFlag[j] == 2 || edges->bdFlag[j] == 3 || edges->bdFlag[j] == 4) // Dirichlet boundary
		{
			nfFlag->val[j] = 1;
			nnf++;
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
* \fn void getFreenodesInfoMINI2d(dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of MINI element in two dimensions
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoMINI2d(dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i, j, k, estride, fstride, nnf;

	int nn = nodes->row;
	// int ne = edges->row;
	int dof = elementDOF->dof;

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
			nnf++;
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
* \fn void getFreenodesInfoMorley2d(EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of Morley element in two dimensions
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoMorley2d(EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i, j, k, estride, fstride, nnf;

	int nn = nodes->row;
	int ne = edges->row;
	int dof = elementDOF->dof;
	int dop = 2;

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
			nnf++;
		}
	}

	for (j = 0; j<ne; j++)
	{
		if (edges->bdFlag[j] == 1 || edges->bdFlag[j] == 2 || edges->bdFlag[j] == 3 || edges->bdFlag[j] == 4) // Dirichlet boundary
		{
			nfFlag->val[nn + j] = 1;
			nnf++;
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
* \fn void getFreenodesInfoDG(ELEMENT_DOF *elementDOF)
* \brief get freenodes information of no non-free nodes
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoDG(ELEMENT_DOF *elementDOF)
{
	int i;

	int dof = elementDOF->dof;
	// int dop = elementDOF->dop;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(dof, nfFlag);
	create_ivector(dof, index);
	create_ivector(0, nfreenodes);
	create_ivector(dof, freenodes);

	for (i = 0; i<dof; i++){
		freenodes->val[i] = i;
		index->val[i] = i;
	}
}

/**
* \fn void assembleMassmatrixLagrange2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double mu)
* \brief assemble mass matrix
* \param *A pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assembleMassmatrixLagrange2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double mu)
{
	int i, j, k, l;

	//	create_dvector(A->row, b);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element, edge, node;

	double phi, phi1[2], phi2[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double val, s, **gradLambda;
	int count;

	int num_qp;
	double lambdas[100][3], weight[100];

	
	/***************************** stiffness matrix A ***************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF->col*elementDOF->col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(elementDOF->col, elementDOF->col, &lA);

	num_qp=getNumQuadPoints(elementDOF->dop * 2, 2); // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		// gradLambda = elements->gradLambda[k];
		// end set parameters

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			for (k2 = 0; k2<elementDOF->col; k2++)
			{
				val = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					lagrange_basis(lambdas[i1], k1, elementDOF->dop, phi1);
					lagrange_basis(lambdas[i1], k2, elementDOF->dop, phi2);
					val += s*weight[i1] * phi1[0] * phi2[0] * 2 * mu;
				}
				lA.val[k1][k2] += val;
			} // k2
		} // k1

		l = elementDOF->col*elementDOF->col * k;
		for (i = 0; i<elementDOF->col; i++)
		{
			for (j = 0; j<elementDOF->col; j++)
			{
				ia[l] = elementDOF->val[k][i];
				ja[l] = elementDOF->val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // i
		} // j
	} // k
	free_dden_matrix(&lA);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	double eps = 0;
	if(eps<1e-20){
		dIJtoCSR(A, ia, ja, va, N, 0, 0);
	free(ia); free(ja); free(va);
	}
	else{
		int nzmax = 0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ) nzmax++;
		}

		int *Ia, *Ja;
		double *Va;
		Ia = (int*)malloc(nzmax * sizeof(int));
		Ja = (int*)malloc(nzmax * sizeof(int));
		Va = (double*)malloc(nzmax * sizeof(double));
		int cur=0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ){
				Ia[cur] = ia[i];
				Ja[cur] = ja[i];
				Va[cur] = va[i];
				cur++;
			}
		}
		free(ia); free(ja); free(va);
		dIJtoCSR(A, Ia, Ja, Va, nzmax, 0, 0);
		free(Ia); free(Ja); free(Va);
	}
}

/**
* \fn void assembleMassmatrixLagrangeDiagBlockRepeat2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double mu)
* \brief assemble mass matrix
* \param *A pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assembleMassmatrixLagrangeDiagBlockRepeat2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double mu)
{
	int i, j, k, l;

	//	create_dvector(A->row, b);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element, edge, node;

	double phi, phi1[2], phi2[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double val, s, **gradLambda;
	int count;

	int num_qp;
	double lambdas[100][3], weight[100];

	
	/***************************** stiffness matrix A ***************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF->col*elementDOF->col*elements->row;
	ia = (int*)malloc(2*N * sizeof(int));
	ja = (int*)malloc(2*N * sizeof(int));
	va = (double*)malloc(2*N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(elementDOF->col, elementDOF->col, &lA);

	num_qp=getNumQuadPoints(elementDOF->dop * 2, 2); // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		// gradLambda = elements->gradLambda[k];
		// end set parameters

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			for (k2 = 0; k2<elementDOF->col; k2++)
			{
				val = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					lagrange_basis(lambdas[i1], k1, elementDOF->dop, phi1);
					lagrange_basis(lambdas[i1], k2, elementDOF->dop, phi2);
					val += s*weight[i1] * phi1[0] * phi2[0] * 2 * mu;
				}
				lA.val[k1][k2] += val;
			} // k2
		} // k1

		l = elementDOF->col*elementDOF->col * k;
		for (i = 0; i<elementDOF->col; i++)
		{
			for (j = 0; j<elementDOF->col; j++)
			{
				ia[l] = elementDOF->val[k][i];
				ja[l] = elementDOF->val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // i
		} // j
	} // k
	free_dden_matrix(&lA);

	for(i=0;i<N;i++){
		ia[i+N] = ia[i] + elementDOF->dof;
		ja[i+N] = ja[i] + elementDOF->dof;
		va[i+N] = va[i];
	}
	N *= 2;
// remove zero elements and transform matrix A from its IJ format to its CSR format
	double eps = 0;
	if(eps<1e-20){
		dIJtoCSR(A, ia, ja, va, N, 0, 0);
	free(ia); free(ja); free(va);
	}
	else{
		int nzmax = 0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ) nzmax++;
		}

		int *Ia, *Ja;
		double *Va;
		Ia = (int*)malloc(nzmax * sizeof(int));
		Ja = (int*)malloc(nzmax * sizeof(int));
		Va = (double*)malloc(nzmax * sizeof(double));
		int cur=0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ){
				Ia[cur] = ia[i];
				Ja[cur] = ja[i];
				Va[cur] = va[i];
				cur++;
			}
		}
		free(ia); free(ja); free(va);
		dIJtoCSR(A, Ia, Ja, Va, nzmax, 0, 0);
		free(Ia); free(Ja); free(Va);
	}
}

/**
* \fn void assembleBiGradLagrange2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double mu)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assembleBiGradLagrange2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double mu)
{
	int i, j, k, l;

	//	create_dvector(A->row, b);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element, edge, node;

	double phi, phi1[2], phi2[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double val, s, **gradLambda;
	int count;

	int num_qp;
	double lambdas[100][3], weight[100];

	
	/***************************** stiffness matrix A ***************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF->col*elementDOF->col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(elementDOF->col, elementDOF->col, &lA);

	num_qp=getNumQuadPoints(elementDOF->dop * 2 - 2, 2); // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		// end set parameters

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			for (k2 = 0; k2<elementDOF->col; k2++)
			{
				val = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					lagrange_basis1(lambdas[i1], gradLambda, k1, elementDOF->dop, phi1);
					lagrange_basis1(lambdas[i1], gradLambda, k2, elementDOF->dop, phi2);
					val += s*weight[i1] * (phi1[0] * phi2[0] + phi1[1] * phi2[1]) * 2 * mu;
				}
				lA.val[k1][k2] += val;
			} // k2
		} // k1

		l = elementDOF->col*elementDOF->col * k;
		for (i = 0; i<elementDOF->col; i++)
		{
			for (j = 0; j<elementDOF->col; j++)
			{
				ia[l] = elementDOF->val[k][i];
				ja[l] = elementDOF->val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // i
		} // j
	} // k
	free_dden_matrix(&lA);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	double eps = 0;
	if(eps<1e-20){
		dIJtoCSR(A, ia, ja, va, N, 0, 0);
	free(ia); free(ja); free(va);
	}
	else{
		int nzmax = 0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ) nzmax++;
		}

		int *Ia, *Ja;
		double *Va;
		Ia = (int*)malloc(nzmax * sizeof(int));
		Ja = (int*)malloc(nzmax * sizeof(int));
		Va = (double*)malloc(nzmax * sizeof(double));
		int cur=0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ){
				Ia[cur] = ia[i];
				Ja[cur] = ja[i];
				Va[cur] = va[i];
				cur++;
			}
		}
		free(ia); free(ja); free(va);
		dIJtoCSR(A, Ia, Ja, Va, nzmax, 0, 0);
		free(Ia); free(Ja); free(Va);
	}
}

/**
 * \fn void assembleRHSLagrange2d(dvector *b, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double (*f)(double *))
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
 * \param *elementdofTran pointer to transpose of elementDOF
 * \return void
 */
void assembleRHSLagrange2d(dvector *b, ELEMENT *elements, ELEMENT_DOF *elementDOF, double (*f)(double *, double *), double *paras)
{
	int i,j,k,k1,i1;
	
	double phi;
	double x[2], **gradLambda, **vertices;
	double s;

	int num_qp;
	double lambdas[100][3], weight[100];
			
	dvector lb;
	create_dvector(elementDOF->col, &lb);
	/************************************************** right hand side b *****************************************************************/
	create_dvector(elementDOF->dof, b);
	num_qp = 49; 
	init_Gauss2d(num_qp, lambdas, weight);

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		vertices = elements->vertices[k];
		// end set parameters

		init_dvector(&lb, 0.0);
		for (i = 0; i<elementDOF->col; i++)
		{
			for (i1 = 0; i1<num_qp; i1++)
			{
				lagrange_basis(lambdas[i1], i, elementDOF->dop, &phi);
				baryToCart2d(lambdas[i1], x, vertices);
				lb.val[i] += s*weight[i1] * f(x, paras)*phi;
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
* \fn void assembleBiGradCrouzeixRaviart2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double mu)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assembleBiGradCrouzeixRaviart2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double mu)
{
	int i, j, k, l;

	//	create_dvector(A->row, b);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element, edge, node;

	double phi, phi1[2], phi2[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double val, s, **gradLambda;
	int count;

	int num_qp;
	double lambdas[100][3], weight[100];

	
	/***************************** stiffness matrix A ***************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF->col*elementDOF->col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(elementDOF->col, elementDOF->col, &lA);

	num_qp=getNumQuadPoints(elementDOF->dop * 2 - 2, 2); // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		// end set parameters

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			for (k2 = 0; k2<elementDOF->col; k2++)
			{
				val = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					cr_basis1(2, gradLambda, k1, phi1);
					cr_basis1(2, gradLambda, k2, phi2);
					val += s*weight[i1] * (phi1[0] * phi2[0] + phi1[1] * phi2[1]) * 2 * mu;
				}
				lA.val[k1][k2] += val;
			} // k2
		} // k1

		l = elementDOF->col*elementDOF->col * k;
		for (i = 0; i<elementDOF->col; i++)
		{
			for (j = 0; j<elementDOF->col; j++)
			{
				ia[l] = elementDOF->val[k][i];
				ja[l] = elementDOF->val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // i
		} // j
	} // k
	free_dden_matrix(&lA);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	double eps = 0;
	if(eps<1e-20){
		dIJtoCSR(A, ia, ja, va, N, 0, 0);
	free(ia); free(ja); free(va);
	}
	else{
		int nzmax = 0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ) nzmax++;
		}

		int *Ia, *Ja;
		double *Va;
		Ia = (int*)malloc(nzmax * sizeof(int));
		Ja = (int*)malloc(nzmax * sizeof(int));
		Va = (double*)malloc(nzmax * sizeof(double));
		int cur=0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ){
				Ia[cur] = ia[i];
				Ja[cur] = ja[i];
				Va[cur] = va[i];
				cur++;
			}
		}
		free(ia); free(ja); free(va);
		dIJtoCSR(A, Ia, Ja, Va, nzmax, 0, 0);
		free(Ia); free(Ja); free(Va);
	}
}

/**
 * \fn void assembleRHSCrouzeixRaviart2d(dvector *b, ELEMENT *elements, ELEMENT_DOF *elementDOF, double (*f)(double *, double *), double *paras)
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
 * \param *elementdofTran pointer to transpose of elementDOF
 * \return void
 */
void assembleRHSCrouzeixRaviart2d(dvector *b, ELEMENT *elements, ELEMENT_DOF *elementDOF, double (*f)(double *, double *), double *paras)
{
	int i,j,k,k1,i1;
	
	double phi;
	double x[2], **gradLambda, **vertices;
	double s;

	int num_qp;
	double lambdas[100][3], weight[100];
			
	dvector lb;
	create_dvector(elementDOF->col, &lb);
	/************************************************** right hand side b *****************************************************************/
	create_dvector(elementDOF->dof, b);
	num_qp = 49; 
	init_Gauss2d(num_qp, lambdas, weight);

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		vertices = elements->vertices[k];
		// end set parameters

		init_dvector(&lb, 0.0);
		for (i = 0; i<elementDOF->col; i++)
		{
			for (i1 = 0; i1<num_qp; i1++)
			{
				cr_basis(2, lambdas[i1], i, &phi);
				baryToCart2d(lambdas[i1], x, vertices);
				lb.val[i] += s*weight[i1] * f(x, paras)*phi;
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
* \fn void assembleBiGradMorley2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \return void
*/
void assembleBiGradMorley2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i, j, k, l;

	//	create_dvector(A->row, b);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element, edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, j1, l1, l2, ej;
	double val, s, **gradLambda, *nve[3];
	int count;

	int num_qp;
	double lambdas[100][3], weight[100];

	
	/***************************** stiffness matrix A ***************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF->col*elementDOF->col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(elementDOF->col, elementDOF->col, &lA);

	num_qp=getNumQuadPoints(elementDOF->dop * 2 - 2, 2); // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		for(i=0;i<3;i++){
			j = elementEdge->val[k][i];
			nve[i] = edges->nvector[j];
		}
		// end set parameters

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF->col; k1++){
			for (k2 = 0; k2<elementDOF->col; k2++){
				val = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					morley_basis1(lambdas[i1], gradLambda, nve, k1, phi1);
					morley_basis1(lambdas[i1], gradLambda, nve, k2, phi2);
					val += s*weight[i1] * (phi1[0] * phi2[0] + phi1[1] * phi2[1]);
				}
				lA.val[k1][k2] += val;
			} // k2
		} // k1

		l = elementDOF->col*elementDOF->col * k;
		for (i = 0; i<elementDOF->col; i++)
		{
			for (j = 0; j<elementDOF->col; j++)
			{
				ia[l] = elementDOF->val[k][i];
				ja[l] = elementDOF->val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // i
		} // j
	} // k
	free_dden_matrix(&lA);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	double eps = 0;
	if(eps<1e-20){
		dIJtoCSR(A, ia, ja, va, N, 0, 0);
	free(ia); free(ja); free(va);
	}
	else{
		int nzmax = 0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ) nzmax++;
		}

		int *Ia, *Ja;
		double *Va;
		Ia = (int*)malloc(nzmax * sizeof(int));
		Ja = (int*)malloc(nzmax * sizeof(int));
		Va = (double*)malloc(nzmax * sizeof(double));
		int cur=0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ){
				Ia[cur] = ia[i];
				Ja[cur] = ja[i];
				Va[cur] = va[i];
				cur++;
			}
		}
		free(ia); free(ja); free(va);
		dIJtoCSR(A, Ia, Ja, Va, nzmax, 0, 0);
		free(Ia); free(Ja); free(Va);
	}
}

/**
* \fn void assembleBiHessMorley2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \return void
*/
void assembleBiHessMorley2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i, j, k, l;

	//	create_dvector(A->row, b);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element, edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, j1, l1, l2, ej;
	double val, s, **gradLambda, *nve[3];
	int count;

	int num_qp;
	double lambdas[100][3], weight[100];

	
	/***************************** stiffness matrix A ***************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF->col*elementDOF->col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(elementDOF->col, elementDOF->col, &lA);

	num_qp=getNumQuadPoints(elementDOF->dop * 2 - 4, 2); // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		for(i=0;i<3;i++){
			j = elementEdge->val[k][i];
			nve[i] = edges->nvector[j];
		}
		// end set parameters

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			for (k2 = 0; k2<elementDOF->col; k2++)
			{
				val = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					morley_basis2(gradLambda, nve, k1, phi1);
					morley_basis2(gradLambda, nve, k2, phi2);
					val += s*weight[i1] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]);
				}
				lA.val[k1][k2] += val;
			} // k2
		} // k1

		l = elementDOF->col*elementDOF->col * k;
		for (i = 0; i<elementDOF->col; i++)
		{
			for (j = 0; j<elementDOF->col; j++)
			{
				ia[l] = elementDOF->val[k][i];
				ja[l] = elementDOF->val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // i
		} // j
	} // k
	free_dden_matrix(&lA);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	double eps = 0;
	if(eps<1e-20){
		dIJtoCSR(A, ia, ja, va, N, 0, 0);
	free(ia); free(ja); free(va);
	}
	else{
		int nzmax = 0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ) nzmax++;
		}

		int *Ia, *Ja;
		double *Va;
		Ia = (int*)malloc(nzmax * sizeof(int));
		Ja = (int*)malloc(nzmax * sizeof(int));
		Va = (double*)malloc(nzmax * sizeof(double));
		int cur=0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ){
				Ia[cur] = ia[i];
				Ja[cur] = ja[i];
				Va[cur] = va[i];
				cur++;
			}
		}
		free(ia); free(ja); free(va);
		dIJtoCSR(A, Ia, Ja, Va, nzmax, 0, 0);
		free(Ia); free(Ja); free(Va);
	}
}

/**
 * \fn void assembleRHSMorley2d(dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, double (*f)(double *, double *), double *paras)
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
 * \param *elementdofTran pointer to transpose of elementDOF
 * \return void
 */
void assembleRHSMorley2d(dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, double (*f)(double *, double *), double *paras)
{
	int i,j,k,k1,i1;
	
	double phi;
	double x[2], **gradLambda, *nve[3], **vertices;
	double s;

	int num_qp;
	double lambdas[100][3], weight[100];
			
	dvector lb;
	create_dvector(elementDOF->col, &lb);
	/************************************************** right hand side b *****************************************************************/
	create_dvector(elementDOF->dof, b);
	num_qp = 49; 
	init_Gauss2d(num_qp, lambdas, weight);

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vertices = elements->vertices[k];
		s = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		for(i=0;i<3;i++){
			j = elementEdge->val[k][i];
			nve[i] = edges->nvector[j];
		}
		// end set parameters
        
		init_dvector(&lb, 0.0);
		for (i = 0; i<elementDOF->col; i++)
		{
			for (i1 = 0; i1<num_qp; i1++)
			{
				morley_basis(lambdas[i1], gradLambda, nve, i, &phi);
				baryToCart2d(lambdas[i1], x, vertices);
				lb.val[i] += s*weight[i1] * f(x, paras)*phi;
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
* \fn void assembleBiHessC0ipdg2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double parapenalty)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \return void
*/
void assembleBiHessC0ipdg2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double parapenalty)
{
	int i, j, k, l, m;

	//	create_dvector(A->row, b);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[4], edge, node;

	int dop = elementDOF->dop;

	double phi, phi1[3], phi2[3], jump[2], avg[2];
	int k1, k2, i1, j1, l1, l2, ei;
	double val, s, **gradLambda, *nve[3];
	int count;

	int num_qp;
	double lambdas[100][3], weight[100];

	
	/***************************** stiffness matrix A ***************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF->col*elementDOF->col*elements->row;
	int patchnodes[200];
	int *index;
	int istart;
	index = (int*)calloc(elementDOF->dof, sizeof(int));
	for (i = 0; i<elementDOF->dof; i++)
		index[i] = -1;
	for(edge=0;edge<edges->row;edge++){
		element[0] = edges->val[edge][2];
		element[1] = edges->val[edge][3];
		// elen = edges->length[edge];

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
		N += count*count;
	}

	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));

	ddenmat lA; // local A
	create_dden_matrix(elementDOF->col, elementDOF->col, &lA);
	// step 3A1: Loop element by element and compute the actual entries storing them in A
	num_qp=getNumQuadPoints(elementDOF->dop * 2 - 4, 2); // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		// for(i=0;i<3;i++){
		// 	j = elementEdge->val[k][i];
		// 	nve[i] = edges->nvector[j];
		// }
		// end set parameters

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			for (k2 = 0; k2<elementDOF->col; k2++)
			{
				val = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					lagrange_basis2(lambdas[i1], gradLambda, k1, dop, phi1);
					lagrange_basis2(lambdas[i1], gradLambda, k2, dop, phi2);
					val += s*weight[i1] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]);
				}
				lA.val[k1][k2] += val;
			} // k2
		} // k1

		l = elementDOF->col*elementDOF->col * k;
		for (i = 0; i<elementDOF->col; i++)
		{
			for (j = 0; j<elementDOF->col; j++)
			{
				ia[l] = elementDOF->val[k][i];
				ja[l] = elementDOF->val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // i
		} // j
	} // k
	free_dden_matrix(&lA);

	
	// additonal terms in IPDG
	// search edge by edge
	double elen;
	int num_qp1;
	double lambdas1[100][2], weight1[100];
	num_qp1=getNumQuadPoints(elementDOF->dop * 2-2, 1);
	if (num_qp1>5) num_qp1 = 5;
	init_Gauss1d(num_qp1, lambdas1, weight1);

	for(edge=0;edge<edges->row;edge++){
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

		create_dden_matrix(count, count, &lA);
		init_dden_matrix(&lA, 0.0);

		for (k1 = 0; k1<count; k1++)
		{
			i = patchnodes[k1];
			for (k2 = 0; k2<count; k2++)
			{
				j = patchnodes[k2];
				for (i1 = 0; i1<num_qp1; i1++)
				{
					averageNormalDerivative2Lagrange2d(lambdas1[i1], edge, elements, elementEdge, edges, elementDOF, i, avg);
					averageNormalDerivative2Lagrange2d(lambdas1[i1], edge, elements, elementEdge, edges, elementDOF, j, avg+1);

					jumpNormalDerivativeLagrange2d(lambdas1[i1], edge, elements, elementEdge, edges, elementDOF, i, jump);
					jumpNormalDerivativeLagrange2d(lambdas1[i1], edge, elements, elementEdge, edges, elementDOF, j, jump+1);

					lA.val[k1][k2] -= elen * weight1[i1] * (avg[0] * jump[1] + avg[1] * jump[0]);
					lA.val[k1][k2] += parapenalty * weight1[i1] * jump[0] * jump[1];

					// averageNormalDerivative2Lagrange2d(lambdas1[i1], edge, elements, elementEdge, edges, elementDOF, i, phi1);
					// jumpNormalDerivativeLagrange2d(lambdas1[i1], edge, elements, elementEdge, edges, elementDOF, j, phi2);
					// lA.val[k1][k2] -= elen * weight1[i1] * phi1[0] * phi2[0];

					// averageNormalDerivative2Lagrange2d(lambdas1[i1], edge, elements, elementEdge, edges, elementDOF, j, phi1);
					// jumpNormalDerivativeLagrange2d(lambdas1[i1], edge, elements, elementEdge, edges, elementDOF, i, phi2);
					// lA.val[k1][k2] -= elen * weight1[i1] * phi1[0] * phi2[0];

					// jumpNormalDerivativeLagrange2d(lambdas1[i1], edge, elements, elementEdge, edges, elementDOF, i, phi1);
					// jumpNormalDerivativeLagrange2d(lambdas1[i1], edge, elements, elementEdge, edges, elementDOF, j, phi2);
					// lA.val[k1][k2] += parapenalty * weight1[i1] * phi1[0] * phi2[0];
				}
			} // k2
		} // k1

		for (i = 0; i<count; i++)
		{
			for (j = 0; j<count; j++)
			{
				ia[l] = patchnodes[i];
				ja[l] = patchnodes[j];
				va[l] = lA.val[i][j];
				l++;
			} // i
		} // j

		free_dden_matrix(&lA);
	}
	free(index);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	double eps = 0;
	if(eps<1e-20){
		dIJtoCSR(A, ia, ja, va, N, 0, 0);
	free(ia); free(ja); free(va);
	}
	else{
		int nzmax = 0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ) nzmax++;
		}

		int *Ia, *Ja;
		double *Va;
		Ia = (int*)malloc(nzmax * sizeof(int));
		Ja = (int*)malloc(nzmax * sizeof(int));
		Va = (double*)malloc(nzmax * sizeof(double));
		int cur=0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ){
				Ia[cur] = ia[i];
				Ja[cur] = ja[i];
				Va[cur] = va[i];
				cur++;
			}
		}
		free(ia); free(ja); free(va);
		dIJtoCSR(A, Ia, Ja, Va, nzmax, 0, 0);
		free(Ia); free(Ja); free(Va);
	}
}

/**
* \fn void assembleBiGradCrouzeixRaviartDiagBlockRepeat2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double mu)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assembleBiGradCrouzeixRaviartDiagBlockRepeat2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double mu)
{
	int i, j, k, l;

	//	create_dvector(A->row, b);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element, edge, node;

	double phi, phi1[2], phi2[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double val, s, **gradLambda;
	int count;

	int num_qp;
	double lambdas[100][3], weight[100];

	
	/***************************** stiffness matrix A ***************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF->col*elementDOF->col*elements->row;
	ia = (int*)malloc(2*N * sizeof(int));
	ja = (int*)malloc(2*N * sizeof(int));
	va = (double*)malloc(2*N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(elementDOF->col, elementDOF->col, &lA);

	num_qp=getNumQuadPoints(elementDOF->dop * 2 - 2, 2); // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		// end set parameters

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			for (k2 = 0; k2<elementDOF->col; k2++)
			{
				val = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					cr_basis1(2, gradLambda, k1, phi1);
					cr_basis1(2, gradLambda, k2, phi2);
					val += s*weight[i1] * (phi1[0] * phi2[0] + phi1[1] * phi2[1]) * 2 * mu;
				}
				lA.val[k1][k2] += val;
			} // k2
		} // k1

		l = elementDOF->col*elementDOF->col * k;
		for (i = 0; i<elementDOF->col; i++)
		{
			for (j = 0; j<elementDOF->col; j++)
			{
				ia[l] = elementDOF->val[k][i];
				ja[l] = elementDOF->val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // i
		} // j
	} // k
	free_dden_matrix(&lA);

	for(i=0;i<N;i++){
		ia[i+N] = ia[i] + elementDOF->dof;
		ja[i+N] = ja[i] + elementDOF->dof;
		va[i+N] = va[i];
	}
	N *= 2;
	// remove zero elements and transform matrix A from its IJ format to its CSR format
	double eps = 0;
	if(eps<1e-20){
		dIJtoCSR(A, ia, ja, va, N, 0, 0);
	free(ia); free(ja); free(va);
	}
	else{
		int nzmax = 0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ) nzmax++;
		}

		int *Ia, *Ja;
		double *Va;
		Ia = (int*)malloc(nzmax * sizeof(int));
		Ja = (int*)malloc(nzmax * sizeof(int));
		Va = (double*)malloc(nzmax * sizeof(double));
		int cur=0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ){
				Ia[cur] = ia[i];
				Ja[cur] = ja[i];
				Va[cur] = va[i];
				cur++;
			}
		}
		free(ia); free(ja); free(va);
		dIJtoCSR(A, Ia, Ja, Va, nzmax, 0, 0);
		free(Ia); free(Ja); free(Va);
	}
}

/**
* \fn void assembleDivCrouzeixRaviartL2poly2d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *BT pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \return void
*/
void assembleDivCrouzeixRaviartL2poly2d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1)
{
	int i, j, k, l;

	int element[2], edge, node;

	double phi, phi1[3], phi2[3], val[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, xs[3], ys[3], s, **gradLambda;

	int num_qp;
	double lambdas[100][3], weight[100], gauss[100][3];

	/***************** matrix A *****************/
	int *ia, *ja;
	double *va;
	int N = 2*elementDOF0->col*elementDOF1->col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA[2]; // local A
	create_dden_matrix(elementDOF0->col, elementDOF1->col, lA);
	create_dden_matrix(elementDOF0->col, elementDOF1->col, lA+1);
	
	// num_qp = getNumQuadPoints_ShunnWilliams(elementDOF->dop * 2, 2); // the number of numerical intergation points
	// init_ShunnWilliams2d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	num_qp = getNumQuadPoints(elementDOF0->dop + elementDOF1->dop - 1, 2); // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		// end set parameters

		init_dden_matrix(&lA[0], 0.0);init_dden_matrix(&lA[1], 0.0);
		for (k1 = 0; k1<elementDOF0->col; k1++){
			for (k2 = 0; k2<elementDOF1->col; k2++){
				val[0] = 0; val[1] = 0;
				for (i1 = 0; i1<num_qp; i1++){
					// huzhang_basisDIV(lambdas[i1], gradLambda, elements->nvector[k], elements->tvector[k], k1, elementDOF[0].dop, phi1);
					cr_basis1(2, gradLambda, k1, phi1);
					val[0] += s*weight[i1] * phi1[0];
					val[1] += s*weight[i1] * phi1[1];
				}
				lA[0].val[k1][k2] += val[0]; lA[1].val[k1][k2] += val[1];
			} // k2
		} // k1

		l = 2*elementDOF0->col*elementDOF1->col * k;
		for (i = 0; i<elementDOF0->col; i++)
		{
			for (j = 0; j<elementDOF1->col; j++)
			{
				ia[l] = elementDOF0->val[k][i];
				ja[l] = elementDOF1->val[k][j];
				va[l] = lA[0].val[i][j];
				l++;
				ia[l] = elementDOF0->val[k][i] + elementDOF0->dof;
				ja[l] = elementDOF1->val[k][j];
				va[l] = lA[1].val[i][j];
				l++;
			} // i
		} // j
	} // k
	free_dden_matrix(&lA[0]);
	free_dden_matrix(&lA[1]);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	double eps = 0;
	if(eps<1e-20){
		dIJtoCSR(A, ia, ja, va, N, 0, 0);
		free(ia); free(ja); free(va);
	}
	else{
		int nzmax = 0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ) nzmax++;
		}

		int *Ia, *Ja;
		double *Va;
		Ia = (int*)malloc(nzmax * sizeof(int));
		Ja = (int*)malloc(nzmax * sizeof(int));
		Va = (double*)malloc(nzmax * sizeof(double));
		int cur=0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ){
				Ia[cur] = ia[i];
				Ja[cur] = ja[i];
				Va[cur] = va[i];
				cur++;
			}
		}
		free(ia); free(ja); free(va);
		dIJtoCSR(A, Ia, Ja, Va, nzmax, 0, 0);
		free(Ia); free(Ja); free(Va);
	}
}

/**
 * \fn void assembleRHScurlMorleyCrouzeixRaviart2d(dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1, dvector *uh)
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
 * \param *elementdofTran pointer to transpose of elementDOF
 * \return void
 */
void assembleRHScurlMorleyCrouzeixRaviart2d(dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1, dvector *uh)
{
	int i,j,k,i1,j1,k1;
	
// write_dvector4Matlab(uh, "output/uh.dat");//////////
// printf("elementDOF0:\n");
// for(i=0;i<elementDOF0->row;i++)
// {
// for(j=0;j<elementDOF0->col;j++)
// printf("%d, ", elementDOF0->val[i][j]);
// printf("\n");
// }
// printf("elementDOF1:\n");
// for(i=0;i<elementDOF1->row;i++)
// {
// for(j=0;j<elementDOF1->col;j++)
// printf("%d, ", elementDOF1->val[i][j]);
// printf("\n");
// }

	double phi, phi1[2], val[2];
	double x[2], **gradLambda, *nve[3], **vertices;
	double s;

	int num_qp;
	double lambdas[100][3], weight[100];
			
	dvector lb;
	create_dvector(2*elementDOF1->col, &lb);
	/************************************************** right hand side b *****************************************************************/
	create_dvector(2*elementDOF1->dof, b);
	num_qp = getNumQuadPoints(elementDOF0->dop+elementDOF1->dop-1, 2); // elementDOF0->dop+elementDOF1->dop-1
	init_Gauss2d(num_qp, lambdas, weight);

// printf("dsdfsfsdfdsfds\n");////////////////////////

	for (k = 0; k<elements->row; k++){
		// set parameters
		vertices = elements->vertices[k];
		s = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		for(i=0;i<3;i++){
			j = elementEdge->val[k][i];
			nve[i] = edges->nvector[j];
		}
		// end set parameters
        
		init_dvector(&lb, 0.0);
		for (i = 0; i<elementDOF1->col; i++){
			for (i1 = 0; i1<num_qp; i1++){
				cr_basis(2, lambdas[i1], i, &phi);
				val[0]=0; val[1]=0;
				for (j = 0; j<elementDOF0->col; j++){
					morley_basis1(lambdas[i1], gradLambda, nve, j, phi1);
					j1 = elementDOF0->val[k][j];
					val[0] += phi1[1] * uh->val[j1];
					val[1] -= phi1[0] * uh->val[j1];
				}
				lb.val[i] += s*weight[i1]*val[0]*phi;
				lb.val[i+elementDOF1->col] += s*weight[i1]*val[1]*phi;
			} // i1
		} // k1

		for (k1 = 0; k1<elementDOF1->col; k1++){
			i = elementDOF1->val[k][k1];
			b->val[i] += lb.val[k1];
			b->val[i+elementDOF1->dof] += lb.val[k1+elementDOF1->col];
		} // k1
	} // k
	free_dvector(&lb);
}

/**
 * \fn void assembleRHSCrouzeixRaviartcurlMorley2d(dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1, dvector *uh)
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
 * \param *elementdofTran pointer to transpose of elementDOF
 * \return void
 */
void assembleRHSCrouzeixRaviartcurlMorley2d(dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1, dvector *uh)
{
	int i,j,k,i1,j1,k1;
	
	double phi, phi1[2], val[2];
	double x[2], **gradLambda, *nve[3], **vertices;
	double s;

	int num_qp;
	double lambdas[100][3], weight[100];
			
	dvector lb;
	create_dvector(elementDOF1->col, &lb);
	/************************************************** right hand side b *****************************************************************/
	create_dvector(elementDOF1->dof, b);
	num_qp = getNumQuadPoints(elementDOF0->dop+elementDOF1->dop-1, 2); // elementDOF0->dop+elementDOF1->dop-1
	init_Gauss2d(num_qp, lambdas, weight);

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vertices = elements->vertices[k];
		s = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		for(i=0;i<3;i++){
			j = elementEdge->val[k][i];
			nve[i] = edges->nvector[j];
		}
		// end set parameters
        
		init_dvector(&lb, 0.0);
		for (i = 0; i<elementDOF1->col; i++)
		{
			for (i1 = 0; i1<num_qp; i1++)
			{
				morley_basis1(lambdas[i1], gradLambda, nve, i, phi1);
				val[0]=0; val[1]=0;
				for (j = 0; j<elementDOF0->col; j++){
					cr_basis(2, lambdas[i1], j, &phi);
					j1 = elementDOF0->val[k][j];
					val[0] += phi * uh->val[j1];
					val[1] += phi * uh->val[j1+elementDOF0->dof];
				}
				lb.val[i] += s*weight[i1] * (val[0]*phi1[1]-val[1]*phi1[0]);
			} // i1
		} // k1

		for (k1 = 0; k1<elementDOF1->col; k1++)
		{
			i = elementDOF1->val[k][k1];
			b->val[i] += lb.val[k1];
		} // k1
	} // k
	free_dvector(&lb);
}

/**
* \fn void assembleMassmatrixP02d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF)
* \brief assemble mass matrix
* \param *M pointer to mass matrix
* \param *elements pointer to the structure of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void assembleMassmatrixP02d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF)
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
* \fn void assembleBiGradMINIsymtensorDiagBlockRepeat2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double mu)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assembleBiGradMINIsymtensorDiagBlockRepeat2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double mu)
{
	int i, j, k, l;

	//	create_dvector(A->row, b);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element, edge, node;

	double phi, phi1[2], phi2[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double val, s, **gradLambda;
	int count;

	int num_qp;
	double lambdas[100][3], weight[100];

	
	/***************************** stiffness matrix A ***************************************/
	int *ia, *ja;
	double *va;
	int N = elementDOF->col*elementDOF->col*elements->row;
	ia = (int*)malloc(3*N * sizeof(int));
	ja = (int*)malloc(3*N * sizeof(int));
	va = (double*)malloc(3*N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(elementDOF->col, elementDOF->col, &lA);

	num_qp=getNumQuadPoints(elementDOF->dop * 2 - 2, 2); // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		// end set parameters

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			for (k2 = 0; k2<elementDOF->col; k2++)
			{
				val = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					mini_basis1(lambdas[i1], gradLambda, k1, phi1);
					mini_basis1(lambdas[i1], gradLambda, k2, phi2);
					val += s*weight[i1] * (phi1[0] * phi2[0] + phi1[1] * phi2[1]) * 2 * mu;
				}
				lA.val[k1][k2] += val;
			} // k2
		} // k1

		l = elementDOF->col*elementDOF->col * k;
		for (i = 0; i<elementDOF->col; i++)
		{
			for (j = 0; j<elementDOF->col; j++)
			{
				ia[l] = elementDOF->val[k][i];
				ja[l] = elementDOF->val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // i
		} // j
	} // k
	free_dden_matrix(&lA);

	for(i=0;i<N;i++){
		ia[i+N] = ia[i] + elementDOF->dof;
		ja[i+N] = ja[i] + elementDOF->dof;
		va[i+N] = va[i];
	}
	for(i=0;i<N;i++){
		ia[i+N*2] = ia[i] + elementDOF->dof*2;
		ja[i+N*2] = ja[i] + elementDOF->dof*2;
		va[i+N*2] = va[i]*2;
	}
	N *= 3;
	// remove zero elements and transform matrix A from its IJ format to its CSR format
	double eps = 0;
	if(eps<1e-20){
		dIJtoCSR(A, ia, ja, va, N, 0, 0);
	free(ia); free(ja); free(va);
	}
	else{
		int nzmax = 0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ) nzmax++;
		}

		int *Ia, *Ja;
		double *Va;
		Ia = (int*)malloc(nzmax * sizeof(int));
		Ja = (int*)malloc(nzmax * sizeof(int));
		Va = (double*)malloc(nzmax * sizeof(double));
		int cur=0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ){
				Ia[cur] = ia[i];
				Ja[cur] = ja[i];
				Va[cur] = va[i];
				cur++;
			}
		}
		free(ia); free(ja); free(va);
		dIJtoCSR(A, Ia, Ja, Va, nzmax, 0, 0);
		free(Ia); free(Ja); free(Va);
	}
}

/**
* \fn void assembleRotSMINILagrange2d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *BT pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \return void
*/
void assembleRotSMINILagrange2d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1)
{
	int i, j, k, l;

	int element[2], edge, node;

	double phi, phi1[3], phi2[3], val[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, xs[3], ys[3], s, **gradLambda;

	int num_qp;
	double lambdas[100][3], weight[100], gauss[100][3];

	/***************** matrix A *****************/
	int *ia, *ja;
	double *va;
	int N = 4*elementDOF0->col*elementDOF1->col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA[2]; // local A
	create_dden_matrix(elementDOF0->col, elementDOF1->col, lA);
	create_dden_matrix(elementDOF0->col, elementDOF1->col, lA+1);
	
	// num_qp = getNumQuadPoints_ShunnWilliams(elementDOF->dop * 2, 2); // the number of numerical intergation points
	// init_ShunnWilliams2d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	num_qp = getNumQuadPoints(elementDOF0->dop + elementDOF1->dop - 1, 2); // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		// end set parameters

		init_dden_matrix(&lA[0], 0.0);init_dden_matrix(&lA[1], 0.0);
		for (k1 = 0; k1<elementDOF0->col; k1++){
			for (k2 = 0; k2<elementDOF1->col; k2++){
				val[0] = 0; val[1] = 0;
				for (i1 = 0; i1<num_qp; i1++){
					mini_basis1(lambdas[i1], gradLambda, k1, phi1);
					lagrange_basis(lambdas[i1], k2, elementDOF1->dop, &phi);
					val[0] -= s*weight[i1] * phi1[1] * phi;
					val[1] += s*weight[i1] * phi1[0] * phi;
				}
				lA[0].val[k1][k2] += val[0]; lA[1].val[k1][k2] += val[1];
			} // k2
		} // k1

		l = 4*elementDOF0->col*elementDOF1->col * k;
		for (i = 0; i<elementDOF0->col; i++)
		{
			for (j = 0; j<elementDOF1->col; j++)
			{
				ia[l] = elementDOF0->val[k][i];
				ja[l] = elementDOF1->val[k][j];
				va[l] = lA[0].val[i][j];
				l++;
				ia[l] = elementDOF0->val[k][i] + elementDOF0->dof;
				ja[l] = elementDOF1->val[k][j] + elementDOF1->dof;
				va[l] = lA[1].val[i][j];
				l++;
				ia[l] = elementDOF0->val[k][i] + elementDOF0->dof*2;
				ja[l] = elementDOF1->val[k][j];
				va[l] = lA[1].val[i][j];
				l++;
				ia[l] = elementDOF0->val[k][i] + elementDOF0->dof*2;
				ja[l] = elementDOF1->val[k][j] + elementDOF1->dof;
				va[l] = lA[0].val[i][j];
				l++;
			} // i
		} // j
	} // k
	free_dden_matrix(&lA[0]);
	free_dden_matrix(&lA[1]);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	double eps = 0;
	if(eps<1e-20){
		dIJtoCSR(A, ia, ja, va, N, 0, 0);
		free(ia); free(ja); free(va);
	}
	else{
		int nzmax = 0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ) nzmax++;
		}

		int *Ia, *Ja;
		double *Va;
		Ia = (int*)malloc(nzmax * sizeof(int));
		Ja = (int*)malloc(nzmax * sizeof(int));
		Va = (double*)malloc(nzmax * sizeof(double));
		int cur=0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ){
				Ia[cur] = ia[i];
				Ja[cur] = ja[i];
				Va[cur] = va[i];
				cur++;
			}
		}
		free(ia); free(ja); free(va);
		dIJtoCSR(A, Ia, Ja, Va, nzmax, 0, 0);
		free(Ia); free(Ja); free(Va);
	}
}

/**
 * \fn void assembleRHShessMorleySMINI2d(dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1, dvector *uh)
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
 * \param *elementdofTran pointer to transpose of elementDOF
 * \return void
 */
void assembleRHShessMorleySMINI2d(dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1, dvector *uh)
{
	int i,j,k,i1,j1,k1;
	
// write_dvector4Matlab(uh, "output/uh.dat");//////////
// printf("elementDOF0:\n");
// for(i=0;i<elementDOF0->row;i++)
// {
// for(j=0;j<elementDOF0->col;j++)
// printf("%d, ", elementDOF0->val[i][j]);
// printf("\n");
// }
// printf("elementDOF1:\n");
// for(i=0;i<elementDOF1->row;i++)
// {
// for(j=0;j<elementDOF1->col;j++)
// printf("%d, ", elementDOF1->val[i][j]);
// printf("\n");
// }

	double phi, phi1[3], val[3];
	double x[2], **gradLambda, *nve[3], **vertices;
	double s;

	int num_qp;
	double lambdas[100][3], weight[100];
			
	dvector lb;
	create_dvector(3*elementDOF1->col, &lb);
	/************************************************** right hand side b *****************************************************************/
	create_dvector(3*elementDOF1->dof, b);
	num_qp = getNumQuadPoints(elementDOF0->dop+elementDOF1->dop-2, 2); 
	init_Gauss2d(num_qp, lambdas, weight);

	for (k = 0; k<elements->row; k++){
		// set parameters
		vertices = elements->vertices[k];
		s = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		for(i=0;i<3;i++){
			j = elementEdge->val[k][i];
			nve[i] = edges->nvector[j];
		}
		// end set parameters
        
		init_dvector(&lb, 0.0);
		for (i = 0; i<elementDOF1->col; i++){
			for (i1 = 0; i1<num_qp; i1++){
				mini_basis(lambdas[i1], i, &phi);
				val[0]=0; val[1]=0;; val[2]=0;
				for (j = 0; j<elementDOF0->col; j++){
					morley_basis2(gradLambda, nve, j, phi1);
					j1 = elementDOF0->val[k][j];
					val[0] += phi1[0] * uh->val[j1];
					val[1] += phi1[1] * uh->val[j1];
					val[2] += phi1[2] * uh->val[j1];
				}
				lb.val[i] += s*weight[i1]*val[0]*phi;
				lb.val[i+elementDOF1->col] += s*weight[i1]*val[1]*phi;
				lb.val[i+elementDOF1->col*2] += s*weight[i1]*val[2]*phi*2;
			} // i1
		} // k1

		for (k1 = 0; k1<elementDOF1->col; k1++){
			i = elementDOF1->val[k][k1];
			b->val[i] += lb.val[k1];
			b->val[i+elementDOF1->dof] += lb.val[k1+elementDOF1->col];
			b->val[i+elementDOF1->dof*2] += lb.val[k1+elementDOF1->col*2];
		} // k1
	} // k
	free_dvector(&lb);
}

/**
 * \fn void assembleRHSSMINIhessMorley2d(dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1, dvector *uh)
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
 * \param *elementdofTran pointer to transpose of elementDOF
 * \return void
 */
void assembleRHSSMINIhessMorley2d(dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1, dvector *uh)
{
	int i,j,k,i1,j1,k1;
	
	double phi, phi1[3], val[3];
	double x[2], **gradLambda, *nve[3], **vertices;
	double s;

	int num_qp;
	double lambdas[100][3], weight[100];
			
	dvector lb;
	create_dvector(elementDOF1->col, &lb);
	/************************************************** right hand side b *****************************************************************/
	create_dvector(elementDOF1->dof, b);
	num_qp = getNumQuadPoints(elementDOF0->dop+elementDOF1->dop-2, 2);
	init_Gauss2d(num_qp, lambdas, weight);

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vertices = elements->vertices[k];
		s = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		for(i=0;i<3;i++){
			j = elementEdge->val[k][i];
			nve[i] = edges->nvector[j];
		}
		// end set parameters
        
		init_dvector(&lb, 0.0);
		for (i = 0; i<elementDOF1->col; i++)
		{
			for (i1 = 0; i1<num_qp; i1++)
			{
				// morley_basis1(lambdas[i1], gradLambda, nve, i, phi1);
				morley_basis2(gradLambda, nve, i, phi1);
				val[0]=0; val[1]=0; val[2]=0;
				for (j = 0; j<elementDOF0->col; j++){
					// cr_basis(2, lambdas[i1], j, &phi);
					mini_basis(lambdas[i1], j, &phi);
					j1 = elementDOF0->val[k][j];
					val[0] += phi * uh->val[j1];
					val[1] += phi * uh->val[j1+elementDOF0->dof];
					val[2] += phi * uh->val[j1+elementDOF0->dof*2];
				}
				lb.val[i] += s*weight[i1] * (val[0]*phi1[0]+val[1]*phi1[1]+val[2]*phi1[2]*2);
			} // i1
		} // k1

		for (k1 = 0; k1<elementDOF1->col; k1++)
		{
			i = elementDOF1->val[k][k1];
			b->val[i] += lb.val[k1];
		} // k1
	} // k
	free_dvector(&lb);
}

/**
* \fn void assembleweightedMassmatrixHuZhang2d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF, double lambda, double mu)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *BT pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assembleweightedMassmatrixHuZhang2d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF, double lambda, double mu)
{
	int i, j, k, l;

	int element[2], edge, node;

	double phi, phi1[3], phi2[3], val1;
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart;
	int count;

	int num_qp;
	double lambdas[100][3], weight[100], gauss[100][3];

	/******************** matrix A ********************/
	int *ia, *ja;
	double *va;
	int N = elementDOF->col*elementDOF->col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(elementDOF->col, elementDOF->col, &lA);

	// num_qp = getNumQuadPoints_ShunnWilliams(elementDOF->dop * 2, 2); // the number of numerical intergation points
	// init_ShunnWilliams2d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	num_qp = getNumQuadPoints(elementDOF->dop * 2, 2); // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		// end set parameters

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF->col; k1++){
			for (k2 = 0; k2<elementDOF->col; k2++){
				val1 = 0;
				for (i1 = 0; i1<num_qp; i1++){
					huzhang_basis(lambdas[i1], elements->nvector[k], elements->tvector[k], k1, elementDOF[0].dop, phi1);
					huzhang_basis(lambdas[i1], elements->nvector[k], elements->tvector[k], k2, elementDOF[0].dop, phi2);
					if (lambda>-0.5)
						val1 += s*weight[i1] * ((phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) - lambda / (2 * lambda + 2 * mu)*(phi1[0] + phi1[1])*(phi2[0] + phi2[1])) / (2 * mu);
					else
						val1 += s*weight[i1] * ((phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) - 1.0 / 2.0*(phi1[0] + phi1[1])*(phi2[0] + phi2[1])) / (2 * mu);
				}
				lA.val[k1][k2] += val1;
			} // k2
		} // k1

		l = elementDOF->col*elementDOF->col * k;
		for (i = 0; i<elementDOF->col; i++)
		{
			for (j = 0; j<elementDOF->col; j++)
			{
				ia[l] = elementDOF->val[k][i];
				ja[l] = elementDOF->val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // i
		} // j
	} // k
	free_dden_matrix(&lA);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	double eps = 0;
	if(eps<1e-20){
		dIJtoCSR(A, ia, ja, va, N, 0, 0);
	free(ia); free(ja); free(va);
	}
	else{
		int nzmax = 0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ) nzmax++;
		}

		int *Ia, *Ja;
		double *Va;
		Ia = (int*)malloc(nzmax * sizeof(int));
		Ja = (int*)malloc(nzmax * sizeof(int));
		Va = (double*)malloc(nzmax * sizeof(double));
		int cur=0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ){
				Ia[cur] = ia[i];
				Ja[cur] = ja[i];
				Va[cur] = va[i];
				cur++;
			}
		}
		free(ia); free(ja); free(va);
		dIJtoCSR(A, Ia, Ja, Va, nzmax, 0, 0);
		free(Ia); free(Ja); free(Va);
	}
}

/**
* \fn void assembleDivHuZhangL2poly2d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *BT pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \return void
*/
void assembleDivHuZhangL2poly2d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF)
{
	int i, j, k, l;

	int element[2], edge, node;

	double phi, phi1[3], phi2[3], val[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, xs[3], ys[3], s, **gradLambda;

	int num_qp;
	double lambdas[100][3], weight[100], gauss[100][3];

	/***************** matrix A *****************/
	int *ia, *ja;
	double *va;
	int N = 2*elementDOF[0].col*elementDOF[1].col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA[2]; // local A
	create_dden_matrix(elementDOF[0].col, elementDOF[1].col, lA);
	create_dden_matrix(elementDOF[0].col, elementDOF[1].col, lA+1);
	
	// num_qp = getNumQuadPoints_ShunnWilliams(elementDOF->dop * 2, 2); // the number of numerical intergation points
	// init_ShunnWilliams2d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	num_qp = getNumQuadPoints(elementDOF[0].dop + elementDOF[1].dop - 1, 2); // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		// end set parameters

		init_dden_matrix(&lA[0], 0.0);init_dden_matrix(&lA[1], 0.0);
		for (k1 = 0; k1<elementDOF[0].col; k1++){
			for (k2 = 0; k2<elementDOF[1].col; k2++){
				val[0] = 0; val[1] = 0;
				for (i1 = 0; i1<num_qp; i1++){
					huzhang_basisDIV(lambdas[i1], gradLambda, elements->nvector[k], elements->tvector[k], k1, elementDOF[0].dop, phi1);
					lagrange_basis(lambdas[i1], k2, elementDOF[1].dop, &phi);
					val[0] += s*weight[i1] * phi1[0] * phi;
					val[1] += s*weight[i1] * phi1[1] * phi;
				}
				lA[0].val[k1][k2] += val[0]; lA[1].val[k1][k2] += val[1];
			} // k2
		} // k1

		l = 2*elementDOF[0].col*elementDOF[1].col * k;
		for (i = 0; i<elementDOF[0].col; i++)
		{
			for (j = 0; j<elementDOF[1].col; j++)
			{
				ia[l] = elementDOF[0].val[k][i];
				ja[l] = elementDOF[1].val[k][j];
				va[l] = lA[0].val[i][j];
				l++;
				ia[l] = elementDOF[0].val[k][i];
				ja[l] = elementDOF[1].val[k][j] + elementDOF[1].dof;
				va[l] = lA[1].val[i][j];
				l++;
			} // i
		} // j
	} // k
	free_dden_matrix(&lA[0]);
	free_dden_matrix(&lA[1]);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	double eps = 0;
	if(eps<1e-20){
		dIJtoCSR(A, ia, ja, va, N, 0, 0);
		free(ia); free(ja); free(va);
	}
	else{
		int nzmax = 0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ) nzmax++;
		}

		int *Ia, *Ja;
		double *Va;
		Ia = (int*)malloc(nzmax * sizeof(int));
		Ja = (int*)malloc(nzmax * sizeof(int));
		Va = (double*)malloc(nzmax * sizeof(double));
		int cur=0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ){
				Ia[cur] = ia[i];
				Ja[cur] = ja[i];
				Va[cur] = va[i];
				cur++;
			}
		}
		free(ia); free(ja); free(va);
		dIJtoCSR(A, Ia, Ja, Va, nzmax, 0, 0);
		free(Ia); free(Ja); free(Va);
	}
}

/**
* \fn void assembleRHSdgPolyVector2d(dvector *b, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, void (*f)(double *, double *, double *), double *paras)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *BT pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \return void
*/
void assembleRHSdgPolyVector2d(dvector *b, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, void (*f)(double *, double *, double *), double *paras)
{
	int i, j, k, l;

	int element[2], edge, node;

	double phi, phi1[3], phi2[3], val[2];
	int k1, k2, i1, j1, l1, l2, ej;
	// double x, y, xs[3], ys[3], s, *eta, *xi;
	double x[2], **vertices, s;
	int rowstart[2], row21[2], taustart;
	int count;

	int num_qp;
	double lambdas[100][3], weight[100], gauss[100][3];

	create_dvector(elementDOF->dof * 2, b);
	dvector lb;
	create_dvector(elementDOF->col * 2, &lb);
	/************************************************** right hand side b *****************************************************************/
	num_qp = 49; // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		vertices = elements->vertices[k];
		// end set parameters

		init_dvector(&lb, 0.0);
		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			for (i1 = 0; i1<num_qp; i1++)
			{
				lagrange_basis(lambdas[i1], k1, elementDOF->dop, &phi);
				baryToCart2d(lambdas[i1], x, vertices);
				// linearElas2d_f(x, val, paras);
				f(x, val, paras);
				lb.val[k1] += s*weight[i1] * val[0]*phi;
				lb.val[k1 + elementDOF->col] += s*weight[i1] * val[1]*phi;
			} // i1
		} // k1 

		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			i = elementDOF->val[k][k1];
			b->val[i] += lb.val[k1];
			b->val[i + elementDOF->dof] += lb.val[k1 + elementDOF->col];
		} // k1 
	} // k
	free_dvector(&lb);
}

/**
* \fn void assembleJumpL2poly2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double parapenalty)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *BT pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assembleJumpL2poly2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double parapenalty)
{
	int i, j, k, l;

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, xs[3], ys[3], s, *eta, *xi;
	int curnode[2];
	int count;

	// int num_qp = getNumQuadPoints(elementDOF->dop * 2, 1); // the number of numerical intergation points
	// double gauss[num_qp][2];
	// init_Gauss1D(num_qp, 1, gauss); // gauss intergation initial

	int num_qp;
	double lambdas[100][2], weight[100];
	

	/***************************** stiffness matrix A ***************************************/
	int *ia, *ja;
	double *va;
	int N = 0;
	int patchnodes[200];
	for(edge=0;edge<edges->row;edge++){
		element[0] = edges->val[edge][2];
		element[1] = edges->val[edge][3];

		count = elementDOF->dop + 1;
		if (element[1] != -1) count*=2;

		N += count*count;
	}

	ia = (int*)calloc(2*N, sizeof(int));
	ja = (int*)calloc(2*N, sizeof(int));
	va = (double*)calloc(2*N, sizeof(double));

	double elen, C11;
	ddenmat lA[2]; // local A
	num_qp = getNumQuadPoints(elementDOF->dop * 2, 1); // the number of numerical intergation points
	init_Gauss1d(num_qp, lambdas, weight);
	l = 0;
	for (edge = 0; edge<edges->row; edge++)
	{
		//		edgeNode[0]=edges->val[edge][0];
		//		edgeNode[1]=edges->val[edge][1];
		element[0] = edges->val[edge][2];
		element[1] = edges->val[edge][3];
		elen = edges->length[edge];

//		C11 = elen;

		count = 0;
		if (elementDOF->dop == 0)
		{
			patchnodes[count] = elementDOF->val[element[0]][0];
			count++;
			if (element[1] != -1)
			{
				patchnodes[count] = elementDOF->val[element[1]][0];
				count++;
			}
		} // if(elementDOF->dop==0)
		else
		{
			for (i = 0; i<3; i++)
			{
				if (elementEdge->val[element[0]][i] == edge)
					break;
			}
			patchnodes[count] = elementDOF->val[element[0]][(i + 1) % 3];
			count++;
			patchnodes[count] = elementDOF->val[element[0]][(i + 2) % 3];
			count++;
			for (j = 0; j<elementDOF->dop - 1; j++)
			{
				patchnodes[count] = elementDOF->val[element[0]][3 + i*(elementDOF->dop - 1) + j];
				count++;
			}

			if (element[1] != -1)
			{
				for (i = 0; i<3; i++)
				{
					if (elementEdge->val[element[1]][i] == edge)
						break;
				}
				patchnodes[count] = elementDOF->val[element[1]][(i + 1) % 3];
				count++;
				patchnodes[count] = elementDOF->val[element[1]][(i + 2) % 3];
				count++;
				for (j = 0; j<elementDOF->dop - 1; j++)
				{
					patchnodes[count] = elementDOF->val[element[1]][3 + i*(elementDOF->dop - 1) + j];
					count++;
				}
			}
		} // if(elementDOF->dop==0) else

		create_dden_matrix(count, count, lA);
		create_dden_matrix(count, count, lA+1);
		for (k1 = 0; k1<count; k1++)
		{
			curnode[0] = patchnodes[k1];
			curnode[1] = curnode[0] + elementDOF->dof;
			for (k2 = 0; k2<count; k2++)
			{
				j = patchnodes[k2];
				// c11
				for (i1 = 0; i1<num_qp; i1++)
				{
					jumpOperatorVector(lambdas[i1][0], lambdas[i1][1], edge, elements, elementEdge, edges, elementDOF, curnode[0], phi1);
					jumpOperatorVector(lambdas[i1][0], lambdas[i1][1], edge, elements, elementEdge, edges, elementDOF, j, phi2);
					lA[0].val[k1][k2] += parapenalty * weight[i1] * phi1[0] * phi2[0];
				}
				// c22
				for (i1 = 0; i1<num_qp; i1++)
				{
					jumpOperatorVector(lambdas[i1][0], lambdas[i1][1], edge, elements, elementEdge, edges, elementDOF, curnode[1], phi1);
					jumpOperatorVector(lambdas[i1][0], lambdas[i1][1], edge, elements, elementEdge, edges, elementDOF, j + elementDOF->dof, phi2);
					lA[1].val[k1][k2] += parapenalty * weight[i1] * phi1[1] * phi2[1];
				}
			} // k2
		} // k1

		for (i = 0; i<count; i++)
		{
			for (j = 0; j<count; j++)
			{
				ia[l] = patchnodes[i];
				ja[l] = patchnodes[j];
				va[l] = lA[0].val[i][j];
				l++;
			} // i
		} // j

		for (i = 0; i<count; i++)
		{
			for (j = 0; j<count; j++)
			{
				ia[l] = patchnodes[i] + elementDOF->dof;
				ja[l] = patchnodes[j] + elementDOF->dof;
				va[l] = lA[1].val[i][j];
				l++;
			} // i
		} // j

		free_dden_matrix(lA);
		free_dden_matrix(lA+1);
	} // edge

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	double eps = 0;
	if(eps<1e-20){
		dIJtoCSR(A, ia, ja, va, 2*N, 0, 0);
		free(ia); free(ja); free(va);
	}
	else{
		int nzmax = 0;
		for(i=0; i<2*N; i++){
			if(fabs(va[i]) > eps ) nzmax++;
		}

		int *Ia, *Ja;
		double *Va;
		Ia = (int*)malloc(nzmax * sizeof(int));
		Ja = (int*)malloc(nzmax * sizeof(int));
		Va = (double*)malloc(nzmax * sizeof(double));
		int cur=0;
		for(i=0; i<2*N; i++){
			if(fabs(va[i]) > eps ){
				Ia[cur] = ia[i];
				Ja[cur] = ja[i];
				Va[cur] = va[i];
				cur++;
			}
		}
		free(ia); free(ja); free(va);
		dIJtoCSR(A, Ia, Ja, Va, nzmax, 2*elementDOF->dof, 2*elementDOF->dof);
		free(Ia); free(Ja); free(Va);
	}
}

/**
* \fn void assembleweightedMassmatrixHuangZhou2d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF, double lambda, double mu)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *BT pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assembleweightedMassmatrixHuangZhou2d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF, double lambda, double mu)
{
	int i, j, k, l;

	int element[2], edge, node;

	double phi, phi1[3], phi2[3], val1;
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, xs[3], ys[3], s, *eta, *xi;
	int count;

	int num_qp;
	double lambdas[100][3], weight[100], gauss[100][3];

	/******************** matrix A ********************/
	int *ia, *ja;
	double *va;
	int N = elementDOF->col*elementDOF->col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(elementDOF->col, elementDOF->col, &lA);

	// num_qp = getNumQuadPoints_ShunnWilliams(elementDOF->dop * 2, 2); // the number of numerical intergation points
	// init_ShunnWilliams2d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	num_qp = getNumQuadPoints(elementDOF->dop * 2, 2); // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		// end set parameters

		init_dden_matrix(&lA, 0.0);
		for (k1 = 0; k1<elementDOF->col; k1++){
			for (k2 = 0; k2<elementDOF->col; k2++){
				val1 = 0;
				for (i1 = 0; i1<num_qp; i1++){
					divS_huangzhou_basis(lambdas[i1], s, elements->nvector[k], elements->tvector[k], elements->tij[k], k1, phi1);
					divS_huangzhou_basis(lambdas[i1], s, elements->nvector[k], elements->tvector[k], elements->tij[k], k2, phi2);
					if (lambda>-0.5)
						val1 += s*weight[i1] * ((phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) - lambda / (2 * lambda + 2 * mu)*(phi1[0] + phi1[1])*(phi2[0] + phi2[1])) / (2 * mu);
					else // lambda = inf
						val1 += s*weight[i1] * ((phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) - 1.0 / 2.0*(phi1[0] + phi1[1])*(phi2[0] + phi2[1])) / (2 * mu);
				}
				lA.val[k1][k2] += val1;
			} // k2
		} // k1

		l = elementDOF->col*elementDOF->col * k;
		for (i = 0; i<elementDOF->col; i++)
		{
			for (j = 0; j<elementDOF->col; j++)
			{
				ia[l] = elementDOF->val[k][i];
				ja[l] = elementDOF->val[k][j];
				va[l] = lA.val[i][j];
				l++;
			} // i
		} // j
	} // k
	free_dden_matrix(&lA);
	
	// remove zero elements and transform matrix A from its IJ format to its CSR format
	double eps = 0;
	if(eps<1e-20){
		dIJtoCSR(A, ia, ja, va, N, 0, 0);
	free(ia); free(ja); free(va);
	}
	else{
		int nzmax = 0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ) nzmax++;
		}

		int *Ia, *Ja;
		double *Va;
		Ia = (int*)malloc(nzmax * sizeof(int));
		Ja = (int*)malloc(nzmax * sizeof(int));
		Va = (double*)malloc(nzmax * sizeof(double));
		int cur=0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ){
				Ia[cur] = ia[i];
				Ja[cur] = ja[i];
				Va[cur] = va[i];
				cur++;
			}
		}
		free(ia); free(ja); free(va);
		dIJtoCSR(A, Ia, Ja, Va, nzmax, 0, 0);
		free(Ia); free(Ja); free(Va);
	}
}

/**
* \fn void assembleDivHuangZhouL2poly2d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *BT pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \return void
*/
void assembleDivHuangZhouL2poly2d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF)
{
	int i, j, k, l;

	int element[2], edge, node;

	double phi, phi1[3], phi2[3], val[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, xs[3], ys[3], s, **gradLambda;

	int num_qp;
	double lambdas[100][3], weight[100], gauss[100][3];

	/***************** matrix A *****************/
	int *ia, *ja;
	double *va;
	int N = 2*elementDOF[0].col*elementDOF[1].col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA[2]; // local A
	create_dden_matrix(elementDOF[0].col, elementDOF[1].col, lA);
	create_dden_matrix(elementDOF[0].col, elementDOF[1].col, lA+1);
	
	// num_qp = getNumQuadPoints_ShunnWilliams(elementDOF->dop * 2, 2); // the number of numerical intergation points
	// init_ShunnWilliams2d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	num_qp = getNumQuadPoints(elementDOF[0].dop + elementDOF[1].dop - 1, 2); // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		// end set parameters

		init_dden_matrix(&lA[0], 0.0); init_dden_matrix(&lA[1], 0.0);
		for (k1 = 0; k1<elementDOF[0].col; k1++){
			for (k2 = 0; k2<elementDOF[1].col; k2++){
				val[0] = 0; val[1] = 0;
				for (i1 = 0; i1<num_qp; i1++){
					divS_huangzhou_basisDIV(lambdas[i1], gradLambda, s, elements->nvector[k], elements->tvector[k], elements->tij[k], k1, phi1);
					lagrange_basis(lambdas[i1], k2, elementDOF[1].dop, &phi);
					val[0] += s*weight[i1] * phi1[0] * phi;
					val[1] += s*weight[i1] * phi1[1] * phi;
				}
				lA[0].val[k1][k2] += val[0]; lA[1].val[k1][k2] += val[1];
			} // k2
		} // k1

		l = 2*elementDOF[0].col*elementDOF[1].col * k;
		for (i = 0; i<elementDOF[0].col; i++)
		{
			for (j = 0; j<elementDOF[1].col; j++)
			{
				ia[l] = elementDOF[0].val[k][i];
				ja[l] = elementDOF[1].val[k][j];
				va[l] = lA[0].val[i][j];
				l++;
				ia[l] = elementDOF[0].val[k][i];
				ja[l] = elementDOF[1].val[k][j] + elementDOF[1].dof;
				va[l] = lA[1].val[i][j];
				l++;
			} // i
		} // j
	} // k
	free_dden_matrix(&lA[0]);
	free_dden_matrix(&lA[1]);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	double eps = 0;
	if(eps<1e-20){
		dIJtoCSR(A, ia, ja, va, N, 0, 0);
		free(ia); free(ja); free(va);
	}
	else{
		int nzmax = 0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ) nzmax++;
		}

		int *Ia, *Ja;
		double *Va;
		Ia = (int*)malloc(nzmax * sizeof(int));
		Ja = (int*)malloc(nzmax * sizeof(int));
		Va = (double*)malloc(nzmax * sizeof(double));
		int cur=0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ){
				Ia[cur] = ia[i];
				Ja[cur] = ja[i];
				Va[cur] = va[i];
				cur++;
			}
		}
		free(ia); free(ja); free(va);
		dIJtoCSR(A, Ia, Ja, Va, nzmax, 0, 0);
		free(Ia); free(Ja); free(Va);
	}
}

/**
 * \fn void assembleStiffmatrixElasLagrange(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param *nodes pointer to the nodes location of the triangulation
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *elementdofTran pointer to transpose of elementDOF
 * \param mu Lame constant or Poisson ratio of plate
 * \return void
 */
void assembleStiffmatrixElasLagrange(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
{
	int i,j,k,l;

	// int nvertices=nodes->row;
	// int nedges=edges->row;
	int element[2], edge, node;
	
	double phi, phi1[2], phi2[2];
	int k1,k2,i1,j1,l1,l2,ej;
	double x, y, xs[3], ys[3], s, **gradLambda;
	int count;
	double val;
	
	int num_qp;
	double lambdas[100][3], weight[100];

	num_qp=getNumQuadPoints(elementDOF->dop*2-2, 2); // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);

	/************************************************** stiffness matrix A *****************************************************************/
	int *ia, *ja;
	double *va;
	int N = 4*elementDOF->col*elementDOF->col*elements->row;
	ia = (int*)malloc(N * sizeof(int));
	ja = (int*)malloc(N * sizeof(int));
	va = (double*)malloc(N * sizeof(double));
	ddenmat lA; // local A
	create_dden_matrix(2*elementDOF->col, 2*elementDOF->col, &lA);

	for(k=0;k<elements->row;k++)
	{
		// set parameters
		s=elements->vol[k];
		// xi=elements->xi[k];
		// eta=elements->eta[k];
		gradLambda=elements->gradLambda[k];
		// end set parameters

		init_dden_matrix(&lA, 0.0);
		for(k1=0;k1<elementDOF->col;k1++)
		{
			for(k2=0;k2<elementDOF->col;k2++)
			{
				// a11
				val = 0;
				for (i1 = 0; i1<num_qp; i1++){
					lagrange_basis1(lambdas[i1], gradLambda, k1, elementDOF->dop, phi1);
					lagrange_basis1(lambdas[i1], gradLambda, k2, elementDOF->dop, phi2);
					val += s*weight[i1] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] / 2) * 2 * mu;
				}
				lA.val[k1][k2] += val;
				// a12
				val = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					lagrange_basis1(lambdas[i1], gradLambda, k1, elementDOF->dop, phi1);
					lagrange_basis1(lambdas[i1], gradLambda, k2, elementDOF->dop, phi2);
					val += s*weight[i1] * phi1[1] * phi2[0]/ 2 * 2 * mu;
				}
				lA.val[k1][k2+elementDOF->col] += val;
				// a21
				val = 0;
				for (i1 = 0; i1<num_qp; i1++)
				{
					lagrange_basis1(lambdas[i1], gradLambda, k1, elementDOF->dop, phi1);
					lagrange_basis1(lambdas[i1], gradLambda, k2, elementDOF->dop, phi2);
					val += s*weight[i1] * phi1[0] * phi2[1] / 2 * 2 * mu;
				}
				lA.val[k1+elementDOF->col][k2] += val;
				// a22
				val = 0;
				for (i1 = 0; i1<num_qp; i1++){
					lagrange_basis1(lambdas[i1], gradLambda, k1, elementDOF->dop, phi1);
					lagrange_basis1(lambdas[i1], gradLambda, k2, elementDOF->dop, phi2);
					val += s*weight[i1] * (phi1[0] * phi2[0] / 2 + phi1[1] * phi2[1]) * 2 * mu;
				}
				lA.val[k1+elementDOF->col][k2+elementDOF->col] += val;
			} // k2
		} // k1

		l = 4*elementDOF->col*elementDOF->col * k;
		for (i = 0; i<elementDOF->col; i++)
		{
			for (j = 0; j<elementDOF->col; j++)
			{
				// a11
				ia[l] = elementDOF->val[k][i];
				ja[l] = elementDOF->val[k][j];
				va[l] = lA.val[i][j];
				l++;
				// a12
				ia[l] = elementDOF->val[k][i];
				ja[l] = elementDOF->val[k][j] + elementDOF->dof;
				va[l] = lA.val[i][j+elementDOF->col];
				l++;
				// a21
				ia[l] = elementDOF->val[k][i] + elementDOF->dof;
				ja[l] = elementDOF->val[k][j];
				va[l] = lA.val[i+elementDOF->col][j];
				l++;
				// a22
				ia[l] = elementDOF->val[k][i] + elementDOF->dof;
				ja[l] = elementDOF->val[k][j] + elementDOF->dof;
				va[l] = lA.val[i+elementDOF->col][j+elementDOF->col];
				l++;
			} // i
		} // j
	} // k
	free_dden_matrix(&lA);

	// remove zero elements and transform matrix A from its IJ format to its CSR format
	double eps = 0;
	if(eps<1e-20){
		dIJtoCSR(A, ia, ja, va, N, 0, 0);
	free(ia); free(ja); free(va);
	}
	else{
		int nzmax = 0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ) nzmax++;
		}

		int *Ia, *Ja;
		double *Va;
		Ia = (int*)malloc(nzmax * sizeof(int));
		Ja = (int*)malloc(nzmax * sizeof(int));
		Va = (double*)malloc(nzmax * sizeof(double));
		int cur=0;
		for(i=0; i<N; i++){
			if(fabs(va[i]) > eps ){
				Ia[cur] = ia[i];
				Ja[cur] = ja[i];
				Va[cur] = va[i];
				cur++;
			}
		}
		free(ia); free(ja); free(va);
		dIJtoCSR(A, Ia, Ja, Va, nzmax, 0, 0);
		free(Ia); free(Ja); free(Va);
	}
}

/**
* \fn void assembleStiffmatrixVecLagrange(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assembleStiffmatrixVecLagrange(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu)
{
	int i, j, k, l;

	A->row = elementDOF->dof * 2;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	//	create_dvector(A->row, b);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[2], phi2[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, xs[3], ys[3], s, **gradLambda;
	int rowstart[2], row21[2], curnode[2];
	int count;


	int num_qp;
	double lambdas[100][3], weight[100];

	num_qp=getNumQuadPoints(elementDOF->dop * 2 - 2, 2); // the number of numerical intergation points
	init_Gauss2d(num_qp, lambdas, weight);

									  /************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)calloc(A->col, sizeof(int));
	for (i = 0; i<A->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (i = 0; i<elementDOF->dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF->col; k++)
			{
				node = elementDOF->val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count * 2;
		A->IA[elementDOF->dof + i + 1] = count * 2;

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
	for (i = 0; i<elementDOF->dof; i++)
	{
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF->col; k++)
			{

				node = elementDOF->val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		row21[0] = (A->IA[i + 1] - A->IA[i]) / 2;
		rowstart[1] = A->IA[i + elementDOF->dof];
		for (j = A->IA[i]; j<A->IA[i] + row21[0]; j++)
		{
			A->JA[j] = istart;
			A->JA[j + row21[0]] = istart + elementDOF->dof;
			A->JA[j - A->IA[i] + rowstart[1]] = A->JA[j];
			A->JA[j - A->IA[i] + rowstart[1] + row21[0]] = A->JA[j + row21[0]];
			istart = index[istart];
			index[A->JA[j]] = -1;
		}
	} // i
	free(index);

	// step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		gradLambda = elements->gradLambda[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			curnode[0] = elementDOF->val[k][k1];
			curnode[1] = curnode[0] + elementDOF->dof;

			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}

			for (k2 = 0; k2<elementDOF->col; k2++)
			{
				j = elementDOF->val[k][k2];
				// a11
				for (j1 = A->IA[curnode[0]]; j1<A->IA[curnode[0]] + row21[0]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp; i1++)
						{
							lagrange_basis1(lambdas[i1], gradLambda, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas[i1], gradLambda, k2, elementDOF->dop, phi2);
							A->val[j1] += s*weight[i1] * (phi1[0] * phi2[0] + phi1[1] * phi2[1]) * 2 * mu;
						}
						break;
					}
				} // j1
				  // a22
				for (j1 = A->IA[curnode[1]] + row21[1]; j1<A->IA[curnode[1] + 1]; j1++)
				{
					if (A->JA[j1] == (j + elementDOF->dof))
					{
						for (i1 = 0; i1<num_qp; i1++)
						{
							lagrange_basis1(lambdas[i1], gradLambda, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas[i1], gradLambda, k2, elementDOF->dop, phi2);
							A->val[j1] += s*weight[i1] * (phi1[0] * phi2[0]  + phi1[1] * phi2[1]) * 2 * mu;
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k
}

/**
 * \fnvoid sumNormalDerivativeEta(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, ddenmat *etas, double *sum)
 * \brief the sum operator for normal derivative of basis function with penalty parameter
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param edge the index of current edge
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *etas pointer to penalty parameter 
 * \param *sum pointer to the result of sum operator
 * \return void
 */
void sumNormalDerivativeEta(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, ddenmat *etas, double *sum)
{
	int i,j,l;
	double lambdas[3];
	double s, **gradLambda, phi[2], nve[2];
	int nodeindex1,nodeindex2;
	int element[2], edgeNode[2];
	int dop=elementDOF->dop;

	*sum=0;
	if(dop<=0)
		return;

	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];
	element[0]=edges->val[edge][2];
	element[1]=edges->val[edge][3];
	nve[0]=edges->nvector[edge][0];
	nve[1]=edges->nvector[edge][1];
	
	for(nodeindex1=0;nodeindex1<elementDOF->col;nodeindex1++)
    {
		if(elementDOF->val[element[0]][nodeindex1]==node)
			break;
	}

	if(nodeindex1<elementDOF->col)
	{
		for(i=0;i<3;i++)
		{
			if(elements->val[element[0]][i]==edgeNode[0])
				break;
		}
	
		for(j=0;j<3;j++)
		{
			if(elements->val[element[0]][j]==edgeNode[1])
				break;
		}
	
		l=3-i-j;
		lambdas[i]=lambda1;
		lambdas[j]=lambda2;
		lambdas[l]=0;
		s=elements->vol[element[0]];
		gradLambda=elements->gradLambda[element[0]];
		lagrange_basis1(lambdas, gradLambda, nodeindex1, dop, phi);
		*sum+=(phi[0]*nve[0]+phi[1]*nve[1])*etas->val[element[0]][l];
	}

	if(element[1]==-1)
			return;

	for(nodeindex2=0;nodeindex2<elementDOF->col;nodeindex2++)
    {
		if(elementDOF->val[element[1]][nodeindex2]==node)
			break;
	}

	if(nodeindex2<elementDOF->col)
	{
		for(i=0;i<3;i++)
		{
			if(elements->val[element[1]][i]==edgeNode[0])
				break;
		}
	
		for(j=0;j<3;j++)
		{
			if(elements->val[element[1]][j]==edgeNode[1])
				break;
		}
	
		l=3-i-j;
		lambdas[i]=lambda1;
		lambdas[j]=lambda2;
		lambdas[l]=0;
		s=elements->vol[element[1]];
		gradLambda=elements->gradLambda[element[1]];
		lagrange_basis1(lambdas, gradLambda, nodeindex2, dop, phi);
		*sum+=(phi[0]*nve[0]+phi[1]*nve[1])*etas->val[element[1]][l];
	}
}

/**
 * \fn void jumpOperatorVectorTensor(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, DOF *dofs, int node, double *jump)
 * \brief the jump operator for vector otimes normal derivative
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param edge the index of current edge
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *jump pointer to the result of jump operator
 * \return void
 */
void jumpOperatorVectorTensor(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i,j,l;
	int nodeindex, edgeindex, ni2;
	int element, edgeNode[2];
	double lambdas[3];
	double *eta, *xi, nv[2];
	double phi;
	double elen=edges->length[edge];

	jump[0]=0;
	jump[1]=0;
	jump[2]=0;

	if(node<0 && node>=elementDOF->dof*2)
		return;

	element=node/(elementDOF->col*2);
	ni2=node%(elementDOF->col*2);

	for(edgeindex=0;edgeindex<3;edgeindex++)
	{
		if(elementEdge->val[element][edgeindex]==edge)
			break;
	}

	if(edgeindex==3)
		return;

	xi=elements->xi[element];
	eta=elements->eta[element];
	nv[0]=-eta[edgeindex]/elen;
	nv[1]=xi[edgeindex]/elen;

	if(elementDOF->dop==0)
	{
		if(ni2<elementDOF->col)
		{
			jump[0]=nv[0];
			jump[1]=0;
			jump[2]=nv[1]/2;
		}
		else
		{
			jump[0]=0;
			jump[1]=nv[1];
			jump[2]=nv[0]/2;
		}
		return;
	}

	nodeindex=node%elementDOF->col;

	if(nodeindex<3 && nodeindex>=0)
	{
		if(nodeindex==edgeindex)
			return;
	}
	else if(nodeindex<3*elementDOF->dop && nodeindex>=3)
	{
		if((nodeindex-3)/(elementDOF->dop-1) != edgeindex)
			return;
	}
	else
		return;

	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];
	for(i=0;i<3;i++)
	{
		if(elements->val[element][i]==edgeNode[0])
			break;
	}
	for(j=0;j<3;j++)
	{
		if(elements->val[element][j]==edgeNode[1])
			break;
	}
	l=3-i-j;
	lambdas[i]=lambda1;
	lambdas[j]=lambda2;
	lambdas[l]=0;

	lagrange_basis(lambdas, nodeindex, elementDOF->dop, &phi);


	if(ni2<elementDOF->col)
	{
		jump[0]=phi*nv[0];
		jump[1]=0;
		jump[2]=phi*nv[1]/2;
	}
	else
	{
		jump[0]=0;
		jump[1]=phi*nv[1];
		jump[2]=phi*nv[0]/2;
	}
}

/**
* \fn void jumpOperatorVector(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, DOF *dofs, int node, double *jump)
* \brief the jump operator for vector
* \param lambda1 the first length coordiante
* \param lambda2 the second length coordiante
* \param edge the index of current edge
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *elementDOF pointer to relation between elements and DOFs
* \param node index of current node variable
* \param *jump pointer to the result of jump operator
* \return void
*/
void jumpOperatorVector(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i, j, l;
	int nodeindex, edgeindex, ni2;
	int element, edgeNode[2];
	double lambdas[3];
	double *eta, *xi, nv[2], nve[2];
	double phi;
	double elen = edges->length[edge];

	jump[0] = 0;
	jump[1] = 0;

	if (node<0 && node >= elementDOF->dof * 2)
		return;

	element = (node % elementDOF->dof) / elementDOF->col;
	int li = node / elementDOF->dof;

	for (edgeindex = 0; edgeindex<3; edgeindex++)
	{
		if (elementEdge->val[element][edgeindex] == edge)
			break;
	}

	if (edgeindex == 3)
		return;

	double sgn = elements->eorien[element][edgeindex];

	if (elementDOF->dop == 0)
	{
		jump[li] = sgn;
		return;
	}

	nodeindex = node%elementDOF->col;

	if (nodeindex<3 && nodeindex >= 0)
	{
		if (nodeindex == edgeindex)
			return;
	}
	else if (nodeindex<3 * elementDOF->dop && nodeindex >= 3)
	{
		if ((nodeindex - 3) / (elementDOF->dop - 1) != edgeindex)
			return;
	}
	else
		return;

	edgeNode[0] = edges->val[edge][0];
	edgeNode[1] = edges->val[edge][1];
	for (i = 0; i<3; i++)
	{
		if (elements->val[element][i] == edgeNode[0])
			break;
	}
	for (j = 0; j<3; j++)
	{
		if (elements->val[element][j] == edgeNode[1])
			break;
	}
	l = 3 - i - j;
	lambdas[i] = lambda1;
	lambdas[j] = lambda2;
	lambdas[l] = 0;

	lagrange_basis(lambdas, nodeindex, elementDOF->dop, &phi);

	jump[li] = phi*sgn;
}


/**
 * \fn void jumpOperator(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
 * \brief the jump operator for normal derivative of basis function
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param edge the index of current edge
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *jump pointer to the result of jump operator
 * \return void
 */
void jumpOperator(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i;
	int nodeindex1,nodeindex2, edgeindex;
	int element[2], edgeNode[2];
	double jumpminus;
	int isInedge=0;
	int nedges=edges->row;
	int dop=elementDOF->dop;

	*jump=0;
	if(dop<=0)
		return;

	double elen=edges->length[edge];
	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];
	element[0]=edges->val[edge][2];
	element[1]=edges->val[edge][3];
	
	for(nodeindex1=0;nodeindex1<elementDOF->col;nodeindex1++)
    {
		if(elementDOF->val[element[0]][nodeindex1]==node)
			break;
	}
	
	if(element[1]==-1)
	{
		if(nodeindex1==elementDOF->col)
			return;
	}
	else
	{
		for(nodeindex2=0;nodeindex2<elementDOF->col;nodeindex2++)
		{
			if(elementDOF->val[element[1]][nodeindex2]==node)
				break;
		}
		
		if(nodeindex1==elementDOF->col && nodeindex2==elementDOF->col)
			return;
	}
	
	for(edgeindex=0;edgeindex<elementEdge->col;edgeindex++)
	{
		if(elementEdge->val[element[0]][edgeindex]==edge)
			break;
	}
	
	for(i=0;i<2;i++)
	{
		if(edgeNode[i]==node)
		{
			isInedge=1;
			break;
		}
	}
	
	if(isInedge==0)
	{
		for(i=0;i<dop-1;i++)
		{
			if(elementDOF->val[element[0]][3+edgeindex*(dop-1)+i]==node)
			{
				isInedge=1;
				break;
			}
		}
	}
	
	if(isInedge==1)
	{
		getinfo4jump(lambda1, lambda2, edgeNode, elen, element[0], elements, dop, nodeindex1, jump);
		if(element[1]!=-1)
		{			
			getinfo4jump(lambda1, lambda2, edgeNode, elen, element[1], elements, dop, nodeindex2, &jumpminus);
				*jump += jumpminus;
		}
	}
	else
	{
		if(nodeindex1<elementDOF->col)
			getinfo4jump(lambda1, lambda2, edgeNode, elen, element[0], elements, dop, nodeindex1, jump);
		else{
				if(element[1]!=-1)
				{
					if(nodeindex2<elementDOF->col)
						getinfo4jump(lambda1, lambda2, edgeNode, elen, element[1], elements, dop, nodeindex2, jump);
				}
			}
	}
}

/**
 * \fn void getinfo4jump(double lambda1, double lambda2, int *edgeNode, double elen, int element, ELEMENT *elements, int dop, int index, double *jump)
 * \brief the normal derivative of basis function
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param *edgeNode pointer to the two indices of current edge
 * \param elen the length of cuurent edge
 * \param element the index of current element
 * \param *elements pointer to the structure of the triangulation
 * \param dop number of degrees of polynomial
 * \param index index of current node variable
 * \param *jump pointer to the normal tensor of gradient of basis function
 * \return void
 */
void getinfo4jump(double lambda1, double lambda2, int *edgeNode, double elen, int element, ELEMENT *elements, int dop, int index, double *jump)
{
	int i,j,l;
	double lambdas[3];
	double s, **gradLambda, phi[2], *nv;
		
	for(i=0;i<3;i++)
	{
		if(elements->val[element][i]==edgeNode[0])
			break;
	}
	
	for(j=0;j<3;j++)
	{
		if(elements->val[element][j]==edgeNode[1])
			break;
	}
	
	l=3-i-j;
	lambdas[i]=lambda1;
	lambdas[j]=lambda2;
	lambdas[l]=0;
	
	s=elements->vol[element];
	gradLambda=elements->gradLambda[element];
	nv=elements->nvector[element][l];
	
	lagrange_basis1(lambdas, gradLambda, index, dop, phi);
		
	*jump=phi[0]*nv[0]+phi[1]*nv[1];
}

/**
 * \fn void TangentDerivative4Edge(double lambda1, double lambda2, int edge, int element, ELEMENT *elements, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *val)
 * \brief the tangential derivative of basis function
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param edge the index of current edge
 * \param *edgeNode pointer to the two indices of current edge
 * \param element the index of current element
 * \param *elements pointer to the structure of the triangulation
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *val pointer to the tangential derivative of basis function on edge
 * \return void
 */
void TangentDerivative4Edge(double lambda1, double lambda2, int edge, int element, ELEMENT *elements, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *val)
{
	int i,j,l;
	double lambdas[3];
	double s, **gradLambda, phi[2], *tve;
	int nodeindex, edgeNode[2];
	int dop=elementDOF->dop;

	*val=0;

	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];

	for(nodeindex=0;nodeindex<elementDOF->col;nodeindex++)
    {
		if(elementDOF->val[element][nodeindex]==node)
			break;
	}

	if(nodeindex==elementDOF->col)
			return;
		
	for(i=0;i<3;i++)
	{
		if(elements->val[element][i]==edgeNode[0])
			break;
	}
	
	for(j=0;j<3;j++)
	{
		if(elements->val[element][j]==edgeNode[1])
			break;
	}
	
	l=3-i-j;
	lambdas[i]=lambda1;
	lambdas[j]=lambda2;
	lambdas[l]=0;
	
	s=elements->vol[element];
	gradLambda=elements->gradLambda[element];
	tve=edges->tvector[edge];
	
	lagrange_basis1(lambdas, gradLambda, nodeindex, dop, phi);
		
	*val=phi[0]*tve[0]+phi[1]*tve[1];
}

/**
 * \fn void averageOperator(double *lambdas1d, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *average)
 * \brief the average of second order normal derivative of basis function
 * \param lambdas1d the barycentric coordiante in one dimension
 * \param edge the index of current edge
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *avg pointer to the result of average
 * \return void
 */
void averageNormalDerivative2Lagrange2d(double *lambdas1d, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *avg)
{
	int i, j, k, l;
	int ni, ei;
	int element, edgeNode[2];
	int dop=elementDOF->dop;

	double lambdas[3], **gradLambda, *nve, phi[3], val;
	double lambda1 = lambdas1d[0];
	double lambda2 = lambdas1d[1];
	
	*avg =0;

	if (node<0 && node >= elementDOF->dof)
		return;

	double elen=edges->length[edge];
	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];
	nve = edges->nvector[edge];

	for (k = 0; k < 2; k++)
	{
		element = edges->val[edge][2 + k];
		if (element == -1)
			continue;

		for (ni = 0; ni<elementDOF->col; ni++)
		{
			if (elementDOF->val[element][ni] == node)
				break;
		}
		if (ni == elementDOF->col)
			continue;

		for (ei = 0; ei<3; ei++)
		{
			if (elementEdge->val[element][ei] == edge)
				break;
		}

		for (i = 0; i<3; i++)
		{
			if (elements->val[element][i] == edgeNode[0])
				break;
		}
		for (j = 0; j<3; j++)
		{
			if (elements->val[element][j] == edgeNode[1])
				break;
		}
		l = 3 - i - j;  //  l == ei

		lambdas[i] = lambda1;
		lambdas[j] = lambda2;
		lambdas[l] = 0;

		gradLambda = elements->gradLambda[element];
	
		lagrange_basis2(lambdas, gradLambda, ni, dop, phi);

		val = phi[0] * nve[0] * nve[0] + phi[1] * nve[1] * nve[1] + 2 * phi[2] * nve[0] * nve[1];
		*avg += val;
	} // k	

	if(edges->val[edge][3]!=-1)
		*avg /= 2.0;
}

/**
 * \fn void jumpNormalDerivativeLagrange2d(double *lambdas1d, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
 * \brief the jump of normal derivative of basis function
 * \param lambdas1d the barycentric coordiante in one dimension
 * \param edge the index of current edge
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *jump pointer to the result of jump
 * \return void
 */
void jumpNormalDerivativeLagrange2d(double *lambdas1d, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i, j, k, l;
	int ni, ei;
	int element, edgeNode[2];
	int dop=elementDOF->dop;

	double lambdas[3], **gradLambda, *nv, phi[3], val;
	double lambda1 = lambdas1d[0];
	double lambda2 = lambdas1d[1];
	
	*jump =0;

	if (node<0 && node >= elementDOF->dof)
		return;

	double elen=edges->length[edge];
	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];
	// nve = edges->nvector[edge];

	for (k = 0; k < 2; k++)
	{
		element = edges->val[edge][2 + k];
		if (element == -1)
			continue;

		for (ni = 0; ni<elementDOF->col; ni++)
		{
			if (elementDOF->val[element][ni] == node)
				break;
		}
		if (ni == elementDOF->col)
			continue;

		for (ei = 0; ei<3; ei++)
		{
			if (elementEdge->val[element][ei] == edge)
				break;
		}

		for (i = 0; i<3; i++)
		{
			if (elements->val[element][i] == edgeNode[0])
				break;
		}
		for (j = 0; j<3; j++)
		{
			if (elements->val[element][j] == edgeNode[1])
				break;
		}
		l = 3 - i - j; //  l == ei
		lambdas[i] = lambda1;
		lambdas[j] = lambda2;
		lambdas[l] = 0;

		gradLambda = elements->gradLambda[element];
		nv = elements->nvector[element][ei];
		
		lagrange_basis1(lambdas, gradLambda, ni, dop, phi);

		val = phi[0] * nv[0] + phi[1] * nv[1];
		*jump += val;
	} // k	
}

/**
 * \fn void jumpOperatorTensor(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
 * \brief the jump operator for tensor
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param edge the index of current edge
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *jump pointer to the result of jump operator
 * \return void
 */
void jumpOperatorTensor(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i,j,l;
	int nodeindex, edgeindex;
	int element, edgeNode[2];
	double lambdas[3];
	double *eta, *xi, nv[2];
	double phi;
	double elen=edges->length[edge];

	jump[0]=0;
	jump[1]=0;

	if(node<0 && node>=elementDOF->dof*3)
		return;

	element=node/(3*elementDOF->col);
	
	for(edgeindex=0;edgeindex<3;edgeindex++)
	{
		if(elementEdge->val[element][edgeindex]==edge)
			break;
	}

	if(edgeindex==3)
		return;

	xi=elements->xi[element];
	eta=elements->eta[element];
	nv[0]=-eta[edgeindex]/elen;
	nv[1]=xi[edgeindex]/elen;

	if(elementDOF->dop==0)
	{
		if(node%3<1)
		{
			jump[0]=nv[0];
			jump[1]=0;
		}
		else if(node%3<2)
		{
			jump[0]=0;
			jump[1]=nv[1];
		}
		else
		{
			jump[0]=nv[1];
			jump[1]=nv[0];
		}
		return;
	}

	nodeindex=node%elementDOF->col;

	if(nodeindex<3 && nodeindex>=0)
	{
		if(nodeindex==edgeindex)
			return;
	}
	else if(nodeindex<3*elementDOF->dop && nodeindex>=3)
	{
		if((nodeindex-3)/(elementDOF->dop-1) != edgeindex)
			return;
	}
	else
		return;

	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];
	for(i=0;i<3;i++)
	{
		if(elements->val[element][i]==edgeNode[0])
			break;
	}
	for(j=0;j<3;j++)
	{
		if(elements->val[element][j]==edgeNode[1])
			break;
	}
	l=3-i-j;
	lambdas[i]=lambda1;
	lambdas[j]=lambda2;
	lambdas[l]=0;

	lagrange_basis(lambdas, nodeindex, elementDOF->dop, &phi);

	int lnode=node%(3*elementDOF->col);
	if(lnode<elementDOF->col)
	{
		jump[0]=phi*nv[0];
		jump[1]=0;
	}
	else if(lnode<elementDOF->col*2)
	{
		jump[0]=0;
		jump[1]=phi*nv[1];
	}
	else
	{
		jump[0]=phi*nv[1];
		jump[1]=phi*nv[0];
	}
}

/**
 * \fn void jumpOperatorTensorNormal(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
 * \brief the jump operator for M_{nn}(\tau)
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param edge the index of current edge
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *jump pointer to the result of jump operator
 * \return void
 */
void jumpOperatorTensorNormal(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i,j,l;
	int nodeindex, edgeindex;
	int element, edgeNode[2];
	double lambdas[3];
	double *eta, *xi, nv[2], nve[2];
	double phi;
	double elen=edges->length[edge];

	*jump=0;

	if(node<0 && node>=elementDOF->dof*3)
		return;

	element=node/(3*elementDOF->col);
	
	for(edgeindex=0;edgeindex<3;edgeindex++)
	{
		if(elementEdge->val[element][edgeindex]==edge)
			break;
	}

	if(edgeindex==3)
		return;

	xi=elements->xi[element];
	eta=elements->eta[element];
	nv[0]=-eta[edgeindex]/elen;
	nv[1]=xi[edgeindex]/elen;
	nve[0]=edges->nvector[edge][0];
	nve[1]=edges->nvector[edge][1];
	
	int lnode=node%(3*elementDOF->col);
	if(elementDOF->dop==0)
	{
		if(lnode<elementDOF->col)
			*jump=nve[0]*nv[0];
		else if(lnode<elementDOF->col*2)
			*jump=nve[1]*nv[1];
		else
			*jump=nve[0]*nv[1]+nve[1]*nv[0];
	}

	nodeindex=node%elementDOF->col;

	if(nodeindex<3 && nodeindex>=0)
	{
		if(nodeindex==edgeindex)
			return;
	}
	else if(nodeindex<3*elementDOF->dop && nodeindex>=3)
	{
		if((nodeindex-3)/(elementDOF->dop-1) != edgeindex)
			return;
	}
	else
		return;

	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];
	for(i=0;i<3;i++)
	{
		if(elements->val[element][i]==edgeNode[0])
			break;
	}
	for(j=0;j<3;j++)
	{
		if(elements->val[element][j]==edgeNode[1])
			break;
	}
	l=3-i-j;
	lambdas[i]=lambda1;
	lambdas[j]=lambda2;
	lambdas[l]=0;

	lagrange_basis(lambdas, nodeindex, elementDOF->dop, &phi);

	if(lnode<elementDOF->col)
		*jump=phi*nve[0]*nv[0];
	else if(lnode<elementDOF->col*2)
		*jump=phi*nve[1]*nv[1];
	else
		*jump=phi*nve[0]*nv[1]+phi*nve[1]*nv[0];
}

/**
 * \fn void jumpOperatorTensoTangent (double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
 * \brief the jump operator for M_{nt}(\tau)
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param edge the index of current edge
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *jump pointer to the result of jump operator
 * \return void
 */
void jumpOperatorTensorTangent(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i,j,l;
	int nodeindex, edgeindex;
	int element, edgeNode[2];
	double lambdas[3];
	double *eta, *xi, nv[2], tve[2];
	double phi;
	double elen=edges->length[edge];

	*jump=0;

	if(node<0 && node>=elementDOF->dof*3)
		return;

	element=node/(3*elementDOF->col);
	
	for(edgeindex=0;edgeindex<3;edgeindex++)
	{
		if(elementEdge->val[element][edgeindex]==edge)
			break;
	}

	if(edgeindex==3)
		return;

	xi=elements->xi[element];
	eta=elements->eta[element];
	nv[0]=-eta[edgeindex]/elen;
	nv[1]=xi[edgeindex]/elen;
	tve[0]=edges->tvector[edge][0];
	tve[1]=edges->tvector[edge][1];

	int lnode=node%(3*elementDOF->col);
	if(elementDOF->dop==0)
	{
		if(lnode<elementDOF->col)
			*jump=tve[0]*nv[0];
		else if(lnode<elementDOF->col*2)
			*jump=tve[1]*nv[1];
		else
			*jump=tve[0]*nv[1]+tve[1]*nv[0];
		return;
	}

	nodeindex=node%elementDOF->col;

	if(nodeindex<3 && nodeindex>=0)
	{
		if(nodeindex==edgeindex)
			return;
	}
	else if(nodeindex<3*elementDOF->dop && nodeindex>=3)
	{
		if((nodeindex-3)/(elementDOF->dop-1) != edgeindex)
			return;
	}
	else
		return;

	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];
	for(i=0;i<3;i++)
	{
		if(elements->val[element][i]==edgeNode[0])
			break;
	}
	for(j=0;j<3;j++)
	{
		if(elements->val[element][j]==edgeNode[1])
			break;
	}
	l=3-i-j;
	lambdas[i]=lambda1;
	lambdas[j]=lambda2;
	lambdas[l]=0;

	lagrange_basis(lambdas, nodeindex, elementDOF->dop, &phi);

	if(lnode<elementDOF->col)
		*jump=phi*tve[0]*nv[0];
	else if(lnode<elementDOF->col*2)
		*jump=phi*tve[1]*nv[1];
	else
		*jump=phi*tve[0]*nv[1]+phi*tve[1]*nv[0];
}

/**
 * \fn void jumpOperatorTensorDIVPt(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
 * \brief the jump operator for divergence of tensor and tangential derivative of M_{nt}(\tau)
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param edge the index of current edge
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *jump pointer to the result of jump operator
 * \return void
 */
void jumpOperatorTensorDIVPt(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i,j,l;
	int nodeindex, edgeindex;
	int element, edgeNode[2];
	double lambdas[3];
	double s, **gradLambda, *nv, *tv;
	double phi[3];
	double elen=edges->length[edge];

	*jump=0;

	if(elementDOF->dop<1)
		return;
	if(node<0 && node>=elementDOF->dof*3)
		return;

	element=node/(3*elementDOF->col);
	
	for(edgeindex=0;edgeindex<3;edgeindex++)
	{
		if(elementEdge->val[element][edgeindex]==edge)
			break;
	}

	if(edgeindex==3)
		return;

	nodeindex=node%elementDOF->col;

	s=elements->vol[element];
	gradLambda=elements->gradLambda[element];
	nv=elements->nvector[element][edgeindex];
	tv=elements->tvector[element][edgeindex];
	
	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];
	for(i=0;i<3;i++)
	{
		if(elements->val[element][i]==edgeNode[0])
			break;
	}
	for(j=0;j<3;j++)
	{
		if(elements->val[element][j]==edgeNode[1])
			break;
	}
	l=3-i-j;
	lambdas[i]=lambda1;
	lambdas[j]=lambda2;
	lambdas[l]=0;

	lagrange_basis1(lambdas, gradLambda, nodeindex, elementDOF->dop, phi);

	int lnode=node%(3*elementDOF->col);
	if(lnode<elementDOF->col)
		*jump+=((1+tv[0]*tv[0])*phi[0] + tv[0]*tv[1]*phi[1])*nv[0];
	else if(lnode<elementDOF->col*2)
		*jump+=(tv[0]*tv[1]*phi[0] + (1+tv[1]*tv[1])*phi[1])*nv[1];
	else
		*jump+=((tv[0]*tv[1]*phi[0] + (1+tv[1]*tv[1])*phi[1])*nv[0] + ((1+tv[0]*tv[0])*phi[0] + tv[0]*tv[1]*phi[1])*nv[1]);

	/*if(lnode<elementDOF->col)
		*jump+=phi[0]*nv[0];
	else if(lnode<elementDOF->col*2)
		*jump+=phi[1]*nv[1];
	else
		*jump+=phi[0]*nv[1]+phi[1]*nv[0];

	int Pt=0;
	if(nodeindex<3 && nodeindex>=0)
	{
		if(nodeindex!=edgeindex)
			Pt=1;
	}
	else if(nodeindex<3*elementDOF->dop && nodeindex>=3)
	{
		if((nodeindex-3)/(elementDOF->dop-1) == edgeindex)
			Pt=1;
	}

	if(Pt==1)
	{
		phi[2]=phi[0]*tv[0]+phi[1]*tv[1];
		if(lnode<elementDOF->col)
			*jump+=phi[2]*tv[0]*nv[0];
		else if(lnode<elementDOF->col*2)
			*jump+=phi[2]*tv[1]*nv[1];
		else
			*jump+=phi[2]*(tv[0]*nv[1]+tv[1]*nv[0]);
	}*/
}

 /**
 * \fn void jumpOperatorTensorDIV(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
 * \brief the jump operator for divergence of tensor
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param edge the index of current edge
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *jump pointer to the result of jump operator
 * \return void
 */
void jumpOperatorTensorDIV(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i,j,l;
	int nodeindex, edgeindex;
	int element, edgeNode[2];
	double lambdas[3];
	double s, **gradLambda, *nv, *tv;
	double phi[3];
	double elen=edges->length[edge];

	*jump=0;

	if(elementDOF->dop<1)
		return;
	if(node<0 && node>=elementDOF->dof*3)
		return;

	element=node/(3*elementDOF->col);
	
	for(edgeindex=0;edgeindex<3;edgeindex++)
	{
		if(elementEdge->val[element][edgeindex]==edge)
			break;
	}

	if(edgeindex==3)
		return;

	nodeindex=node%elementDOF->col;

	s=elements->vol[element];
	gradLambda=elements->gradLambda[element];
	nv=elements->nvector[element][edgeindex];
	
	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];
	for(i=0;i<3;i++)
	{
		if(elements->val[element][i]==edgeNode[0])
			break;
	}
	for(j=0;j<3;j++)
	{
		if(elements->val[element][j]==edgeNode[1])
			break;
	}
	l=3-i-j;
	lambdas[i]=lambda1;
	lambdas[j]=lambda2;
	lambdas[l]=0;

	lagrange_basis1(lambdas, gradLambda, nodeindex, elementDOF->dop, phi);

	int lnode=node%(3*elementDOF->col);
	if(lnode<elementDOF->col)
		*jump+=phi[0]*nv[0];
	else if(lnode<elementDOF->col*2)
		*jump+=phi[1]*nv[1];
	else
		*jump+=phi[0]*nv[1]+phi[1]*nv[0];
}

/**
* \fn void jumpOperatorATensorTangent2(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double lambda, double mu, double *jump)
* \brief the jump of M_{tt}(Asigmah) for Hu-Zhang element
* \param lambda1 the first length coordiante
* \param lambda2 the second length coordiante
* \param edge the index of current edge
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *elementDOF pointer to relation between elements and DOFs
* \param node index of current node variable
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \param *jump pointer to the result of jump operator
* \return void
*/
void jumpOperatorATensorTangent2(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double lambda, double mu, double *jump)
{
	int i, j, k, l;
	int nodeindex, edgeindex;
	int element, edgeNode[2];
	double lambdas[3];
	double *eta, *xi, *nv, *nve, *tv;
	double phi[3], sgn, value[2];
	
	*jump = 0;
	
	if (node<0 && node >= elementDOF->dof)
		return;

	double elen = edges->length[edge];
	edgeNode[0] = edges->val[edge][0];
	edgeNode[1] = edges->val[edge][1];
	nve = edges->nvector[edge];

	for (k = 0; k < 2; k++)
	{
		element = edges->val[edge][2 + k];
		if (element == -1)
			continue;

		for (nodeindex = 0; nodeindex<elementDOF->col; nodeindex++)
		{
			if (elementDOF->val[element][nodeindex] == node)
				break;
		}
		if (nodeindex == elementDOF->col)
			continue;

		for (edgeindex = 0; edgeindex<3; edgeindex++)
		{
			if (elementEdge->val[element][edgeindex] == edge)
				break;
		}

		for (i = 0; i<3; i++)
		{
			if (elements->val[element][i] == edgeNode[0])
				break;
		}
		for (j = 0; j<3; j++)
		{
			if (elements->val[element][j] == edgeNode[1])
				break;
		}
		l = 3 - i - j;
		lambdas[i] = lambda1;
		lambdas[j] = lambda2;
		lambdas[l] = 0;

		double sgn = elements->eorien[element][edgeindex];
		nv=elements->nvector[element][edgeindex];
		tv=elements->tvector[element][edgeindex];
	

		huzhang_basis(lambdas, elements->nvector[element], elements->tvector[element], nodeindex, elementDOF->dop, phi);

		value[0] = phi[0] * tv[0] * tv[0] + phi[1] * tv[1] * tv[1] + 2 * phi[2] * tv[0] * tv[1];
		value[1] = phi[0] + phi[1];

		*jump += (value[0] - value[1] * lambda / (2 * lambda + 2 * mu)) / (2 * mu)*sgn;
	} // k	
}

/**
* \fn void jumpOperatorRotATensorTangentPt(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double lambda, double mu, double *jump)
* \brief the jump of rot(Asigmah) t - \partial_t(M_{nt}(Asigmah)) for Hu-Zhang element
* \param lambda1 the first length coordiante
* \param lambda2 the second length coordiante
* \param edge the index of current edge
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *elementDOF pointer to relation between elements and DOFs
* \param node index of current node variable
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \param *jump pointer to the result of jump operator
* \return void
*/
void jumpOperatorRotATensorTangentPt(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double lambda, double mu, double *jump)
{
	int i, j, k, l;
	int nodeindex, edgeindex;
	int element, edgeNode[2];
	double lambdas[3];
	double s, **gradLambda, *nv, nve[2], *tv;
	double phi1[2], phi2[2], phi3[3][2], sgn, value[3];

	*jump = 0;

	if (node<0 && node >= elementDOF->dof)
		return;

	double elen = edges->length[edge];
	edgeNode[0] = edges->val[edge][0];
	edgeNode[1] = edges->val[edge][1];

	for (k = 0; k < 2; k++)
	{
		element = edges->val[edge][2 + k];
		if (element == -1)
			continue;

		for (nodeindex = 0; nodeindex<elementDOF->col; nodeindex++)
		{
			if (elementDOF->val[element][nodeindex] == node)
				break;
		}
		if (nodeindex == elementDOF->col)
			continue;

		for (edgeindex = 0; edgeindex<3; edgeindex++)
		{
			if (elementEdge->val[element][edgeindex] == edge)
				break;
		}

		for (i = 0; i<3; i++)
		{
			if (elements->val[element][i] == edgeNode[0])
				break;
		}
		for (j = 0; j<3; j++)
		{
			if (elements->val[element][j] == edgeNode[1])
				break;
		}
		l = 3 - i - j;
		lambdas[i] = lambda1;
		lambdas[j] = lambda2;
		lambdas[l] = 0;

		s = elements->vol[element];
		gradLambda = elements->gradLambda[element];
		nv = elements->nvector[element][edgeindex];
		tv = elements->tvector[element][edgeindex];
//		sgn = nve[0] * nv[0] + nve[1] * nv[1];

		huzhang_basisROT(lambdas, gradLambda, elements->nvector[element], elements->tvector[element], nodeindex, elementDOF->dop, phi1);
		huzhang_basisCurlTrace(lambdas, gradLambda, elements->nvector[element], elements->tvector[element], nodeindex, elementDOF->dop, phi2);

		value[0] = phi1[0] * tv[0] + phi1[1] * tv[1];
		value[1] = phi2[0] * tv[0] + phi2[1] * tv[1];

		*jump += (value[0] - value[1] * lambda / (2 * lambda + 2 * mu)) / (2 * mu);

		huzhang_basis1(lambdas, gradLambda, elements->nvector[element], elements->tvector[element], nodeindex, elementDOF->dop, phi3);
		value[0] = phi3[0][0] * tv[0] + phi3[0][1] * tv[1];
		value[1] = phi3[1][0] * tv[0] + phi3[1][1] * tv[1];
		value[2] = phi3[2][0] * tv[0] + phi3[2][1] * tv[1];

		*jump -= (value[0] * nv[0] * tv[0] + value[1] * nv[1] * tv[1] + value[2] * (nv[0] * tv[1] + nv[1] * tv[0])) / (2 * mu);
	} // k	
}

/**
 * \fn void getElementEdgeGeoInfo(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes)
 * \brief compute the geometric information of elements and edges
 * \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
 * \return void
 */
void getElementEdgeGeoInfo(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes)
{
	int i,j;
	int vertex, edge;
	int element1, element2;
	double tri[3][2], xi[3], eta[3];
/*	if(!elements->h)
		free(elements->h);
	if(!edges->h)
		free(edges->h);
	elements->h=(double*)calloc(elements->row, sizeof(double));
	edges->h=(double*)calloc(edges->row, sizeof(double)); */
	for(i=0;i<elements->row;i++)
	{
		for(j=0;j<3;j++)
		{
			vertex=elements->val[i][j];
			copy_array(2, nodes->val[vertex], tri[j]);
			copy_array(2, nodes->val[vertex], elements->vertices[i][j]); 
		}

		elements->vol[i] = area(elements->vertices[i]);

		for(j=0;j<2;j++)
			elements->barycenter[i][j] = (tri[0][j]+tri[1][j]+tri[2][j]) / 3.0;

		// elements->xi[i][0] = nodes->val[vertex2][0]-nodes->val[vertex3][0];
		// elements->xi[i][1] = nodes->val[vertex3][0]-nodes->val[vertex1][0];
		// elements->xi[i][2] = nodes->val[vertex1][0]-nodes->val[vertex2][0];
		// elements->eta[i][0] = nodes->val[vertex2][1]-nodes->val[vertex3][1];
		// elements->eta[i][1] = nodes->val[vertex3][1]-nodes->val[vertex1][1];
		// elements->eta[i][2] = nodes->val[vertex1][1]-nodes->val[vertex2][1];

		xi[0] = tri[1][0] - tri[2][0];
		xi[1] = tri[2][0] - tri[0][0];
		xi[2] = tri[0][0] - tri[1][0];
		eta[0] = tri[1][1] - tri[2][1];
		eta[1] = tri[2][1] - tri[0][1];
		eta[2] = tri[0][1] - tri[1][1];

		elements->xi[i][0] = xi[0];
		elements->xi[i][1] = xi[1];
		elements->xi[i][2] = xi[2];
		elements->eta[i][0] = eta[0];
		elements->eta[i][1] = eta[1];
		elements->eta[i][2] = eta[2];

		axpyz_array(2, -1.0, tri[1], tri[2], elements->tij[i][0]);
		axpyz_array(2, -1.0, tri[2], tri[0], elements->tij[i][1]);
		axpyz_array(2, -1.0, tri[0], tri[1], elements->tij[i][2]);
		
		for(j=0;j<3;j++){
			elements->edgeslength[i][j] = sqrt(xi[j]*xi[j]+eta[j]*eta[j]);
			elements->nvector[i][j][0] = -eta[j]/elements->edgeslength[i][j];
		    elements->nvector[i][j][1] = xi[j]/elements->edgeslength[i][j];
			elements->tvector[i][j][0] = -elements->nvector[i][j][1];
			elements->tvector[i][j][1] = elements->nvector[i][j][0];
			elements->gradLambda[i][j][0] = eta[j]/(2*elements->vol[i]);
			elements->gradLambda[i][j][1] = -xi[j]/(2*elements->vol[i]);
		}

		// for(j=0;j<3;j++){
		// 	printf("(%le, %le), (%le, %le), (%le, %le)\n", elements->tij[i][j][0]-elements->tvector[i][j][0]*elements->edgeslength[i][j], elements->tij[i][j][1]-elements->tvector[i][j][1]*elements->edgeslength[i][j], elements->tij[i][j][0], elements->tij[i][j][1], elements->tvector[i][j][0]*elements->edgeslength[i][j], elements->tvector[i][j][1]*elements->edgeslength[i][j]);
		// }

		for(j=0;j<3;j++){
			edge = elementEdge->val[i][j];
			if(elements->val[i][(j+1)%3]== edges->val[edge][0]){
				elements->eperm[i][j][0] = 0;
				elements->eperm[i][j][1] = 1;
				elements->eorien[i][j] = 1;
			}
			else{
				elements->eperm[i][j][0] = 1;
				elements->eperm[i][j][1] = 0;
				elements->eorien[i][j] = -1;
			}
		}
	}

	for(i=0;i<edges->row;i++)
	{
		for(j=0;j<2;j++)
		{
			vertex=edges->val[i][j];
			copy_array(2, nodes->val[vertex], tri[j]);
			// for(k=0;k<3;k++)
			// 	tet[j][k]=nodes->val[vertex][k];
		}
		// vertex1=edges->val[i][0];
		// vertex2=edges->val[i][1];
		// edges->xi[i]=nodes->val[vertex1][0]-nodes->val[vertex2][0];
		// edges->eta[i]=nodes->val[vertex1][1]-nodes->val[vertex2][1];
		xi[0] = tri[0][0] - tri[1][0];
		eta[0] = tri[0][1] - tri[1][1];
		edges->length[i]=sqrt(xi[0]*xi[0]+eta[0]*eta[0]);
		edges->nvector[i][0]=-eta[0]/edges->length[i];
	    edges->nvector[i][1]=xi[0]/edges->length[i];
		edges->tvector[i][0]=-edges->nvector[i][1];
	    edges->tvector[i][1]=edges->nvector[i][0];

		element1=edges->val[i][2];
		element2=edges->val[i][3];
		if(element2>-1)
			edges->h[i]=sqrt((elements->vol[element1]+elements->vol[element2])/2.0);
		else
			edges->h[i]=sqrt(elements->vol[element1]);
	}
}

/**
 * \fn void getBoundaryInfo(EDGE *edges, dennode *nodes, int dof, int dop, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index)
 * \brief get information of triangulation on the boundary
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
 * \param dof number of degrees of freedom
 * \param dop number of degrees of polynomial
 * \param *isInNode pointer to boundary information of nodes: if the node is interior node, it will be 0; if the node is on the boundary, it will be -1
 * \param *dirichlet pointer to the index of dirichlet nodes
 * \param *nondirichlet pointer to the index of nondirichlet nodes
 * \param *index pointer to the transpose of dirichlet and nondirichlet nodes
 * \return void
 */
void getBoundaryInfo(EDGE *edges, dennode *nodes, int dof, int dop, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index)
{
	int nn=nodes->row;
	int i,j,k,dr;

	create_ivector(dof, isInNode);
	create_ivector(dof, index);
	
	for(i=0;i<nn;i++)
		isInNode->val[i]=1-nodes->bdFlag[i]-1;
		// isInNode->val[i]=nodes->isInNode[i];
	for(j=0;j<edges->row;j++)
	{
		if(edges->val[j][3]==-1)
		{
			for(i=0;i<dop-1;i++)
				isInNode->val[nn+j*(dop-1)+i]=-1;
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
* \fn void getBoundaryInfoVector2d(EDGE *edges, dennode *nodes, int dof, int dop, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index)
* \brief get information of triangulation on the boundary in 2d-vector version
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param dof number of degrees of freedom
* \param dop number of degrees of polynomial
* \param *isInNode pointer to boundary information of nodes: if the node is interior node, it will be 0; if the node is on the boundary, it will be -1
* \param *dirichlet pointer to the index of dirichlet nodes
* \param *nondirichlet pointer to the index of nondirichlet nodes
* \param *index pointer to the transpose of dirichlet and nondirichlet nodes
* \return void
*/
void getBoundaryInfoVector2d(EDGE *edges, dennode *nodes, int dof, int dop, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index)
{
	int nn = nodes->row;
	int i, j, k, dr;

	create_ivector(dof * 2, isInNode);
	create_ivector(dof * 2, index);

	// for (i = 0; i<nn; i++)
	// {
	// 	isInNode->val[i] = nodes->isInNode[i];
	// 	isInNode->val[i + dof] = nodes->isInNode[i];
	// }
	for (j = 0; j<edges->row; j++)
	{
		if (edges->val[j][3] == -1)
		{
			for (i = 0; i<dop - 1; i++)
			{
				isInNode->val[nn + j*(dop - 1) + i] = -1;
				isInNode->val[nn + j*(dop - 1) + i + dof] = -1;
			}
		}
	}

	dr = 0;
	for (i = 0; i<isInNode->row; i++)
	{
		if (isInNode->val[i] == -1)
			dr++;
	}
	create_ivector(dr, dirichlet);
	create_ivector(isInNode->row - dr, nondirichlet);

	j = 0; k = 0;
	for (i = 0; i<isInNode->row; i++)
	{
		if (isInNode->val[i] == -1) //  Dirichlet boundary node
		{
			dirichlet->val[k] = i;
			index->val[i] = k;
			k++;
		}
		else // free variable
		{
			nondirichlet->val[j] = i;
			index->val[i] = j;
			j++;
		}
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
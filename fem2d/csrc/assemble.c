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
 * \fn int getmesh(ELEMENT *ptr_elements, idenmat *ptr_elementEdge, EDGE *ptr_edges, dennode *ptr_nodes, iCSRmat *edgesTran, ivector *nodeCEdge, int levelNum)
 * \brief generate mesh information
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
int getmesh(ELEMENT *ptr_elements, idenmat *ptr_elementEdge, EDGE *ptr_edges, dennode *ptr_nodes, iCSRmat *edgesTran, ivector *nodeCEdge, int levelNum)
{
	ddenmat nodes; // the first column stores the x coordinate of points, the second column stores the y coordinate of points
	idenmat elements[levelNum]; // triangulation: store 3 points corresponding to the element in each row
	idenmat edges[levelNum]; // the first two columns store the two vertice, the third column stores the affiliated element
	iCSRmat elementsTran[levelNum]; // the transpose of elements. 
	ivector isInNode; // if the node is interior node, it will be 0; if the node is on the boundary, it will be -1
	
	int IsExist=getCoarseInfo(&nodes, &elements[0], &edges[0], &elementsTran[0], &edgesTran[0], &isInNode, nodeCEdge);
	if(IsExist==0)
	{
		printf("Constructing coarse grid fails!\n");
		return 0;
	}
	int i,j,k,l;
	int nvertices;

	for(l=0;l<levelNum-1;l++)
	{
		nvertices = nodes.row;
		create_dennode(nvertices, 2, ptr_nodes+l);
		for (i = 0; i < nvertices; i++)
		{
			ptr_nodes[l].val[i][0] = nodes.val[i][0];
			ptr_nodes[l].val[i][1] = nodes.val[i][1];
			ptr_nodes[l].isInNode[i] = isInNode.val[i];
		}
		
		refine(&nodes, &elements[l], &edges[l], &elementsTran[l], &elements[l+1], &edges[l+1], &elementsTran[l+1], &edgesTran[l+1], &isInNode, nodeCEdge);
	}
	
	nvertices = nodes.row;
	create_dennode(nvertices, 2, ptr_nodes + l);
	for (i = 0; i < nvertices; i++)
	{
		ptr_nodes[l].val[i][0] = nodes.val[i][0];
		ptr_nodes[l].val[i][1] = nodes.val[i][1];
		ptr_nodes[l].isInNode[i] = isInNode.val[i];
	}
	
	// generate ptr_elements, ptr_edges, ptr_elementEdge	
	int point1, point2, element1, element2;
	
	for (l = 0; l < levelNum; l++)
	{
		create_ELEMENT(elements[l].row, 3, ptr_elements + l);
		create_EDGE(edges[l].row, edges[l].col, ptr_edges + l);
		create_iden_matrix(ptr_elements[l].row, 3, ptr_elementEdge + l);

		for (i = 0; i<ptr_elements[l].row; i++)
		{
			ptr_elements[l].val[i][0] = elements[l].val[i][0];
			ptr_elements[l].val[i][1] = elements[l].val[i][1];
			ptr_elements[l].val[i][2] = elements[l].val[i][2];
		}

		for (i = 0; i<ptr_edges[l].row; i++)
		{
			for (j = 0; j<ptr_edges[l].col; j++)
				ptr_edges[l].val[i][j] = edges[l].val[i][j];
		}

		for (i = 0; i<edges[l].row; i++)
		{
			point1 = edges[l].val[i][0];
			point2 = edges[l].val[i][1];
			element1 = edges[l].val[i][2];
			element2 = edges[l].val[i][3];

			for (j = 0; j<3; j++)
			{
				if (ptr_elements[l].val[element1][j] == point1)
					break;
			}
			for (k = 0; k<3; k++)
			{
				if (ptr_elements[l].val[element1][k] == point2)
					break;
			}
			ptr_elementEdge[l].val[element1][3 - j - k] = i;

			if (element2>-1)
			{
				for (j = 0; j<3; j++)
				{
					if (ptr_elements[l].val[element2][j] == point1)
						break;
				}
				for (k = 0; k<3; k++)
				{
					if (ptr_elements[l].val[element2][k] == point2)
						break;
				}
				ptr_elementEdge[l].val[element2][3 - j - k] = i;
			}
		}
	}
	// end generate ptr_elements, ptr_edges, ptr_elementEdge
	
	
	for(i=0;i<nodes.row;i++)
		free(nodes.val[i]);
	free(nodes.val);
	
	for(l=0;l<levelNum;l++)
	{
		for(i=0;i<elements[l].row;i++)
			free(elements[l].val[i]);
		free(elements[l].val);
		
		for(i=0;i<edges[l].row;i++)
			free(edges[l].val[i]);
		free(edges[l].val);
		
		free(elementsTran[l].IA);
		free(elementsTran[l].JA);
	}
	
	free(isInNode.val);
	
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
 * \fn void getElementDOF_Lagrange(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, int nedges, int nvertices, int dop)
 * \brief get the degrees of freedom of Lagrange element
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param nvertices number of vertices
 * \param dop degree of polynomial
 */
void getElementDOF_Lagrange(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
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
* \fn void assemble(dCSRmat *ptr_A, dvector *ptr_b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu)
* \brief assemble stiffness matrix *C and righ hand side *ptr_b
* \param *ptr_A pointer to stiffness matrix
* \param *ptr_b pointer to right hand side
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
void assemble(dCSRmat *ptr_A, dvector *ptr_b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu)
{
	/**	A[0]=A, A[1]=B^T, A{2]=-B
	Ax1 + B^Tx2 = 0
	-Bx1 +       = b
	where x1: sigma_h, x2: u_h
	**/

	if (elementDOF[0].dop > 2)
	{
		assembleStiffmatrixHuZhang2d(&ptr_A[0], &ptr_A[1], &ptr_b[1], elements, elementEdge, edges, nodes, elementDOF, elementdofTran, lambda, mu);
		ptr_A[3].row = 0;
		ptr_A[3].col = 0;
		ptr_A[3].IA = NULL;
		ptr_A[3].JA = NULL;
		ptr_A[3].val = NULL;
	}
	else
		assembleStiffmatrixHuZhangIP2d(&ptr_A[0], &ptr_A[1], &ptr_A[3], &ptr_b[1], elements, elementEdge, edges, nodes, elementDOF, elementdofTran, lambda, mu);


	getTransposeOfSparse(&ptr_A[1], &ptr_A[2]);

	int i;
//	for (i = 0; i<ptr_A[2].nnz; i++)
//		ptr_A[2].val[i] *= -1.0;
	if(ptr_A[3].row>0)
	{
		for (i = 0; i<ptr_A[3].nnz; i++)
			ptr_A[3].val[i] *= -1.0;
	}

	for (i = 0; i<ptr_b[1].row; i++)
		ptr_b[1].val[i] *= -1.0;

	create_dvector(ptr_A[0].row, &ptr_b[0]);

//	free_csr_matrix(&ptr_A[0]);
//	assembleStiffmatrixHuZhangA11_2d(&ptr_A[0], elements, elementDOF, elementdofTran, lambda, mu);
}

/**
* \fn void assembleStiffmatrixHuZhang2d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu)
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
void assembleStiffmatrixHuZhang2d(dCSRmat *A, dCSRmat *BT, dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu)
{
	int i, j, k, l;

	A->row = elementDOF[0].dof;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	BT->row = elementDOF[0].dof;
	BT->col = elementDOF[1].dof * 2;
	BT->IA = (int*)calloc(BT->row + 1, sizeof(int));
	BT->JA = NULL;
	BT->val = NULL;

	create_dvector(BT->col, b);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart;
	int count;

	int num_qp21 = getNumQuadPoints(elementDOF[0].dop * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(elementDOF[0].dop + elementDOF[1].dop - 1, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

	int num_qp11 = getNumQuadPoints(elementDOF[1].dop * 2, 1); // the number of numerical intergation points
	double gauss11[num_qp11][2];
	init_Gauss1D(num_qp11, 1, gauss11); // gauss intergation initial

										/************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)calloc(A->col, sizeof(int));
	for (i = 0; i<A->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count;

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
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{

				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		for (j = A->IA[i]; j<A->IA[i + 1]; j++)
		{
			A->JA[j] = istart;
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
		//	xi=elements->xi[k];
		//	eta=elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = A->IA[i]; j1<A->IA[i + 1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							huzhang_basis(lambdas, elements->nvector[k], elements->tvector[k], k1, elementDOF[0].dop, phi1);
							huzhang_basis(lambdas, elements->nvector[k], elements->tvector[k], k2, elementDOF[0].dop, phi2);
							if (lambda>-0.5)
								A->val[j1] += 2 * s*gauss21[i1][2] * ((phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) - lambda / (2 * lambda + 2 * mu)*(phi1[0] + phi1[1])*(phi2[0] + phi2[1])) / (2 * mu);
							else
								A->val[j1] += 2 * s*gauss21[i1][2] * ((phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) - 1.0 / 2.0*(phi1[0] + phi1[1])*(phi2[0] + phi2[1])) / (2 * mu);
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k

	  /************************************************** stiffness matrix BT *****************************************************************/
	  // step 1BT: Find first the structure IA of the stiffness matrix BT
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		BT->IA[i + 1] += 2 * elementDOF[1].col*(elementdofTran->IA[i + 1] - elementdofTran->IA[i]);
	} // i

	for (i = 0; i<BT->row; i++)
		BT->IA[i + 1] += BT->IA[i];

	BT->nnz = BT->IA[BT->row];

	// step 2BT: Find the structure JA of the stiffness matrix BT
	BT->JA = (int*)calloc(BT->nnz, sizeof(int));
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		rowstart[0] = BT->IA[i];
		row21[0] = (BT->IA[i + 1] - BT->IA[i]) / 2;
		count = 0;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];
			for (i1 = 0; i1<elementDOF[1].col; i1++)
			{
				BT->JA[rowstart[0] + count] = elementDOF[1].val[element[0]][i1];
				BT->JA[rowstart[0] + count + row21[0]] = elementDOF[1].val[element[0]][i1] + elementDOF[1].dof;
				count++;
			}
		} // j
	} // i

	  // step 3BT: Loop element by element and compute the actual entries storing them in BT
	BT->val = (double*)calloc(BT->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			row21[0] = (BT->IA[i + 1] - BT->IA[i]) / 2;
			for (k2 = 0; k2<elementDOF[1].col; k2++)
			{
				j = elementDOF[1].val[k][k2];
				// b11
				for (j1 = BT->IA[i]; j1<BT->IA[i] + row21[0]; j1++)
				{
					if (BT->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							huzhang_basisDIV(lambdas, s, eta, xi, elements->nvector[k], elements->tvector[k], k1, elementDOF[0].dop, phi1);
							lagrange_basis(lambdas, k2, elementDOF[1].dop, &phi);
							BT->val[j1] += 2 * s*gauss22[i1][2] * phi1[0] * phi;
						}
						break;
					}
					if (j1 == (BT->IA[i] + row21[0]))
					{
						printf("There is something wrong in constructing b11 of BT\n");
						exit(1);
					}
				}
				// b12
				for (j1 = BT->IA[i] + row21[0]; j1<BT->IA[i + 1]; j1++)
				{
					if (BT->JA[j1] == (j + elementDOF[1].dof))
					{
						for (i1 = 0; i1<num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							huzhang_basisDIV(lambdas, s, eta, xi, elements->nvector[k], elements->tvector[k], k1, elementDOF[0].dop, phi1);
							lagrange_basis(lambdas, k2, elementDOF[1].dop, &phi);
							BT->val[j1] += 2 * s*gauss22[i1][2] * phi1[1] * phi;
						}
						break;
					}
					if (j1 == (BT->IA[i + 1]))
					{
						printf("There is something wrong in constructing b12 of BT\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k


	  /************************************************** right hand side b *****************************************************************/
	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		s = elements->vol[k];

		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				lagrange_basis(lambdas, k1, elementDOF[1].dop, &phi);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b->val[i] += 2 * s*gauss0[i1][2] * f1(x, y, lambda, mu)*phi;
				b->val[i + elementDOF[1].dof] += 2 * s*gauss0[i1][2] * f2(x, y, lambda, mu)*phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assembleStiffmatrixHuZhang2d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu)
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
void assembleStiffmatrixHuZhangA11_2d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu)
{
	int i, j, k, l;

	A->row = elementDOF[0].dof;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart;
	int count;

	int num_qp21 = getNumQuadPoints(elementDOF[0].dop * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

										/************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)calloc(A->col, sizeof(int));
	for (i = 0; i<A->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count;

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
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{

				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		for (j = A->IA[i]; j<A->IA[i + 1]; j++)
		{
			A->JA[j] = istart;
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
		//	xi=elements->xi[k];
		//	eta=elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = A->IA[i]; j1<A->IA[i + 1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							huzhang_basis(lambdas, elements->nvector[k], elements->tvector[k], k1, elementDOF[0].dop, phi1);
							huzhang_basis(lambdas, elements->nvector[k], elements->tvector[k], k2, elementDOF[0].dop, phi2);
							if (lambda>-0.5)
								A->val[j1] += 2 * s*gauss21[i1][2] * ((phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) - lambda / (2 * lambda + 2 * mu)*(phi1[0] + phi1[1])*(phi2[0] + phi2[1])) / (2 * mu);
							else
								A->val[j1] += 2 * s*gauss21[i1][2] * ((phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) - 1.0 / 2.0*(phi1[0] + phi1[1])*(phi2[0] + phi2[1])) / (2 * mu);
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k
}

/**
* \fn void assembleStiffmatrixHuZhangIP2d(dCSRmat *A, dCSRmat *BT, dCSRmat *C, dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu)
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
void assembleStiffmatrixHuZhangIP2d(dCSRmat *A, dCSRmat *BT, dCSRmat *C, dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu)
{
	int i, j, k, l;

	A->row = elementDOF[0].dof;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	BT->row = elementDOF[0].dof;
	BT->col = elementDOF[1].dof * 2;
	BT->IA = (int*)calloc(BT->row + 1, sizeof(int));
	BT->JA = NULL;
	BT->val = NULL;

	C->row = elementDOF[1].dof * 2;
	C->col = C->row;
	C->IA = (int*)calloc(C->row + 1, sizeof(int));
	C->JA = NULL;
	C->val = NULL;

	create_dvector(C->row, b);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart;
	int count;

	int num_qp21 = getNumQuadPoints(elementDOF[0].dop * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(elementDOF[0].dop + elementDOF[1].dop - 1, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

	int num_qp11 = getNumQuadPoints(elementDOF[1].dop * 2, 1); // the number of numerical intergation points
	double gauss11[num_qp11][2];
	init_Gauss1D(num_qp11, 1, gauss11); // gauss intergation initial

										/************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)calloc(A->col, sizeof(int));
	for (i = 0; i<A->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count;

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
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{

				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		for (j = A->IA[i]; j<A->IA[i + 1]; j++)
		{
			A->JA[j] = istart;
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
		//	xi=elements->xi[k];
		//	eta=elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = A->IA[i]; j1<A->IA[i + 1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							huzhang_basis(lambdas, elements->nvector[k], elements->tvector[k], k1, elementDOF[0].dop, phi1);
							huzhang_basis(lambdas, elements->nvector[k], elements->tvector[k], k2, elementDOF[0].dop, phi2);
							if (lambda>-0.5)
								A->val[j1] += 2 * s*gauss21[i1][2] * ((phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) - lambda / (2 * lambda + 2 * mu)*(phi1[0] + phi1[1])*(phi2[0] + phi2[1])) / (2 * mu);
							else
								A->val[j1] += 2 * s*gauss21[i1][2] * ((phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) - 1.0 / 2.0*(phi1[0] + phi1[1])*(phi2[0] + phi2[1])) / (2 * mu);
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k

	  /************************************************** stiffness matrix BT *****************************************************************/
	  // step 1BT: Find first the structure IA of the stiffness matrix BT
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		BT->IA[i + 1] += 2 * elementDOF[1].col*(elementdofTran->IA[i + 1] - elementdofTran->IA[i]);
	} // i

	for (i = 0; i<BT->row; i++)
		BT->IA[i + 1] += BT->IA[i];

	BT->nnz = BT->IA[BT->row];

	// step 2BT: Find the structure JA of the stiffness matrix BT
	BT->JA = (int*)calloc(BT->nnz, sizeof(int));
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		rowstart[0] = BT->IA[i];
		row21[0] = (BT->IA[i + 1] - BT->IA[i]) / 2;
		count = 0;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];
			for (i1 = 0; i1<elementDOF[1].col; i1++)
			{
				BT->JA[rowstart[0] + count] = elementDOF[1].val[element[0]][i1];
				BT->JA[rowstart[0] + count + row21[0]] = elementDOF[1].val[element[0]][i1] + elementDOF[1].dof;
				count++;
			}
		} // j
	} // i

	  // step 3BT: Loop element by element and compute the actual entries storing them in BT
	BT->val = (double*)calloc(BT->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			row21[0] = (BT->IA[i + 1] - BT->IA[i]) / 2;
			for (k2 = 0; k2<elementDOF[1].col; k2++)
			{
				j = elementDOF[1].val[k][k2];
				// b11
				for (j1 = BT->IA[i]; j1<BT->IA[i] + row21[0]; j1++)
				{
					if (BT->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							huzhang_basisDIV(lambdas, s, eta, xi, elements->nvector[k], elements->tvector[k], k1, elementDOF[0].dop, phi1);
							lagrange_basis(lambdas, k2, elementDOF[1].dop, &phi);
							BT->val[j1] += 2 * s*gauss22[i1][2] * phi1[0] * phi;
						}
						break;
					}
					if (j1 == (BT->IA[i] + row21[0]))
					{
						printf("There is something wrong in constructing b11 of BT\n");
						exit(1);
					}
				}
				// b12
				for (j1 = BT->IA[i] + row21[0]; j1<BT->IA[i + 1]; j1++)
				{
					if (BT->JA[j1] == (j + elementDOF[1].dof))
					{
						for (i1 = 0; i1<num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							huzhang_basisDIV(lambdas, s, eta, xi, elements->nvector[k], elements->tvector[k], k1, elementDOF[0].dop, phi1);
							lagrange_basis(lambdas, k2, elementDOF[1].dop, &phi);
							BT->val[j1] += 2 * s*gauss22[i1][2] * phi1[1] * phi;
						}
						break;
					}
					if (j1 == (BT->IA[i + 1]))
					{
						printf("There is something wrong in constructing b12 of BT\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k


	  /************************************************** stiffness matrix C *****************************************************************/
	int curnode[2];
	// step 1C: Find first the structure IA of the stiffness matrix C
	for (k = 0; k<elements->row; k++)
	{
		if (elementDOF[1].dop == 0)
		{
			curnode[0] = elementDOF[1].val[k][0];
			curnode[1] = curnode[0] + elementDOF[1].dof;

			C->IA[curnode[0] + 1] += 2;
			C->IA[curnode[1] + 1] += 2;
			for (i = 0; i<3; i++)
			{
				edge = elementEdge->val[k][i];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					C->IA[curnode[0] + 1] += 2;
					C->IA[curnode[1] + 1] += 2;
				}
			}
			continue;
		}

		for (i = 0; i<3 * elementDOF[1].dop; i++) //  for each node
		{
			curnode[0] = elementDOF[1].val[k][i];
			curnode[1] = curnode[0] + elementDOF[1].dof;

			if (i<3)
			{
				C->IA[curnode[0] + 1] += 2 * (elementDOF[1].dop * 2 + 1);
				C->IA[curnode[1] + 1] += 2 * (elementDOF[1].dop * 2 + 1);
				edge = elementEdge->val[k][(i + 1) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					C->IA[curnode[0] + 1] += 2 * (elementDOF[1].dop + 1);
					C->IA[curnode[1] + 1] += 2 * (elementDOF[1].dop + 1);
				}
				edge = elementEdge->val[k][(i + 2) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					C->IA[curnode[0] + 1] += 2 * (elementDOF[1].dop + 1);
					C->IA[curnode[1] + 1] += 2 * (elementDOF[1].dop + 1);
				}
			}
			else
			{
				C->IA[curnode[0] + 1] += 2 * (elementDOF[1].dop + 1);
				C->IA[curnode[1] + 1] += 2 * (elementDOF[1].dop + 1);
				edge = elementEdge->val[k][(i - 3) / (elementDOF[1].dop - 1)];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					C->IA[curnode[0] + 1] += 2 * (elementDOF[1].dop + 1);
					C->IA[curnode[1] + 1] += 2 * (elementDOF[1].dop + 1);
				}
			}
		} // i
	} // k

	for (i = 0; i<C->row; i++)
		C->IA[i + 1] += C->IA[i];

	C->nnz = C->IA[C->row];

	// step 2C: Find the structure JA of the stiffness matrix C
	C->JA = (int*)calloc(C->nnz, sizeof(int));
	for (k = 0; k<elements->row; k++)
	{
		if (elementDOF[1].dop == 0)
		{
			curnode[0] = elementDOF[1].val[k][0];
			curnode[1] = curnode[0] + elementDOF[1].dof;

			for (j = 0; j<2; j++)
			{
				rowstart[j] = C->IA[curnode[j]];
				row21[j] = (C->IA[curnode[j] + 1] - C->IA[curnode[j]]) / 2;
			}
			count = 0;
			for (j = 0; j<2; j++)
			{
				C->JA[rowstart[j] + count] = curnode[0];
				C->JA[rowstart[j] + count + row21[j]] = curnode[1];
			}
			count++;
			for (i = 0; i<3; i++)
			{
				edge = elementEdge->val[k][i];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					for (j = 0; j<2; j++)
					{
						C->JA[rowstart[j] + count] = elementDOF[1].val[element[0]][0];
						C->JA[rowstart[j] + count + row21[j]] = elementDOF[1].val[element[0]][0] + elementDOF[1].dof;
					}
					count++;
				}
			}
			continue;
		}

		// elementDOF[1].dop > 0
		for (i = 0; i<3 * elementDOF[1].dop; i++) //  for each node
		{
			curnode[0] = elementDOF[1].val[k][i];
			curnode[1] = curnode[0] + elementDOF[1].dof;

			for (j = 0; j<2; j++)
			{
				rowstart[j] = C->IA[curnode[j]];
				row21[j] = (C->IA[curnode[j] + 1] - C->IA[curnode[j]]) / 2;
			}
			count = 0;
			if (i<3)
			{
				for (i1 = 0; i1<3; i1++)
				{
					node = elementDOF[1].val[k][i1];
					for (j = 0; j<2; j++)
					{
						C->JA[rowstart[j] + count] = node;
						C->JA[rowstart[j] + count + row21[j]] = node + elementDOF[1].dof;
					}
					count++;
				}
				l = (i + 1) % 3;
				for (i1 = 0; i1<elementDOF[1].dop - 1; i1++)
				{
					node = elementDOF[1].val[k][3 + l*(elementDOF[1].dop - 1) + i1];
					for (j = 0; j<2; j++)
					{
						C->JA[rowstart[j] + count] = node;
						C->JA[rowstart[j] + count + row21[j]] = node + elementDOF[1].dof;
					}
					count++;
				}
				l = (i + 2) % 3;
				for (i1 = 0; i1<elementDOF[1].dop - 1; i1++)
				{
					node = elementDOF[1].val[k][3 + l*(elementDOF[1].dop - 1) + i1];
					for (j = 0; j<2; j++)
					{
						C->JA[rowstart[j] + count] = node;
						C->JA[rowstart[j] + count + row21[j]] = node + elementDOF[1].dof;
					}
					count++;
				}

				edge = elementEdge->val[k][(i + 1) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
					count = getEdgeDOFsVector(C, count, element[0], edge, elementEdge, &elementDOF[1], rowstart, row21);

				edge = elementEdge->val[k][(i + 2) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
					count = getEdgeDOFsVector(C, count, element[0], edge, elementEdge, &elementDOF[1], rowstart, row21);
			} // if(i<3)
			else
			{
				l = (i - 3) / (elementDOF[1].dop - 1);
				node = elementDOF[1].val[k][(l + 1) % 3];
				for (j = 0; j<2; j++)
				{
					C->JA[rowstart[j] + count] = node;
					C->JA[rowstart[j] + count + row21[j]] = node + elementDOF[1].dof;
				}
				count++;
				node = elementDOF[1].val[k][(l + 2) % 3];
				for (j = 0; j<2; j++)
				{
					C->JA[rowstart[j] + count] = node;
					C->JA[rowstart[j] + count + row21[j]] = node + elementDOF[1].dof;
				}
				count++;
				for (i1 = 0; i1<elementDOF[1].dop - 1; i1++)
				{
					node = elementDOF[1].val[k][3 + l*(elementDOF[1].dop - 1) + i1];
					for (j = 0; j<2; j++)
					{
						C->JA[rowstart[j] + count] = node;
						C->JA[rowstart[j] + count + row21[j]] = node + elementDOF[1].dof;
					}
					count++;
				}

				edge = elementEdge->val[k][l];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
					count = getEdgeDOFsVector(C, count, element[0], edge, elementEdge, &elementDOF[1], rowstart, row21);
			}  // if(i<3) else
		} // i
	} // k

	  // step 3C: Loop edge by edge and compute the actual entries storing them in C
	double elen, C11;
	int patchnodes[100];
	C->val = (double*)calloc(C->nnz, sizeof(double));
	for (edge = 0; edge<edges->row; edge++)
	{
		//		edgeNode[0]=edges->val[edge][0];
		//		edgeNode[1]=edges->val[edge][1];
		element[0] = edges->val[edge][2];
		element[1] = edges->val[edge][3];
		elen = edges->length[edge];

//		C11 = elen;

		count = 0;
		if (elementDOF[1].dop == 0)
		{
			patchnodes[count] = elementDOF[1].val[element[0]][0];
			count++;
			if (element[1] != -1)
			{
				patchnodes[count] = elementDOF[1].val[element[1]][0];
				count++;
			}
		} // if(elementDOF[1].dop==0)
		else
		{
			for (i = 0; i<3; i++)
			{
				if (elementEdge->val[element[0]][i] == edge)
					break;
			}
			patchnodes[count] = elementDOF[1].val[element[0]][(i + 1) % 3];
			count++;
			patchnodes[count] = elementDOF[1].val[element[0]][(i + 2) % 3];
			count++;
			for (j = 0; j<elementDOF[1].dop - 1; j++)
			{
				patchnodes[count] = elementDOF[1].val[element[0]][3 + i*(elementDOF[1].dop - 1) + j];
				count++;
			}

			if (element[1] != -1)
			{
				for (i = 0; i<3; i++)
				{
					if (elementEdge->val[element[1]][i] == edge)
						break;
				}
				patchnodes[count] = elementDOF[1].val[element[1]][(i + 1) % 3];
				count++;
				patchnodes[count] = elementDOF[1].val[element[1]][(i + 2) % 3];
				count++;
				for (j = 0; j<elementDOF[1].dop - 1; j++)
				{
					patchnodes[count] = elementDOF[1].val[element[1]][3 + i*(elementDOF[1].dop - 1) + j];
					count++;
				}
			}
		} // if(elementDOF[1].dop==0) else

		for (k1 = 0; k1<count; k1++)
		{
			curnode[0] = patchnodes[k1];
			curnode[1] = curnode[0] + elementDOF[1].dof;
			for (j = 0; j<2; j++)
			{
				rowstart[j] = C->IA[curnode[j]];
				row21[j] = (C->IA[curnode[j] + 1] - C->IA[curnode[j]]) / 2;
			}
			for (k2 = 0; k2<count; k2++)
			{
				j = patchnodes[k2];
				// c11
				for (j1 = C->IA[curnode[0]]; j1<C->IA[curnode[0]] + row21[0]; j1++)
				{
					if (C->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp11; i1++)
						{
							lambdas[0] = gauss11[i1][0];
							lambdas[1] = 1 - lambdas[0];
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, &elementDOF[1], curnode[0], phi1);
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, &elementDOF[1], j, phi2);
							C->val[j1] += gauss11[i1][1] * phi1[0] * phi2[0];
						}
						break;
					}
				}
				// c22
				for (j1 = C->IA[curnode[1]] + row21[1]; j1<C->IA[curnode[1] + 1]; j1++)
				{
					if (C->JA[j1] == (j + elementDOF[1].dof))
					{
						for (i1 = 0; i1<num_qp11; i1++)
						{
							lambdas[0] = gauss11[i1][0];
							lambdas[1] = 1 - lambdas[0];
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, &elementDOF[1], curnode[1], phi1);
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, &elementDOF[1], j + elementDOF[1].dof, phi2);
							C->val[j1] += gauss11[i1][1] * phi1[1] * phi2[1];
						}
						break;
					}
				}
			} // k2
		} // k1
	} // edge



	  /************************************************** right hand side b *****************************************************************/
	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		s = elements->vol[k];

		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				lagrange_basis(lambdas, k1, elementDOF[1].dop, &phi);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b->val[i] += 2 * s*gauss0[i1][2] * f1(x, y, lambda, mu)*phi;
				b->val[i + elementDOF[1].dof] += 2 * s*gauss0[i1][2] * f2(x, y, lambda, mu)*phi;
			} // i1
		} // k1 
	} // k
}

/**
 * \fn void assembleStiffmatrixElasLagrange(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu)
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
void assembleStiffmatrixElasLagrange(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu)
{
	int i,j,k,l;

	A->row = elementDOF->dof * 2;
	A->col = A->row;
	A->IA=(int*)calloc(A->row+1, sizeof(int));
	A->JA=NULL;
	A->val=NULL;
	
//	create_dvector(A->row, b);

	int nvertices=nodes->row;
	int nedges=edges->row;
	int element[2], edge, node;
	
	double phi, phi1[2], phi2[2];
	int k1,k2,i1,j1,l1,l2,ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], curnode[2];
	int count;
	
	int num_qp21=getNumQuadPoints(elementDOF->dop*2-2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial
	
	/************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index=(int*)calloc(A->col, sizeof(int));
	for(i=0;i<A->col;i++)
		index[i]=-1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for(i=0;i<elementDOF->dof;i++)
	{
		count=0;
		istart=-2;
		for(j=elementdofTran->IA[i];j<elementdofTran->IA[i+1];j++)
		{
			element[0]=elementdofTran->JA[j];

			for(k=0;k<elementDOF->col;k++)
			{
				node=elementDOF->val[element[0]][k];
				if(index[node]==-1)
				{
					index[node]=istart;
					istart=node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count * 2;
		A->IA[elementDOF->dof + i + 1] = count * 2;

		for(j=0;j<count;j++)
		{
			l=istart;
			istart=index[l];
			index[l]=-1;
		}		
	} // i
	
	for(i=0;i<A->row;i++)
		A->IA[i+1]+=A->IA[i];
	
	A->nnz=A->IA[A->row];
	
	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA=(int*)calloc(A->nnz,sizeof(int));
	for(i=0;i<elementDOF->dof;i++)
	{
		istart=-2;
		for(j=elementdofTran->IA[i];j<elementdofTran->IA[i+1];j++)
		{
			element[0]=elementdofTran->JA[j];

			for(k=0;k<elementDOF->col;k++)
			{

				node=elementDOF->val[element[0]][k];
				if(index[node]==-1)
				{
					index[node]=istart;
					istart=node;
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
			istart=index[istart];
			index[A->JA[j]]=-1;
		}
	} // i
	free(index);
	
	// step 3A: Loop element by element and compute the actual entries storing them in A
	A->val=(double*)calloc(A->nnz, sizeof(double));
	for(k=0;k<elements->row;k++)
	{
		// set parameters
		s=elements->vol[k];
		xi=elements->xi[k];
		eta=elements->eta[k];
		// end set parameters

		for(k1=0;k1<elementDOF->col;k1++)
		{
			curnode[0] = elementDOF->val[k][k1];
			curnode[1] = curnode[0] + elementDOF->dof;

			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}

			for(k2=0;k2<elementDOF->col;k2++)
			{
				j=elementDOF->val[k][k2];
				// a11
				for (j1 = A->IA[curnode[0]]; j1<A->IA[curnode[0]] + row21[0]; j1++)
				{
					if(A->JA[j1]==j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] / 2) * 2 * mu;
						}
						break;
					}
				} // j1
				// a12
				for (j1 = A->IA[curnode[0]] + row21[0]; j1<A->IA[curnode[0] + 1]; j1++)
				{
					if (A->JA[j1] == j + elementDOF->dof)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * phi1[1] * phi2[0]/ 2 * 2 * mu;
						}
						break;
					}
				} // j1
				// a21
				for (j1 = A->IA[curnode[1]]; j1<A->IA[curnode[1]] + row21[1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * phi1[0] * phi2[1] / 2 * 2 * mu;
						}
						break;
					}
				} // j1
				// a22
				for (j1 = A->IA[curnode[1]] + row21[1]; j1<A->IA[curnode[1] + 1]; j1++)
				{
					if (A->JA[j1] == (j + elementDOF->dof))
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] / 2 + phi1[1] * phi2[1]) * 2 * mu;
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k

	
	/************************************************** right hand side b *****************************************************************/
/*	int num_qp0=49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for(k=0;k<elements->row;k++)
	{
		for(i=0;i<3;i++)
		{
			j=elements->val[k][i];
			xs[i]=nodes->val[j][0];
			ys[i]=nodes->val[j][1];
		}
		s=elements->vol[k];
		
		for(k1=0;k1<elementDOF->col;k1++)
		{
			i=elementDOF->val[k][k1];
			for(i1=0;i1<num_qp0;i1++)
			{
				lambdas[0]=gauss0[i1][0];
				lambdas[1]=gauss0[i1][1];
				lambdas[2]=1-lambdas[0]-lambdas[1];
				lagrange_basis(lambdas, k1, elementDOF->dop, &phi);
				x=xs[0]*lambdas[0]+xs[1]*lambdas[1]+xs[2]*lambdas[2];
				y=ys[0]*lambdas[0]+ys[1]*lambdas[1]+ys[2]*lambdas[2];
				b->val[i]+=2*s*gauss0[i1][2]*f1(x,y, lambda, mu)*phi;
				b->val[i+elementDOF->dof]+=2*s*gauss0[i1][2]*f2(x,y, lambda, mu)*phi;
			} // i1
		} // k1 
	} // k  */
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
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], curnode[2];
	int count;

	int num_qp21 = getNumQuadPoints(elementDOF->dop * 2 - 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

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
		xi = elements->xi[k];
		eta = elements->eta[k];
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
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1]) * 2 * mu;
						}
						break;
					}
				} // j1
				  // a22
				for (j1 = A->IA[curnode[1]] + row21[1]; j1<A->IA[curnode[1] + 1]; j1++)
				{
					if (A->JA[j1] == (j + elementDOF->dof))
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0]  + phi1[1] * phi2[1]) * 2 * mu;
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k


	  /************************************************** right hand side b *****************************************************************/
	  /*	int num_qp0=49; // the number of numerical intergation points
	  double gauss0[num_qp0][3];
	  init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	  for(k=0;k<elements->row;k++)
	  {
	  for(i=0;i<3;i++)
	  {
	  j=elements->val[k][i];
	  xs[i]=nodes->val[j][0];
	  ys[i]=nodes->val[j][1];
	  }
	  s=elements->vol[k];

	  for(k1=0;k1<elementDOF->col;k1++)
	  {
	  i=elementDOF->val[k][k1];
	  for(i1=0;i1<num_qp0;i1++)
	  {
	  lambdas[0]=gauss0[i1][0];
	  lambdas[1]=gauss0[i1][1];
	  lambdas[2]=1-lambdas[0]-lambdas[1];
	  lagrange_basis(lambdas, k1, elementDOF->dop, &phi);
	  x=xs[0]*lambdas[0]+xs[1]*lambdas[1]+xs[2]*lambdas[2];
	  y=ys[0]*lambdas[0]+ys[1]*lambdas[1]+ys[2]*lambdas[2];
	  b->val[i]+=2*s*gauss0[i1][2]*f1(x,y, lambda, mu)*phi;
	  b->val[i+elementDOF->dof]+=2*s*gauss0[i1][2]*f2(x,y, lambda, mu)*phi;
	  } // i1
	  } // k1
	  } // k  */
}

/**
* \fn void assembleStiffmatrixLagrange(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu)
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
void assembleStiffmatrixLagrange(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu)
{
	int i, j, k, l;

	A->row = elementDOF->dof;
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
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int count;

	int num_qp21 = getNumQuadPoints(elementDOF->dop * 2 - 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

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
		A->IA[i + 1] = count;

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

		for (j = A->IA[i]; j<A->IA[i + 1]; j++)
		{
			A->JA[j] = istart;
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
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			i = elementDOF->val[k][k1];
			for (k2 = 0; k2<elementDOF->col; k2++)
			{
				j = elementDOF->val[k][k2];
				for (j1 = A->IA[i]; j1<A->IA[i+1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1]) * 2 * mu;
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k


	  /************************************************** right hand side b *****************************************************************/
	  /*	int num_qp0=49; // the number of numerical intergation points
	  double gauss0[num_qp0][3];
	  init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	  for(k=0;k<elements->row;k++)
	  {
	  for(i=0;i<3;i++)
	  {
	  j=elements->val[k][i];
	  xs[i]=nodes->val[j][0];
	  ys[i]=nodes->val[j][1];
	  }
	  s=elements->vol[k];

	  for(k1=0;k1<elementDOF->col;k1++)
	  {
	  i=elementDOF->val[k][k1];
	  for(i1=0;i1<num_qp0;i1++)
	  {
	  lambdas[0]=gauss0[i1][0];
	  lambdas[1]=gauss0[i1][1];
	  lambdas[2]=1-lambdas[0]-lambdas[1];
	  lagrange_basis(lambdas, k1, elementDOF->dop, &phi);
	  x=xs[0]*lambdas[0]+xs[1]*lambdas[1]+xs[2]*lambdas[2];
	  y=ys[0]*lambdas[0]+ys[1]*lambdas[1]+ys[2]*lambdas[2];
	  b->val[i]+=2*s*gauss0[i1][2]*f(x,y)*phi;
	  } // i1
	  } // k1
	  } // k  */
}


/**
* \fn void assembleStiffmatrixElasIPDG0(dCSRmat *A, dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void assembleStiffmatrixElasIPDG0(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int dop = elementDOF->dop;
	int dof = elementDOF->dof;
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
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], curnode[2];
	int count;

	int num_qp21 = getNumQuadPoints(elementDOF->dop * 2 - 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp11 = getNumQuadPoints(elementDOF->dop * 2, 1); // the number of numerical intergation points
	double gauss11[num_qp11][2];
	init_Gauss1D(num_qp11, 1, gauss11); // gauss intergation initial


	/************************************************** stiffness matrix A *****************************************************************/
	
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (k = 0; k<elements->row; k++)
	{
		if (dop == 0)
		{
			curnode[0] = elementDOF->val[k][0];
			curnode[1] = curnode[0] + dof;

			A->IA[curnode[0] + 1] += 2;
			A->IA[curnode[1] + 1] += 2;
			for (i = 0; i<3; i++)
			{
				edge = elementEdge->val[k][i];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					A->IA[curnode[0] + 1] += 2;
					A->IA[curnode[1] + 1] += 2;
				}
			}
			continue;
		}

		for (i = 0; i<elementDOF->col; i++) //  for each node
		{
			curnode[0] = elementDOF->val[k][i];
			curnode[1] = curnode[0] + dof;

			A->IA[curnode[0] + 1] += 2 * elementDOF->col;
			A->IA[curnode[1] + 1] += 2 * elementDOF->col;
			if (i<3)
			{
				edge = elementEdge->val[k][(i + 1) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					A->IA[curnode[0] + 1] += 2 * (dop + 1);
					A->IA[curnode[1] + 1] += 2 * (dop + 1);
				}
				edge = elementEdge->val[k][(i + 2) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					A->IA[curnode[0] + 1] += 2 * (dop + 1);
					A->IA[curnode[1] + 1] += 2 * (dop + 1);
				}
			}
			else if (i<3 * dop)
			{
				edge = elementEdge->val[k][(i - 3) / (dop - 1)];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					A->IA[curnode[0] + 1] += 2 * (dop + 1);
					A->IA[curnode[1] + 1] += 2 * (dop + 1);
				}
			}
		} // i
	} // k

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
//	count = 0;
	for (k = 0; k<elements->row; k++)
	{
		if (dop == 0)
		{
			curnode[0] = elementDOF->val[k][0];
			curnode[1] = curnode[0] + dof;

			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}
			count = 0;
			for (j = 0; j<2; j++)
			{
				A->JA[rowstart[j] + count] = curnode[0];
				A->JA[rowstart[j] + count + row21[j]] = curnode[1];
			}
			count++;
			for (i = 0; i<3; i++)
			{
				edge = elementEdge->val[k][i];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					for (j = 0; j<2; j++)
					{
						A->JA[rowstart[j] + count] = elementDOF->val[element[0]][0];
						A->JA[rowstart[j] + count + row21[j]] = elementDOF->val[element[0]][0] + dof;
					}
					count++;
				}
			}
			continue;
		}

		for (i = 0; i<elementDOF->col; i++) //  for each node
		{
			curnode[0] = elementDOF->val[k][i];
			curnode[1] = curnode[0] + dof;

			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}
			count = 0;
			for (i1 = 0; i1<elementDOF->col; i1++)
			{
				node = elementDOF->val[k][i1];
				for (j = 0; j<2; j++)
				{
					A->JA[rowstart[j] + count] = node;
					A->JA[rowstart[j] + count + row21[j]] = node + dof;
				}
				count++;
			}

			if (i<3)
			{
				edge = elementEdge->val[k][(i + 1) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
					count = getEdgeDOFsVector(A, count, element[0], edge, elementEdge, elementDOF, rowstart, row21);

				edge = elementEdge->val[k][(i + 2) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
					count = getEdgeDOFsVector(A, count, element[0], edge, elementEdge, elementDOF, rowstart, row21);
			}
			else if (i<3 * dop)
			{
				edge = elementEdge->val[k][(i - 3) / (dop - 1)];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
					count = getEdgeDOFsVector(A, count, element[0], edge, elementEdge, elementDOF, rowstart, row21);
			}
		} // i
	} // k
	
	// step 3A1: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		if (dop == 0)
			continue;

		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
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
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] / 2);
						}
						break;
					}
				} // j1
				 // a12
				for (j1 = A->IA[curnode[0]] + row21[0]; j1<A->IA[curnode[0] + 1]; j1++)
				{
					if (A->JA[j1] == j + elementDOF->dof)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * phi1[1] * phi2[0] / 2;
						}
						break;
					}
				} // j1
				  // a21
				for (j1 = A->IA[curnode[1]]; j1<A->IA[curnode[1]] + row21[1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * phi1[0] * phi2[1] / 2;
						}
						break;
					}
				} // j1
				  // a22
				for (j1 = A->IA[curnode[1]] + row21[1]; j1<A->IA[curnode[1] + 1]; j1++)
				{
					if (A->JA[j1] == (j + elementDOF->dof))
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] / 2 + phi1[1] * phi2[1]);
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k

	 // step 3A2: Loop edge by edge and compute the actual entries storing them in A
	double elen;
	int patchnodes[100];
	for (edge = 0; edge<edges->row; edge++)
	{
		element[0] = edges->val[edge][2];
		element[1] = edges->val[edge][3];
		elen = edges->length[edge];

		count = 0;
		if (dop == 0)
		{
			patchnodes[count] = elementDOF->val[element[0]][0];
			count++;
			if (element[1] != -1)
			{
				patchnodes[count] = elementDOF->val[element[1]][0];
				count++;
			}
		} // if(dop==0)
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
			for (j = 0; j<dop - 1; j++)
			{
				patchnodes[count] = elementDOF->val[element[0]][3 + i*(dop - 1) + j];
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
				for (j = 0; j<dop - 1; j++)
				{
					patchnodes[count] = elementDOF->val[element[1]][3 + i*(dop - 1) + j];
					count++;
				}
			}
		} // if(dop==0) else

		for (k1 = 0; k1<count; k1++)
		{
			curnode[0] = patchnodes[k1];
			curnode[1] = curnode[0] + elementDOF->dof;
			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}
			for (k2 = 0; k2<count; k2++)
			{
				j = patchnodes[k2];
				// a11
				for (j1 = A->IA[curnode[0]]; j1<A->IA[curnode[0]] + row21[0]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp11; i1++)
						{
							lambdas[0] = gauss11[i1][0];
							lambdas[1] = 1 - lambdas[0];
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, elementDOF, curnode[0], phi1);
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, elementDOF, j, phi2);
							A->val[j1] += gauss11[i1][1]*phi1[0] * phi2[0];
						}
						break;
					}
				}
				// a22
				for (j1 = A->IA[curnode[1]] + row21[1]; j1<A->IA[curnode[1] + 1]; j1++)
				{
					if (A->JA[j1] == (j + elementDOF->dof))
					{
						for (i1 = 0; i1<num_qp11; i1++)
						{
							lambdas[0] = gauss11[i1][0];
							lambdas[1] = 1 - lambdas[0];
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, elementDOF, curnode[1], phi1);
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, elementDOF, j + elementDOF->dof, phi2);
							A->val[j1] += gauss11[i1][1]*phi1[1] * phi2[1];
						}
						break;
					}
				}
			} // k2
		} // k1
	} // edge
	
	  /************************************************** right hand side b *****************************************************************/
/*	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		s = elements->vol[k];

		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			i = elementDOF->val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				lagrange_basis(lambdas, k1, elementDOF->dop, &phi);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b->val[i] += 2 * s*gauss0[i1][2] * f1(x, y, lambda, mu)*phi;
				b->val[i + elementDOF->dof] += 2 * s*gauss0[i1][2] * f2(x, y, lambda, mu)*phi;
			} // i1
		} // k1 
	} // k  */
}

/**
* \fn void assembleStiffmatrixVecIPDG0(dCSRmat *A, dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void assembleStiffmatrixVecIPDG0(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int dop = elementDOF->dop;
	int dof = elementDOF->dof;
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
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], curnode[2];
	int count;

	int num_qp21 = getNumQuadPoints(elementDOF->dop * 2 - 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp11 = getNumQuadPoints(elementDOF->dop * 2, 1); // the number of numerical intergation points
	double gauss11[num_qp11][2];
	init_Gauss1D(num_qp11, 1, gauss11); // gauss intergation initial


										/************************************************** stiffness matrix A *****************************************************************/

										// step 1A: Find first the structure IA of the stiffness matrix A
	for (k = 0; k<elements->row; k++)
	{
		if (dop == 0)
		{
			curnode[0] = elementDOF->val[k][0];
			curnode[1] = curnode[0] + dof;

			A->IA[curnode[0] + 1] += 2;
			A->IA[curnode[1] + 1] += 2;
			for (i = 0; i<3; i++)
			{
				edge = elementEdge->val[k][i];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					A->IA[curnode[0] + 1] += 2;
					A->IA[curnode[1] + 1] += 2;
				}
			}
			continue;
		}

		for (i = 0; i<elementDOF->col; i++) //  for each node
		{
			curnode[0] = elementDOF->val[k][i];
			curnode[1] = curnode[0] + dof;

			A->IA[curnode[0] + 1] += 2 * elementDOF->col;
			A->IA[curnode[1] + 1] += 2 * elementDOF->col;
			if (i<3)
			{
				edge = elementEdge->val[k][(i + 1) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					A->IA[curnode[0] + 1] += 2 * (dop + 1);
					A->IA[curnode[1] + 1] += 2 * (dop + 1);
				}
				edge = elementEdge->val[k][(i + 2) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					A->IA[curnode[0] + 1] += 2 * (dop + 1);
					A->IA[curnode[1] + 1] += 2 * (dop + 1);
				}
			}
			else if (i<3 * dop)
			{
				edge = elementEdge->val[k][(i - 3) / (dop - 1)];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					A->IA[curnode[0] + 1] += 2 * (dop + 1);
					A->IA[curnode[1] + 1] += 2 * (dop + 1);
				}
			}
		} // i
	} // k

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
	//	count = 0;
	for (k = 0; k<elements->row; k++)
	{
		if (dop == 0)
		{
			curnode[0] = elementDOF->val[k][0];
			curnode[1] = curnode[0] + dof;

			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}
			count = 0;
			for (j = 0; j<2; j++)
			{
				A->JA[rowstart[j] + count] = curnode[0];
				A->JA[rowstart[j] + count + row21[j]] = curnode[1];
			}
			count++;
			for (i = 0; i<3; i++)
			{
				edge = elementEdge->val[k][i];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					for (j = 0; j<2; j++)
					{
						A->JA[rowstart[j] + count] = elementDOF->val[element[0]][0];
						A->JA[rowstart[j] + count + row21[j]] = elementDOF->val[element[0]][0] + dof;
					}
					count++;
				}
			}
			continue;
		}

		for (i = 0; i<elementDOF->col; i++) //  for each node
		{
			curnode[0] = elementDOF->val[k][i];
			curnode[1] = curnode[0] + dof;

			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}
			count = 0;
			for (i1 = 0; i1<elementDOF->col; i1++)
			{
				node = elementDOF->val[k][i1];
				for (j = 0; j<2; j++)
				{
					A->JA[rowstart[j] + count] = node;
					A->JA[rowstart[j] + count + row21[j]] = node + dof;
				}
				count++;
			}

			if (i<3)
			{
				edge = elementEdge->val[k][(i + 1) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
					count = getEdgeDOFsVector(A, count, element[0], edge, elementEdge, elementDOF, rowstart, row21);

				edge = elementEdge->val[k][(i + 2) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
					count = getEdgeDOFsVector(A, count, element[0], edge, elementEdge, elementDOF, rowstart, row21);
			}
			else if (i<3 * dop)
			{
				edge = elementEdge->val[k][(i - 3) / (dop - 1)];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
					count = getEdgeDOFsVector(A, count, element[0], edge, elementEdge, elementDOF, rowstart, row21);
			}
		} // i
	} // k

	  // step 3A1: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		if (dop == 0)
			continue;

		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
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
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1]);
						}
						break;
					}
				} // j1
				  // a22
				for (j1 = A->IA[curnode[1]] + row21[1]; j1<A->IA[curnode[1] + 1]; j1++)
				{
					if (A->JA[j1] == (j + elementDOF->dof))
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1]);
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k

	  // step 3A2: Loop edge by edge and compute the actual entries storing them in A
	double elen;
	int patchnodes[100];
	for (edge = 0; edge<edges->row; edge++)
	{
		element[0] = edges->val[edge][2];
		element[1] = edges->val[edge][3];
		elen = edges->length[edge];

		count = 0;
		if (dop == 0)
		{
			patchnodes[count] = elementDOF->val[element[0]][0];
			count++;
			if (element[1] != -1)
			{
				patchnodes[count] = elementDOF->val[element[1]][0];
				count++;
			}
		} // if(dop==0)
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
			for (j = 0; j<dop - 1; j++)
			{
				patchnodes[count] = elementDOF->val[element[0]][3 + i*(dop - 1) + j];
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
				for (j = 0; j<dop - 1; j++)
				{
					patchnodes[count] = elementDOF->val[element[1]][3 + i*(dop - 1) + j];
					count++;
				}
			}
		} // if(dop==0) else

		for (k1 = 0; k1<count; k1++)
		{
			curnode[0] = patchnodes[k1];
			curnode[1] = curnode[0] + elementDOF->dof;
			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}
			for (k2 = 0; k2<count; k2++)
			{
				j = patchnodes[k2];
				// a11
				for (j1 = A->IA[curnode[0]]; j1<A->IA[curnode[0]] + row21[0]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp11; i1++)
						{
							lambdas[0] = gauss11[i1][0];
							lambdas[1] = 1 - lambdas[0];
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, elementDOF, curnode[0], phi1);
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, elementDOF, j, phi2);
							A->val[j1] += gauss11[i1][1] * phi1[0] * phi2[0];
						}
						break;
					}
				}
				// a22
				for (j1 = A->IA[curnode[1]] + row21[1]; j1<A->IA[curnode[1] + 1]; j1++)
				{
					if (A->JA[j1] == (j + elementDOF->dof))
					{
						for (i1 = 0; i1<num_qp11; i1++)
						{
							lambdas[0] = gauss11[i1][0];
							lambdas[1] = 1 - lambdas[0];
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, elementDOF, curnode[1], phi1);
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, elementDOF, j + elementDOF->dof, phi2);
							A->val[j1] += gauss11[i1][1] * phi1[1] * phi2[1];
						}
						break;
					}
				}
			} // k2
		} // k1
	} // edge

	  /************************************************** right hand side b *****************************************************************/
	  /*	int num_qp0 = 49; // the number of numerical intergation points
	  double gauss0[num_qp0][3];
	  init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	  for (k = 0; k<elements->row; k++)
	  {
	  for (i = 0; i<3; i++)
	  {
	  j = elements->val[k][i];
	  xs[i] = nodes->val[j][0];
	  ys[i] = nodes->val[j][1];
	  }
	  s = elements->vol[k];

	  for (k1 = 0; k1<elementDOF->col; k1++)
	  {
	  i = elementDOF->val[k][k1];
	  for (i1 = 0; i1<num_qp0; i1++)
	  {
	  lambdas[0] = gauss0[i1][0];
	  lambdas[1] = gauss0[i1][1];
	  lambdas[2] = 1 - lambdas[0] - lambdas[1];
	  lagrange_basis(lambdas, k1, elementDOF->dop, &phi);
	  x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
	  y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
	  b->val[i] += 2 * s*gauss0[i1][2] * f1(x, y, lambda, mu)*phi;
	  b->val[i + elementDOF->dof] += 2 * s*gauss0[i1][2] * f2(x, y, lambda, mu)*phi;
	  } // i1
	  } // k1
	  } // k  */
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
	double s, *eta, *xi, phi[2], nve[2];
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
		xi=elements->xi[element[0]];
		eta=elements->eta[element[0]];
		lagrange_basis1(lambdas, s, eta, xi, nodeindex1, dop, phi);
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
		xi=elements->xi[element[1]];
		eta=elements->eta[element[1]];
		lagrange_basis1(lambdas, s, eta, xi, nodeindex2, dop, phi);
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

	xi = elements->xi[element];
	eta = elements->eta[element];
	nv[0] = -eta[edgeindex] / elen;
	nv[1] = xi[edgeindex] / elen;
	nve[0] = edges->nvector[edge][0];
	nve[1] = edges->nvector[edge][1];
	double sgn = nve[0] * nv[0] + nve[1] * nv[1];

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
	double s, *eta, *xi, phi[2], nv[2];
		
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
	xi=elements->xi[element];
	eta=elements->eta[element];
	nv[0]=-eta[l]/elen;
	nv[1]=xi[l]/elen;
	
	lagrange_basis1(lambdas, s, eta, xi, index, dop, phi);
		
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
	double s, *eta, *xi, phi[2], tve[2];
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
	xi=elements->xi[element];
	eta=elements->eta[element];
	tve[0]=edges->tvector[edge][0];
	tve[1]=edges->tvector[edge][1];
	
	lagrange_basis1(lambdas, s, eta, xi, nodeindex, dop, phi);
		
	*val=phi[0]*tve[0]+phi[1]*tve[1];
}

/**
 * \fn void averageOperator(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *average)
 * \brief the average operator for gradient of basis function
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
void averageOperator(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *average)
{
	int i;
	int nodeindex1,nodeindex2, edgeindex;
	int element[2], edgeNode[2];
	double averageplus[2];
	int isInedge=0;
	int nedges=edges->row;
	int dop=elementDOF->dop;

	average[0]=0;
	average[1]=0;
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
		getinfo4average(lambda1, lambda2, edgeNode, element[0], elements, dop, nodeindex1, average);
		if(element[1]!=-1)
		{			
			getinfo4average(lambda1, lambda2, edgeNode, element[1], elements, dop, nodeindex2, averageplus);
			average[0] += averageplus[0];
			average[1] += averageplus[1];
		}
	}
	else
	{
		if(nodeindex1<elementDOF->col)
			getinfo4average(lambda1, lambda2, edgeNode, element[0], elements, dop, nodeindex1, average);
		else{
				if(element[1]!=-1)
				{
					if(nodeindex2<elementDOF->col)
						getinfo4average(lambda1, lambda2, edgeNode, element[1], elements, dop, nodeindex2, average);
				}
			}
	}

	if(element[1]!=-1)
	{
		average[0]/=2.0;
		average[1]/=2.0;
	}
}

/**
 * \fn void getinfo4average(double lambda1, double lambda2, int *edgeNode, int element, ELEMENT *elements, int dop, int index, double *average)
 * \brief the gradient of basis function
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param *edgeNode pointer to the two indices of current edge
 * \param element the index of current element
 * \param *elements pointer to the structure of the triangulation
 * \param dop number of degrees of polynomial
 * \param index index of current node variable
 * \param *average pointer to gradient of basis function
 * \return void
 */
void getinfo4average(double lambda1, double lambda2, int *edgeNode, int element, ELEMENT *elements, int dop, int index, double *average)
{
	int i,j,l;
	double lambdas[3];
	double s, *eta, *xi;
		
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
	xi=elements->xi[element];
	eta=elements->eta[element];
	
	lagrange_basis1(lambdas, s, eta, xi, index, dop, average);
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
	double s, *eta, *xi, nv[2], tv[2];
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
	xi=elements->xi[element];
	eta=elements->eta[element];
	nv[0]=-eta[edgeindex]/elen;
	nv[1]=xi[edgeindex]/elen;
	tv[0]=-nv[1];
	tv[1]=nv[0];
	
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

	lagrange_basis1(lambdas, s, eta, xi, nodeindex, elementDOF->dop, phi);

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
	double s, *eta, *xi, nv[2], tv[2];
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
	xi=elements->xi[element];
	eta=elements->eta[element];
	nv[0]=-eta[edgeindex]/elen;
	nv[1]=xi[edgeindex]/elen;
	
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

	lagrange_basis1(lambdas, s, eta, xi, nodeindex, elementDOF->dop, phi);

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
	double *eta, *xi, nv[2], nve[2], tv[2];
	double phi[3], sgn, value[2];
	
	*jump = 0;
	
	if (node<0 && node >= elementDOF->dof)
		return;

	double elen = edges->length[edge];
	edgeNode[0] = edges->val[edge][0];
	edgeNode[1] = edges->val[edge][1];
	nve[0] = edges->nvector[edge][0];
	nve[1] = edges->nvector[edge][1];

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

		xi = elements->xi[element];
		eta = elements->eta[element];
		nv[0] = -eta[edgeindex] / elen;
		nv[1] = xi[edgeindex] / elen;
		tv[0] = -nv[1];
		tv[1] = nv[0];
		sgn = nve[0] * nv[0] + nve[1] * nv[1];

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
	double s, *eta, *xi, nv[2], nve[2], tv[2];
	double phi1[2], phi2[2], phi3[3][2], sgn, value[3];

	*jump = 0;

	if (node<0 && node >= elementDOF->dof)
		return;

	double elen = edges->length[edge];
	edgeNode[0] = edges->val[edge][0];
	edgeNode[1] = edges->val[edge][1];
//	nve[0] = edges->nvector[edge][0];
//	nve[1] = edges->nvector[edge][1];

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
		xi = elements->xi[element];
		eta = elements->eta[element];
		nv[0] = -eta[edgeindex] / elen;
		nv[1] = xi[edgeindex] / elen;
		tv[0] = -nv[1];
		tv[1] = nv[0];
//		sgn = nve[0] * nv[0] + nve[1] * nv[1];

		huzhang_basisROT(lambdas, s, eta, xi, elements->nvector[element], elements->tvector[element], nodeindex, elementDOF->dop, phi1);
		huzhang_basisCurlTrace(lambdas, s, eta, xi, elements->nvector[element], elements->tvector[element], nodeindex, elementDOF->dop, phi2);

		value[0] = phi1[0] * tv[0] + phi1[1] * tv[1];
		value[1] = phi2[0] * tv[0] + phi2[1] * tv[1];

		*jump += (value[0] - value[1] * lambda / (2 * lambda + 2 * mu)) / (2 * mu);

		huzhang_basis1(lambdas, s, eta, xi, elements->nvector[element], elements->tvector[element], nodeindex, elementDOF->dop, phi3);
		value[0] = phi3[0][0] * tv[0] + phi3[0][1] * tv[1];
		value[1] = phi3[1][0] * tv[0] + phi3[1][1] * tv[1];
		value[2] = phi3[2][0] * tv[0] + phi3[2][1] * tv[1];

		*jump -= (value[0] * nv[0] * tv[0] + value[1] * nv[1] * tv[1] + value[2] * (nv[0] * tv[1] + nv[1] * tv[0])) / (2 * mu);
	} // k	
}

/**
 * \fn void getElementEdgeGeoInfo(ELEMENT *elements, EDGE *edges, dennode *nodes)
 * \brief compute the geometric information of elements and edges
 * \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
 * \return void
 */
void getElementEdgeGeoInfo(ELEMENT *elements, EDGE *edges, dennode *nodes)
{
	int i,j;
	int vertex1, vertex2, vertex3;
	int element1, element2;
/*	if(!elements->h)
		free(elements->h);
	if(!edges->h)
		free(edges->h);
	elements->h=(double*)calloc(elements->row, sizeof(double));
	edges->h=(double*)calloc(edges->row, sizeof(double)); */
	for(i=0;i<elements->row;i++)
	{
		vertex1 = elements->val[i][0];
		vertex2 = elements->val[i][1];
		vertex3 = elements->val[i][2];
		elements->vol[i] = area(nodes->val[vertex1][0],nodes->val[vertex2][0],nodes->val[vertex3][0],nodes->val[vertex1][1],nodes->val[vertex2][1],nodes->val[vertex3][1]);
		elements->xi[i][0] = nodes->val[vertex2][0]-nodes->val[vertex3][0];
		elements->xi[i][1] = nodes->val[vertex3][0]-nodes->val[vertex1][0];
		elements->xi[i][2] = nodes->val[vertex1][0]-nodes->val[vertex2][0];
		elements->eta[i][0] = nodes->val[vertex2][1]-nodes->val[vertex3][1];
		elements->eta[i][1] = nodes->val[vertex3][1]-nodes->val[vertex1][1];
		elements->eta[i][2] = nodes->val[vertex1][1]-nodes->val[vertex2][1];
		for(j=0;j<3;j++)
		{
			elements->edgeslength[i][j]=sqrt(elements->xi[i][j]*elements->xi[i][j]+elements->eta[i][j]*elements->eta[i][j]);
			elements->nvector[i][j][0]=-elements->eta[i][j]/elements->edgeslength[i][j];
		    elements->nvector[i][j][1]=elements->xi[i][j]/elements->edgeslength[i][j];
			elements->tvector[i][j][0]=-elements->nvector[i][j][1];
			elements->tvector[i][j][1]=elements->nvector[i][j][0];
		}
	}

	for(i=0;i<edges->row;i++)
	{
		vertex1=edges->val[i][0];
		vertex2=edges->val[i][1];
		edges->xi[i]=nodes->val[vertex1][0]-nodes->val[vertex2][0];
		edges->eta[i]=nodes->val[vertex1][1]-nodes->val[vertex2][1];
		edges->length[i]=sqrt(edges->xi[i]*edges->xi[i]+edges->eta[i]*edges->eta[i]);
		edges->nvector[i][0]=-edges->eta[i]/edges->length[i];
	    edges->nvector[i][1]=edges->xi[i]/edges->length[i];
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
		isInNode->val[i]=nodes->isInNode[i];
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

	for (i = 0; i<nn; i++)
	{
		isInNode->val[i] = nodes->isInNode[i];
		isInNode->val[i + dof] = nodes->isInNode[i];
	}
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
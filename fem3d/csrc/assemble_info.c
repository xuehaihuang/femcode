/*
 *  assemble_info.c
 *
 *  Created by Xuehai Huang on 3/29/09.
 *  Copyright 2009 PSU. All rights reserved.
 *
 */

/*! \file assemble_info.c
 *  \brief Subroutines for assembling purpose
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "header.h"
#include "matvec.h"

static void addFace(int i, int point1, int point2, int element, int iSTART, FACE *faces);
//static int vv2edge3d(int v1, int v2);
static void getEdgeElement3(ELEMENT *elements, iCSRmat *elementsTran, EDGE *edges, iCSRmat *edgeElement);

/**
 * \fn int getCoarseInfo(int domain_num, dennode *nodes, ELEMENT *elements, FACE *faces, EDGE *edges, iCSRmat *elementsTran, idenmat *elementEdge)
 * \brief generate the coarse grid information and store it into nodes, elements, faces respectively
 * \param domain_num number of domain
 * \param *nodes the first column stores the x coordinate of points, the second column stores the y coordinate of points, 
                 the third column stores the z coordinate of points
 * \param *elements stores 4 nodes corresponding to the element, the 5th column stores the relation with coarse grid elements( store -1 if itself is the coarset grid),
                    the 6th column stores the type of element
 * \param *faces the first three columns store the three vertices corresponding to the face; 
 *				 the 4rd and 5th columns store the elements which the face belongs to;
 *				 if the face is a boundary, the 5th column will stores -1;
 *				 the first column is in ascend order.
 * \param *elementsTran pointer to the transpose of *elements
 * \param *facesTran the relation between nodes and faces. JA stores face index, A stores another two vertices
 * \param *isInNode if the node is interior node, it will be 0; if the node is on the boundary, it will be -1
 * \return 1 if succeed 0 if fail
 */
int getCoarseInfo(int domain_num, dennode *nodes, ELEMENT *elements, FACE *faces, EDGE *edges, iCSRmat *elementsTran, idenmat *elementEdge)
{
	// get data from inputFile
	char *str1 = "data/testdata.dat";
	char *str2 = "data/testdata.dat";
	char *str3 = "data/testdata.dat";

	char *filename;
	switch (domain_num) {
	case 1:filename = str1; break;
	case 2:filename = str2; break;
	case 3:filename = str3; break;
	default:filename = str1;
	}

	FILE *inputFile;
	inputFile=fopen(filename, "r");
	if(inputFile==NULL)
	{
		printf("Opening file %s fails!\n", filename);
		return 0;
	}
	
	int NumNode, Ncoor, NumElem, Nnode;
	int i,j;
	
	// get the nodes' coordinates
	fscanf(inputFile, "%d %d", &NumNode, &Ncoor);
	create_dennode(NumNode, Ncoor, nodes);	
	for(i=0;i<nodes->row;i++)
	{
		for(j=0;j<nodes->col;j++)
			fscanf(inputFile, "%lf", &nodes->val[i][j]);
		fscanf(inputFile, "%d", &nodes->bdFlag[i]);
	}
	
	// get triangular grid
	fscanf(inputFile, "%d %d", &NumElem, &Nnode);
	create_ELEMENT(NumElem, Nnode, elements);

	for(i=0;i<elements->row;i++)
	{
		for(j=0;j<elements->col;j++)
		{
			fscanf(inputFile, "%d", &elements->val[i][j]);
			elements->val[i][j]--; // the data from matlab start with 1 while from c start with 0
		}
		elements->parent[i] = -1;
		fscanf(inputFile, "%d", &elements->type[i]);
	}	
	
	fclose(inputFile);
	
	getTransposeOfELEMENT(elements, elementsTran, elements->col, nodes->row);
	
	// get face and edge information
	getFaceInfo(elements, elementsTran, faces, nodes);
	getEdgeInfo3(elements, elementsTran, faces, elementEdge, edges, nodes);
	
	// get nodeCEdge
/*	nodeCEdge->row = nodes->row;
	nodeCEdge->val = (int*)calloc(nodeCEdge->row, sizeof(int));
	for (i = 0; i<nodeCEdge->row; i++)
		nodeCEdge->val[i] = -1;
	*/
	return 1;
}

/**
 * \fn void uniformrefine(dennode *Cnodes, ELEMENT *Celements, FACE *Cfaces, EDGE *Cedges, iCSRmat *CelementsTran, dennode *Fnodes, ELEMENT *Felements, FACE *Ffaces, EDGE *Fedges, iCSRmat *FelementsTran, idenmat *FelementEdge)
 * \brief generate fine grid using regular section
 * \param *nodes pointer to the nodes location of the triangulation
 * \param *Celements pointer to the structure of the triangulation on coarse grid
 * \param *Cfaces pointer to the face information of the triangulation on coarse grid
 * \param *CelementsTran pointer to the transpose of *Celements
 * \param *CfacesTran the relation between nodes and faces on coarse grid. JA stores face index, A stores another two vertices
 * \param *Felements pointer to the structure of the triangulation on fine grid
 * \param *Ffaces pointer to the face information of the triangulation on fine grid
 * \param *FelementsTran pointer to the transpose of *Felements
 * \param *FfacesTran the relation between nodes and faces on fine grid. JA stores face index, A stores another two vertices
 * \return void
 */
void uniformrefine(dennode *Cnodes, ELEMENT *Celements, FACE *Cfaces, EDGE *Cedges, iCSRmat *CelementsTran, dennode *Fnodes, ELEMENT *Felements, FACE *Ffaces, EDGE *Fedges, iCSRmat *FelementsTran, idenmat *FelementEdge)
{
	int i,j,k,l;
	int element;
	int NumCNodes;
	NumCNodes=Cnodes->row;

	iCSRmat CedgeElement;
	getEdgeElement3(Celements, CelementsTran, Cedges, &CedgeElement);
	
	// generate fine grid's nodes information
	//	int midElements[Celements->row][6]; // store the 6 middle point of each tetrahedron
	int **midElements; // store the 6 middle point of each tetrahedron
	midElements = (int**)malloc(Celements->row * sizeof(int *));
	for (i = 0; i<Celements->row; i++)
		midElements[i] = (int*)malloc(6 * sizeof(int));

	create_dennode(Cnodes->row + Cedges->row, Cnodes->col, Fnodes); 
	// copy Cnodes into Fnodes
	for (i = 0; i < Cnodes->row; i++)
	{
		for (j = 0; j < Cnodes->col; j++)
			Fnodes->val[i][j] = Cnodes->val[i][j];
		Fnodes->bdFlag[i] = Cnodes->bdFlag[i];
	}

	int col=Celements->col;
	int point1, point2;
	for(i=0;i<Cedges->row;i++)
	{
		point1=Cedges->val[i][0];
		point2=Cedges->val[i][1];
		Fnodes->val[NumCNodes+i][0]=(Fnodes->val[point1][0]+Fnodes->val[point2][0])/2;
		Fnodes->val[NumCNodes+i][1]=(Fnodes->val[point1][1]+Fnodes->val[point2][1])/2;
		Fnodes->val[NumCNodes+i][2]=(Fnodes->val[point1][2]+Fnodes->val[point2][2])/2;
		Fnodes->bdFlag[NumCNodes + i] = Cedges->bdFlag[i];

		// get the relation between tetrahedron and edge
		for(j=CedgeElement.IA[i];j<CedgeElement.IA[i+1];j++)
		{
			element=CedgeElement.JA[j];
			l=CedgeElement.val[j];
			midElements[element][l]=i;
		}
	}
	
	free_icsr_matrix(&CedgeElement);

	// generate fine grid tetrahedrons information
	int type;
	create_ELEMENT(Celements->row * 8, Celements->col, Felements);
	
	for(i=0;i<Celements->row;i++)
	{
		// regular section start
		type=Celements->type[i];
		Felements->val[i*8][0] = Celements->val[i][0];
		Felements->val[i*8][1] = midElements[i][vv2edge3d(0,1)]+NumCNodes;
		Felements->val[i*8][2] = midElements[i][vv2edge3d(0,2)]+NumCNodes;
		Felements->val[i*8][3] = midElements[i][vv2edge3d(0,3)]+NumCNodes;
		Felements->parent[i*8] = i;
		Felements->type[i*8] = type;
		Felements->val[i*8+1][0] = midElements[i][vv2edge3d(1,0)]+NumCNodes;
		Felements->val[i*8+1][1] = Celements->val[i][1];
		Felements->val[i*8+1][2] = midElements[i][vv2edge3d(1,2)]+NumCNodes;
		Felements->val[i*8+1][3] = midElements[i][vv2edge3d(1,3)]+NumCNodes;
		Felements->parent[i*8+1] = i;
		Felements->type[i*8+1] = type;
		Felements->val[i*8+2][0] = midElements[i][vv2edge3d(2,0)]+NumCNodes;
		Felements->val[i*8+2][1] = midElements[i][vv2edge3d(2,1)]+NumCNodes;
		Felements->val[i*8+2][2] = Celements->val[i][2];
		Felements->val[i*8+2][3] = midElements[i][vv2edge3d(2,3)]+NumCNodes;
		Felements->parent[i*8+2] = i;
		Felements->type[i*8+2] = type;
		Felements->val[i*8+3][0] = midElements[i][vv2edge3d(3,0)]+NumCNodes;
		Felements->val[i*8+3][1] = midElements[i][vv2edge3d(3,1)]+NumCNodes;
		Felements->val[i*8+3][2] = midElements[i][vv2edge3d(3,2)]+NumCNodes;
		Felements->val[i*8+3][3] = Celements->val[i][3];
		Felements->parent[i*8+3] = i;
		Felements->type[i*8+3] = type;
		if(type==1) // type = 1
		{
			Felements->val[i*8+4][0] = midElements[i][vv2edge3d(2,0)]+NumCNodes;
			Felements->val[i*8+4][1] = midElements[i][vv2edge3d(2,1)]+NumCNodes;
			Felements->val[i*8+4][2] = midElements[i][vv2edge3d(2,3)]+NumCNodes;
			Felements->val[i*8+4][3] = midElements[i][vv2edge3d(3,0)]+NumCNodes;
			Felements->parent[i*8+4] = i;
			Felements->type[i*8+4] = 2;
			Felements->val[i*8+5][0] = midElements[i][vv2edge3d(2,0)]+NumCNodes;
			Felements->val[i*8+5][1] = midElements[i][vv2edge3d(2,1)]+NumCNodes;
			Felements->val[i*8+5][2] = midElements[i][vv2edge3d(3,0)]+NumCNodes;
			Felements->val[i*8+5][3] = midElements[i][vv2edge3d(0,1)]+NumCNodes;
			Felements->parent[i*8+5] = i;
			Felements->type[i*8+5] = 1;
			Felements->val[i*8+6][0] = midElements[i][vv2edge3d(3,0)]+NumCNodes;
			Felements->val[i*8+6][1] = midElements[i][vv2edge3d(3,1)]+NumCNodes;
			Felements->val[i*8+6][2] = midElements[i][vv2edge3d(0,1)]+NumCNodes;
			Felements->val[i*8+6][3] = midElements[i][vv2edge3d(1,2)]+NumCNodes;
			Felements->parent[i*8+6] = i;
			Felements->type[i*8+6] = 1;
			Felements->val[i*8+7][0] = midElements[i][vv2edge3d(3,0)]+NumCNodes;
			Felements->val[i*8+7][1] = midElements[i][vv2edge3d(3,1)]+NumCNodes;
			Felements->val[i*8+7][2] = midElements[i][vv2edge3d(1,2)]+NumCNodes;
			Felements->val[i*8+7][3] = midElements[i][vv2edge3d(2,3)]+NumCNodes;
			Felements->parent[i*8+7] = i;
			Felements->type[i*8+7] = 2;
		}
		else // type = 2
		{
			Felements->val[i*8+4][0] = midElements[i][vv2edge3d(2,0)]+NumCNodes;
			Felements->val[i*8+4][1] = midElements[i][vv2edge3d(2,1)]+NumCNodes;
			Felements->val[i*8+4][2] = midElements[i][vv2edge3d(2,3)]+NumCNodes;
			Felements->val[i*8+4][3] = midElements[i][vv2edge3d(3,1)]+NumCNodes;
			Felements->parent[i*8+4] = i;
			Felements->type[i*8+4] = 1;
			Felements->val[i*8+5][0] = midElements[i][vv2edge3d(2,0)]+NumCNodes;
			Felements->val[i*8+5][1] = midElements[i][vv2edge3d(2,1)]+NumCNodes;
			Felements->val[i*8+5][2] = midElements[i][vv2edge3d(3,1)]+NumCNodes;
			Felements->val[i*8+5][3] = midElements[i][vv2edge3d(1,0)]+NumCNodes;
			Felements->parent[i*8+5] = i;
			Felements->type[i*8+5] = 2;
			Felements->val[i*8+6][0] = midElements[i][vv2edge3d(3,0)]+NumCNodes;
			Felements->val[i*8+6][1] = midElements[i][vv2edge3d(3,1)]+NumCNodes;
			Felements->val[i*8+6][2] = midElements[i][vv2edge3d(1,0)]+NumCNodes;
			Felements->val[i*8+6][3] = midElements[i][vv2edge3d(0,2)]+NumCNodes;
			Felements->parent[i*8+6] = i;
			Felements->type[i*8+6] = 2;
			Felements->val[i*8+7][0] = midElements[i][vv2edge3d(3,0)]+NumCNodes;
			Felements->val[i*8+7][1] = midElements[i][vv2edge3d(3,1)]+NumCNodes;
			Felements->val[i*8+7][2] = midElements[i][vv2edge3d(0,2)]+NumCNodes;
			Felements->val[i*8+7][3] = midElements[i][vv2edge3d(2,3)]+NumCNodes;
			Felements->parent[i*8+7] = i;
			Felements->type[i*8+7] = 1;
		}
		// regular section end
	}
	for (i = 0; i<Celements->row; i++)
		free(midElements[i]);
	free(midElements);

	getTransposeOfELEMENT(Felements, FelementsTran, Felements->col, Fnodes->row);

	getFaceInfo(Felements, FelementsTran, Ffaces, Fnodes);
	getEdgeInfo3(Felements, FelementsTran, Ffaces, FelementEdge, Fedges, Fnodes);
}

/**
 * \fn void extractNondirichletMatrix(dCSRmat *A, dCSRmat *A1, ivector *dirichlet, ivector *nondirichlet)
 * \brief extract stiffness matrix by removing the corresponding dirichlet boundary condition  
 * \param *A pointer to the stiffness matrix with dirichelt boundary condition(without removed)
 * \param *A1 pointer to the stiffness matrix without dirichelt boundary condition(removed)
 * \param *dirichlet pointer to the indicator of the dirichlet boundary
 * \param *nondirichlet pointer to the indicator of the node which is not in the dirichlet boundary
 * \return void
 */
void extractNondirichletMatrix(dCSRmat *A, dCSRmat *A1, ivector *dirichlet, ivector *nondirichlet)
{
	// achiveve A1 due to dirichlet boundary condition
	int i,k,l;
	int count;
	A1->row=A->row-dirichlet->row;
	A1->col=A->col;
	A1->IA=(int*)calloc(A1->row+1, sizeof(int));
	A1->JA=NULL;
	A1->val=NULL;
	
	// form A1->IA
	for(i=0;i<A1->row;i++)
	{
		l=nondirichlet->val[i];
		A1->IA[i+1]=A->IA[l+1]-A->IA[l];
	}
	
	for(i=0;i<A1->row;i++)
		A1->IA[i+1]+=A1->IA[i];

	A1->nnz=A1->IA[A1->row];
	
	// form A1->JA, A1->val
	A1->JA=(int*)calloc(A1->nnz, sizeof(int));
	A1->val=(double*)calloc(A1->nnz, sizeof(double));
	count=0;
	for(i=0;i<A1->row;i++)
	{
		l=nondirichlet->val[i];
		for(k=A->IA[l];k<A->IA[l+1];k++)
		{
			A1->JA[count]=A->JA[k];
			A1->val[count]=A->val[k];
			count++;
		}
	}	
}

/**
 * \fn void extractNondirichletVector(dvector *b, dvector *b1, ivector *dirichlet, ivector *nondirichlet)
 * \brief extract vector by removing the corresponding dirichlet boundary condition  
 * \param *b pointer to the vector with dirichelt boundary condition(without removed)
 * \param *b1 pointer to the vector without dirichelt boundary condition(removed)
 * \param *dirichlet pointer to the indicator of the dirichlet boundary
 * \param *nondirichlet pointer to the indicator of the node which is not in the dirichlet boundary
 * \return void
 */
void extractNondirichletVector(dvector *b, dvector *b1, ivector *dirichlet, ivector *nondirichlet)
{
	int i,i1;
	b1->row=b->row-dirichlet->row;
	b1->val=(double*)calloc(b1->row, sizeof(double));	
		
	// achiveve b1 due to dirichlet boundary condition
	for(i1=0;i1<b1->row;i1++)
	{
		i=nondirichlet->val[i1];
		b1->val[i1]=b->val[i];
	}		
}

/**
* \fn void extractFreenodesVector2StokesDirichlet(dCSRmat *A, dCSRmat *B, dvector *b, dvector *b1, ELEMENT_DOF *elementDOF, dvector *uh)
* \brief extract vector according to freenodes
* \param *A pointer to left-top stiffness submatrix
* \param *B pointer to left-bottom stiffness submatrix
* \param *b pointer to the orginal vector
* \param *b1 pointer to the vector related to freenodes
* \param *elementDOF pointer to relation between elements and DOFs
* \param *uh pointer to the initial solution
* \return void
*/
void extractFreenodesVector2StokesDirichlet(dCSRmat *A, dCSRmat *B, dvector *b, dvector *b1, ELEMENT_DOF *elementDOF, dvector *uh)
{
	int i, j, k, i1, j1;

	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *nfFlag = &elementDOF->nfFlag;

	create_dvector(freenodes->row, &b1[0]);

	if (uh == NULL)
		return;

	// achieve b1[0] due to freenodes
	for (i1 = 0; i1<b1[0].row; i1++)
	{
		i = freenodes->val[i1];
		b1[0].val[i1] = b->val[i];

		for (j1 = A->IA[i]; j1 < A->IA[i + 1]; j1++)
		{
			j = A->JA[j1];
			if (nfFlag->val[j] == 1)
			{
				b1[0].val[i1] -= A->val[j1] * uh->val[j];
			}
		}
	}

	// achieve b1[1] due to freenodes
	for (i = 0; i<b1[1].row; i++)
	{
		for (j1 = B->IA[i]; j1 < B->IA[i + 1]; j1++)
		{
			j = B->JA[j1];
			if (nfFlag->val[j] == 1)
			{
				b1[1].val[i] -= B->val[j1] * uh->val[j];
			}
		}
	}
}

/**
* \fn void extractFreenodesVector(dCSRmat *A, dvector *b, dvector *b1, ELEMENT_DOF *elementDOF, dvector *uh)
* \brief extract vector according to freenodes
* \param *A pointer to left-top stiffness submatrix
* \param *b pointer to the orginal vector
* \param *b1 pointer to the vector related to freenodes
* \param *elementDOF pointer to relation between elements and DOFs
* \param *uh pointer to the initial solution
* \return void
*/
void extractFreenodesVector(dCSRmat *A, dvector *b, dvector *b1, ELEMENT_DOF *elementDOF, dvector *uh)
{
	int i, j, k, i1, j1;

	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;

	create_dvector(freenodes->row, b1);

	// achieve b1 due to freenodes
	for (i1 = 0; i1<b1->row; i1++)
	{
		i = freenodes->val[i1];
		b1->val[i1] = b->val[i];

		if (uh != NULL)
		{
			for (j1 = 0; j1 < nfreenodes->row; j1++)
			{
				j = nfreenodes->val[j1];
				for (k = A->IA[i]; k < A->IA[i + 1]; k++)
				{
					if (A->JA[k] == j)
					{
						b1->val[i1] -= A->val[k] * uh->val[j];
						break;
					}
				}
			}
		}
	}
}

/**
* \fn void updateFreenodesRHS(dCSRmat *A, dvector *b, dvector *x, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc, int flag)
* \brief update right hand side b in Ax = b when some elements of x are given
* \param *A pointer to matrix
* \param *b pointer to the right hand side
* \param *elementDOF pointer to relation between elements and DOFs
* \param *x pointer to the solution
* \return void
*/
void updateFreenodesRHS(dCSRmat *A, dvector *b, dvector *x, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc, int flag)
{
	if (x == NULL) return;

	int i, j, k, i1, j1;

	ivector *freenodesr = &elementDOFr->freenodes;
	ivector *nfreenodesr = &elementDOFr->nfreenodes;
	ivector *freenodesc = &elementDOFc->freenodes;
	ivector *nfreenodesc = &elementDOFc->nfreenodes;
	ivector *nfFlagc = &elementDOFc->nfFlag;
	
	if(flag==1){
		if(nfreenodesr->row != nfreenodesc->row) return;
	}

	// set b = x for non-freenodes
	if(flag==1){
    	for (i1 = 0; i1<nfreenodesr->row; i1++){
			i = nfreenodesr->val[i1];
			j = nfreenodesc->val[i1];
			b->val[i] = x->val[j];
		}
	}

	// update b for freenodes
	for (i1 = 0; i1<freenodesr->row; i1++){
		i = freenodesr->val[i1];
		for (j1 = A->IA[i]; j1 < A->IA[i + 1]; j1++){
			j = A->JA[j1];
			if (nfFlagc->val[j] == 1){
				b->val[i] -= A->val[j1] * x->val[j];
				break;
			}
		}
	}
	// for (i1 = 0; i1<freenodesr->row; i1++){
	// 	i = freenodesr->val[i1];
	// 	for (j1 = 0; j1 < nfreenodesc->row; j1++){
	// 		j = nfreenodesc->val[j1];
	// 		for (k = A->IA[i]; k < A->IA[i + 1]; k++){
	// 			if (A->JA[k] == j){
	// 				b->val[i] -= A->val[k] * x->val[j];
	// 				break;
	// 			}
	// 		}
	// 	}
	// }
}

/**
* \fn void updateFreenodes2bRHS(dCSRmat *A, dvector *b, dvector *x, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1)
* \brief update right hand side b in Ax = b when some elements of x are given
* \param *A pointer to matrix
* \param *b pointer to the right hand side
* \param *elementDOF pointer to relation between elements and DOFs
* \param *x pointer to the solution
* \return void
*/
void updateFreenodes2bRHS(dCSRmat *A, dvector *b, dvector *x, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1)
{

	updateFreenodesRHS(&A[0], &b[0], &x[0], elementDOF0, elementDOF0, 1);
	updateFreenodesRHS(&A[1], &b[0], &x[1], elementDOF0, elementDOF1, 0);
	updateFreenodesRHS(&A[2], &b[1], &x[0], elementDOF1, elementDOF0, 0);
	updateFreenodesRHS(&A[3], &b[1], &x[1], elementDOF1, elementDOF1, 1);
}

/**
* \fn void updateFreenodesRHS0(dCSRmat *A, dvector *b, dvector *x, ELEMENT_DOF *elementDOF)
* \brief update right hand side b in Ax = b when some elements of x are given
* \param *A pointer to matrix
* \param *b pointer to the right hand side
* \param *elementDOF pointer to relation between elements and DOFs
* \param *x pointer to the solution
* \return void
*/
void updateFreenodesRHS0(dCSRmat *A, dvector *b, dvector *x, ELEMENT_DOF *elementDOF)
{
	if (x == NULL) return;

	int i, j, k, i1, j1;

	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;

	// set b = x for non-freenodes
    for (i1 = 0; i1<nfreenodes->row; i1++){
		i = nfreenodes->val[i1];
		b->val[i] = x->val[i];
	}

	// update b for freenodes
	for (i1 = 0; i1<freenodes->row; i1++){
		i = freenodes->val[i1];
		for (j1 = 0; j1 < nfreenodes->row; j1++){
			j = nfreenodes->val[j1];
			for (k = A->IA[i]; k < A->IA[i + 1]; k++){
				if (A->JA[k] == j){
					b->val[i] -= A->val[k] * x->val[j];
					break;
				}
			}
		}
	}
}

/**
* \fn void extractFreenodesMatrix11(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc)
* \brief extract stiffness matrix A11 according to freenodes
* \param *A pointer to the original matrix
* \param *A11 pointer to the extracted matrix according to freenodes
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void extractFreenodesMatrix11(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc)
{
	// achiveve A11 due to dirichlet boundary condition
	int i, j, k, l;
	int count;

	//	ivector *nfFlagr = &elementDOFr->nfFlag;
	ivector *freenodesr = &elementDOFr->freenodes;
	//	ivector *nfreenodesr = &elementDOFr->nfreenodes;
	//	ivector *indexr = &elementDOFr->index;

	ivector *nfFlagc = &elementDOFc->nfFlag;
	ivector *freenodesc = &elementDOFc->freenodes;
	//	ivector *nfreenodesc = &elementDOFc->nfreenodes;
	ivector *indexc = &elementDOFc->index;

	A11->row = freenodesr->row;
	A11->col = freenodesc->row;
	A11->IA = (int*)calloc(A11->row + 1, sizeof(int));
	A11->JA = NULL;
	A11->val = NULL;

	// form A11->IA
	for (i = 0; i<A11->row; i++)
	{
		l = freenodesr->val[i];
		for (k = A->IA[l]; k<A->IA[l + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlagc->val[j] == 0)
				A11->IA[i + 1]++;
		}
	}

	for (i = 0; i<A11->row; i++)
		A11->IA[i + 1] += A11->IA[i];

	A11->nnz = A11->IA[A11->row];

	// form A11->JA, A11->val
	A11->JA = (int*)calloc(A11->nnz, sizeof(int));
	A11->val = (double*)calloc(A11->nnz, sizeof(double));
	count = 0;
	for (i = 0; i<A11->row; i++)
	{
		l = freenodesr->val[i];
		for (k = A->IA[l]; k<A->IA[l + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlagc->val[j] == 0)
			{
				A11->JA[count] = indexc->val[j];
				A11->val[count] = A->val[k];
				count++;
			}
		}
	}
}

/**
* \fn void updateFreenodesMatrix11(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc, int flag)
* \brief extract stiffness matrix A11 according to freenodes
* \param *A pointer to the original matrix
* \param *A11 pointer to the extracted matrix according to freenodes
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void updateFreenodesMatrix11(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc, int flag)
{
	// achiveve A11 due to dirichlet boundary condition
	int i, j, k, l;
	int count;

	ivector *nfFlagr = &elementDOFr->nfFlag;
	ivector *freenodesr = &elementDOFr->freenodes;
	ivector *nfreenodesr = &elementDOFr->nfreenodes;
	//	ivector *indexr = &elementDOFr->index;

	ivector *nfFlagc = &elementDOFc->nfFlag;
	// ivector *freenodesc = &elementDOFc->freenodes;
	//	ivector *nfreenodesc = &elementDOFc->nfreenodes;
	// ivector *indexc = &elementDOFc->index;

	// A11->row = freenodesr->row;
	// A11->col = freenodesc->row;
	A11->row = A->row;
	A11->col = A->col;
	A11->IA = (int*)calloc(A11->row + 1, sizeof(int));
	A11->JA = NULL;
	A11->val = NULL;

	// form A11->IA
	if(flag==1){
		for (i = 0; i < nfreenodesr->row; i++){
			l = nfreenodesr->val[i];
			A11->IA[l + 1] = 1;
		}
	}
	for (i = 0; i < freenodesr->row; i++){
		l = freenodesr->val[i];
		for (k = A->IA[l]; k<A->IA[l + 1]; k++){
			j = A->JA[k];
			if (nfFlagc->val[j] == 0)
				A11->IA[l + 1]++;
		}
	}

	for (i = 0; i<A11->row; i++)
		A11->IA[i + 1] += A11->IA[i];

	A11->nnz = A11->IA[A11->row];

	// form A11->JA, A11->val
	A11->JA = (int*)calloc(A11->nnz, sizeof(int));
	A11->val = (double*)calloc(A11->nnz, sizeof(double));
	count = 0;
	for (i = 0; i<A11->row; i++){
		if(nfFlagr->val[i] == 1){
			if(flag==1){
				A11->JA[count] = i;
				A11->val[count] = 1.0;
				count++;
			}
		}
		else{
			for (k = A->IA[i]; k<A->IA[i + 1]; k++)
			{
				j = A->JA[k];
				if (nfFlagc->val[j] == 0)
				{
					A11->JA[count] = j;
					A11->val[count] = A->val[k];
					count++;
				}
			}		
		}
	}
}

/**
* \fn void updateFreenodes2bMatrix11(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1)
* \brief extract stiffness matrix A11 according to freenodes
* \param *A pointer to the original matrix
* \param *A11 pointer to the extracted matrix according to freenodes
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void updateFreenodes2bMatrix11(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1)
{
	updateFreenodesMatrix11(&A[0], &A11[0], elementDOF0, elementDOF0, 1);
	updateFreenodesMatrix11(&A[1], &A11[1], elementDOF0, elementDOF1, 0);
	updateFreenodesMatrix11(&A[2], &A11[2], elementDOF1, elementDOF0, 0);
	updateFreenodesMatrix11(&A[3], &A11[3], elementDOF1, elementDOF1, 1);
}

/**
* \fn void extractFreenodesMatrix1r(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF)
* \brief extract stiffness matrix A1 according to freenodes
* \param *A pointer to the original matrix
* \param *A1 pointer to the extracted matrix according to freenodes
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void extractFreenodesMatrix1r(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF)
{
	// achiveve A1 due to dirichlet boundary condition
	int i, j, k, l;
	int count;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	A1->row = freenodes->row;
	A1->col = A->col;
	A1->IA = (int*)calloc(A1->row + 1, sizeof(int));
	A1->JA = NULL;
	A1->val = NULL;

	// form A1->IA
	for (i = 0; i<A1->row; i++)
	{
		l = freenodes->val[i];
		A1->IA[i + 1] = A->IA[l + 1] - A->IA[l];
	}

	for (i = 0; i<A1->row; i++)
		A1->IA[i + 1] += A1->IA[i];

	A1->nnz = A1->IA[A1->row];

	// form A1->JA, A1->val
	A1->JA = (int*)calloc(A1->nnz, sizeof(int));
	A1->val = (double*)calloc(A1->nnz, sizeof(double));
	count = 0;
	for (i = 0; i<A1->row; i++)
	{
		l = freenodes->val[i];
		for (k = A->IA[l]; k<A->IA[l + 1]; k++)
		{
			A1->JA[count] = A->JA[k];
			A1->val[count] = A->val[k];
			count++;
		}
	}
}

/**
* \fn void extractFreenodesMatrix1c(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF)
* \brief extract stiffness matrix A1 according to freenodes
* \param *A pointer to the original matrix
* \param *A1 pointer to the extracted matrix according to freenodes
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void extractFreenodesMatrix1c(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF)
{
	// achiveve A1 due to dirichlet boundary condition
	int i, j, k;
	int count;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	A1->row = A->row;
	A1->col = freenodes->row;
	A1->IA = (int*)calloc(A1->row + 1, sizeof(int));
	A1->JA = NULL;
	A1->val = NULL;

	// form A1->IA
	for (i = 0; i<A1->row; i++)
	{
		for (k = A->IA[i]; k<A->IA[i + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlag->val[j] == 0)
				A1->IA[i + 1]++;
		}
	}

	for (i = 0; i<A1->row; i++)
		A1->IA[i + 1] += A1->IA[i];

	A1->nnz = A1->IA[A1->row];

	// form A1->JA, A1->val
	A1->JA = (int*)calloc(A1->nnz, sizeof(int));
	A1->val = (double*)calloc(A1->nnz, sizeof(double));
	count = 0;
	for (i = 0; i<A1->row; i++)
	{
		for (k = A->IA[i]; k<A->IA[i + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlag->val[j] == 0)
			{
				A1->JA[count] = index->val[j];
				A1->val[count] = A->val[k];
				count++;
			}
		}
	}
}

/**
* \fn void extractFreenodesMatrix1cBlock(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF)
* \brief extract stiffness matrix A1 according to freenodes
* \param *A pointer to the original matrix
* \param *A1 pointer to the extracted matrix according to freenodes
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void extractFreenodesMatrix1cBlock(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF)
{
	// achiveve A1 according to freenodes
	int i, j, k;
	int count;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	A1->row = A->row;
	A1->col = freenodes->row * 2;
	A1->IA = (int*)calloc(A1->row + 1, sizeof(int));
	A1->JA = NULL;
	A1->val = NULL;

	int col = A->col / 2;

	// form A1->IA
	for (i = 0; i<A1->row; i++)
	{
		for (k = A->IA[i]; k<A->IA[i + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlag->val[j%col] == 0)
				A1->IA[i + 1]++;
		}
	}

	for (i = 0; i<A1->row; i++)
		A1->IA[i + 1] += A1->IA[i];

	A1->nnz = A1->IA[A1->row];

	// form A1->JA, A1->val
	A1->JA = (int*)calloc(A1->nnz, sizeof(int));
	A1->val = (double*)calloc(A1->nnz, sizeof(double));
	count = 0;
	for (i = 0; i<A1->row; i++)
	{
		for (k = A->IA[i]; k<A->IA[i + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlag->val[j%col] == 0)
			{
				A1->JA[count] = index->val[j%col] + freenodes->row*(j / col);
				A1->val[count] = A->val[k];
				count++;
			}
		}
	}
}

/**
* \fn void extractFreenodesMatrix11cBlock(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc)
* \brief extract stiffness matrix A11 according to freenodes
* \param *A pointer to the original matrix
* \param *A1 pointer to the extracted matrix according to freenodes
* \param *elementDOFr pointer to relation between elements and DOFs in rows
* \param *elementDOFc pointer to relation between elements and DOFs in columns, blockwise
* \return void
*/
void extractFreenodesMatrix11cBlock(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc)
{
	// achiveve A1 according to freenodes
	int i, j, k, l;
	int count;

	//	ivector *nfFlagr = &elementDOFr->nfFlag;
	ivector *freenodesr = &elementDOFr->freenodes;
	//	ivector *nfreenodesr = &elementDOFr->nfreenodes;
	//	ivector *indexr = &elementDOFr->index;

	ivector *nfFlagc = &elementDOFc->nfFlag;
	ivector *freenodesc = &elementDOFc->freenodes;
	//	ivector *nfreenodesc = &elementDOFc->nfreenodes;
	ivector *indexc = &elementDOFc->index;

	A11->row = freenodesr->row;
	A11->col = freenodesc->row * 2;
	A11->IA = (int*)calloc(A11->row + 1, sizeof(int));
	A11->JA = NULL;
	A11->val = NULL;

	int col = A->col / 2;

	// form A11->IA
	for (i = 0; i<A11->row; i++)
	{
		l = freenodesr->val[i];
		for (k = A->IA[l]; k<A->IA[l + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlagc->val[j%col] == 0)
				A11->IA[i + 1]++;
		}
	}

	for (i = 0; i<A11->row; i++)
		A11->IA[i + 1] += A11->IA[i];

	A11->nnz = A11->IA[A11->row];

	// form A11->JA, A11->val
	A11->JA = (int*)calloc(A11->nnz, sizeof(int));
	A11->val = (double*)calloc(A11->nnz, sizeof(double));
	count = 0;
	for (i = 0; i<A11->row; i++)
	{
		l = freenodesr->val[i];
		for (k = A->IA[l]; k<A->IA[l + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlagc->val[j%col] == 0)
			{
				A11->JA[count] = indexc->val[j%col] + freenodesc->row*(j / col);
				A11->val[count] = A->val[k];
				count++;
			}
		}
	}
}

/**
* \fn void extractFreenodesMatrix11cBlock3d(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc)
* \brief extract stiffness matrix A11 according to freenodes
* \param *A pointer to the original matrix
* \param *A1 pointer to the extracted matrix according to freenodes
* \param *elementDOFr pointer to relation between elements and DOFs in rows
* \param *elementDOFc pointer to relation between elements and DOFs in columns, blockwise
* \return void
*/
void extractFreenodesMatrix11cBlock3d(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc)
{
	// achiveve A1 according to freenodes
	int i, j, k, l;
	int count;

	//	ivector *nfFlagr = &elementDOFr->nfFlag;
	ivector *freenodesr = &elementDOFr->freenodes;
	//	ivector *nfreenodesr = &elementDOFr->nfreenodes;
	//	ivector *indexr = &elementDOFr->index;

	ivector *nfFlagc = &elementDOFc->nfFlag;
	ivector *freenodesc = &elementDOFc->freenodes;
	//	ivector *nfreenodesc = &elementDOFc->nfreenodes;
	ivector *indexc = &elementDOFc->index;

	A11->row = freenodesr->row;
	A11->col = freenodesc->row * 3;
	A11->IA = (int*)calloc(A11->row + 1, sizeof(int));
	A11->JA = NULL;
	A11->val = NULL;

	int col = A->col / 3;

	// form A11->IA
	for (i = 0; i<A11->row; i++)
	{
		l = freenodesr->val[i];
		for (k = A->IA[l]; k<A->IA[l + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlagc->val[j%col] == 0)
				A11->IA[i + 1]++;
		}
	}

	for (i = 0; i<A11->row; i++)
		A11->IA[i + 1] += A11->IA[i];

	A11->nnz = A11->IA[A11->row];

	// form A11->JA, A11->val
	A11->JA = (int*)calloc(A11->nnz, sizeof(int));
	A11->val = (double*)calloc(A11->nnz, sizeof(double));
	count = 0;
	for (i = 0; i<A11->row; i++)
	{
		l = freenodesr->val[i];
		for (k = A->IA[l]; k<A->IA[l + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlagc->val[j%col] == 0)
			{
				A11->JA[count] = indexc->val[j%col] + freenodesc->row*(j / col);
				A11->val[count] = A->val[k];
				count++;
			}
		}
	}
}

/**
 * \fn int getEdgeDOFsTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row31, int *row32)
 * \brief generate JA of A from the DOFs of edge in element 
 * \param *A pointer to stiffness matrix
 * \param count current index of A->JA
 * \param element index of current element
 * \param edge index of current edge
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *rowstart the starting location of i-th row, i.e. A->IA[i] 
 * \param *row31 1/3 of A->IA[i+1] - A->IA[i] 
 * \param *row32 2/3 of A->IA[i+1] - A->IA[i] 
 * \return count, the current index of A->JA
 */
int getEdgeDOFsTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row31, int *row32)
{
	int i,j,l;
	int node;
	
	for(l=0;l<elementEdge->col;l++)
	{
		if(elementEdge->val[element][l]==edge)
			break;
	}
	
	node=elementDOF->val[element][(l+1)%3];
	for(j=0;j<3;j++)
	{
		A->JA[rowstart[j]+count]=node;
		A->JA[rowstart[j]+count+row31[j]]=node+elementDOF->dof;
		A->JA[rowstart[j]+count+row32[j]]=node+elementDOF->dof*2;
	}
	count++;
	
	node=elementDOF->val[element][(l+2)%3];
	for(j=0;j<3;j++)
	{
		A->JA[rowstart[j]+count]=node;
		A->JA[rowstart[j]+count+row31[j]]=node+elementDOF->dof;
		A->JA[rowstart[j]+count+row32[j]]=node+elementDOF->dof*2;
	}
	count++;
	
	for(i=0;i<elementDOF->dop-1;i++)
	{
		node=elementDOF->val[element][3+l*(elementDOF->dop-1)+i];
		for(j=0;j<3;j++)
		{
			A->JA[rowstart[j]+count]=node;
			A->JA[rowstart[j]+count+row31[j]]=node+elementDOF->dof;
			A->JA[rowstart[j]+count+row32[j]]=node+elementDOF->dof*2;
		}
		count++;
	}
	
	return count;
}

/**
 * \fn int getEdgeDOFsScalarTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int rowstart, int row31, int row32)
 * \brief generate JA of A from the DOFs of edge in element 
 * \param *A pointer to stiffness matrix
 * \param count current index of A->JA
 * \param element index of current element
 * \param edge index of current edge
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *rowstart the starting location of i-th row, i.e. A->IA[i] 
 * \param *row31 1/3 of A->IA[i+1] - A->IA[i] 
 * \param *row32 2/3 of A->IA[i+1] - A->IA[i] 
 * \return count, the current index of A->JA
 */
int getEdgeDOFsScalarTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int rowstart, int row31, int row32)
{
	int i,l;
	int node, taustart;

	if(elementDOF->dop==0)
	{
		taustart=3*element;
		A->JA[rowstart+count]=taustart;
		A->JA[rowstart+count+row31]=taustart+1;
		A->JA[rowstart+count+row32]=taustart+2;
		count++;
		return count;
	}
	
	for(l=0;l<elementEdge->col;l++)
	{
		if(elementEdge->val[element][l]==edge)
			break;
	}
	
	taustart=3*elementDOF->col*element;
	node=(l+1)%3;
	A->JA[rowstart+count]=taustart+node;
	A->JA[rowstart+count+row31]=taustart+node+elementDOF->col;
	A->JA[rowstart+count+row32]=taustart+node+elementDOF->col*2;
	count++;
	
	node=(l+2)%3;
	A->JA[rowstart+count]=taustart+node;
	A->JA[rowstart+count+row31]=taustart+node+elementDOF->col;
	A->JA[rowstart+count+row32]=taustart+node+elementDOF->col*2;
	count++;
	
	for(i=0;i<elementDOF->dop-1;i++)
	{
		node=3+l*(elementDOF->dop-1)+i;
		A->JA[rowstart+count]=taustart+node;
		A->JA[rowstart+count+row31]=taustart+node+elementDOF->col;
		A->JA[rowstart+count+row32]=taustart+node+elementDOF->col*2;
		count++;
	}
	
	return count;
}

/**
 * \fn int getEdgeDOFsVectorTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row31, int *row32)
 * \brief generate JA of A from the DOFs of edge in element 
 * \param *A pointer to stiffness matrix
 * \param count current index of A->JA
 * \param element index of current element
 * \param edge index of current edge
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *rowstart the starting location of i-th row, i.e. A->IA[i] 
 * \param *row31 1/3 of A->IA[i+1] - A->IA[i] 
 * \param *row32 2/3 of A->IA[i+1] - A->IA[i] 
 * \return count, the current index of A->JA
 */
int getEdgeDOFsVectorTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row31, int *row32)
{
	int i,j,l;
	int node;

	if(elementDOF->dop==0)
	{
		node=elementDOF->val[element][0];
		for(j=0;j<2;j++)
		{
			A->JA[rowstart[j]+count]=node;
			A->JA[rowstart[j]+count+row31[j]]=node+elementDOF->dof;
			A->JA[rowstart[j]+count+row32[j]]=node+elementDOF->dof*2;
		}
		count++;
		return count;
	}
	
	for(l=0;l<elementEdge->col;l++)
	{
		if(elementEdge->val[element][l]==edge)
			break;
	}
	
	node=elementDOF->val[element][(l+1)%3];
	for(j=0;j<2;j++)
	{
		A->JA[rowstart[j]+count]=node;
		A->JA[rowstart[j]+count+row31[j]]=node+elementDOF->dof;
		A->JA[rowstart[j]+count+row32[j]]=node+elementDOF->dof*2;
	}
	count++;
	
	node=elementDOF->val[element][(l+2)%3];
	for(j=0;j<2;j++)
	{
		A->JA[rowstart[j]+count]=node;
		A->JA[rowstart[j]+count+row31[j]]=node+elementDOF->dof;
		A->JA[rowstart[j]+count+row32[j]]=node+elementDOF->dof*2;
	}
	count++;
	
	for(i=0;i<elementDOF->dop-1;i++)
	{
		node=elementDOF->val[element][3+l*(elementDOF->dop-1)+i];
		for(j=0;j<2;j++)
		{
			A->JA[rowstart[j]+count]=node;
			A->JA[rowstart[j]+count+row31[j]]=node+elementDOF->dof;
			A->JA[rowstart[j]+count+row32[j]]=node+elementDOF->dof*2;
		}
		count++;
	}
	
	return count;
}

/**
 * \fn int getEdgeDOFsVector(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row21)
 * \brief generate JA of A from the DOFs of edge in element 
 * \param *A pointer to stiffness matrix
 * \param count current index of A->JA
 * \param element index of current element
 * \param edge index of current edge
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *rowstart the starting location of i-th row, i.e. A->IA[i] 
 * \param *row21 1/2 of A->IA[i+1] - A->IA[i] 
 * \return count, the current index of A->JA
 */
int getEdgeDOFsVector(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row21)
{
	int i,j,l;
	int node;
	
	for(l=0;l<elementEdge->col;l++)
	{
		if(elementEdge->val[element][l]==edge)
			break;
	}
	
	node=2*element*elementDOF->col + (l+1)%3;
	for(j=0;j<2;j++)
	{
		A->JA[rowstart[j]+count]=node;
		A->JA[rowstart[j]+count+row21[j]]=node+elementDOF->col;
	}
	count++;
	
	node=2*element*elementDOF->col + (l+2)%3;
	for(j=0;j<2;j++)
	{
		A->JA[rowstart[j]+count]=node;
		A->JA[rowstart[j]+count+row21[j]]=node+elementDOF->col;
	}
	count++;
	
	for(i=0;i<elementDOF->dop-1;i++)
	{
		node=2*element*elementDOF->col + 3+l*(elementDOF->dop-1)+i;
		for(j=0;j<2;j++)
		{
			A->JA[rowstart[j]+count]=node;
			A->JA[rowstart[j]+count+row21[j]]=node+elementDOF->col;
		}
		count++;
	}
	
	return count;
}

/**
 * \fn int getFaceDOFsVector(dCSRmat *A, int count, int element, int face, idenmat *elementFace, ELEMENT_DOF *elementDOF, int *rowstart, int *row31)
 * \brief generate JA of A from the DOFs of face in element 
 * \param *A pointer to stiffness matrix
 * \param count current index of A->JA
 * \param element index of current element
 * \param face index of current face
 * \param *elementFace pointer to relation between tetrahedrons and faces: each row stores 4 faces index
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *rowstart the starting location of i-th row, i.e. A->IA[i] 
 * \param *row31 1/3 of A->IA[i+1] - A->IA[i] 
 * \return count, the current index of A->JA
 */
int getFaceDOFsVector(dCSRmat *A, int count, int element, int face, idenmat *elementFace, ELEMENT_DOF *elementDOF, int *rowstart, int *row31)
{
	int l;
	
	for(l=0;l<elementFace->col;l++)
	{
		if(elementFace->val[element][l]==face)
			break;
	}

	count=getFaceDOFsVectorlocalFace(A, count, element, l, elementDOF, rowstart, row31);
	
	return count;
}

/**
* \fn int getFaceDOFsVectorOld(dCSRmat *A, int count, int element, int face, idenmat *elementFace, ELEMENT_DOF *elementDOF, int *rowstart, int *row31)
* \brief generate JA of A from the DOFs of face in element
* \param *A pointer to stiffness matrix
* \param count current index of A->JA
* \param element index of current element
* \param face index of current face
* \param *elementFace pointer to relation between tetrahedrons and faces: each row stores 4 faces index
* \param *elementDOF pointer to relation between elements and DOFs
* \param *rowstart the starting location of i-th row, i.e. A->IA[i]
* \param *row31 1/3 of A->IA[i+1] - A->IA[i]
* \return count, the current index of A->JA
*/
int getFaceDOFsVectorOld(dCSRmat *A, int count, int element, int face, idenmat *elementFace, ELEMENT_DOF *elementDOF, int *rowstart, int *row31)
{
	int l;

	for (l = 0; l<elementFace->col; l++)
	{
		if (elementFace->val[element][l] == face)
			break;
	}

	count = getFaceDOFsVectorlocalFaceOld(A, count, element, l, elementDOF, rowstart, row31);

	return count;
}

/**
 * \fn int getFaceDOFsVectorlocalFace(dCSRmat *A, int count, int element, int lface, ELEMENT_DOF *elementDOF, int *rowstart, int *row31)
 * \brief generate JA of A from the DOFs of face in element 
 * \param *A pointer to stiffness matrix
 * \param count current index of A->JA
 * \param element index of current element
 * \param lface local index of current face
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *rowstart the starting location of i-th row, i.e. A->IA[i] 
 * \param *row31 1/3 of A->IA[i+1] - A->IA[i] 
 * \return count, the current index of A->JA
 */
int getFaceDOFsVectorlocalFace(dCSRmat *A, int count, int element, int lface, ELEMENT_DOF *elementDOF, int *rowstart, int *row31)
{
	int i,j,l;
	int node, edge, vertices[3];

	l=lface;

	for(i=0;i<3;i++)
		vertices[i]=(l+1+i)%4;
	
	for(i=0;i<3;i++) // 3 vertices
	{
		node = elementDOF->val[element][vertices[i]];
		for(j=0;j<3;j++)
		{
			A->JA[rowstart[j]+count]=node;
			A->JA[rowstart[j]+count+row31[j]]=node+elementDOF->dof;
			A->JA[rowstart[j]+count+row31[j]*2]=node+elementDOF->dof*2;
		}
		count++;
	}

	// interior points in 3 edges
	edge=vv2edge3d(vertices[0], vertices[1]);
	for(i=0;i<elementDOF->dop-1;i++)
	{
		node = elementDOF->val[element][4 + edge*(elementDOF->dop - 1) + i];
		for(j=0;j<3;j++)
		{
			A->JA[rowstart[j]+count]=node;
			A->JA[rowstart[j]+count+row31[j]]=node+elementDOF->dof;
			A->JA[rowstart[j]+count+row31[j]*2]=node+elementDOF->dof * 2;
		}
		count++;
	}
	edge=vv2edge3d(vertices[0], vertices[2]);
	for(i=0;i<elementDOF->dop-1;i++)
	{
		node = elementDOF->val[element][4 + edge*(elementDOF->dop - 1) + i];
		for(j=0;j<3;j++)
		{
			A->JA[rowstart[j]+count]=node;
			A->JA[rowstart[j]+count+row31[j]]=node+elementDOF->dof;
			A->JA[rowstart[j]+count+row31[j]*2]=node+elementDOF->dof*2;
		}
		count++;
	}
	edge=vv2edge3d(vertices[1], vertices[2]);
	for(i=0;i<elementDOF->dop-1;i++)
	{
		node = elementDOF->val[element][4 + edge*(elementDOF->dop - 1) + i];
		for(j=0;j<3;j++)
		{
			A->JA[rowstart[j]+count]=node;
			A->JA[rowstart[j]+count+row31[j]]=node+elementDOF->dof;
			A->JA[rowstart[j]+count+row31[j]*2]=node+elementDOF->dof*2;
		}
		count++;
	}

	// interior points in face
	for(i=0;i<(elementDOF->dop-1)*(elementDOF->dop-2)/2;i++)
	{
		node = elementDOF->val[element][6 * elementDOF->dop - 2 + l*(elementDOF->dop - 1)*(elementDOF->dop - 2) / 2 + i];
		for(j=0;j<3;j++)
		{
			A->JA[rowstart[j]+count]=node;
			A->JA[rowstart[j]+count+row31[j]]=node+elementDOF->dof;
			A->JA[rowstart[j]+count+row31[j]*2]=node+elementDOF->dof*2;
		}
		count++;
	}
	
	return count;
}

/**
* \fn int getFaceDOFsVectorlocalFaceOld(dCSRmat *A, int count, int element, int lface, ELEMENT_DOF *elementDOF, int *rowstart, int *row31)
* \brief generate JA of A from the DOFs of face in element
* \param *A pointer to stiffness matrix
* \param count current index of A->JA
* \param element index of current element
* \param lface local index of current face
* \param *elementDOF pointer to relation between elements and DOFs
* \param *rowstart the starting location of i-th row, i.e. A->IA[i]
* \param *row31 1/3 of A->IA[i+1] - A->IA[i]
* \return count, the current index of A->JA
*/
int getFaceDOFsVectorlocalFaceOld(dCSRmat *A, int count, int element, int lface, ELEMENT_DOF *elementDOF, int *rowstart, int *row31)
{
	int i, j, l;
	int node, edge, vertices[3];

	l = lface;

	for (i = 0; i<3; i++)
		vertices[i] = (l + 1 + i) % 4;

	for (i = 0; i<3; i++) // 3 vertices
	{
		node = 3 * element*elementDOF->col + vertices[i];
		for (j = 0; j<3; j++)
		{
			A->JA[rowstart[j] + count] = node;
			A->JA[rowstart[j] + count + row31[j]] = node + elementDOF->col;
			A->JA[rowstart[j] + count + row31[j] * 2] = node + elementDOF->col * 2;
		}
		count++;
	}

	// interior points in 3 edges
	edge = vv2edge3d(vertices[0], vertices[1]);
	for (i = 0; i<elementDOF->dop - 1; i++)
	{
		node = 3 * element*elementDOF->col + 4 + edge*(elementDOF->dop - 1) + i;
		for (j = 0; j<3; j++)
		{
			A->JA[rowstart[j] + count] = node;
			A->JA[rowstart[j] + count + row31[j]] = node + elementDOF->col;
			A->JA[rowstart[j] + count + row31[j] * 2] = node + elementDOF->col * 2;
		}
		count++;
	}
	edge = vv2edge3d(vertices[0], vertices[2]);
	for (i = 0; i<elementDOF->dop - 1; i++)
	{
		node = 3 * element*elementDOF->col + 4 + edge*(elementDOF->dop - 1) + i;
		for (j = 0; j<3; j++)
		{
			A->JA[rowstart[j] + count] = node;
			A->JA[rowstart[j] + count + row31[j]] = node + elementDOF->col;
			A->JA[rowstart[j] + count + row31[j] * 2] = node + elementDOF->col * 2;
		}
		count++;
	}
	edge = vv2edge3d(vertices[1], vertices[2]);
	for (i = 0; i<elementDOF->dop - 1; i++)
	{
		node = 3 * element*elementDOF->col + 4 + edge*(elementDOF->dop - 1) + i;
		for (j = 0; j<3; j++)
		{
			A->JA[rowstart[j] + count] = node;
			A->JA[rowstart[j] + count + row31[j]] = node + elementDOF->col;
			A->JA[rowstart[j] + count + row31[j] * 2] = node + elementDOF->col * 2;
		}
		count++;
	}

	// interior points in face
	for (i = 0; i<(elementDOF->dop - 1)*(elementDOF->dop - 2) / 2; i++)
	{
		node = 3 * element*elementDOF->col + 6 * elementDOF->dop - 2 + l*(elementDOF->dop - 1)*(elementDOF->dop - 2) / 2 + i;
		for (j = 0; j<3; j++)
		{
			A->JA[rowstart[j] + count] = node;
			A->JA[rowstart[j] + count + row31[j]] = node + elementDOF->col;
			A->JA[rowstart[j] + count + row31[j] * 2] = node + elementDOF->col * 2;
		}
		count++;
	}

	return count;
}

/**
 * \fn int getFaceDOFs(int *A, int count, int element, int face, idenmat *elementFace, ELEMENT_DOF *elementDOF)
 * \brief generate A from the DOFs of face in element 
 * \param *A pointer to patchnodes related to the face
 * \param count current index of A
 * \param element index of current element
 * \param face index of current face
 * \param *elementFace pointer to relation between tetrahedrons and faces: each row stores 4 faces index
 * \param *elementDOF pointer to relation between elements and DOFs
 * \return count, the current index of A
 */
int getFaceDOFs(int *A, int count, int element, int face, idenmat *elementFace, ELEMENT_DOF *elementDOF)
{
	int i,l;
	int edge, vertices[3];

	for(l=0;l<elementFace->col;l++)
	{
		if(elementFace->val[element][l]==face)
			break;
	}

	for(i=0;i<3;i++)
		vertices[i]=(l+1+i)%4;
	
	for(i=0;i<3;i++) // 3 vertices
	{
		A[count]=elementDOF->val[element][vertices[i]];
		count++;
	}

	// interior points in 3 edges
	edge=vv2edge3d(vertices[0], vertices[1]);
	for(i=0;i<elementDOF->dop-1;i++)
	{
		A[count]=elementDOF->val[element][4+edge*(elementDOF->dop-1)+i];
		count++;
	}
	edge=vv2edge3d(vertices[0], vertices[2]);
	for(i=0;i<elementDOF->dop-1;i++)
	{
		A[count]=elementDOF->val[element][4+edge*(elementDOF->dop-1)+i];
		count++;
	}
	edge=vv2edge3d(vertices[1], vertices[2]);
	for(i=0;i<elementDOF->dop-1;i++)
	{
		A[count]=elementDOF->val[element][4+edge*(elementDOF->dop-1)+i];
		count++;
	}

	// interior points in face
	for(i=0;i<(elementDOF->dop-1)*(elementDOF->dop-2)/2;i++)
	{
		A[count]=elementDOF->val[element][6*elementDOF->dop-2 +l*(elementDOF->dop-1)*(elementDOF->dop-2)/2+i];
		count++;
	}
	
	return count;
}


/**
 * \fn void patchtogether22(dCSRmat *A11, dCSRmat *A12, dCSRmat *A21, dCSRmat *A22, dCSRmat *A)
 * \brief concatenate submatrices A11, A12, A21, A22 to matrix A by A=[A11 A12 
																	   A21 A22]
 * \param *A11 pointer to the submatrix located at the top-left corner
 * \param *A12 pointer to the submatrix located at the top-right corner
 * \param *A21 pointer to the submatrix located at the bottom-left corner
 * \param *A22 pointer to the submatrix located at the bottom-right corner
 * \param *A pointer to the resulted matrix
 * \return void
 */
void patchtogether22(dCSRmat *A11, dCSRmat *A12, dCSRmat *A21, dCSRmat *A22, dCSRmat *A)
{
	if((A11->row!=A12->row) || (A21->row!=A22->row) || (A11->col!=A21->col) || (A12->col!=A22->col))
	{
		printf("The dimensions of four submatrices do not match in function patchtogether!\n");
		exit(0);
	}

	int i,k,count;
	create_csr_matrix(A11->row+A21->row, A11->col+A12->col, A11->nnz+A12->nnz+A21->nnz+A22->nnz, A);
	
	// form A->IA
	for(i=0;i<A11->row;i++)
		A->IA[i+1]=A11->IA[i+1]-A11->IA[i]+A12->IA[i+1]-A12->IA[i];

	for(i=0;i<A21->row;i++)
		A->IA[A11->row+i+1]=A21->IA[i+1]-A21->IA[i]+A22->IA[i+1]-A22->IA[i];
	
	for(i=0;i<A->row;i++)
		A->IA[i+1]+=A->IA[i];

	// form A->JA, A->val
	count=0;
	for(i=0;i<A11->row;i++)
	{
		for(k=A11->IA[i];k<A11->IA[i+1];k++)
		{
			A->JA[count]=A11->JA[k];
			A->val[count]=A11->val[k];
			count++;
		}
		for(k=A12->IA[i];k<A12->IA[i+1];k++)
		{
			A->JA[count]=A12->JA[k]+A11->col;
			A->val[count]=A12->val[k];
			count++;
		}
	}

	for(i=0;i<A21->row;i++)
	{
		for(k=A21->IA[i];k<A21->IA[i+1];k++)
		{
			A->JA[count]=A21->JA[k];
			A->val[count]=A21->val[k];
			count++;
		}
		for(k=A22->IA[i];k<A22->IA[i+1];k++)
		{
			A->JA[count]=A22->JA[k]+A11->col;
			A->val[count]=A22->val[k];
			count++;
		}
	}	
}

/**
 * \fn void patchtogether33(dCSRmat *A11, dCSRmat *A12, dCSRmat *A13, dCSRmat *A21, dCSRmat *A22, dCSRmat *A23, dCSRmat *A31, dCSRmat *A32, dCSRmat *A33, dCSRmat *A)
 * \brief concatenate submatrices A11, A12, A13, A21, A22, A23, A31, A32, A33 to matrix A by 
 *                                                 A=[A11 A12 A13
 *  												  A21 A22 A23
 *											          A31 A32 A33]
 * \param *A11 pointer to the submatrix located at the top-left corner
 * \param *A12 pointer to the submatrix located at the top-middle corner
 * \param *A13 pointer to the submatrix located at the top-right corner
 * \param *A21 pointer to the submatrix located at the middle-left corner
 * \param *A22 pointer to the submatrix located at the middle-middle corner
 * \param *A23 pointer to the submatrix located at the middle-right corner
 * \param *A31 pointer to the submatrix located at the bottom-left corner
 * \param *A32 pointer to the submatrix located at the bottom-middle corner
 * \param *A33 pointer to the submatrix located at the bottom-right corner
 * \param *A pointer to the resulted matrix
 * \return void
 */
void patchtogether33(dCSRmat *A11, dCSRmat *A12, dCSRmat *A13, dCSRmat *A21, dCSRmat *A22, dCSRmat *A23, dCSRmat *A31, dCSRmat *A32, dCSRmat *A33, dCSRmat *A)
{
	if((A11->row!=A12->row) || (A11->row!=A13->row) || (A21->row!=A22->row) || (A21->row!=A23->row) || (A31->row!=A32->row) || (A31->row!=A33->row))
	{
		printf("The row dimensions of nine submatrices do not match in function patchtogether!\n");
		exit(0);
	}
	
	if((A11->col!=A21->col) || (A11->col!=A31->col) || (A12->col!=A22->col) || (A12->col!=A32->col) || (A13->col!=A23->col) || (A13->col!=A33->col))
	{
		printf("The column dimensions of nine submatrices do not match in function patchtogether!\n");
		exit(0);
	}

	int i,k,count;
	create_csr_matrix(A11->row+A21->row+A31->row, A11->col+A12->col+A13->col, A11->nnz+A12->nnz+A13->nnz+A21->nnz+A22->nnz+A23->nnz+A31->nnz+A32->nnz+A33->nnz, A);
	
	// form A->IA
	for(i=0;i<A11->row;i++)
		A->IA[i+1]=A11->IA[i+1]-A11->IA[i]+A12->IA[i+1]-A12->IA[i]+A13->IA[i+1]-A13->IA[i];

	for(i=0;i<A21->row;i++)
		A->IA[A11->row+i+1]=A21->IA[i+1]-A21->IA[i]+A22->IA[i+1]-A22->IA[i]+A23->IA[i+1]-A23->IA[i];
	
	for(i=0;i<A31->row;i++)
		A->IA[A11->row+A21->row+i+1]=A31->IA[i+1]-A31->IA[i]+A32->IA[i+1]-A32->IA[i]+A33->IA[i+1]-A33->IA[i];

	for(i=0;i<A->row;i++)
		A->IA[i+1]+=A->IA[i];

	// form A->JA, A->val
	count=0;
	for(i=0;i<A11->row;i++)
	{
		for(k=A11->IA[i];k<A11->IA[i+1];k++)
		{
			A->JA[count]=A11->JA[k];
			A->val[count]=A11->val[k];
			count++;
		}
		for(k=A12->IA[i];k<A12->IA[i+1];k++)
		{
			A->JA[count]=A12->JA[k]+A11->col;
			A->val[count]=A12->val[k];
			count++;
		}
		for(k=A13->IA[i];k<A13->IA[i+1];k++)
		{
			A->JA[count]=A13->JA[k]+A11->col+A12->col;
			A->val[count]=A13->val[k];
			count++;
		}
	}

	for(i=0;i<A21->row;i++)
	{
		for(k=A21->IA[i];k<A21->IA[i+1];k++)
		{
			A->JA[count]=A21->JA[k];
			A->val[count]=A21->val[k];
			count++;
		}
		for(k=A22->IA[i];k<A22->IA[i+1];k++)
		{
			A->JA[count]=A22->JA[k]+A11->col;
			A->val[count]=A22->val[k];
			count++;
		}
		for(k=A23->IA[i];k<A23->IA[i+1];k++)
		{
			A->JA[count]=A23->JA[k]+A11->col+A12->col;
			A->val[count]=A23->val[k];
			count++;
		}
	}

	for(i=0;i<A31->row;i++)
	{
		for(k=A31->IA[i];k<A31->IA[i+1];k++)
		{
			A->JA[count]=A31->JA[k];
			A->val[count]=A31->val[k];
			count++;
		}
		for(k=A32->IA[i];k<A32->IA[i+1];k++)
		{
			A->JA[count]=A32->JA[k]+A11->col;
			A->val[count]=A32->val[k];
			count++;
		}
		for(k=A33->IA[i];k<A33->IA[i+1];k++)
		{
			A->JA[count]=A33->JA[k]+A11->col+A12->col;
			A->val[count]=A33->val[k];
			count++;
		}
	}
}


/**
 * \fn void getPenaltyParameters(ddenmat *etas, idenmat *elementEdge, EDGE *edges, double *alpha, double *beta)
 * \brief get penalty parameters for HDG method 
 * \param *etas pointer to penalty parameter
  * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param alpha power of penalty parameter
 * \param beta coefficient of penalty parameter
 * \return void
 */
void getPenaltyParameters(ddenmat *etas, idenmat *elementEdge, EDGE *edges, double *alpha, double *beta)
{
	int i,j,k,l;
	
	for(i=0;i<3;i++)
		create_dden_matrix(elementEdge->row, elementEdge->col, etas+i);

	for(k=0;k<elementEdge->row;k++)
	{
		for(j=0;j<elementEdge->col;j++)
		{
			l=elementEdge->val[k][j];
			for(i=0;i<3;i++)
				etas[i].val[k][j] = pow(edges->length[l], alpha[i])*beta[i];
		}
	}	
}

/**
 * \fn void getFaceInfo(idenmat *elements, iCSRmat *elementsTran, idenmat *faces, iCSRmat *facesTran, idenmat *elementFace)
 * \brief get faces information from elements
 *			 ALgorithm: node-->element-->face
 * \param *elements stores 4 nodes corresponding to the element, the 5th column stores the relation with coarse grid elements( store -1 if itself is the coarset grid),
                    the 6th column stores the type of element
 * \param *elementsTran pointer to the transpose of *elements
 * \param *faces the first three columns store the three vertices corresponding to the face; 
 *				 the 4rd and 5th columns store the elements which the face belongs to;
 *				 if the face is a boundary, the 5th column will stores -1;
 *				 the first column is in ascend order.
 * \param *facesTran the relation between nodes and faces. JA stores face index, A stores another two vertices
 * \param *elementFace the relation between elements and faces
 * \return void
 */
void getFaceInfo(ELEMENT *elements, iCSRmat *elementsTran, FACE *faces, dennode *nodes)
{
	int i,j,k,l;
	int element, point1, point2, point3, col;	
	
	int iSTART=0;
	faces->row=0;
	faces->col=5;
	faces->val=NULL;
	col=elements->col;
	for(i=0;i<elementsTran->row;i++)
	{
		for(j=elementsTran->IA[i];j<elementsTran->IA[i+1];j++)
		{
			element=elementsTran->JA[j]; // achieve the element
			for(k=0;k<col;k++)
			{
				if(elements->val[element][k]==i)
				{
					point1=elements->val[element][(k+1)%col];
					point2=elements->val[element][(k+2)%col];
					point3=elements->val[element][(k+3)%col];
					break;
				}
			} // k
			
			if(i<point1 && i<point2) // the case i>point1 or i>point2 already exists in faces
				addFace(i, point1, point2, element, iSTART, faces);
								
			if(i<point1 && i<point3) // the case i>point1 or i>point3 already exists in faces
				addFace(i, point1, point3, element, iSTART, faces);

			if(i<point2 && i<point3) // the case i>point2 or i>point3 already exists in faces
				addFace(i, point2, point3, element, iSTART, faces);
		} // j
		
		iSTART=faces->row;
	} // i

	faces->nvector = (double**)calloc(faces->row, sizeof(double *));
	faces->t1vector = (double**)calloc(faces->row, sizeof(double *));
	faces->t2vector = (double**)calloc(faces->row, sizeof(double *));
	faces->t01 = (double**)calloc(faces->row, sizeof(double *));
	faces->t02 = (double**)calloc(faces->row, sizeof(double *));
	faces->t12 = (double**)calloc(faces->row, sizeof(double *));
	faces->barycenter = (double**)calloc(faces->row, sizeof(double *));
	for (i = 0; i<faces->row; i++)
	{
		faces->nvector[i] = (double*)calloc(3, sizeof(double));
		faces->t1vector[i] = (double*)calloc(3, sizeof(double));
		faces->t2vector[i] = (double*)calloc(3, sizeof(double));
		faces->t01[i] = (double*)calloc(3, sizeof(double));
		faces->t02[i] = (double*)calloc(3, sizeof(double));
		faces->t12[i] = (double*)calloc(3, sizeof(double));
		faces->barycenter[i] = (double*)calloc(3, sizeof(double));
	}
	faces->area = (double*)calloc(faces->row, sizeof(double));
	faces->bdFlag = (int*)calloc(faces->row, sizeof(int));
	
	getFaceBdflag(faces, nodes);
}

/**
* \fn void getFaceBdflag(FACE *faces, dennode *nodes)
* \brief get faces boundary information from vertices boundary information
* \param *faces the first two columns store the two vertice corresponding to the edge;
*				 the 3rd and 4th columns store the triangles which the edge belongs to;
*				 if the edge is a boundary, the 4th column will stores -1;
*				 the first column is in ascend order.
* \param *nodes pointer to the nodes location of the triangulation
* \return void
************* boundary type of vertex *************
0 : non-boundary, i.e., an interior vertex.
1 : first type, i.e., a Dirichlet boundary vertex.
2 : second type, i.e., a Neumann boundary vertex.
3 : third type, i.e., a Robin boundary vertex.
12: a Dirichlet-Neumann boundary vertex.
22: a Neumann-Neumann boundary vertex.
************* boundary type of face *************
* 0 : non - boundary, i.e., an interior edge or face.
* 1 : first type, i.e., a Dirichlet boundary edge or face.
* 2 : second type, i.e., a Neumann boundary edge or face.
* 3 : third type, i.e., a Robin boundary edge or face.
*/
void getFaceBdflag(FACE *faces, dennode *nodes)
{
	int *vbdflag = nodes->bdFlag;
	int *fbdflag = faces->bdFlag;
	int i, j, k, temp;
	int ev[2];
	int lflag[2], rflag[2];

	for (i = 0; i < faces->row; i++)
	{
		if (faces->val[i][4] == -1)
		{
			fbdflag[i] = 1;
			// will be discussed further for complicated boundary condition
		}
		else
			fbdflag[i] = 0;
	}

}

/**
 * \fn void getEdgeInfo3(ELEMENT *elements, iCSRmat *elementsTran, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes)
 * \brief get edges information from elements
 *			 ALgorithm: node-->element-->edge
 * \param *elements stores 4 nodes corresponding to the element, the 5th column stores the relation with coarse grid elements( store -1 if itself is the coarset grid),
                    the 6th column stores the type of element
 * \param *elementsTran pointer to the transpose of *elements
 * \param *elementEdge the relation between elements and edges
 * \param *edges stores the two vertice corresponding to the edge
 *				 the first column is in ascend order.
 * \return void
 */
void getEdgeInfo3(ELEMENT *elements, iCSRmat *elementsTran, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes)
{
	int i,j,k,l;
	int element, face, point[3], col;	
	
	int iSTART=0;
	edges->row=0;
	edges->col=2;
	edges->val=NULL;
	col = elements->col;
	for(i=0;i<elementsTran->row;i++)
	{
		for(j=elementsTran->IA[i];j<elementsTran->IA[i+1];j++)
		{
			element=elementsTran->JA[j]; // achieve the element
			for(k=0;k<col;k++)
			{
				if(elements->val[element][k]==i)
				{
					point[0]=elements->val[element][(k+1)%col];
					point[1]=elements->val[element][(k+2)%col];
					point[2]=elements->val[element][(k+3)%col];
					break;
				}
			} // k
			
			for(k=0;k<3;k++)
			{
				if(i<point[k]) // the case i>point[k] already exists in edges
				{
					for(l=iSTART;l<edges->row;l++)
					{					
						if((edges->val[l][0]==i) && (edges->val[l][1]==point[k]))///////////
							break;
					}
					if(l==edges->row) // the edge(i, point[k]) is new 
					{
						edges->val=(int**)realloc(edges->val, sizeof(int *)*(edges->row+1));
						edges->val[edges->row]=(int*)calloc(edges->col, sizeof(int));
						edges->row++;
						edges->val[edges->row-1][0]=i;
						edges->val[edges->row-1][1]=point[k];
					}
				}
			} // k
		} // j
		
		iSTART=edges->row;
	} // i
	
	edges->n1vector=(double**)calloc(edges->row, sizeof(double *));
	edges->n2vector=(double**)calloc(edges->row, sizeof(double *));
	edges->tvector=(double**)calloc(edges->row, sizeof(double *));
	for(i=0;i<edges->row;i++)
	{
		edges->n1vector[i]=(double*)calloc(3, sizeof(double));
		edges->n2vector[i]=(double*)calloc(3, sizeof(double));
		edges->tvector[i]=(double*)calloc(3, sizeof(double));
	}
	edges->length=(double*)calloc(edges->row, sizeof(double));
	edges->bdFlag = (int*)calloc(edges->row, sizeof(int));

	// generate elementEdge
	create_iden_matrix(elements->row, 6, elementEdge);

	int edge;
	for(edge=0;edge<edges->row;edge++)
	{
		point[0]=edges->val[edge][0];
		point[1]=edges->val[edge][1];
		for(k=elementsTran->IA[point[0]];k<elementsTran->IA[point[0]+1];k++)
		{
			element=elementsTran->JA[k];
			for(j=0;j<col;j++)
			{
				if(elements->val[element][j]==point[1])
					break;
			} // j

			if(j<col)
			{
				for(i=0;i<col;i++)
				{
					if(elements->val[element][i]==point[0])
						break;
				} // i

				elementEdge->val[element][vv2edge3d(i,j)]=edge;
			}
		} // k
	} // edge

	getEdgeBdflag3(edges, elements, faces, elementEdge, nodes);
}

/**
* \fn void getEdgeBdflag3(EDGE *edges, ELEMENT *elements, FACE *faces, idenmat *elementEdge, dennode *nodes)
* \brief get edges boundary information from vertices boundary information
* \param *faces the first two columns store the two vertice corresponding to the edge;
*				 the 3rd and 4th columns store the triangles which the edge belongs to;
*				 if the edge is a boundary, the 4th column will stores -1;
*				 the first column is in ascend order.
* \param *nodes pointer to the nodes location of the triangulation
* \return void
************* boundary type of vertex *************
0 : non-boundary, i.e., an interior vertex.
1 : first type, i.e., a Dirichlet boundary vertex.
2 : second type, i.e., a Neumann boundary vertex.
3 : third type, i.e., a Robin boundary vertex.
12: a Dirichlet-Neumann boundary vertex.
22: a Neumann-Neumann boundary vertex.
************* boundary type of face *************
* 0 : non - boundary, i.e., an interior edge or face.
* 1 : first type, i.e., a Dirichlet boundary edge or face.
* 2 : second type, i.e., a Neumann boundary edge or face.
* 3 : third type, i.e., a Robin boundary edge or face.
*/
void getEdgeBdflag3(EDGE *edges, ELEMENT *elements, FACE *faces, idenmat *elementEdge, dennode *nodes)
{
//	int *vbdflag = nodes->bdFlag;
	int *ebdflag = edges->bdFlag;
	int i, j, element, face, edge, vertex, temp;
	int vi[3];
	
	for (i = 0; i < edges->row; i++)
		ebdflag[i] = 0;

	for (face = 0; face < faces->row; face++)
	{
		if (faces->val[face][4] > -1)
			continue;

		element = faces->val[face][3];
		for (i = 0; i < 3; i++)
		{
			vertex = faces->val[face][i];
			for (vi[i] = 0; vi[i] < elements->col; vi[i]++)
			{
				if (elements->val[element][vi[i]] == vertex)
					break;
			}
		}

		edge = elementEdge->val[element][vv2edge3d(vi[0], vi[1])];
		ebdflag[edge] = 1;
		edge = elementEdge->val[element][vv2edge3d(vi[1], vi[2])];
		ebdflag[edge] = 1;
		edge = elementEdge->val[element][vv2edge3d(vi[2], vi[0])];
		ebdflag[edge] = 1;
	}
}

/**
 * \fn void getPermutation(int *a, int *b, int *perm, int n)
 * \brief The array b is a permutation of the array a. Get the permutation relation betweeen a and b
 * \param *a pointer to the first array
 * \param *b pointer to the second array, which is the permutation of the first array a
 * \param *perm pointer to the relation betweeen a and b, i.e. the permutation of 0, 1, ..., n-1 
 * \param n the length of array.
 * \return void
 */
void getPermutation(int *a, int *b, int *perm, int n)
{
	int i, j;
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			if(b[j]==a[i])
			{
				perm[i]=j;
				break;
			}
		}
	}
}

/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/

static void addFace(int i, int point1, int point2, int element, int iSTART, FACE *faces)
{
	int pointmin, pointmax;
	int l;

	if(point1<point2)
	{
		pointmin=point1;
		pointmax=point2;
	}
	else
	{
		pointmin=point2;
		pointmax=point1;
	}

	for(l=iSTART;l<faces->row;l++)
	{					
		if((faces->val[l][0]==i) && (faces->val[l][1]==pointmin) && (faces->val[l][2]==pointmax))
			break;
	}
	if(l==faces->row) // the face(i, pointmin, pointmax) is new 
	{
		faces->val=(int**)realloc(faces->val, sizeof(int *)*(faces->row+1));
		faces->val[faces->row]=(int*)calloc(faces->col, sizeof(int));
		faces->row++;
		faces->val[faces->row-1][0]=i;
		faces->val[faces->row-1][1]=pointmin;
		faces->val[faces->row-1][2]=pointmax;
		faces->val[faces->row-1][3]=element;
		faces->val[faces->row-1][4]=-1;
	}
	else // the face(i, pointmin, pointmax) already exits
	{
		faces->val[l][4]=element;
	}
}

/******************************************************
 *********************mapping**************************
	      vertex1    vertex2   ---->   midpoint
		     0          1      ---->       0
			 1          2      ---->       1
			 2          3      ---->       2
			 3          0      ---->       3
			 0          2      ---->       4
			 1          3      ---->       5
 *******************************************************/
int vv2edge3d(int v1, int v2)
{
	if((v1==0&&v2==1)||(v1==1&&v2==0))
		return 0;

	if((v1==1&&v2==2)||(v1==2&&v2==1))
		return 1;

	if((v1==2&&v2==3)||(v1==3&&v2==2))
		return 2;

	if((v1==3&&v2==0)||(v1==0&&v2==3))
		return 3;

	if((v1==0&&v2==2)||(v1==2&&v2==0))
		return 4;

	if((v1==1&&v2==3)||(v1==3&&v2==1))
		return 5;

	return -1;
}

/******************************************************
 *********************mapping**************************
  	   edge   ---->    vertex1    vertex2
		0     ---->       0          1      
		1     ---->       1          2      
		2     ---->       2          3      
		3     ---->   	  3          0     
		4     ---->   	  0          2     
		5     ---->   	  1          3     
 *******************************************************/
void edge2vv3d(int edge, int *v)
{
	v[0] = edge%4;
	v[1] = (v[0]+edge/4 +1)%4;
}

/******************************************************
 *********************mapping**************************
  	   face   ---->    vertex1    vertex2    vertex3
		0     ---->       1          2          3      
		1     ---->       0          3          2      
		2     ---->       0          1          3      
		3     ---->   	  0          2          1
 *******************************************************/
void face2vertices3d(int face, int *v)
{
	switch(face)
	{
	case 0: v[0]=1;v[1]=2;v[2]=3;break;
	case 1: v[0]=0;v[1]=3;v[2]=2;break;
	case 2: v[0]=0;v[1]=1;v[2]=3;break;
	case 3: v[0]=0;v[1]=2;v[2]=1;break;
	default: v[0]=-1;v[1]=-1;v[2]=-1;
	}
}


/**
 * \fn static void getEdgeElement3(ELEMENT *elements, iCSRmat *elementsTran, EDGE *edges, iCSRmat *edgeElement)
 * \brief get edges information from elements
 *			 ALgorithm: node-->element-->edge
 * \param *elements stores 4 nodes corresponding to the element, the 5th column stores the relation with coarse grid elements( store -1 if itself is the coarset grid),
                    the 6th column stores the type of element
 * \param *elementsTran pointer to the transpose of *elements
 * \param *faces the first three columns store the three vertices corresponding to the face; 
 *				 the 4rd and 5th columns store the elements which the face belongs to;
 *				 if the face is a boundary, the 5th column will stores -1;
 *				 the first column is in ascend order.
 * \param *facesTran the relation between nodes and faces. JA stores face index, A stores another two vertices
 * \param *edges stores the two vertice corresponding to the edge
 *				 the first column is in ascend order.
 * \param *edgeElement the relation between edges and elements
 * \param *edgeFace the relation between edges and faces
 * \return void
 */
static void getEdgeElement3(ELEMENT *elements, iCSRmat *elementsTran, EDGE *edges, iCSRmat *edgeElement)
{
	int i,j,k,l;
	int element, point[3], col;	
	
	col = elements->col;
	
	// generate edgeElement
	int edge;
	int index=0;
	edgeElement->row=edges->row;
	edgeElement->col=elements->row;
	edgeElement->IA=(int*)calloc(edgeElement->row+1, sizeof(int));
	edgeElement->nnz = elements->row * 6;
	edgeElement->JA=(int*)calloc(edgeElement->nnz, sizeof(int));
	edgeElement->val=(int*)calloc(edgeElement->nnz, sizeof(int));
	for(edge=0;edge<edges->row;edge++)
	{
		point[0]=edges->val[edge][0];
		point[1]=edges->val[edge][1];
		for(k=elementsTran->IA[point[0]];k<elementsTran->IA[point[0]+1];k++)
		{
			element=elementsTran->JA[k];
			for(j=0;j<col;j++)
			{
				if(elements->val[element][j]==point[1])
					break;
			} // j

			if(j<col)
			{
				for(i=0;i<col;i++)
				{
					if(elements->val[element][i]==point[0])
						break;
				} // i
				edgeElement->JA[index]=element;
				edgeElement->val[index]=vv2edge3d(i,j);
				index++;
			}
		} // k
		edgeElement->IA[edge+1]=index;
	} // edge
}
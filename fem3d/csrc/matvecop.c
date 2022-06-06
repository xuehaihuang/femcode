/*
 *  matvecop.c
 *
 *------------------------------------------------------
 *
 *		Created by Chensong Zhang on 03/31/2009.
 *		Copyright 2009 PSU. All rights reserved. 
 *
 *------------------------------------------------------
 */

/*! \file matvecop.c
 *  \brief Matrix-vector operations
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "header.h"
#include "matvec.h"

// Some function for matrix/vector memory management ...

/**
 * \fn int create_csr_matrix(int m, int n, int nnz, dCSRmat *A)
 * \brief Create CSR sparse matrix data memory space
 * \param m integer, number of rows
 * \param n integer, number of columns
 * \param nnz integer, number of nonzeros
 * \param *A pointer to the dCSRmat matrix
 * \return 1 if succeed
 */
int create_csr_matrix(int m, int n, int nnz, dCSRmat *A)
{		
	A->row=m;
	A->col=n;
	A->nnz=nnz;
	A->IA=(int*)calloc(m+1, sizeof(int));
	A->JA=(int*)calloc(nnz, sizeof(int));
	A->val=(double*)calloc(nnz, sizeof(double));
	return 1;
}

/**
 * \fn int free_csr_matrix(dCSRmat *A)
 * \brief Free CSR sparse matrix data memeory space
 * \param *A pointer to the dCSRmat matrix
 * \return 1 if succeed
 */
int free_csr_matrix(dCSRmat *A)
// Free CSR sparse matrix data
{		
	A->row=0;
	A->col=0;
	A->nnz=0;
	free(A->IA);
	free(A->JA);
	free(A->val);
	return 1;
}

/**
 * \fn int free_icsr_matrix(iCSRmat *A)
 * \brief Free int CSR sparse matrix data memeory space
 * \param *A pointer to the iCSRmat matrix
 * \return 1 if succeed
 */
int free_icsr_matrix(iCSRmat *A)
// Free CSR sparse matrix data
{		
	A->row=0;
	A->col=0;
	A->nnz=0;
	free(A->IA);
	free(A->JA);
	free(A->val);
	return 1;
}

/**
 * \fn int create_iden_matrix(int m, int n, idenmat *A)
 * \brief Create int dense matrix data memory space
 * \param m integer, number of rows
 * \param n integer, number of columns
 * \param *A pointer to the idenmat matrix
 * \return 1 if succeed
 */
int create_iden_matrix(int m, int n, idenmat *A)
{		
	A->row=m;
	A->col=n;
	A->val=(int**)calloc(A->row, sizeof(int *));
	int i;
	for(i=0;i<A->row;i++)
		A->val[i]=(int*)calloc(A->col, sizeof(int));
	return 1;
}

/**
 * \fn int free_iden_matrix(idenmat *A)
 * \brief Free int dense matrix data memeory space
 * \param *A pointer to the idenmat matrix
 * \return 1 if succeed
 */
int free_iden_matrix(idenmat *A)
{	
	int i;
	for(i=0;i<A->row;i++)
		free(A->val[i]);
	free(A->val);
	A->row=0;
	A->col=0;
	return 1;
}


/**
 * \fn int create_ELEMENT(int m, int n, ELEMENT *A)
 * \brief Create ELEMENT
 * \param m integer, number of rows
 * \param n integer, number of columns
 * \param *A pointer to ELEMENT
 * \return 1 if succeed
 */
int create_ELEMENT(int m, int n, ELEMENT *A)
{	
	A->row=m;
	A->col=n;
	A->val=(int**)calloc(A->row, sizeof(int *));
	A->parent = (int*)calloc(A->row, sizeof(int));
	A->type = (int*)calloc(A->row, sizeof(int));
	A->vertices=(double***)calloc(A->row, sizeof(double **));
	A->barycenter=(double**)calloc(A->row, sizeof(double *));
	A->bcFace=(double***)calloc(A->row, sizeof(double **));
	A->lambdaConst=(double**)calloc(A->row, sizeof(double *));
	A->gradLambda=(double***)calloc(A->row, sizeof(double **));
	A->nvector=(double***)calloc(A->row, sizeof(double **));
	A->fperm = (int***)malloc(A->row * sizeof(int **));
	A->fpermi = (int***)malloc(A->row * sizeof(int **));
	A->forien = (short**)malloc(A->row * sizeof(short *));
	A->eperm = (int***)malloc(A->row * sizeof(int **));
	A->eorien = (short**)malloc(A->row * sizeof(short *));
	int i,j;
	for(i=0;i<A->row;i++)
	{
		A->val[i]=(int*)calloc(A->col, sizeof(int));
		A->vertices[i]=(double**)calloc(A->col, sizeof(double *));
		A->barycenter[i]=(double*)calloc(A->col-1, sizeof(double));
		A->bcFace[i]=(double**)calloc(A->col, sizeof(double *));
		A->lambdaConst[i]=(double*)calloc(A->col, sizeof(double));
		A->gradLambda[i]=(double**)calloc(A->col, sizeof(double *));
		A->nvector[i]=(double**)calloc(A->col, sizeof(double *));
		A->fperm[i] = (int**)malloc(A->col * sizeof(int *));
		A->fpermi[i] = (int**)malloc(A->col * sizeof(int *));
		A->forien[i] = (short*)malloc(A->col * sizeof(short));
		A->eperm[i] = (int**)malloc(A->col*(A->col - 1) / 2 * sizeof(int *));
		A->eorien[i] = (short*)malloc(A->col*(A->col - 1) / 2 * sizeof(short));
		for(j=0;j<A->col;j++)
		{
			A->vertices[i][j]=(double*)calloc(A->col-1, sizeof(double));
			A->bcFace[i][j]=(double*)calloc(A->col-1, sizeof(double));
			A->gradLambda[i][j]=(double*)calloc(A->col-1, sizeof(double));
			A->nvector[i][j]=(double*)calloc(A->col-1, sizeof(double));
			A->fperm[i][j] = (int*)malloc((A->col - 1) * sizeof(int));
			A->fpermi[i][j] = (int*)malloc((A->col - 1) * sizeof(int));
		}
		for (j = 0; j < A->col*(A->col - 1) / 2; j++)
			A->eperm[i][j] = (int*)malloc(2 * sizeof(int));
	}
	A->vol=(double*)calloc(A->row, sizeof(double));
	return 1;
}

/**
 * \fn int free_ELEMENT(ELEMENT *A)
 * \brief Free ELEMENT
 * \param *A pointer to ELEMENT
 * \return 1 if succeed
 */
int free_ELEMENT(ELEMENT *A)
{	
	int i,j;
	for(i=0;i<A->row;i++)
	{
		free(A->val[i]);
		free(A->barycenter[i]);
		free(A->lambdaConst[i]);
		for(j=0;j<A->col;j++)
		{
			free(A->vertices[i][j]);
			free(A->bcFace[i][j]);
			free(A->gradLambda[i][j]);
			free(A->nvector[i][j]);
			free(A->fperm[i][j]);
			free(A->fpermi[i][j]);
		}
		free(A->vertices[i]);
		free(A->bcFace[i]);
		free(A->gradLambda[i]);
		free(A->nvector[i]);
		free(A->fperm[i]);
		free(A->fpermi[i]);
		free(A->forien[i]);
		for (j = 0; j < A->col*(A->col - 1) / 2; j++)
			free(A->eperm[i][j]);
		free(A->eperm[i]);
		free(A->eorien[i]);
	}
	free(A->parent);
	free(A->type);
	free(A->val);
	free(A->vol);
	free(A->vertices);
	free(A->barycenter);
	free(A->bcFace);
	free(A->lambdaConst);
	free(A->gradLambda);
	free(A->nvector);
	free(A->fperm);
	free(A->fpermi);
	free(A->forien);
	free(A->eperm);
	free(A->eorien);
	A->row=0;
	A->col=0;
	return 1;
}

/**
 * \fn int create_FACE(int m, int n, FACGE *A)
 * \brief Create FACE
 * \param m integer, number of rows
 * \param n integer, number of columns
 * \param *A pointer to FACE
 * \return 1 if succeed
 */
int create_FACE(int m, int n, FACE *A)
{	
	A->row=m;
	A->col=n;
	A->val=(int**)calloc(A->row, sizeof(int *));
	A->nvector=(double**)calloc(A->row, sizeof(double *));
	A->t1vector=(double**)calloc(A->row, sizeof(double *));
	A->t2vector=(double**)calloc(A->row, sizeof(double *));
	A->t01=(double**)calloc(A->row, sizeof(double *));
	A->t02=(double**)calloc(A->row, sizeof(double *));
	A->t12=(double**)calloc(A->row, sizeof(double *));
	A->barycenter=(double**)calloc(A->row, sizeof(double *));
	int i;
	for(i=0;i<A->row;i++)
	{
		A->val[i]=(int*)calloc(A->col, sizeof(int));
		A->nvector[i]=(double*)calloc(3, sizeof(double));
		A->t1vector[i]=(double*)calloc(3, sizeof(double));
		A->t2vector[i]=(double*)calloc(3, sizeof(double));
		A->t01[i]=(double*)calloc(3, sizeof(double));
		A->t02[i]=(double*)calloc(3, sizeof(double));
		A->t12[i]=(double*)calloc(3, sizeof(double));
		A->barycenter[i]=(double*)calloc(3, sizeof(double));
	}
	A->area=(double*)calloc(A->row, sizeof(double));
	A->bdFlag = (int*)calloc(A->row, sizeof(int));
	return 1;
}

/**
 * \fn int free_FACE(FACE *A)
 * \brief Free FACE
 * \param *A pointer to FACE
 * \return 1 if succeed
 */
int free_FACE(FACE *A)
{	
	int i;
	for(i=0;i<A->row;i++)
	{
		free(A->val[i]);
		free(A->nvector[i]);
		free(A->t1vector[i]);
		free(A->t2vector[i]);
		free(A->t01[i]);
		free(A->t02[i]);
		free(A->t12[i]);
		free(A->barycenter[i]);
	}
	free(A->val);
	free(A->nvector);
	free(A->t1vector);
	free(A->t2vector);
	free(A->t01);
	free(A->t02);
	free(A->t12);
	free(A->area);
	free(A->barycenter);
	free(A->bdFlag);
	A->row=0;
	A->col=0;
	return 1;
}

/**
 * \fn int create_EDGE(int m, int n, EDGE *A)
 * \brief Create EDGE
 * \param m integer, number of rows
 * \param n integer, number of columns
 * \param *A pointer to EDGE
 * \return 1 if succeed
 */
int create_EDGE(int m, int n, EDGE *A)
{		
	A->row=m;
	A->col=n;
	A->val=(int**)calloc(A->row, sizeof(int *));
	A->n1vector=(double**)calloc(A->row, sizeof(double *));
	A->n2vector=(double**)calloc(A->row, sizeof(double *));
	A->tvector=(double**)calloc(A->row, sizeof(double *));
	int i;
	for(i=0;i<A->row;i++)
	{
		A->val[i]=(int*)calloc(A->col, sizeof(int));
		A->n1vector[i]=(double*)calloc(3, sizeof(double));
		A->n2vector[i]=(double*)calloc(3, sizeof(double));
		A->tvector[i]=(double*)calloc(3, sizeof(double));
	}
	A->length=(double*)calloc(A->row, sizeof(double));
	A->bdFlag = (int*)calloc(A->row, sizeof(int));
	return 1;
}

/**
 * \fn int free_EDGE(EDGE *A)
 * \brief Free EDGE
 * \param *A pointer to EDGE
 * \return 1 if succeed
 */
int free_EDGE(EDGE *A)
{	
	int i;
	for(i=0;i<A->row;i++)
	{
		free(A->val[i]);
		free(A->n1vector[i]);
		free(A->n2vector[i]);
		free(A->tvector[i]);
	}
	free(A->val);
	free(A->n1vector);
	free(A->n2vector);
	free(A->tvector);
	free(A->length);
	free(A->bdFlag);
	A->row=0;
	A->col=0;
	return 1;
}

/**
 * \fn int create_dden_matrix(int m, int n, ddenmat *A)
 * \brief Create double dense matrix data memory space
 * \param m integer, number of rows
 * \param n integer, number of columns
 * \param *A pointer to the ddenmat matrix
 * \return 1 if succeed
 */
int create_dden_matrix(int m, int n, ddenmat *A)
{		
	A->row=m;
	A->col=n;
	A->val=(double**)calloc(A->row, sizeof(double *));
	int i;
	for(i=0;i<A->row;i++)
		A->val[i]=(double*)calloc(A->col, sizeof(double));
	return 1;
}

/**
 * \fn int free_dden_matrix(ddenmat *A)
 * \brief Free double dense matrix data memeory space
 * \param *A pointer to the ddenmat matrix
 * \return 1 if succeed
 */
int free_dden_matrix(ddenmat *A)
{	
	int i;
	for(i=0;i<A->row;i++)
		free(A->val[i]);
	free(A->val);
	A->row=0;
	A->col=0;
	return 1;
}

/**
 * \fn void init_dden_matrix(ddenmat *A, double val)
 * \brief Initial double dense matrix
 * \param *A pointer to the ddenmat matrix
 * \param val initial value
 */
void init_dden_matrix(ddenmat *A, double val)
{		
	int i, j;
	for(i=0;i<A->row;i++)
	{
		for(j=0;j<A->col;j++)
			A->val[i][j] = val;
	}
}

/**
 * \fn int create_dden_matrix3(int l, int m, int n, ddenmat3 *A)
 * \brief Create three-dimensinal double dense matrix data memory space
 * \param l integer, number of pages
 * \param m integer, number of rows
 * \param n integer, number of columns
 * \param *A pointer to the ddenmat3 matrix
 * \return 1 if succeed
 */
int create_dden_matrix3(int l, int m, int n, ddenmat3 *A)
{		
	A->pag=l;
	A->row=m;
	A->col=n;
	A->val=(double***)calloc(A->pag, sizeof(double **));
	int i,j;
	for(i=0;i<A->pag;i++)
	{
		A->val[i]=(double**)calloc(A->row, sizeof(double*));
		for(j=0;j<A->row;j++)
			A->val[i][j]=(double*)calloc(A->col, sizeof(double));
	}
	return 1;
}

/**
 * \fn int free_dden_matrix3(ddenmat3 *A)
 * \brief Free three-dimensinal double dense matrix data memeory space
 * \param *A pointer to the three-dimensinal ddenmat matrix
 * \return 1 if succeed
 */
int free_dden_matrix3(ddenmat3 *A)
{	
	int i,j;
	for(i=0;i<A->pag;i++)
	{
		for(j=0;j<A->row;j++)
			free(A->val[i][j]);
		free(A->val[i]);
	}
    free(A->val);
	A->pag=0;
	A->row=0;
	A->col=0;
	return 1;
}

/**
 * \fn int create_dbd_matrix(int m, int n, int nb, dBDmat *A)
 * \brief Create block diagonal dense matrix data memory space
 * \param m integer, number of rows
 * \param n integer, number of columns
 * \param nb integer, number of blocks
 * \param *A pointer to the dBDmat matrix
 * \return 1 if succeed
 */
int create_dbd_matrix(int m, int n, int nb, dBDmat *A)
{		
	A->row=m;
	A->col=n;
	A->nb=nb;
	A->blk=(ddenmat*)calloc(nb, sizeof(ddenmat));
	A->blk->row=0;
	A->blk->col=0;
	A->blk->val=NULL;
	return 1;
}

/**
 * \fn int free_dbd_matrix(dBDmat *A)
 * \brief Free block diagonal dense matrix data memeory space
 * \param *A pointer to the dBDmat matrix
 * \return 1 if succeed
 */
int free_dbd_matrix(dBDmat *A)
{		
	int i;
	for(i=0;i<A->nb;i++)
		free_dden_matrix(A->blk+i);
	A->blk=NULL;
	A->row=0;
	A->col=0;
	A->nb=0;
	return 1;
}

/**
* \fn int create_dobd_matrix(int m, int n, int nb, dOBDmat *A)
* \brief Create overlapped block diagonal dense matrix data memory space
* \param m integer, number of rows
* \param n integer, number of columns
* \param nb integer, number of blocks
* \param *A pointer to the dBDmat matrix
* \return 1 if succeed
*/
int create_dobd_matrix(int m, int n, int nb, dOBDmat *A)
{
	A->row = m;
	A->col = n;
	A->nb = nb;
	A->blk = (ddenmat*)malloc(nb * sizeof(ddenmat));
	A->blk->row = 0;
	A->blk->col = 0;
	A->blk->val = NULL;
	A->rindices = (ivector*)malloc(nb * sizeof(ivector));
	A->rindices->row = 0;
	A->rindices->val = NULL;
	A->cindices = (ivector*)malloc(nb * sizeof(ivector));
	A->cindices->row = 0;
	A->cindices->val = NULL;
	return 1;
}

/**
* \fn int free_dobd_matrix(dOBDmat *A)
* \brief Free overlapped block diagonal dense matrix data memeory space
* \param *A pointer to the dBDmat matrix
* \return 1 if succeed
*/
int free_dobd_matrix(dOBDmat *A)
{
	int i;
	for (i = 0; i < A->nb; i++)
	{
		free_dden_matrix(A->blk + i);
		free_ivector(A->rindices + i);
		free_ivector(A->cindices + i);
	}
	free(A->blk);
	A->blk = NULL;
	free(A->rindices);
	A->rindices = NULL;
	free(A->cindices);
	A->cindices = NULL;
	A->row = 0;
	A->col = 0;
	A->nb = 0;
	return 1;
}

/**
 * \fn int create_dennode(int m, int n, dennode *A)
 * \brief Create dennode data memory space
 * \param m integer, number of rows
 * \param n integer, number of columns
 * \param *A pointer to the dennode data
 * \return 1 if succeed
 */
int create_dennode(int m, int n, dennode *A)
{		
	A->row=m;
	A->col=n;
	A->bdFlag = (int*)calloc(A->row, sizeof(int));
	A->val=(double**)calloc(A->row, sizeof(double *));
	int i;
	for (i = 0; i < A->row; i++)
		A->val[i] = (double*)calloc(A->col, sizeof(double));
	return 1;
}

/**
 * \fn int free_dennode(dennode *A)
 * \brief Free dennode data memeory space
 * \param *A pointer to the dennode data
 * \return 1 if succeed
 */
int free_dennode(dennode *A)
{	
	int i;
	for(i=0;i<A->row;i++)
		free(A->val[i]);
	free(A->val);
	free(A->bdFlag);
	A->row=0;
	A->col=0;
	return 1;
}

/**
 * \fn int create_elementDOF(int dop, int dof, int row, int col, ELEMENT_DOF *A)
 * \brief Create ELEMENT_DOF data memory space
 * \param dop integer, number of degree of polynomial
 * \param dof integer, number of degree of freedom
 * \param row integer, number of elements
 * \param col integer, number of columns
 * \param *A pointer to the DOF data
 * \return 1 if succeed
 */
int create_elementDOF(int dop, int dof, int row, int col, ELEMENT_DOF *A)
{		
	A->dop=dop;
	A->dof=dof;
	A->row=row;
	A->col=col;
	A->val=(int**)calloc(A->row, sizeof(int *));
	int i;
	for(i=0;i<A->row;i++)
		A->val[i]=(int*)calloc(A->col, sizeof(int));
	A->nfFlag.row = 0; A->nfFlag.val = NULL;
	A->freenodes.row = 0; A->freenodes.val = NULL;
	A->nfreenodes.row = 0; A->nfreenodes.val = NULL;
	A->index.row = 0; A->index.val = NULL;
	return 1;
}

/**
 * \fn int free_elementDOF(ELEMENT_DOF *A);
 * \brief Free ELEMENT_DOF data memeory space
 * \param *A pointer to the ELEMENT_DOF data
 * \return 1 if succeed
 */
int free_elementDOF(ELEMENT_DOF *A)
{	
	int i;
	for(i=0;i<A->row;i++)
		free(A->val[i]);
	free(A->val);
	A->row=0;
	A->col=0;
	A->dop=0;
	A->dof=0;
	if (A->nfFlag.row>0)
		free_ivector(&A->nfFlag);
	if (A->freenodes.row>0)
		free_ivector(&A->freenodes);
	if (A->nfreenodes.row>0)
		free_ivector(&A->nfreenodes);
	if (A->index.row>0)
		free_ivector(&A->index);
	return 1;
}


/**
 * \fn int create_dvector(int m, dvector *u)
 * \brief Create dvector data space of double type
 * \param m integer, number of rows
 * \param *u pointer to the dvector
 * \return 1 if succeed
 */
int create_dvector(int m, dvector *u)
{		
	u->row=m;
	u->val=(double*)calloc(m, sizeof(double)); 	
	return 1;
}


/**
 * \fn int create_ivector(int m, ivector *u)
 * \brief Create vector data space of int type
 * \param m integer, number of rows
 * \param *u pointer to the vector
 * \return 1 if succeed
 */
int create_ivector(int m, ivector *u)
{		
	u->row=m;
	u->val=(int*)calloc(m, sizeof(int)); 	
	return 1;
}

/**
 * \fn int free_vector(vector *u)
 * \brief Free vector data space of double type
 * \param *u pointer to the vector
 * \return 1 if succeed
 */
int free_dvector(dvector *u)
{		
	u->row=0;
	free(u->val);
	return 1;
}

/**
 * \fn int free_ivector(ivector *u)
 * \brief Free vector data space of int type
 * \param *u pointer to the vector
 * \return 1 if succeed
 */
int free_ivector(ivector *u)
{		
	// this function is same as free_dvector except input pointer type
	u->row=0;
	free(u->val);
	return 1;
}



// Some simple array Blas operations follow ...

/**
 * \fn int copy_arrayint(int n, int *x, int *y) 
 * \brief Copy an int array to the other y=x
 * \param n number of variables
 * \param *x pointer to the original vector
 * \param *y pointer to the destination vector
 * \return 1 if succeed, 0 if fail
 */
int copy_arrayint(int n, int *x, int *y) 
{
	int i;
	for (i=0; i<n; i++) y[i]=x[i];
	return 1;
}

/**
 * \fn int init_array(int n, double *x, double val)
 * \brief Initialize an array to be zero x=0
 * \param n number of variables
 * \param *x pointer to the vector
 * \param val initial value for the double array
 * \return 1 if succeed, 0 if fail
 */
int init_array(int n, double *x, double val)
// x = 0
{
	int i;
	for (i=0; i<n; i++) x[i]=val;
	return 1;
}

/**
* \fn int copy_iarray(int n, int *x, int *y)
* \brief Copy an array to the other y=x
* \param n number of variables
* \param *x pointer to the original vector
* \param *y pointer to the destination vector
* \return 1 if succeed, 0 if fail
*/
int copy_iarray(int n, int *x, int *y)
{
	int i;
	for (i = 0; i<n; i++) y[i] = x[i];
	return 1;
}

/**
 * \fn int copy_array(int n, double *x, double *y) 
 * \brief Copy an array to the other y=x
 * \param n number of variables
 * \param *x pointer to the original vector
 * \param *y pointer to the destination vector
 * \return 1 if succeed, 0 if fail
 */
int copy_array(int n, double *x, double *y) 
{
	int i;
	for (i=0; i<n; i++) y[i]=x[i];
	return 1;
}

/**
 * \fn int ax_array(int n, double a, double *x)
 * \brief x = a*x
 * \param n number of variables
 * \param a a real number
 * \param *x pointer to the original vector
 * \return 1 if succeed, 0 if fail
 */
int ax_array(int n, double a, double *x) 
{
	int i;
	for (i=0; i<n; i++) x[i] *= a;
	return 1;
}

/**
 * \fn int ax_array(int n, double a, double *x, double *y)
 * \brief y = a*x
 * \param n number of variables
 * \param a a real number
 * \param *x pointer to the original vector
 * \return 1 if succeed, 0 if fail
 */
int axy_array(int n, double a, double *x, double *y) 
{
	int i;
	for (i=0; i<n; i++) y[i] = a * x[i];
	return 1;
}

/**
 * \fn int axpy_array(int n, double a, double *x, double *y)
 * \brief y = a*x + y
 * \param n number of variables
 * \param a a real number
 * \param *x pointer to the original vector
 * \param *y pointer to the destination vector
 * \return 1 if succeed, 0 if fail
 */
int axpy_array(int n, double a, double *x, double *y) 
{
	int i;
	for (i=0; i<n; i++) y[i] += a*x[i];
	return 1;
}

/**
 * \fn int axpyz_array(int n, double a, double *x, double *y, double *z)
 * \brief z = a*x + y, z is the third vector
 * \param n number of variables
 * \param a a real number
 * \param *x pointer to the original vector 1
 * \param *y pointer to the original vector 2
 * \param *z pointer to the destination vector
 * \return 1 if succeed, 0 if fail
 */
int axpyz_array(int n, double a, double *x, double *y, double *z) 
// z = a*x+y
{
	int i;
	for (i=0; i<n; i++) z[i] = a*x[i]+y[i];
	return 1;
}

/**
 * \fn int axpby_array(int n, double a, double *x, double b, double *y)
 * \brief y = a*x + b*y
 * \param n number of variables
 * \param a real number
 * \param b real number
 * \param *x pointer to the origianl vector
 * \param *y pointer to the destination vector
 * \return 1 if succeed, 0 if fail
 */
int axpby_array(int n, double a, double *x, double b, double *y) 
{
	int i;
	for(i=0; i<n; i++) y[i] = a*x[i]+b*y[i];
	return 1;
}

/**
 * \fn int axpbyz_array(int n, double a, double *x, double b, double *y, double *z)
 * \brief z = a*x + b*y
 * \param n number of variables
 * \param a real number
 * \param b real number
 * \param *x pointer to the original vector 1
 * \param *y pointer to the original vector 2
 * \param *z pointer to the destination vector
 * \return 1 if succeed, 0 if fail
 */
int axpbyz_array(int n, double a, double *x, double b, double *y, double *z) 
{
	int i;
	for(i=0; i<n; i++) z[i] = a*x[i]+b*y[i];
	return 1;
}

/**
 * \fn double dot_array(int n, double *x, double *y) 
 * \brief Inner product of two arraies (x,y)
 * \param n number of variables
 * \param *x pointer to vector 1
 * \param *y pointer to vector 2
 * \return inner product
 */
double dot_array(int n, double *x, double *y) 
{
	int i;
	double value;
	for (value=i=0;i<n;i++) value+=x[i]*y[i];
	return value;
}

/**
 * \fn void cross_array(double *x, double *y, double *z) 
 * \brief Cross product of two arraies x and y
 * \param *x pointer to vector 1
 * \param *y pointer to vector 2
 * \return cross product
 */
void cross_array(double *x, double *y, double *z) 
{
	z[0]=x[1]*y[2]-x[2]*y[1];
	z[1]=x[2]*y[0]-x[0]*y[2];
	z[2]=x[0]*y[1]-x[1]*y[0];
}

/**
 * \fn void orthocomplement_array(double *x, double *y, double *z) 
 * \brief find a set of two normal orthogonal bases of the orthogonal complement space of vector x in three dimensions
 * \param *x pointer to vector 1
 * \param *y pointer to orthogonal basis 1
 * \param *z pointer to orthogonal basis 2
 * \return a set of two normal orthogonal bases of the orthogonal complement space of vector x
 */
void orthocomplement_array(double *x, double *y, double *z) 
{
	int i, mini;
	double temp;
	for(i=0;i<3;i++)
		z[i]=0;

	mini=0;
	for(i=1;i<3;i++)
	{
		if(fabs(x[i])<fabs(x[mini]))
			mini=i;
	}
	z[mini]=1;

	cross_array(x, z, y);
	temp=sqrt(dot_array(3, y, y));
	for(i=0;i<3;i++)
		y[i]/=temp;

	cross_array(x, y, z);
}

/**
 * \fn void tensorproduct3d_array(double *a, double *b, double *tau) 
 * \brief tensor product of two vector in three dimensions
 * \param *a pointer to vector 1
 * \param *b pointer to vector 2
 * \param *tau pointer to tensor
 * \return tensor product of two vector a and b
 */
void tensorproduct3d_array(double *a, double *b, double *tau) 
{
	// double subscript of tensor: 00 11 22 12 20 01
	tau[0]=a[0]*b[0];
	tau[1]=a[1]*b[1];
	tau[2]=a[2]*b[2];
	tau[3]=(a[1]*b[2]+a[2]*b[1])/2;
	tau[4]=(a[2]*b[0]+a[0]*b[2])/2;
	tau[5]=(a[0]*b[1]+a[1]*b[0])/2;
}


// some simple vector Blas operations follows, same operations as above, but for vector

/**
 * \fn int init_dvector(int n, dvector *x, double val) 
 * \brief Initialize dvector x=0
 * \param n number of variables
 * \param *x pointer to dvector
 * \param val initial value for the dvector
 * \return 1 if succeed
 */
int init_dvector(dvector *x, double val)
{
	int i;
	for (i = 0; i<x->row; i++) x->val[i] = val;
	return 1;
}

/**
* \fn int init_dvector2b(dvector *x, double val)
* \brief Initialize dvector x=0
* \param n number of variables
* \param *x pointer to dvector
* \param val initial value for the dvector
* \return 1 if succeed
*/
int init_dvector2b(dvector *x, double val)
{
	int i;
	for (i = 0; i<x[0].row; i++) x[0].val[i] = val;
	for (i = 0; i<x[1].row; i++) x[1].val[i] = val;
	return 1;
}

/**
 * \fn int copy_dvector(dvector *x, dvector *y) 
 * \brief Copy dvector x to dvector y
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 * \return 1 if succeed
 */
int copy_dvector(dvector *x, dvector *y) 
{
	int i;
	y->row=x->row;
	for (i=0; i<y->row; i++) y->val[i]=x->val[i];
	return 1;
}

/**
* \fn int copy_dvector2b(dvector *x, dvector *y)
* \brief Copy dvector x to dvector y in 2 blocks
* \param *x pointer to dvector
* \param *y pointer to dvector
* \return 1 if succeed
*/
int copy_dvector2b(dvector *x, dvector *y)
{
	int i;
	y[0].row = x[0].row;
	for (i = 0; i<y[0].row; i++)
		y[0].val[i] = x[0].val[i];

	y[1].row = x[1].row;
	for (i = 0; i<y[1].row; i++)
		y[1].val[i] = x[1].val[i];

	return 1;
}

/**
 * \fn int axpy_dvector(double a, dvector *x, dvector *y) 
 * \brief y = a*x + y
 * \param a real number
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 * \return 1 if succeed, 0 if fail
 */
int axpy_dvector(double a, dvector *x, dvector *y) 
{
	int i;
	if (y->row == x->row) {
		for (i=0; i<y->row; i++) y->val[i] += a*x->val[i];
		return 1;
	}
	else {
		return 0;
	}
}

/**
* \fn int axpy_dvector(double a, dvector *x, dvector *y)
* \brief y = a*x + y in 2 blocks
* \param a real number
* \param *x pointer to dvector
* \param *y pointer to dvector
* \return 1 if succeed, 0 if fail
*/
int axpy_dvector2b(double a, dvector *x, dvector *y)
{
	int i;
	for (i = 0; i<y[0].row; i++) 
		y[0].val[i] += a*x[0].val[i];

	for (i = 0; i<y[1].row; i++)
		y[1].val[i] += a*x[1].val[i];

	return 1;
}

/**
 * \fn int axpyz_dvector(double a, dvector *x, dvector *y, dvector *z) 
 * \brief z = a*x + y, z is a third vector (z is cleared)
 * \param a real number
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 * \param *z pointer to dvector
 * \return 1 if succeed, 0 if fail
 */
int axpyz_dvector(double a, dvector *x, dvector *y, dvector *z) 
{
	int i;
	if (x->row == y->row) {
		z->row = x->row;
		for (i=0; i<z->row; i++) z->val[i] = a*x->val[i]+y->val[i];
		return 1;
	}
	else {
		return 0;
	}
}

/**
* \fn int axy_dvector(double a, dvector *x, dvector *y)
* \brief y = a*x
* \param a real number
* \param *x pointer to dvector
* \param *y pointer to dvector
* \return 1 if succeed
*/
int axy_dvector(double a, dvector *x, dvector *y)
{
	int i;
	for (i = 0; i < y->row; i++)
		y->val[i] = a*x->val[i];

	return 1;
}

/**
* \fn int axy_dvector2b(double a, dvector *x, dvector *y)
* \brief y = a*x in 2 blocks
* \param a real number
* \param *x pointer to dvector
* \param *y pointer to dvector
* \return 1 if succeed
*/
int axy_dvector2b(double a, dvector *x, dvector *y)
{
	int i;
	for (i = 0; i < y[0].row; i++)
		y[0].val[i] = a*x[0].val[i];

	for (i = 0; i < y[1].row; i++)
		y[1].val[i] = a*x[1].val[i];

	return 1;
}

/**
 * \fn double dot_dvector(dvector *x, dvector *y) 
 * \brief Inner product of two vectors (x,y)
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 * \return Inner product
 */
double dot_dvector(dvector *x, dvector *y) 
{
	int i;
	double value;
	for (value=i=0;i<x->row;i++) value+=x->val[i]*y->val[i];
	return value;
}

/**
* \fn double dot_dvector2b(dvector *x, dvector *y)
* \brief Inner product of two vectors (x,y) in 2 block
* \param *x pointer to dvector
* \param *y pointer to dvector
* \return Inner product
*/
double dot_dvector2b(dvector *x, dvector *y)
{
	int i;
	double value=0;
	for (i = 0; i<x[0].row; i++) 
		value += x[0].val[i] * y[0].val[i];
	
	for (i = 0; i<x[1].row; i++)
		value += x[1].val[i] * y[1].val[i];
	return value;
}

/**
 * \fn double maxnorm_dvector(dvector *x) 
 * \brief Maximal norm of dvector x
 * \param *x pointer to dvector
 * \return maximal norm of x
 */
double maxnorm_dvector(dvector *x) 
{
	int i;
	double maxnorm, value;

	for (maxnorm=i=0;i<x->row;i++) {
		value=fabs(x->val[i]);
		if (maxnorm<value) maxnorm=value;
	}		
	return maxnorm;
}

/**
 * \fn double maxdiff_dvector(dvector *a, dvector *b) 
 * \brief Maximal difference of two dvector a and b
 * \param *a pointer to dvector
 * \param *b pointer to dvector
 * \return maximal norm of a-b
 */
double maxdiff_dvector(dvector *a, dvector *b)
{
    int i;
    double Linf = 0.0, diffi;

    for(i=0;i<a->row;i++)
       if ((diffi = fabs(a->val[i] - b->val[i])) > Linf)
          Linf = diffi;

    return Linf;
}

/**
 * \fn double onenorm_dvector(dvector *x) 
 * \brief L1 norm of dvector x
 * \param *x pointer to dvector
 * \return L1 norm of x
 */
double onenorm_dvector(dvector *x) 
{
	int i;
	double onenorm;
	for (onenorm=i=0;i<x->row;i++) onenorm+=fabs(x->val[i]);
	return onenorm;
}

/**
 * \fn double twonorm_dvector(dvector *x) 
 * \brief L2 norm of dvector x
 * \param *x pointer to dvector
 * \return L2 norm of x
 */
double twonorm_dvector(dvector *x) 
{
	int i;
	double twonorm;
	for (twonorm=i=0;i<x->row;i++) twonorm+=pow(x->val[i],2.0);
	return sqrt(twonorm);
}

/**
 * \fn int sparse_mv0(double alpha, dCSRmat *A, double *x, double *y) 
 * \brief Matrix-vector multiplication y = alpha*A*x
 * \param alpha real number
 * \param *A pointer to dCSRmat CSR matrix
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 * \return 1 if succeed
 */
int sparse_mv0(double alpha, dCSRmat *A, double *x, double *y) 
{
	int i,k,j;	
	double temp;
	
	for(i=0;i<A->row;i++) {
		temp=0.0;
		for(k=A->IA[i];k<A->IA[i+1];k++) {
			j=A->JA[k];
			temp+=A->val[k]*x[j];
		}
		y[i]=temp*alpha;
	}
	return 1;
}


/**
 * \fn int sparse_mv(double alpha, dCSRmat *A, double *x, double *y) 
 * \brief Matrix-vector multiplication y = alpha*A*x + y
 * \param alpha real number
 * \param *A pointer to dCSRmat CSR matrix
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 * \return 1 if succeed
 */
int sparse_mv(double alpha, dCSRmat *A, double *x, double *y) 
{
	int i,k,j;	
	double temp;
	
	for(i=0;i<A->row;i++) {
		temp=0.0;
		for(k=A->IA[i];k<A->IA[i+1];k++) {
			j=A->JA[k];
			temp+=A->val[k]*x[j];
		}
		y[i]+=temp*alpha;
	}
	return 1;
}

/**
* \fn int sparse_mv2b(double alpha, dCSRmat *A, dvector *x, dvector *y)
* \brief Matrix-vector multiplication y = alpha*A*x + y in 2 blocks
* \param alpha real number
* \param *A pointer to dCSRmat CSR matrix
* \param *x pointer to dvector
* \param *y pointer to dvector
* \return 1 if succeed
*/
int sparse_mv2b(double alpha, dCSRmat *A, dvector *x, dvector *y)
{
	int i, k, j;
	double temp;

	for (i = 0; i<y[0].row; i++) {
		temp = 0.0;
		for (k = A[0].IA[i]; k<A[0].IA[i + 1]; k++) {
			j = A[0].JA[k];
			temp += A[0].val[k] * x[0].val[j];
		}
		for (k = A[1].IA[i]; k<A[1].IA[i + 1]; k++) {
			j = A[1].JA[k];
			temp += A[1].val[k] * x[1].val[j];
		}
		y[0].val[i] += temp*alpha;
	}
	
	for (i = 0; i<y[1].row; i++) {
		temp = 0.0;
		for (k = A[2].IA[i]; k<A[2].IA[i + 1]; k++) {
			j = A[2].JA[k];
			temp += A[2].val[k] * x[0].val[j];
		}
		if (A[3].row > 0)
		{
			for (k = A[3].IA[i]; k < A[3].IA[i + 1]; k++) {
				j = A[3].JA[k];
				temp += A[3].val[k] * x[1].val[j];
			}
		}
		y[1].val[i] += temp*alpha;
	}
	
	return 1;
}

/**
 * \fn int denmat_mv(double alpha, dCSRmat *A, double *x, double *y) 
 * \brief Matrix-vector multiplication y = alpha*A*x + y
 * \param alpha real number
 * \param *A pointer to ddenmat matrix
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 * \return 1 if succeed
 */
int denmat_mv(double alpha, ddenmat *A, double *x, double *y) 
{
	int i,j;	
	double temp;
	
	for(i=0;i<A->row;i++) {
		temp=0.0;
		for(j=0;j<A->col;j++)
			temp+=A->val[i][j]*x[j];
		
		y[i]+=temp*alpha;
	}
	return 1;
}

/**
* \fn int dBDmat_mv(double alpha, dBDmat *A, dvector *x, dvector *y)
* \brief Matrix-vector multiplication y = alpha*A*x + y
* \param alpha real number
* \param *A pointer to dBDmat matrix
* \param *x pointer to dvector
* \param *y pointer to dvector
* \return 1 if succeed
*/
int dBDmat_mv(double alpha, dBDmat *A, dvector *x, dvector *y)
{
	if (A->row != y->row || A->col != x->row)
		return 0;

	int i, j, i1, j1, nblk;
	double temp;
	ddenmat *blk;

	int rcount = 0;
	int lcount = 0;
	for (nblk = 0; nblk<A->nb; nblk++)
	{
		blk = A->blk + nblk;
		for (i1 = 0; i1<blk->row; i1++)
		{
			i = rcount + i1;
			temp = 0.0;
			for (j1 = 0; j1<blk->col; j1++)
			{
				j = lcount + j1;
				temp += blk->val[i1][j1] * x->val[j];
			}
			y->val[i] += temp*alpha;
		}

		rcount += blk->row;
		lcount += blk->col;
	}

	return 1;
}

/**
* \fn int dBDmat_mv0(double alpha, dBDmat *A, dvector *x, dvector *y)
* \brief Matrix-vector multiplication y = alpha*A*x
* \param alpha real number
* \param *A pointer to dBDmat matrix
* \param *x pointer to dvector
* \param *y pointer to dvector
* \return 1 if succeed
*/
int dBDmat_mv0(double alpha, dBDmat *A, dvector *x, dvector *y)
{
	if (A->row != y->row || A->col != x->row)
		return 0;

	int i, j, i1, j1, nblk;
	double temp;
	ddenmat *blk;

	int rcount = 0;
	int lcount = 0;

	for (nblk = 0; nblk<A->nb; nblk++)
	{
		blk = A->blk + nblk;
		for (i1 = 0; i1<blk->row; i1++)
		{
			i = rcount + i1;
			temp = 0.0;
			for (j1 = 0; j1<blk->col; j1++)
			{
				j = lcount + j1;
				temp += blk->val[i1][j1] * x->val[j];
			}
			y->val[i] = temp*alpha;
		}

		rcount += blk->row;
		lcount += blk->col;
	}

	return 1;
}

/**
 * \fn int dBDMultiplydvector(double alpha, dBDmat *A, dvector *x, dvector *y) 
 * \brief Matrix-vector multiplication y = alpha*A*x + y
 * \param alpha real number
 * \param *A pointer to dBDmat matrix
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 * \return 1 if succeed
 */
int dBDMultiplydvector(double alpha, dBDmat *A, dvector *x, dvector *y) 
{
	if(A->row!=y->row || A->col!=x->row)
		return 0;

	int i,j,i1,j1, nblk;	
	double temp;
	ddenmat *blk;

	int rcount=0;
	int lcount=0;
	for(nblk=0;nblk<A->nb;nblk++)
	{
		blk=A->blk+nblk;
		for(i1=0;i1<blk->row;i1++)
		{
			i = rcount+i1;
			temp=0.0;
			for(j1=0;j1<blk->col;j1++)
			{
				j = lcount+j1;
				temp+=blk->val[i1][j1]*x->val[j];
			}
			y->val[i]+=temp*alpha;
		}

		rcount += blk->row;
		lcount += blk->col;
	}

	return 1;
}

/**
 * \fn int getdiag(int n, dCSRmat *A, dvector *diag)
 * \brief Get first n diagonal entries of a CSR matrix A
 * \param n an interger (if n=0, then get all diagonal entries)
 * \param *A pointer to dCSRmat CSR matrix
 * \param *diag pointer to the diagonal as a dvector
 * \return 1 if succeed, 0 if fail
 */
int getdiag(int n, dCSRmat *A, dvector *diag) 
{
	int i,k,j,ret;	
	
	if (n==0) n=min(A->row,A->col);
	ret=create_dvector(n,diag);
	
	for(i=0;i<n;i++) {
		for(k=A->IA[i];k<A->IA[i+1];k++) {
			j=A->JA[k];
			if (j==i) diag->val[i] = A->val[k];
		}
	}

	return ret;
}

/**
 * \fn void Axy_ddenmat(double alpha, ddenmat *A, double *x, double *y)
 * \brief Matrix-vector multiplication y = alpha*A*x
 * \param alpha real number
 * \param *A pointer to ddenmat matrix
 * \param *B pointer to ddenmat matrix
 * \param *C pointer to ddenmat matrix
 */
void Axy_ddenmat(double alpha, ddenmat *A, double *x, double *y)
{
	int i,j;	
	double temp;
	
	for(i=0;i<A->row;i++) {
		temp=0.0;
		for(j=0;j<A->col;j++)
			temp += A->val[i][j]*x[j];
		y[i] = temp*alpha;
	}
}

/**
 * \fn void Atxy_ddenmat(double alpha, ddenmat *A, double *x, double *y)
 * \brief Matrix-vector multiplication y = alpha*A^T*x
 * \param alpha real number
 * \param *A pointer to ddenmat matrix
 * \param *B pointer to ddenmat matrix
 * \param *C pointer to ddenmat matrix
 */
void Atxy_ddenmat(double alpha, ddenmat *A, double *x, double *y)
{
	int i,j;	
	double temp;
	
	for(i=0;i<A->col;i++) {
		temp=0.0;
		for(j=0;j<A->row;j++)
			temp += A->val[j][i]*x[j];
		y[i] = temp*alpha;
	}
}

/**
 * \fn void AB_ddenmat(double alpha, ddenmat *A, ddenmat *B, ddenmat *C)
 * \brief Matrix-vector multiplication C = alpha*A*B
 * \param alpha real number
 * \param *A pointer to ddenmat matrix
 * \param *B pointer to ddenmat matrix
 * \param *C pointer to ddenmat matrix
 */
void AB_ddenmat(double alpha, ddenmat *A, ddenmat *B, ddenmat *C)
{
	int i,j,k;	
	double temp;
	if(A->col!=B->row) {
		printf("AB_ddenmat: A->col!=B->row\n");
		exit(1);
	}
	
	for(i=0;i<C->row;i++) {
		for(j=0;j<C->col;j++) {
			temp=0.0;
			for(k=0;k<A->col;k++)
				temp += A->val[i][k]*B->val[k][j];
			C->val[i][j] = temp*alpha;
		}
	}
}

/**
 * \fn void ABt_ddenmat(double alpha, ddenmat *A, ddenmat *B, ddenmat *C)
 * \brief Matrix-vector multiplication C = alpha*A*B^T
 * \param alpha real number
 * \param *A pointer to ddenmat matrix
 * \param *B pointer to ddenmat matrix
 * \param *C pointer to ddenmat matrix
 */
void ABt_ddenmat(double alpha, ddenmat *A, ddenmat *B, ddenmat *C)
{
	int i,j,k;	
	double temp;
	if(A->col!=B->col) {
		printf("ABt_ddenmat: A->col!=B->col\n");
		exit(1);
	}
	
	for(i=0;i<C->row;i++) {
		for(j=0;j<C->col;j++) {
			temp=0.0;
			for(k=0;k<A->col;k++)
				temp += A->val[i][k]*B->val[j][k];
			C->val[i][j] = temp*alpha;
		}
	}
}

/**
 * \fn void AtB_ddenmat(double alpha, ddenmat *A, ddenmat *B, ddenmat *C)
 * \brief Matrix-vector multiplication C = alpha*A^T*B
 * \param alpha real number
 * \param *A pointer to ddenmat matrix
 * \param *B pointer to ddenmat matrix
 * \param *C pointer to ddenmat matrix
 */
void AtB_ddenmat(double alpha, ddenmat *A, ddenmat *B, ddenmat *C)
{
	int i,j,k;	
	double temp;
	if(A->row!=B->row) {
		printf("AtB_ddenmat: A->row!=B->row\n");
		exit(1);
	}
	
	for(i=0;i<C->row;i++) {
		for(j=0;j<C->col;j++) {
			temp=0.0;
			for(k=0;k<A->row;k++)
				temp += A->val[k][i]*B->val[k][j];
			C->val[i][j] = temp*alpha;
		}
	}
}

/**
 * \fn void ABAt_ddenmat(double alpha, ddenmat *A, ddenmat *B, ddenmat *C) 
 * \brief Matrix-vector multiplication C = alpha*A*B*A^T
 * \param alpha real number
 * \param *A pointer to ddenmat matrix
 * \param *B pointer to ddenmat matrix
 * \param *C pointer to ddenmat matrix
 * \return 1 if succeed
 */
void ABAt_ddenmat(double alpha, ddenmat *A, ddenmat *B, ddenmat *C)
{
	int i,j,l,k;	
	double temp;
	ddenmat AB;
	create_dden_matrix(A->row, B->col, &AB);
	AB_ddenmat(alpha, A, B, &AB);
	ABt_ddenmat(1.0, &AB, A, C);
	free_dden_matrix(&AB);
}

/**
* \fn void compress_dcsr(dCSRmat *A, dCSRmat *B, double EPS)
* \brief Remove the zero terms of a CSR matrix A
* \param *A pointer to the original dCSRmat CSR matrix
* \param *B pointer to the compressed dCSRmat CSR matrix
* \param EPS used to characterize zero
* \return void
*/
void compress_dcsr(dCSRmat *A, dCSRmat *B, double EPS)
{
	int i, j, k;
	int count;
	B->row = A->row;
	B->col = A->col;
	B->IA = (int*)calloc(B->row + 1, sizeof(int));
	B->JA = NULL;
	B->val = NULL;

	for (i = 0; i < A->row; i++)
	{
		for (j = A->IA[i]; j < A->IA[i + 1]; j++)
		{
			if (abs(A->val[j]) > EPS)
				B->IA[i + 1]++;
		}
	}

	for (i = 0; i<B->row; i++)
		B->IA[i + 1] += B->IA[i];
	B->nnz = B->IA[B->row];

	B->JA = (int*)calloc(B->nnz, sizeof(int));
	B->val = (double*)calloc(B->nnz, sizeof(double));

	for (i = 0; i < A->row; i++)
	{
		count = 0;
		for (j = A->IA[i]; j < A->IA[i + 1]; j++)
		{
			if (abs(A->val[j]) > EPS)
			{
				B->JA[B->IA[i] + count] = A->JA[j];
				B->val[B->IA[i] + count] = A->val[j];
				count++;
			}
		}
	}
}

/**
 * \fn int inverse_dBDmat(dBDmat *A, dBDmat *Ainv)
 * \brief Get inverse matrix of dBDmat A
 * \param *A pointer to dBDmat matrix
 * \param *Ainv pointer to inverse matrix of A
 * \return 1 if succeed, 0 if fail
 */
int inverse_dBDmat(dBDmat *A, dBDmat *Ainv)
{
	int i,j,k;
	create_dbd_matrix(A->row, A->col, A->nb, Ainv);
	for(i=0;i<Ainv->nb;i++)
		create_dden_matrix(A->blk[i].row, A->blk[i].col, Ainv->blk+i);

	ddenmat block;
	
	for(k=0;k<A->nb;k++)
	{
		create_dden_matrix(A->blk[k].row, A->blk[k].col, &block);
		for(i=0;i<A->blk[k].row;i++)
		{
			for(j=0;j<A->blk[k].col;j++)
				block.val[i][j]=A->blk[k].val[i][j];
			Ainv->blk[k].val[i][i]=1;
		}
		
		AxBrref(&block, Ainv->blk+k);

		free_dden_matrix(&block);
	}
	
	return 1;
}


/**
 * \fn int inverse_dBDmat2(dBDmat *A, dBDmat *Ainv, int col1, int col2)
 * \brief Get inverse matrix of dBDmat A where each block in A also include two subblocks
 * \param *A pointer to dBDmat matrix
 * \param *Ainv pointer to inverse matrix of A
 * \param col1 the column of the first subblock in each block of A
 * \param col2 the column of the second subblock in each block of A
 * \return 1 if succeed, 0 if fail
 */
int inverse_dBDmat2(dBDmat *A, dBDmat *Ainv, int col1, int col2)
{
	if(col1<1 || col2<1 || !((col1+col2)==A->blk[0].col))
		return 0;

	int i,j,k;
	create_dbd_matrix(A->row, A->col, A->nb, Ainv);
	for(i=0;i<Ainv->nb;i++)
		create_dden_matrix(A->blk[i].row, A->blk[i].col, Ainv->blk+i);

	ddenmat A1, A2, A1inv, A2inv;
	create_dden_matrix(col1, col1, &A1);
	create_dden_matrix(col2, col2, &A2);
	create_dden_matrix(col1, col1, &A1inv);
	create_dden_matrix(col2, col2, &A2inv);
	for(k=0;k<A->nb;k++)
	{
		for(i=0;i<col1;i++)
		{
			for(j=0;j<col1;j++)
			{
				A1.val[i][j]=A->blk[k].val[i][j];
				A1inv.val[i][j]=0;
			}
			A1inv.val[i][i]=1;
		}
		
		for(i=0;i<col2;i++)
		{
			for(j=0;j<col2;j++)
			{
				A2.val[i][j]=A->blk[k].val[i+col1][j+col1];
				A2inv.val[i][j]=0;
			}
			A2inv.val[i][i]=1;
		}

		AxBrref(&A1, &A1inv);
		AxBrref(&A2, &A2inv);
		
		for(i=0;i<col1;i++)
		{
			for(j=0;j<col1;j++)
				Ainv->blk[k].val[i][j]=A1inv.val[i][j];
		}
		for(i=0;i<col2;i++)
		{
			for(j=0;j<col2;j++)
				Ainv->blk[k].val[i+col1][j+col1]=A2inv.val[i][j];
		}
	}
	free_dden_matrix(&A1);
	free_dden_matrix(&A2);
	free_dden_matrix(&A1inv);
	free_dden_matrix(&A2inv);

	return 1;
}

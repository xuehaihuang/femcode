/*
 *  output.c
 *
 *  Created by Xuehai Huang on 08/31/2010.
 *  Copyright 2010 WZU. All rights reserved.
 *
 */

/*! \file output.c
 *  \brief Output to files
 */

#include <stdio.h>
#include <stdlib.h>
#include "header.h"


/**
 * \fn int getElementsNodes4Matlab(ELEMENT *elements, dennode *nodes, dvector *uh)
 * \brief output elements and nodes data used for drawing mesh in Matlab
 * \param *elements pointer to the structure of the triangulation
 * \param *nodes pointer to the nodes location of the triangulation
 * \param *uh pointer to the value of solution
 * \return 1: success 0: false
 */
int getElementsNodes4Matlab(ELEMENT *elements, dennode *nodes, dvector *uh)
{
	int i;
	int Nt=elements->row;
	char *fnodes="output/nodes.dat";
	char *felements="output/elements.dat";
	char *fuh="output/uh.dat";

	FILE *outputFile;
	outputFile=fopen(fnodes, "w");
	if(outputFile==NULL)
	{
		printf("Opening file %s fails!\n", fnodes);
		return 0;
	}
	for(i=0;i<nodes->row;i++)
		fprintf(outputFile, "%e %e\n", nodes->val[i][0], nodes->val[i][1]);
	fclose(outputFile);
	
	outputFile=fopen(felements, "w");
	if(outputFile==NULL)
	{
		printf("Opening file %s fails!\n", felements);
		return 0;
	}
	for(i=0;i<elements->row;i++)
		fprintf(outputFile, "%d %d %d\n", elements->val[i][0]+1, elements->val[i][1]+1, elements->val[i][2]+1);
	fclose(outputFile);

	if(uh==NULL)
		return 1;

	outputFile=fopen(fuh, "w");
	if(outputFile==NULL)
	{
		printf("Opening file %s fails!\n", fuh);
		return 0;
	}
//	for(i=0;i<nodes->row;i++)
	for (i = 0; i<uh->row; i++)
		fprintf(outputFile, "%e\n", uh->val[i]);
	fclose(outputFile);

	return 1;
}

/**
 * \fn int writeElementsNodes(ELEMENT *elements, dennode *nodes, char *fname1, char *fname2)
 * \brief Write elements and nodes to disk files
 */
int writeElementsNodes(ELEMENT *elements, dennode *nodes, char *fname1, char *fname2)
{     
	int i, j;
	FILE *outputFile1, *outputFile2;
	
	outputFile1=fopen(fname1, "w");
	if(outputFile1==NULL)
	{
		printf("writeElementsNodes: opening file %s fails!\n", fname1);
		return 0;
	}
	int col=elements->col;
	for(i=0;i<elements->row;i++)
	{
		for(j=0;j<col;j++)
			fprintf(outputFile1,"%d ",elements->val[i][j]+1);
		// fprintf(outputFile1,"%d\n",elements->type[i]);
	}
	fclose(outputFile1);

	outputFile2=fopen(fname2, "w");
	if(outputFile2==NULL)
	{
		printf("writeElementsNodes: opening file %s fails!\n", fname2);
		return 0;
	}
	for(i=0;i<nodes->row;i++)
	{
		for(j=0;j<nodes->col;j++)
			fprintf(outputFile2,"%lf ",nodes->val[i][j]);
		fprintf(outputFile2,"%d\n",nodes->bdFlag[i]);
	}
	fclose(outputFile2);
	
  return 1;
}

/**
* \fn int write_IJ_dCSRmat(dCSRmat *A, char *fname)
* \brief Write dCSRmat *A to disk files
*/
int write_IJ_dCSRmat(dCSRmat *A, char *fname)
{
	int i, j;
	FILE *outputFile;

	outputFile = fopen(fname, "w");
	if (outputFile == NULL)
	{
		printf("write_IJ_dCSRmat: opening file %s fails!\n", fname);
		return 0;
	}

	fprintf(outputFile, "%d\t%d\t%d\n", A[0].row, A[0].col, A[0].nnz);
	for (i = 0; i<A[0].row; i++)
	{
		for (j = A[0].IA[i]; j<A[0].IA[i + 1]; j++)
			fprintf(outputFile, "%d\t%d\t%le\n", i + 1, A[0].JA[j] + 1, A[0].val[j]);
	}
	fclose(outputFile);

	return 1;
}

/**
* \fn void write_IJ_dCSRmat4Matlab(dCSRmat *A, char *fname)
* \brief Write dCSRmat *A to disk files
*/
void write_IJ_dCSRmat4Matlab(dCSRmat *A, char *fname)
{
	int i, j;
	FILE *outputFile;

	outputFile = fopen(fname, "w");
	if (outputFile == NULL)
	{
		printf("write_IJ_dCSRmat4Matlab: opening file %s fails!\n", fname);
		return;
	}

	// fprintf(outputFile, "%d\t%d\t%d\n", A[0].row, A[0].col, A[0].nnz);
	for (i = 0; i<A[0].row; i++)
	{
		for (j = A[0].IA[i]; j<A[0].IA[i + 1]; j++)
			fprintf(outputFile, "%d\t%d\t%le\n", i + 1, A[0].JA[j] + 1, A[0].val[j]);
	}
	fclose(outputFile);
}

/**
* \fn void write_dvector4Matlab(dvector *vec, char *fname)
* \brief Write dvector *vec to disk files
*/
void write_dvector4Matlab(dvector *vec, char *fname)
{
	int i;
	FILE *outputFile;

	outputFile = fopen(fname, "w");
	if (outputFile == NULL)
	{
		printf("write_dvector4Matlab: opening file %s fails!\n", fname);
		return;
	}

	for (i = 0; i<vec[0].row; i++){
		fprintf(outputFile, "%le\n", vec[0].val[i]);
	}
	fclose(outputFile);
}

/**
* \fn void read_dvector4Matlab(dvector *vec, char *fname)
* \brief Read dvector *vec from disk files
*/
void read_dvector4Matlab(dvector *vec, char *fname)
{
	int i, j;
	FILE *inputFile;

	inputFile = fopen(fname, "r");
	if (inputFile == NULL)
	{
		printf("read_dvector4Matlab: opening file %s fails!\n", fname);
		return;
	}

	for (i = 0; i<vec[0].row; i++){
		fscanf(inputFile, "%lf", &vec[0].val[i]);
	}
	fclose(inputFile);
}
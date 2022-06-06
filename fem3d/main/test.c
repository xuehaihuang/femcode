/**
 *		Test function for solving a sparse SPD linear system using PCG and AMG. 
 *
 *		Matrix-Solver Community Project
 *
 *------------------------------------------------------
 *
 *		Created by Xuehai Huang on 07/11/2013.
 *		Revised by Xuehai Huang on 13/11/2015.
 *		Copyright 2015 WZU. All rights reserved. 
 *
 *------------------------------------------------------
 *
 */

/*! \file test.c
 *  \brief Test Function for Solvers
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "header.h"
#include "matvec.h"
#include "checkmat.h"

/**
 * \fn int main (int argc, const char * argv[])
 *
 * This is the main function for test purpose. It contains five steps:
 */
int main (int argc, const char * argv[]) 
{	
	ELEMENT *elements;
	FACE *faces;
	EDGE *edges;
	idenmat *elementEdge, *elementFace;
	dennode *nodes;
	int i;
		
	/** Step 0. Read input parameters */
	char *inputfile="ini/input.dat";
	Input_data Input;
	read_Input_data(inputfile,&Input);
	
	int glevelNum = Input.glevelNum;
	int domain_num = Input.domain_num;
	int problem_num = Input.problem_num;
		
	elements = (ELEMENT*)calloc(glevelNum, sizeof(ELEMENT));
	faces = (FACE*)calloc(glevelNum, sizeof(FACE));
	edges = (EDGE*)calloc(glevelNum, sizeof(EDGE));
	nodes = (dennode*)calloc(glevelNum, sizeof(dennode));
	elementFace = (idenmat*)calloc(glevelNum, sizeof(idenmat));
	elementEdge = (idenmat*)calloc(glevelNum, sizeof(idenmat));
	
	/** Step 1. generate mesh */
	if(getmesh(domain_num, elements, elementFace, faces, elementEdge, edges, nodes, glevelNum)==0)
	{
		printf("It's fail to generate mesh.\n");
		return 0;
	}
	
	getElementFaceEdgeGeoInfo(&elements[glevelNum - 1], &elementFace[glevelNum - 1], &faces[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1]);

	
	// int j,k;
	// double *val1, *val2; 
	// for(i=0;i<elements[glevelNum - 1].row;i++)
	// {
	// 	for(j=0;j<4;j++)
	// 	{
	// 		val1=elements[glevelNum - 1].bcFace[i][j];
	// 		k=elementFace[glevelNum - 1].val[i][j];
	// 		val2=faces[glevelNum - 1].barycenter[k];
	// 		printf("(%lf, %lf, %lf), (%lf, %lf, %lf), (%lf, %lf, %lf)\n", val1[0],val1[1],val1[2],val2[0],val2[1],val2[2],val1[0]-val2[0],val1[1]-val2[1],val1[2]-val2[2]);
	// 	}
	// }

	printf("h=%f\n", edges[glevelNum - 1].length[1]);

	/** Step 2. discete method */
	if(problem_num == 1) // Poisson equation
	{
		poissonfem3d(&elements[glevelNum - 1], &elementFace[glevelNum - 1], &faces[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], &Input);
	}
	else if(problem_num == 2) // Maxwell equation
	{
		maxwellfem(&elements[glevelNum - 1], &elementFace[glevelNum - 1], &faces[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], &Input);
	}
	else if(problem_num == 3) // Quad-curl equation
	{
		quadcurlfem(&elements[glevelNum - 1], &elementFace[glevelNum - 1], &faces[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], &Input);
	}
	else if(problem_num == 4) // Quad-curl perturbation equation
	{
		quadcurlperturbfem(&elements[glevelNum - 1], &elementFace[glevelNum - 1], &faces[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], &Input);
	}
	else
	{
		printf("Please set problem_num = 1 or 2!\n");
		exit(0);
	}


	/***************************************output mesh******************************************************/
	FILE *meshFile;
	meshFile=fopen("output/mesh.dat", "w");
	fprintf(meshFile,"%d %d\n",nodes[glevelNum - 1].row,3);
	for (i = 0; i < nodes[glevelNum - 1].row; i++)
		fprintf(meshFile, "%f %f %f %d\n", nodes[glevelNum - 1].val[i][0], nodes[glevelNum - 1].val[i][1], nodes[glevelNum - 1].val[i][2], nodes[glevelNum - 1].bdFlag[i]);
	fprintf(meshFile,"%d %d\n",elements[glevelNum - 1].row,4);
	for (i = 0; i < elements[glevelNum - 1].row; i++)
		fprintf(meshFile, "%d %d %d %d %d\n", elements[glevelNum - 1].val[i][0] + 1, elements[glevelNum - 1].val[i][1] + 1, elements[glevelNum - 1].val[i][2] + 1, elements[glevelNum - 1].val[i][3] + 1, elements[glevelNum - 1].type[i]);
	fclose(meshFile);
	/********************************************************************************************/

	/** For drawing mesh in Matlab */
	writeElementsNodes(&elements[glevelNum - 1], &nodes[glevelNum - 1], "output/elements.dat", "output/nodes.dat");
	
	for (i = 0; i < glevelNum; i++)
	{
		free_ELEMENT(&elements[i]);
		free_FACE(&faces[i]);
		free_EDGE(&edges[i]);
		free_dennode(&nodes[i]);
		free_iden_matrix(&elementFace[i]);
		free_iden_matrix(&elementEdge[i]);
	}

	return 1;
}

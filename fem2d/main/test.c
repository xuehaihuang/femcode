/**
 *		Test function for solving a sparse SPD linear system using PCG and AMG. 
 *
 *		Matrix-Solver Community Project
 *
 *------------------------------------------------------
 *
 *		Created by Xuehai Huang on 07/11/2013.
 *		Copyright 2013 WZU. All rights reserved. 
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
#include <time.h>


#include "header.h"
#include "matvec.h"
#include "checkmat.h"

/**
 * \fn int main (int argc, const char * argv[])
 *
 * This is the main function for test purpose. It contains five steps:
 */
int main(int argc, const char * argv[])
{
	ELEMENT *elements;
	idenmat *elementEdge;
	EDGE *edges;
	dennode *nodes;
	iCSRmat *edgesTran;
	ivector nodeCEdge;


	/** Step 0. Read input parameters */
	char *inputfile = "ini/input.dat";
	Input_data Input;
	read_Input_data(inputfile, &Input);

	int glevelNum = Input.glevelNum;
	int domain_num = Input.domain_num;
	int problem_num = Input.problem_num;

	elements = (ELEMENT*)calloc(glevelNum, sizeof(ELEMENT));
	elementEdge = (idenmat*)calloc(glevelNum, sizeof(idenmat));
	edges = (EDGE*)calloc(glevelNum, sizeof(EDGE));
	nodes = (dennode*)calloc(glevelNum, sizeof(dennode));
	edgesTran = (iCSRmat*)calloc(glevelNum, sizeof(iCSRmat));


	/** Step 1. generate mesh */
	if (getmesh(domain_num, elements, elementEdge, edges, nodes, edgesTran, &nodeCEdge, glevelNum) == 0)
	{
		printf("It's fail to generate mesh.\n");
		return 0;
	}

	int i;

	getElementEdgeGeoInfo(&elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1]);
//	for (i = 0; i < glevelNum; i++)
//		getElementEdgeGeoInfo(&elements[i], &edges[i], &nodes[i]);


	printf("h = %f\n", edges[glevelNum - 1].length[1]);

	/** Step 2. discete method */
	if(problem_num == 1) // Poisson equation
	{
		poissonfem2d(&elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], &Input);
	}
	else if(problem_num == 2) // Biharmonic equation
	{
		biharmonicfem2d(&elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], &Input);
	}
	else if(problem_num == 3) // Linear Elasticity
	{
		linearElas2dfem(&elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], &edgesTran[glevelNum - 1], &nodeCEdge, &Input);
	}
	else if(problem_num == 4) // Stokes equation
	{
		// quadcurlperturbfem(&elements[glevelNum - 1], &elementFace[glevelNum - 1], &faces[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], &Input);
		exit(0);
	}
	else
	{
		printf("Please set problem_num = 1 or 2!\n");
		exit(0);
	}
	

	getElementsNodes4Matlab(&elements[glevelNum - 1], &nodes[glevelNum - 1], NULL);

	for (i = 0; i < glevelNum; i++)
	{
		free_ELEMENT(&elements[i]);
		free_iden_matrix(&elementEdge[i]);
		free_EDGE(&edges[i]);
		free_dennode(&nodes[i]);
		free_icsr_matrix(&edgesTran[i]);
	}
	free_ivector(&nodeCEdge);
	
	free(elements);
	free(elementEdge);
	free(edges);
	free(nodes);
	free(edgesTran);
	
	return 1;
}

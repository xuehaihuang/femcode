/*
 *  quadcurlfem.c
 *
 *  Created by Xuehai Huang on Aug 8, 2020.
 *  Copyright 2020 SUFE. All rights reserved.
 *
 */

/*! \file quadcurlfem.c
 *  \brief Assembling for stiffness matrix and solve it
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "header.h"
#include "matvec.h"
#include "checkmat.h"

 /**
 * \fn void quadcurlperturbfem(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
 * \brief finite element methods for Qual-curl perturbation equation in three dimensions
 * \param *uh pointer to solution
 * \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices, the last 3 columns store the indexes of edge's midpoints
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
 the fourth column stores -1 if the edge is on boundary
 * \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
 the third column stores the Dirichlet status of the points(0: nondirichlet, -1: dirichlet)
 * \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
 * \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
 * \param Input input data structure
 * \return void
 */
void quadcurlperturbfem(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i;
	int variationalform_type = Input->variationalform_type;
	int fem_num = Input->fem_num;
	short nitsche = Input->nitsche;
	if (variationalform_type == 1) // primal formulation
	{
		if(nitsche == 0)
			printf("Discretize quad-curl perturbation equation in primal formulation.\n");
		else
			printf("Discretize quad-curl perturbation equation in primal formulation with the Nitsche’s technique.\n");

		if(fem_num == 1)
		{
			printf("Discretized by the Huang element to be developed.\n");
			exit(0);
			// quadcurlHuang3d(elements, elementFace, faces, elementEdge, edges, nodes, Input);
		}
		
		if(fem_num == 2)
		{
			printf("Discretized by the Huang-Zhang element.\n");
			quadcurlperturbHuangZhang3d(elements, elementFace, faces, elementEdge, edges, nodes, Input);
		}
	}
	else if (variationalform_type == 2) // decoupled mixed formulation
	{
		printf("Discretize quad-curl equation in Maxwell-Stokes-Maxwell form.\n");
		if(fem_num == 1)
		{
			printf("Discretized by the Huang element to be developed.\n");
			exit(0);
			// quadcurlHuang3d_msm(elements, elementFace, faces, elementEdge, edges, nodes, Input);
		}

		if(fem_num == 2)
		{
			printf("Discretized by the Huang-Zhang element to be developed.\n");
			exit(0);
		}
	}
	else
	{
		printf("Please set variationalform_type = 1 or 2!\n");
		exit(0);
	}
}

/**
* \fn void quadcurlperturbHuangZhang3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief Huang-Zhang element method for quad-curl equation in three dimensions
* \param *uh pointer to solution
* \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices, the last 3 columns store the indexes of edge's midpoints
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								  the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
the third column stores the Dirichlet status of the points(0: nondirichlet, -1: dirichlet)
* \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
* \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
* \param Input input data structure
* \return void
*/
void quadcurlperturbHuangZhang3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i,j;
	dCSRmat A;
	dvector b, uh, _uh;
	ELEMENT_DOF elementDOF[2];
	iCSRmat elementdofTran[2];
		
	/** Step 0. Read input parameters */
	int print_level = Input->print_level;
	short nitsche = Input->nitsche;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;
//	int dop1 = Input->dop1;
	int dop1 = 2;
	int dop2 = 2;
	double paraeps = Input->paraeps;
		
	/** Step 1. generate degrees of freedom */
	getElementDOF_HuangZhang3d(elementDOF, elements, elementFace, faces, elementEdge, edges);
	if(nitsche == 0)
		getFreenodesInfoHuangZhang3d(faces, edges, nodes, 1, elementDOF);
	else
		getFreenodesInfoHuangZhang3d(faces, edges, nodes, 2, elementDOF);
	getTransposeOfelementDoF(elementDOF, elementdofTran, 0);
	getElementDOF_Lagrange3d(elementDOF+1, elements, elementFace, faces, elementEdge, edges, nodes->row, dop2);
	getFreenodesInfoLagrange3d(faces, edges, nodes, elementDOF+1);
	getTransposeOfelementDoF(elementDOF+1, elementdofTran+1, 0);

	/***************************Generate coefficient of basis functions**************************/
	ddenmat3 basisCoeffs;
	create_dden_matrix3(elements->row, elementDOF->col, elementDOF->col, &basisCoeffs);
	getHuangZhangBasisCoeffs(&basisCoeffs, elements, elementFace, faces, elementEdge, edges);
	/********************************************************************************************/

	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_quadcurlperturbHuangZhang3d(&A, &b, &uh, paraeps, nitsche, &basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, elementdofTran);
	for(i=0;i<2;i++)
	{
		free(elementdofTran[i].IA);
		free(elementdofTran[i].JA);
	}

//	print_dcsr_matrix(&A);///////////
//	print_dvector(0, &b);///////
	//print_darray(20, A.val);///
	
	/** Step 3. Check matrix properties */

//	print_dcsr_matrix(&A);
//	print_dvector(b.row, &b);

	/* output A, b to a diskfile */
	// char *outputfileA="output/A.dat";
	// char *outputfileb="output/b.dat";
	// write_IJ_dCSRmat4Matlab(&A, outputfileA); 
	// write_IJ_dvector4Matlab(&b, outputfileb); 
	
	check_symm(&A);
	check_diagpos(&A);
	check_diagdom(&A);

//	exit(1);/////////


	/** Step 4. Solve the system */

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", Input->itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", Input->itsolver_tol);
	}

	create_dvector(b.row, &_uh);
	init_dvector(&_uh, 0.0);

	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A.row, A.col, A.nnz);

	AMG_param amgparam; /* parameters for AMG */

	amgparam.print_level = Input->print_level;
	amgparam.max_levels = Input->AMG_levels;
	amgparam.max_iter = Input->itsolver_maxit;
	amgparam.tol = Input->itsolver_tol;
	amgparam.AMG_max_iter = Input->MG_maxit;
	amgparam.AMG_tol = Input->MG_tol;

	amgparam.smoother = Input->MG_smoother;
	amgparam.presmooth_iter = Input->MG_smooth_iter;
	amgparam.postsmooth_iter = Input->MG_smooth_iter;

	amgparam.coarsening_type = Input->AMG_coarsening_type;
	amgparam.interpolation_type = Input->AMG_interpolation_type;
	amgparam.coarse_dof = Input->AMG_coarse_dof;
	amgparam.strong_threshold = Input->AMG_strong_threshold;
	amgparam.truncation_threshold = Input->AMG_truncation_threshold;
	amgparam.max_row_sum = Input->AMG_max_row_sum;

	/* Solve Ax = b by Matlab */
	if (itsolver_type == 0) {
		/* output A, b to a diskfile */
		char *outputfileA="output/A.dat";
		char *outputfileb="output/b.dat";
		char *outputfileuh="output/uh.dat";
		write_IJ_dCSRmat4Matlab(&A, outputfileA); 
		write_dvector4Matlab(&b, outputfileb);
		printf("Wait for Matlab solver\n");
		char c;
		while ( (c = getchar()) != '\n' && c != EOF ) ;

		read_dvector4Matlab(&_uh, outputfileuh);
	}

	/* AMG solver */
	if (itsolver_type == 1) {
		printf("AMG preconditioned CG solver\n");
		classicAMG_PCG(&A, &b, &_uh, &amgparam, print_level);
	}

	
	/* PCG+diag */
	else if (itsolver_type == 2) {
		printf("Diagonal preconditioned CG solver\n");
		diag_PCG(&A, &b, &_uh, itsolver_maxit, itsolver_tol, print_level);
	}

	/* PCG+AMG */
	else if (itsolver_type == 3) {
		printf("AMG iterative solver\n");
		classicAMG(&A, &b, &_uh, &amgparam);
	}

	/* CG */
	else if (itsolver_type == 4) {
		printf("Classical CG solver\n");
		standard_CG(&A, &b, &_uh, itsolver_maxit, itsolver_tol, print_level);
	}

	/* GMRES+AMG */
	else if (itsolver_type == 5) {
		printf("AMG preconditioned GMRES solver\n");
		classicAMG_GMRES(&A, &b, &_uh, &amgparam, print_level);
	}

	for (i = 0; i < _uh.row; i++)
		uh.val[elementDOF[0].freenodes.val[i]] = _uh.val[i];

	free_dvector(&_uh);

	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */

	/** Step 5. Compute the error between numerical solution and exact solution */
	double errors[4];

	geterrorsQuadcurlperturbHuangZhang3d(errors, paraeps, nitsche, &uh, &basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);

	printf("\nA Priori Errors:\n");
	printf("L2 norm of u-u_h = %e\n", errors[0]);
	printf("L2 norm of curl(u-u_h) = %e\n", errors[1]);
	printf("H1 semi-norm of curl(u-u_h) = %e\n", errors[2]);
	printf("Energy norm of u-u_h = %e\n", errors[3]);

	/*********************************************************************************************/
	FILE *outputFile;
	outputFile = fopen("output/error.dat", "w");
	fprintf(outputFile, "%e\t%e\t%e\t%e\n", errors[0], errors[1], errors[2], errors[3]);
	fclose(outputFile);
	/********************************************************************************************/

	
	free_dden_matrix3(&basisCoeffs);	
	free_csr_matrix(&A);
	free_dvector(&b);
	free_dvector(&uh);
	free_elementDOF(elementDOF);
	free_elementDOF(elementDOF+1);
}

/**
 * \fn void assemble_quadcurlperturbHuangZhang3d(dCSRmat *A, dvector *b, dvector *uh, double paraeps, short nitsche, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran)
 * \brief assemble stiffness matrix *C and righ hand side *ptr_b
 * \param *ptr_A pointer to stiffness matrix
 * \param *ptr_b pointer to right hand side
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
 * \param lambda Lame constant
 * \param mu Lame constant or Poisson ratio of plate
 * \return void
 */
void assemble_quadcurlperturbHuangZhang3d(dCSRmat *A, dvector *b, dvector *uh, double paraeps, short nitsche, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Dx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	dCSRmat AA, BB, CC, tempA;
	dCSRmat B, BT, C;
	dvector bb, D;

	assembleBiGradcurlperturbHuangZhang3d(&AA, paraeps, nitsche, basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, elementdofTran);
	assembleHuangZhangGradLagrange3d(&BB, basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, elementdofTran);
	assembleMassmatrixLagrange3d(&CC, elements, elementDOF+1, elementdofTran+1);
    // assembleRHSHuangZhang3d(&bb, basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, elementdofTran, quadcurl3d_f);
	assembleRHSHuangZhang3d(&bb, basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, elementdofTran, maxwell3d_f);
    // initial solution
	create_dvector(bb.row, uh);
	
	// extract
	extractFreenodesVector(&AA, &bb, b, elementDOF, uh);
	free_dvector(&bb);
	extractFreenodesMatrix11(&AA, &tempA, elementDOF, elementDOF);
	free_csr_matrix(&AA);
	extractFreenodesMatrix11(&BB, &B, elementDOF+1, elementDOF);
	free_csr_matrix(&BB);
	extractFreenodesMatrix11(&CC, &C, elementDOF+1, elementDOF+1);
	free_csr_matrix(&CC);
	getdiag(C.row, &C, &D);
	free_csr_matrix(&C);
	getTransposeOfSparse(&B, &BT);
	
	// reduce: A = tempA + BT Dinv B
	dDiagVectorInvMultiplydCSR(&D, &B, &BB);
	free_dvector(&D);
	free_csr_matrix(&B);
	sparseMultiplication(&BT, &BB, &AA);
	free_csr_matrix(&BT);
	free_csr_matrix(&BB);
	sparseAddition(&tempA, &AA, A);
	free_csr_matrix(&tempA);
	free_csr_matrix(&AA);

/**********************************
	double EPS = 1e-300;
	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A.row, A.col, A.nnz);
	compress_dcsr(&A, ptr_A, EPS);
	free_csr_matrix(&A);
	printf("A.row=%d, A.col=%d, A.nnz=%d\n", ptr_A->row, ptr_A->col, ptr_A->nnz);
	*******************************************/
}


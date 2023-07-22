/*
 *  poissonfem.c
 *
 *  Created by Xuehai Huang on Apr 30, 2022.
 *  Copyright 2022 SUFE. All rights reserved.
 *
 */

/*! \file poissonfem.c
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
 * \fn void poissonfem2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
 * \brief finite element methods for Poisson equation in two dimensions
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
void poissonfem2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i;
	int variationalform_type = Input->variationalform_type;
	int fem_num = Input->fem_num;
	if (variationalform_type == 1){ // primal formulation
		printf("Discretize Poisson equation in primal formulation.\n");
		if(fem_num == 1){
			printf("Discretized by Lagrange element.\n");
			poissonLagrange2d(elements, elementEdge, edges, nodes, Input);
		}
		
		if(fem_num == 2){
			printf("Discretized by Crouzeix-Raviart element.\n");
			poissonCrouzeixRaviart2d(elements, elementEdge, edges, nodes, Input);
		}
	}
	else if (variationalform_type == 2){ // mixed formulation
		printf("Discretize Poission equation in mixed formulation.\n");
		if(fem_num == 1){
			printf("Discretized by the Raviart-Thomas element with k = %d.\n", Input->dop1);
			poissonRaviartThomas2d(elements, elementEdge, edges, nodes, Input);
		}

		if(fem_num == 2){
			printf("Discretized by the Brezzi-Douglas-Marini element with k = %d.\n", Input->dop1);
			poissonBrezziDouglasMarini2d(elements, elementEdge, edges, nodes, Input);
		}
	}
	else{
		printf("Please set variationalform_type = 1 or 2!\n");
		exit(0);
	}
}

/**
* \fn void poissonLagrange2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief Lagrange element method for Poisson equation in two dimensions
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
void poissonLagrange2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i,j;
	dCSRmat A;
	dvector b, uh;
	ELEMENT_DOF elementDOF;
	// iCSRmat elementdofTran;
		
	/** Step 0. Read input parameters */
	int print_level = Input->print_level;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;
	double lambda = Input->lambda;
	double mu = Input->mu;
	int glevelNum = Input->glevelNum;
	int CglevelNum = Input->CglevelNum;
	int FglevelNum = Input->FglevelNum;
	int dop = Input->dop1;
		
	printf("k = %d\n", dop);
	
	/** Step 1. generate degrees of freedom */
	getElementDOF_Lagrange2d(&elementDOF, elements, elementEdge, edges, nodes->row, dop);
	getFreenodesInfoLagrange2d(edges, nodes, &elementDOF);
	// getTransposeOfelementDoF(&elementDOF, &elementdofTran, 0);
		
	/** Step 2. assemble stiffmatrix and right hand side term */
	// create_dvector(b.row, &uh);
	// init_dvector(&uh, 0.0);
	assemble_poissonLagrange2d(&A, &b, &uh, elements, elementEdge, edges, nodes, &elementDOF);
	// free(elementdofTran.IA);
	// free(elementdofTran.JA);

//	print_dcsr_matrix(&A);///////////
//	print_dvector(0, &b);///////
	//print_darray(20, A.val);///
	
	/** Step 3. Check matrix properties */

	check_symm(&A);
	check_diagpos(&A);
	check_diagdom(&A);


	/** Step 4. Solve the system */

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", Input->itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", Input->itsolver_tol);
	}

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

	/* PCG+AMG */
	if (itsolver_type == 1) {
		printf("AMG preconditioned CG solver\n");
		classicAMG_PCG(&A, &b, &uh, &amgparam, print_level);
	}

	/* AMG solver */
	else if (itsolver_type == 2) {
		printf("AMG iterative solver\n");
		classicAMG(&A, &b, &uh, &amgparam);
	}

	/* PCG+diag */
	else if (itsolver_type == 3) {
		printf("Diagonal preconditioned CG solver\n");
		diag_PCG(&A, &b, &uh, itsolver_maxit, itsolver_tol, print_level);
	}

	/* CG */
	else if (itsolver_type == 4) {
		printf("Classical CG solver\n");
		standard_CG(&A, &b, &uh, itsolver_maxit, itsolver_tol, print_level);
	}

	/* GMRES+AMG */
	else if (itsolver_type == 5) {
		printf("AMG preconditioned GMRES solver\n");
		classicAMG_GMRES(&A, &b, &uh, &amgparam, print_level);
	}

	// for (i = 0; i < _uh.row; i++)
	// 	uh.val[elementDOF.freenodes.val[i]] = _uh.val[i];

	// free_dvector(&_uh);

	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */

	/** Step 5. Compute the error between numerical solution and exact solution */
	double errors[3];

	geterrorsPoissonLagrange2d(errors, &uh, elements, elementEdge, edges, nodes, &elementDOF);

	printf("\nA Priori Errors:\n");
	printf("L2 norm of u-u_h = %e\n", errors[0]);
	printf("H1 semi-norm of u-u_h = %e\n", errors[1]);
	printf("H1 norm of u-u_h = %e\n", errors[2]);
	
	/*********************************************************************************************/
	FILE *outputFile;
	outputFile = fopen("output/error.dat", "w");
	fprintf(outputFile, "%e\t%e\t%e\n", errors[0], errors[1], errors[2]);
	fclose(outputFile);
	/********************************************************************************************/

	
	free_csr_matrix(&A);
	free_dvector(&b);
	free_dvector(&uh);
	free_elementDOF(&elementDOF);
}

/**
 * \fn void assemble_poissonLagrange2d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
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
 * \return void
 */
void assemble_poissonLagrange2d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	/**	
	Ax = b
	**/
	dCSRmat AA;

	assembleBiGradLagrange2d(&AA, elements, elementEdge, edges, nodes, elementDOF, 0.5);
	assembleRHSLagrange2d(b, elements, elementDOF, poisson2d_f, NULL);
    // initial solution
	create_dvector(b->row, uh);
	
	// Apply boundary condition
	updateFreenodesRHS(&AA, b, uh, elementDOF, elementDOF, 1);
	updateFreenodesMatrix11(&AA, A, elementDOF, elementDOF, 1);

	// /* output A, b to a diskfile */
	// char *outputfileAA="output/AA.dat";
	// char *outputfileA="output/A.dat";
	// char *outputfileb="output/b.dat";
	// write_IJ_dCSRmat4Matlab(&AA, outputfileAA); 
	// write_IJ_dCSRmat4Matlab(A, outputfileA); 
	// write_dvector4Matlab(b, outputfileb); 

// printf("AA\n");
// print_dcsr_matrix(&AA);///////////
// printf("A, row=%d, col=%d, nnz=%d\n", A->row, A->col, A->nnz);
// print_dcsr_matrix(A);///////////

	free_csr_matrix(&AA);
}

/**
* \fn void poissonCrouzeixRaviart2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief Crouzeix-Raviart element method for Poisson equation in two dimensions
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
void poissonCrouzeixRaviart2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i,j;
	dCSRmat A;
	dvector b, uh;
	ELEMENT_DOF elementDOF;
	// iCSRmat elementdofTran;
		
	/** Step 0. Read input parameters */
	int print_level = Input->print_level;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;
	double lambda = Input->lambda;
	double mu = Input->mu;
	int glevelNum = Input->glevelNum;
	int CglevelNum = Input->CglevelNum;
	int FglevelNum = Input->FglevelNum;
	int dop = Input->dop1;
		
	printf("k = %d\n", dop);
	
	/** Step 1. generate degrees of freedom */
	getElementDOF_CrouzeixRaviart2d(&elementDOF, elements, elementEdge, edges);
	getFreenodesInfoCrouzeixRaviart2d(edges, &elementDOF);
	// getTransposeOfelementDoF(&elementDOF, &elementdofTran, 0);
		
	/** Step 2. assemble stiffmatrix and right hand side term */
	// create_dvector(b.row, &uh);
	// init_dvector(&uh, 0.0);
	assemble_poissonCrouzeixRaviart2d(&A, &b, &uh, elements, elementEdge, edges, nodes, &elementDOF);
	// free(elementdofTran.IA);
	// free(elementdofTran.JA);

//	print_dcsr_matrix(&A);///////////
//	print_dvector(0, &b);///////
	//print_darray(20, A.val);///
	
	/** Step 3. Check matrix properties */

	check_symm(&A);
	check_diagpos(&A);
	check_diagdom(&A);


	/** Step 4. Solve the system */

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", Input->itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", Input->itsolver_tol);
	}

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

	/* PCG+AMG */
	if (itsolver_type == 1) {
		printf("AMG iterative solver\n");
		classicAMG_PCG(&A, &b, &uh, &amgparam, print_level);
	}

	/* AMG solver */
	else if (itsolver_type == 2) {
		printf("AMG preconditioned CG solver\n");
		classicAMG(&A, &b, &uh, &amgparam);
	}

	/* PCG+diag */
	else if (itsolver_type == 3) {
		printf("Diagonal preconditioned CG solver\n");
		diag_PCG(&A, &b, &uh, itsolver_maxit, itsolver_tol, print_level);
	}

	/* CG */
	else if (itsolver_type == 4) {
		printf("Classical CG solver\n");
		standard_CG(&A, &b, &uh, itsolver_maxit, itsolver_tol, print_level);
	}

	/* GMRES+AMG */
	else if (itsolver_type == 5) {
		printf("AMG preconditioned GMRES solver\n");
		classicAMG_GMRES(&A, &b, &uh, &amgparam, print_level);
	}

	// for (i = 0; i < _uh.row; i++)
	// 	uh.val[elementDOF.freenodes.val[i]] = _uh.val[i];

	// free_dvector(&_uh);

	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */

	/** Step 5. Compute the error between numerical solution and exact solution */
	double errors[3];

	geterrorsPoissonCrouzeixRaviart2d(errors, &uh, elements, elementEdge, edges, nodes, &elementDOF);

	printf("\nA Priori Errors:\n");
	printf("L2 norm of u-u_h = %e\n", errors[0]);
	printf("H1 semi-norm of u-u_h = %e\n", errors[1]);
	printf("H1 norm of u-u_h = %e\n", errors[2]);
	
	/*********************************************************************************************/
	FILE *outputFile;
	outputFile = fopen("output/error.dat", "w");
	fprintf(outputFile, "%e\t%e\t%e\n", errors[0], errors[1], errors[2]);
	fclose(outputFile);
	/********************************************************************************************/

	
	free_csr_matrix(&A);
	free_dvector(&b);
	free_dvector(&uh);
	free_elementDOF(&elementDOF);
}

/**
 * \fn void assemble_poissonCrouzeixRaviart2d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
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
 * \return void
 */
void assemble_poissonCrouzeixRaviart2d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	/**	
	Ax = b
	**/
	dCSRmat AA;

	assembleBiGradCrouzeixRaviart2d(&AA, elements, elementEdge, edges, nodes, elementDOF, 0.5);
	assembleRHSCrouzeixRaviart2d(b, elements, elementDOF, poisson2d_f, NULL);
    // initial solution
	create_dvector(b->row, uh);
	
	// Apply boundary condition
	updateFreenodesRHS(&AA, b, uh, elementDOF, elementDOF, 1);
	updateFreenodesMatrix11(&AA, A, elementDOF, elementDOF, 1);

	// /* output A, b to a diskfile */
	// char *outputfileAA="output/AA.dat";
	// char *outputfileA="output/A.dat";
	// char *outputfileb="output/b.dat";
	// write_IJ_dCSRmat4Matlab(&AA, outputfileAA); 
	// write_IJ_dCSRmat4Matlab(A, outputfileA); 
	// write_dvector4Matlab(b, outputfileb); 

// printf("AA\n");
// print_dcsr_matrix(&AA);///////////
// printf("A, row=%d, col=%d, nnz=%d\n", A->row, A->col, A->nnz);
// print_dcsr_matrix(A);///////////

	free_csr_matrix(&AA);
}

/**
* \fn void poissonRaviartThomas2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief Raviart-Thomas element method for Poisson equation in two dimensions
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
void poissonRaviartThomas2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i,j;
	dCSRmat A[4];
	dvector b[2], uh[2];
	ELEMENT_DOF elementDOF[2];
	// iCSRmat elementdofTran;
		
	/** Step 0. Read input parameters */
	int print_level = Input->print_level;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;
	double lambda = Input->lambda;
	double mu = Input->mu;
	int glevelNum = Input->glevelNum;
	int CglevelNum = Input->CglevelNum;
	int FglevelNum = Input->FglevelNum;
	int dop1 = Input->dop1;
	if(dop1<1) dop1 = 1;
	int dop2 = dop1-1;
		
	printf("k = %d\n", dop1);
	
	/** Step 1. generate degrees of freedom */
	getElementDOF_RaviartThomas2d(elementDOF, elements, elementEdge, edges, dop1);
	getElementDOFdg(elementDOF+1, elements->row, (dop2+1)*(dop2+2)/2, dop2);
	// getFreenodesInfoLagrange2d(edges, nodes, &elementDOF);
	// getTransposeOfelementDoF(&elementDOF, &elementdofTran, 0);
		
	/** Step 2. assemble stiffmatrix and right hand side term */
	// create_dvector(b.row, &uh);
	// init_dvector(&uh, 0.0);
	assemble_poissonRaviartThomas2d(A, b, uh, elements, elementEdge, edges, nodes, elementDOF);
	// free(elementdofTran.IA);
	// free(elementdofTran.JA);

//////////////////////

//	print_dcsr_matrix(&A);///////////
//	print_dvector(0, &b);///////
	//print_darray(20, A.val);///
	
	/** Step 3. Check matrix properties */
	for (i = 0; i<3; i++){
		check_symm(&A[i]);
		check_diagpos(&A[i]);
		check_diagdom(&A[i]);
	}

	/** Step 4. Solve the system */

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", Input->itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", Input->itsolver_tol);
	}

	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A[0].row, A[0].col, A[0].nnz);
	printf("B.row=%d, B.col=%d, B.nnz=%d\n", A[2].row, A[2].col, A[2].nnz);

	AMG_param amgparam; /* parameters for AMG */

	amgparam.print_level = Input->print_level;
	amgparam.restart = Input->restart;
	amgparam.max_iter = Input->itsolver_maxit;
	amgparam.tol = Input->itsolver_tol;
	amgparam.AMG_max_iter = Input->MG_maxit;
	amgparam.AMG_tol = Input->MG_tol;

	amgparam.max_levels = Input->AMG_levels;
	amgparam.smoother = Input->MG_smoother;
	amgparam.presmooth_iter = Input->MG_smooth_iter;
	amgparam.postsmooth_iter = Input->MG_smooth_iter;

	amgparam.coarsening_type = Input->AMG_coarsening_type;
	amgparam.interpolation_type = Input->AMG_interpolation_type;
	amgparam.coarse_dof = Input->AMG_coarse_dof;
	amgparam.strong_threshold = Input->AMG_strong_threshold;
	amgparam.truncation_threshold = Input->AMG_truncation_threshold;
	amgparam.max_row_sum = Input->AMG_max_row_sum;

	ASP_param aspparam; /* parameters for AMG */

	aspparam.print_level = Input->print_level;
	aspparam.max_iter = Input->itsolver_maxit;
	aspparam.tol = Input->itsolver_tol;
	aspparam.restart = Input->restart;

	aspparam.precond_type = Input->precond_type;
	aspparam.precond_scale = Input->precond_scale;
	aspparam.smoother = Input->smoother;
	aspparam.smooth_iter = Input->smooth_iter;
	
	aspparam.levelNum = glevelNum;
	aspparam.mg_max_iter = Input->MG_maxit;
	aspparam.mg_tol = Input->MG_tol;
	aspparam.mg_smoother = Input->MG_smoother;
	aspparam.mg_smooth_iter = Input->MG_smooth_iter;

	aspparam.elements = elements;
	aspparam.elementEdge = elementEdge;
	aspparam.edges = edges;
	aspparam.nodes = nodes;
	// aspparam.edgesTran = edgesTran;
	// aspparam.nodeCEdge = nodeCEdge;
	aspparam.elementDOF = elementDOF;
	// aspparam.elementdofTran = &elementdofTran;

	aspparam.lambda = Input->lambda;
	aspparam.mu = Input->mu;

	aspparam.max_levels = Input->AMG_levels;
	aspparam.AMG_coarsening_type = Input->AMG_coarsening_type;
	aspparam.AMG_interpolation_type = Input->AMG_interpolation_type;
	aspparam.AMG_coarse_dof = Input->AMG_coarse_dof;
	aspparam.AMG_strong_threshold = Input->AMG_strong_threshold;
	aspparam.AMG_truncation_threshold = Input->AMG_truncation_threshold;
	aspparam.AMG_max_row_sum = Input->AMG_max_row_sum;

	if (itsolver_type == 1){
		printf("Block diagonal preconditioned MINRES solver with auxiliary space method\n");
		DiagAsP1PoissonDG_MINRES(A, b, uh, &aspparam, print_level);
	}
	else{
		printf("Block triangular preconditioned GMRES solver with auxiliary space method\n");
		TriAsP1PoissonDG_GMRES(A, b, uh, &aspparam, print_level);
	}

	// for (i = 0; i < _uh.row; i++)
	// 	uh.val[elementDOF.freenodes.val[i]] = _uh.val[i];

	// free_dvector(&_uh);

	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */

	/** Postprocessing */
	dvector Qhu, uhstar;
	create_dvector(uh[1].row, &Qhu);
	printf("Postprocessing...\n");
	projL2Pk2d(&Qhu, elements, nodes, &elementDOF[1], poisson2d_u, NULL);
	printf("Postprocessing finish!\n");


	/** Step 5. Compute the error between numerical solution and exact solution */
	double errors[6];

	geterrorsPoissonRaviartThomas2d(errors, uh, &Qhu, elements, elementEdge, edges, nodes, elementDOF);

	printf("\nA Priori Errors:\n");
	printf("L2 norm of sigma - sigma_h = %e\n", errors[0]);
	printf("L2 norm of div(sigma - sigma_h) = %e\n", errors[1]);
	printf("L2 norm of u-u_h = %e\n", errors[2]);
	printf("H1 semi-norm of u-u_h = %e\n", errors[3]);
	printf("L2 norm of Qhu-u_h = %e\n", errors[4]);
	printf("H1 semi-norm of Qhu-u_h = %e\n", errors[5]);
	
	/*********************************************************************************************/
	FILE *outputFile;
	outputFile = fopen("output/error.dat", "w");
	fprintf(outputFile, "%e\t%e\t%e\t%e\t%e\t%e\n", errors[0], errors[1], errors[2], errors[3], errors[4], errors[5]);
	fclose(outputFile);
	/********************************************************************************************/

	for(i=0;i<3;i++)
		free_csr_matrix(A+i);
	for(i=0;i<2;i++){
		free_dvector(b+i);
		free_dvector(uh+i);
		free_elementDOF(elementDOF+i);
	}
	free_dvector(&Qhu);
}

/**
 * \fn void assemble_poissonRaviartThomas2d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
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
 * \return void
 */
void assemble_poissonRaviartThomas2d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	/**	A[0]=A, A[1]=B^T, A[2]=B
	Ax1 + B^Tx2 = 0
	Bx1 +       = b
	where x1: sigma_h, x2: u_h
	**/

	assembleMassmatrixRaviartThomas2d(A, elements, elementEdge, edges, elementDOF);
	assembleDivRaviartThomasL2poly2d(A+1, elements, elementEdge, edges, elementDOF, elementDOF+1);
	getTransposeOfSparse(A+1, A+2);
	A[3].row = 0;
	A[3].col = 0;
	A[3].IA = NULL;
	A[3].JA = NULL;
	A[3].val = NULL;
	create_dvector(A->row, b);
	assembleRHSdgPoly2d(b+1, elements, elementDOF+1, poisson2d_f, NULL);
	ax_array(b[1].row, -1, b[1].val);

    // initial solution
	create_dvector(b[0].row, &uh[0]);
	create_dvector(b[1].row, &uh[1]);
	
	// Apply boundary condition
	// updateFreenodesRHS(&AA, b, uh, elementDOF, elementDOF, 1);
	// updateFreenodesMatrix11(&AA, A, elementDOF, elementDOF, 1);

	// /* output A, b to a diskfile */
	// char *outputfileAA="output/AA.dat";
	// char *outputfileA="output/A.dat";
	// char *outputfileb="output/b.dat";
	// write_IJ_dCSRmat4Matlab(&AA, outputfileAA); 
	// write_IJ_dCSRmat4Matlab(A, outputfileA); 
	// write_dvector4Matlab(b, outputfileb); 

// printf("AA\n");
// print_dcsr_matrix(&AA);///////////
// printf("A, row=%d, col=%d, nnz=%d\n", A->row, A->col, A->nnz);
// print_dcsr_matrix(A);///////////

	// free_csr_matrix(&AA);
}

/**
* \fn void poissonBrezziDouglasMarini2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief Brezzi–Douglas–Marini element method for Poisson equation in two dimensions
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
void poissonBrezziDouglasMarini2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i,j;
	dCSRmat A[4];
	dvector b[2], uh[2];
	ELEMENT_DOF elementDOF[2];
	// iCSRmat elementdofTran;
		
	/** Step 0. Read input parameters */
	int print_level = Input->print_level;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;
	double lambda = Input->lambda;
	double mu = Input->mu;
	int glevelNum = Input->glevelNum;
	int CglevelNum = Input->CglevelNum;
	int FglevelNum = Input->FglevelNum;
	int dop1 = Input->dop1;
	if(dop1<1) dop1 = 1;
	int dop2 = dop1-1;
		
	printf("k = %d\n", dop1);
	
	/** Step 1. generate degrees of freedom */
	getElementDOF_BrezziDouglasMarini2d(elementDOF, elements, elementEdge, edges, dop1);
	getElementDOFdg(elementDOF+1, elements->row, (dop2+1)*(dop2+2)/2, dop2);
	// getFreenodesInfoLagrange2d(edges, nodes, &elementDOF);
	// getTransposeOfelementDoF(&elementDOF, &elementdofTran, 0);
		
	/** Step 2. assemble stiffmatrix and right hand side term */
	// create_dvector(b.row, &uh);
	// init_dvector(&uh, 0.0);
	assemble_poissonBrezziDouglasMarini2d(A, b, uh, elements, elementEdge, edges, nodes, elementDOF);
	// free(elementdofTran.IA);
	// free(elementdofTran.JA);

//////////////////////

//	print_dcsr_matrix(&A);///////////
//	print_dvector(0, &b);///////
	//print_darray(20, A.val);///
	
	/** Step 3. Check matrix properties */
	for (i = 0; i<3; i++){
		check_symm(&A[i]);
		check_diagpos(&A[i]);
		check_diagdom(&A[i]);
	}

	/** Step 4. Solve the system */

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", Input->itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", Input->itsolver_tol);
	}

	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A[0].row, A[0].col, A[0].nnz);
	printf("B.row=%d, B.col=%d, B.nnz=%d\n", A[2].row, A[2].col, A[2].nnz);

	AMG_param amgparam; /* parameters for AMG */

	amgparam.print_level = Input->print_level;
	amgparam.restart = Input->restart;
	amgparam.max_iter = Input->itsolver_maxit;
	amgparam.tol = Input->itsolver_tol;
	amgparam.AMG_max_iter = Input->MG_maxit;
	amgparam.AMG_tol = Input->MG_tol;

	amgparam.max_levels = Input->AMG_levels;
	amgparam.smoother = Input->MG_smoother;
	amgparam.presmooth_iter = Input->MG_smooth_iter;
	amgparam.postsmooth_iter = Input->MG_smooth_iter;

	amgparam.coarsening_type = Input->AMG_coarsening_type;
	amgparam.interpolation_type = Input->AMG_interpolation_type;
	amgparam.coarse_dof = Input->AMG_coarse_dof;
	amgparam.strong_threshold = Input->AMG_strong_threshold;
	amgparam.truncation_threshold = Input->AMG_truncation_threshold;
	amgparam.max_row_sum = Input->AMG_max_row_sum;

	ASP_param aspparam; /* parameters for AMG */

	aspparam.print_level = Input->print_level;
	aspparam.max_iter = Input->itsolver_maxit;
	aspparam.tol = Input->itsolver_tol;
	aspparam.restart = Input->restart;

	aspparam.precond_type = Input->precond_type;
	aspparam.precond_scale = Input->precond_scale;
	aspparam.smoother = Input->smoother;
	aspparam.smooth_iter = Input->smooth_iter;
	
	aspparam.levelNum = glevelNum;
	aspparam.mg_max_iter = Input->MG_maxit;
	aspparam.mg_tol = Input->MG_tol;
	aspparam.mg_smoother = Input->MG_smoother;
	aspparam.mg_smooth_iter = Input->MG_smooth_iter;

	aspparam.elements = elements;
	aspparam.elementEdge = elementEdge;
	aspparam.edges = edges;
	aspparam.nodes = nodes;
	// aspparam.edgesTran = edgesTran;
	// aspparam.nodeCEdge = nodeCEdge;
	aspparam.elementDOF = elementDOF;
	// aspparam.elementdofTran = &elementdofTran;

	aspparam.lambda = Input->lambda;
	aspparam.mu = Input->mu;

	aspparam.max_levels = Input->AMG_levels;
	aspparam.AMG_coarsening_type = Input->AMG_coarsening_type;
	aspparam.AMG_interpolation_type = Input->AMG_interpolation_type;
	aspparam.AMG_coarse_dof = Input->AMG_coarse_dof;
	aspparam.AMG_strong_threshold = Input->AMG_strong_threshold;
	aspparam.AMG_truncation_threshold = Input->AMG_truncation_threshold;
	aspparam.AMG_max_row_sum = Input->AMG_max_row_sum;

	if (itsolver_type == 1){
		printf("Block diagonal preconditioned MINRES solver with auxiliary space method\n");
		DiagAsP1PoissonDG_MINRES(A, b, uh, &aspparam, print_level);
	}
	else{
		printf("Block triangular preconditioned GMRES solver with auxiliary space method\n");
		TriAsP1PoissonDG_GMRES(A, b, uh, &aspparam, print_level);
	}

	// for (i = 0; i < _uh.row; i++)
	// 	uh.val[elementDOF.freenodes.val[i]] = _uh.val[i];

	// free_dvector(&_uh);

	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */

	/** Postprocessing */
	dvector Qhu, uhstar;
	create_dvector(uh[1].row, &Qhu);
	printf("Postprocessing...\n");
	projL2Pk2d(&Qhu, elements, nodes, &elementDOF[1], poisson2d_u, NULL);
	printf("Postprocessing finish!\n");


	/** Step 5. Compute the error between numerical solution and exact solution */
	double errors[6];

	geterrorsPoissonBrezziDouglasMarini2d(errors, uh, &Qhu, elements, elementEdge, edges, nodes, elementDOF);

	printf("\nA Priori Errors:\n");
	printf("L2 norm of sigma - sigma_h = %e\n", errors[0]);
	printf("L2 norm of div(sigma - sigma_h) = %e\n", errors[1]);
	printf("L2 norm of u-u_h = %e\n", errors[2]);
	printf("H1 semi-norm of u-u_h = %e\n", errors[3]);
	printf("L2 norm of Qhu-u_h = %e\n", errors[4]);
	printf("H1 semi-norm of Qhu-u_h = %e\n", errors[5]);
	
	/*********************************************************************************************/
	FILE *outputFile;
	outputFile = fopen("output/error.dat", "w");
	fprintf(outputFile, "%e\t%e\t%e\t%e\t%e\t%e\n", errors[0], errors[1], errors[2], errors[3], errors[4], errors[5]);
	fclose(outputFile);
	/********************************************************************************************/

	for(i=0;i<3;i++)
		free_csr_matrix(A+i);
	for(i=0;i<2;i++){
		free_dvector(b+i);
		free_dvector(uh+i);
		free_elementDOF(elementDOF+i);
	}
	free_dvector(&Qhu);
}

/**
 * \fn void assemble_poissonBrezziDouglasMarini2d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
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
 * \return void
 */
void assemble_poissonBrezziDouglasMarini2d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	/**	A[0]=A, A[1]=B^T, A[2]=B
	Ax1 + B^Tx2 = 0
	Bx1 +       = b
	where x1: sigma_h, x2: u_h
	**/

	assembleMassmatrixBrezziDouglasMarini2d(A, elements, elementEdge, edges, elementDOF);
	assembleDivBrezziDouglasMariniL2poly2d(A+1, elements, elementEdge, edges, elementDOF, elementDOF+1);
	getTransposeOfSparse(A+1, A+2);
	A[3].row = 0;
	A[3].col = 0;
	A[3].IA = NULL;
	A[3].JA = NULL;
	A[3].val = NULL;
	create_dvector(A->row, b);
	assembleRHSdgPoly2d(b+1, elements, elementDOF+1, poisson2d_f, NULL);
	ax_array(b[1].row, -1, b[1].val);

    // initial solution
	create_dvector(b[0].row, &uh[0]);
	create_dvector(b[1].row, &uh[1]);
	
	// Apply boundary condition
	// updateFreenodesRHS(&AA, b, uh, elementDOF, elementDOF, 1);
	// updateFreenodesMatrix11(&AA, A, elementDOF, elementDOF, 1);

	// /* output A, b to a diskfile */
	// char *outputfileAA="output/AA.dat";
	// char *outputfileA="output/A.dat";
	// char *outputfileb="output/b.dat";
	// write_IJ_dCSRmat4Matlab(&AA, outputfileAA); 
	// write_IJ_dCSRmat4Matlab(A, outputfileA); 
	// write_dvector4Matlab(b, outputfileb); 

// printf("AA\n");
// print_dcsr_matrix(&AA);///////////
// printf("A, row=%d, col=%d, nnz=%d\n", A->row, A->col, A->nnz);
// print_dcsr_matrix(A);///////////

	// free_csr_matrix(&AA);
}

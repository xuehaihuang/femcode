/*
 *  biharmonicfem.c
 *
 *  Created by Xuehai Huang on Apr 30, 2022.
 *  Copyright 2022 SUFE. All rights reserved.
 *
 */

/*! \file biharmonicfem.c
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
 * \fn void triharmonicfem2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
 * \brief finite element methods for triharmonic equation in two dimensions
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
void triharmonicfem2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i;
	int variationalform_type = Input->variationalform_type;
	int fem_num = Input->fem_num;
	int glevelNum = Input->glevelNum;
	short extrap = Input->extrap; // for extrapolation

	if (variationalform_type == 1) // primal formulation
	{
		printf("Discretize triharmonic equation in primal formulation.\n");
		printf("To be developed!\n");
		exit(1);
		// if(fem_num == 1)
		// {
		// 	printf("Discretized by Morley element.\n");
		// 	// biharmonicMorley2d(elements, elementEdge, edges, nodes, Input);
		// 	biharmonicMorley2d(&elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], Input);
		// }
		
	}
	else if (variationalform_type == 2) // decoupled mixed formulation
	{
		printf("Discretize triharmonic equation in decoupled Biharominc-symTensorStokes-Biharominc formulation.\n");
		if(fem_num == 1)
		{
			printf("Discretized by Morley : Enriched C0P1-C0P1 : Morley element.\n");
			// biharmonicPSP_MorleyNcfmP1P0Morley2d(&elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], Input);
			triharmonicBSB_MorleyAnHuangMorley2d(&elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], Input);
		}
	}
	else
	{
		printf("Please set variationalform_type = 1 or 2!\n");
		exit(1);
	}
}

/**
* \fn void triharmonicBSB_MorleyAnHuangMorley2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief Morley : AnHuang : Morley element element method for the decoupled Biharominc-symTensorStokes-Biharominc formulation of triharmonic equation in two dimensions
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
void triharmonicBSB_MorleyAnHuangMorley2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i,j;
	dvector uh, sigmah;
	ELEMENT_DOF elementDOF[5];
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
	// int dop = Input->dop1;
		
	/** Step 1. generate degrees of freedom */
	getElementDOF_Morley2d(elementDOF, elements, elementEdge, edges, nodes->row);
	getFreenodesInfoMorley2d(edges, nodes, elementDOF);
	getElementDOF_MINI2d(elementDOF+2, elements, nodes->row);
	getFreenodesInfoMINI2d(nodes, elementDOF+2);
	repeatElementDoF(elementDOF+2, elementDOF+1, 3);
	getElementDOF_Lagrange2d(elementDOF+3, elements, elementEdge, edges, nodes->row, 1);
	getFreenodesInfoDG(elementDOF+3);
	repeatElementDoF(elementDOF+3, elementDOF+4, 2);
	
	/** Step 2. assemble stiffmatrix and right hand side term */
	solve_triharmonicBSB_MorleyAnHuangMorley2d(&uh, &sigmah, elementDOF, elements, elementEdge, edges, nodes, Input);

	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */

	/** Step 3. Compute the error between numerical solution and exact solution */
	double errors[7];

	geterrorsTriharmonicMorleyMINIMorley2d(errors, &uh, &sigmah, elements, elementEdge, edges, nodes, elementDOF, elementDOF+2);
	// geterrorsBiharmonicPSP_MorleyNcfmP1P0Morley2d(errors, &uh, elements, elementEdge, edges, nodes, elementDOF);

	printf("\nA Priori Errors:\n");
	printf("L2 norm of u-u_h = %e\n", errors[0]);
	printf("H1 semi-norm of u-u_h = %e\n", errors[1]);
	printf("H2 semi-norm of u-u_h = %e\n", errors[2]);
	printf("H2 norm of u-u_h = %e\n", errors[3]);
	printf("L2 norm of sigma-sigma_h = %e\n", errors[4]);
	printf("H1 semi-norm of sigma-sigma_h = %e\n", errors[5]);
	printf("H1 norm of sigma-sigma_h = %e\n", errors[6]);
	
	/*********************************************************************************************/
	FILE *outputFile;
	outputFile = fopen("output/error.dat", "w");
	fprintf(outputFile, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", errors[0], errors[1], errors[2], errors[3], errors[4], errors[5], errors[6]);
	fclose(outputFile);
	/********************************************************************************************/

	
	// free_csr_matrix(&A);
	// free_dvector(&b);
	free_dvector(&uh);
	free_dvector(&sigmah);
	for(i=0;i<5;i++) free_elementDOF(elementDOF+i);
}

/**
* \fn void solve_triharmonicBSB_MorleyAnHuangMorley2d(dvector *uh, dvector *sigmah, ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief Morley element method for triharmonic equation in two dimensions
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
void solve_triharmonicBSB_MorleyAnHuangMorley2d(dvector *uh, dvector *sigmah, ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i,j;
	// dCSRmat A;
	// dvector b;
	// ELEMENT_DOF elementDOF;
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
	// int dop = Input->dop1;
		
	if (print_level>0) {
		printf("Maximal iteration number = %d\n", Input->itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", Input->itsolver_tol);
	}

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
	aspparam.elementDOF = elementDOF+2;
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

	dCSRmat AAm, AA[4], Am, A[4];
	dvector bb, b[2], wh[2];
	/************ Eqn1: biharmonic equation disctetized by Morley element ************/
	assembleBiHessMorley2d(&AAm, elements, elementEdge, edges, nodes, elementDOF);
	assembleRHSMorley2d(b, elements, elementEdge, edges, elementDOF, triharmonic2d_f, NULL);
    // initial solution
	create_dvector(b->row, wh);
	// Apply boundary condition
	updateFreenodesRHS(&AAm, b, wh, elementDOF, elementDOF, 1);
	updateFreenodesMatrix11(&AAm, &Am, elementDOF, elementDOF, 1);
	// free_csr_matrix(&AA[0]);
	printf("Am.row=%d, Am.col=%d, Am.nnz=%d\n", Am.row, Am.col, Am.nnz);
	printf("Solve...\n");
	printf("Solving Morley by AMG preconditioned CG solver\n");
	classicAMG_PCG(&Am, b, wh, &amgparam, print_level);
	free_dvector(&b[0]);
	
	/************ Eqn2: symmetric tensor Stokes equation disctetized by AnHuang(MINI) element ************/
	assembleBiGradMINIsymtensorDiagBlockRepeat2d(&AA[0], elements, elementEdge, edges, nodes, elementDOF+2, 0.5);
	assembleRotSMINILagrange2d(&AA[1], elements, elementDOF+2, elementDOF+3);
	getTransposeOfSparse(&AA[1], &AA[2]);
	create_csr_matrix(AA[1].col, AA[1].col, 0, &AA[3]);
	assembleRHShessMorleySMINI2d(b, elements, elementEdge, edges, elementDOF, elementDOF+2, wh);
	free_dvector(wh);
	create_dvector(AA[1].col, b+1);
	create_dvector(AA[1].row, wh);
	create_dvector(AA[1].col, wh+1);
	// Apply boundary condition
	updateFreenodes2bRHS(AA, b, wh, elementDOF+1,elementDOF+4);
	updateFreenodes2bMatrix11(AA, A, elementDOF+1,elementDOF+4);
	for(i=0;i<4;i++) free_csr_matrix(AA+i);
	printf("As.row=%d, As.col=%d, As.nnz=%d\n", A[0].row, A[0].col, A[0].nnz);
	printf("Bs.row=%d, Bs.col=%d, Bs.nnz=%d\n", A[2].row, A[2].col, A[2].nnz);

// write_dvector4Matlab(b, "output/b0.dat");//////////
// write_dvector4Matlab(b+1, "output/b1.dat");//////////
// write_IJ_dCSRmat4Matlab(&A[0], "output/A0.dat");
// write_IJ_dCSRmat4Matlab(&A[1], "output/A1.dat");
// write_IJ_dCSRmat4Matlab(&A[2], "output/A2.dat");
// write_IJ_dCSRmat4Matlab(&A[3], "output/A3.dat");
// write_ivector4Matlab(&elementDOF[1].freenodes, "output/ff0.dat");//////////
// write_ivector4Matlab(&elementDOF[1].nfreenodes, "output/nf0.dat");//////////
// write_ivector4Matlab(&elementDOF[3].freenodes, "output/ff1.dat");//////////
// write_ivector4Matlab(&elementDOF[3].nfreenodes, "output/nf1.dat");//////////

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}
	if (itsolver_type == 1){
		printf("\nASP Approximate Block Factorization preconditioned GMRES solver with auxiliary space method.\n\n");
		printf("Auxiliary space: P1 Lagrangian element.\n\n");
		AbfpAsP1symStokesMINI_GMRES(A, b, wh, &aspparam, print_level);
	//	Abfp2AsP1StokesNcP1P0_GMRES(A, b, wh, &aspparam, print_level);
	}
	else{
		printf("\nASP Block diagonal preconditioned MINRES solver with auxiliary space method\n");
		printf("Auxiliary space: P1 Lagrangian element.\n\n");
		DiagAsP1symStokesMINI_MINRES(A, b, wh, &aspparam, print_level);
		// DiagAsP1symStokesMINI_MINRES1(A, b, wh, &aspparam, print_level);
	//	DiagAMGStokesNcP1P0_MINRES(A, b, wh, &aspparam, print_level); // sometimes better
	}

	for(i=0;i<4;i++) free_csr_matrix(A+i);
	free_dvector(&b[0]); free_dvector(&b[1]);

	create_dvector(wh[0].row, sigmah);
	copy_dvector(wh, sigmah);
	free_dvector(&wh[0]); free_dvector(&wh[1]);
	
	
	/************ Eqn3: biharmonic equation disctetized by Morley element ************/
	assembleRHSSMINIhessMorley2d(b, elements, elementEdge, edges, elementDOF+2, elementDOF, sigmah);
	// initial solution
	create_dvector(b->row, uh);
	// Apply boundary condition
	updateFreenodesRHS(&AAm, b, uh, elementDOF, elementDOF, 1);
	free_csr_matrix(&AAm);
	printf("Am.row=%d, Am.col=%d, Am.nnz=%d\n", Am.row, Am.col, Am.nnz);
	printf("Solve...\n");
	printf("Solving Morley by AMG preconditioned CG solver\n");
	classicAMG_PCG(&Am, b, uh, &amgparam, print_level);
	free_dvector(&b[0]);
	free_csr_matrix(&Am);


	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */

}

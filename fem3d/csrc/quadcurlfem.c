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
 * \fn void quadcurlfem(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
 * \brief finite element methods for Qual-curl equation in three dimensions
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
void quadcurlfem(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i;
	int variationalform_type = Input->variationalform_type;
	int fem_num = Input->fem_num;
	if (variationalform_type == 1) // primal formulation
	{
		printf("Discretize quad-curl equation in primal formulation.\n");
		if(fem_num == 1)
		{
			printf("Discretized by the Huang element.\n");
			quadcurlHuang3d(elements, elementFace, faces, elementEdge, edges, nodes, Input);
		}
		
		if(fem_num == 2)
		{
			printf("Discretized by the Huang-Zhang element.\n");
			quadcurlHuangZhang3d(elements, elementFace, faces, elementEdge, edges, nodes, Input);
		}
	}
	else if (variationalform_type == 2) // decoupled mixed formulation
	{
		printf("Discretize quad-curl equation in Maxwell-Stokes-Maxwell form.\n");
		if(fem_num == 1)
		{
			printf("Discretized by the Huang element.\n");
			quadcurlHuang3d_msm(elements, elementFace, faces, elementEdge, edges, nodes, Input);
		}

		if(fem_num == 2)
		{
			printf("Discretized by the Huang-Zhang element to be developed.\n");
			exit(0);
		}
	}
	else if (variationalform_type == 3) // distributional mixed formulation
	{
		printf("Discretize quad-curl equation in hybridized distributional mixed formulation.\n");
		int dop1 = Input->dop1;
		int dop2 = Input->dop2;
		if(((dop1+1)!=dop2) && ((dop1+2)!=dop2)){
			printf("dop2 should equal dop1+1 or dop1+2: dop1 = %d, dop2 = %d.\n", dop1, dop2);
			exit(1);
		}
		else
			printf("dop1 = %d, dop2 = %d.\n", dop1, dop2);

		quadcurlDistribMixedFEM3d(elements, elementFace, faces, elementEdge, edges, nodes, Input);
	}
	else
	{
		printf("Please set variationalform_type = 1, 2 or 3!\n");
		exit(0);
	}
}

/**
* \fn void quadcurlHuang3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief Huang element method for quad-curl equation in three dimensions
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
void quadcurlHuang3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i,j;
	dCSRmat A;
	dvector b, uh;
	ELEMENT_DOF elementDOF[2];
		
	/** Step 0. Read input parameters */
	int print_level = Input->print_level;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;
//	int dop1 = Input->dop1;
	int dop1 = 2;
	int dop2 = 1;
		
	/** Step 1. generate degrees of freedom */
	getElementDOF_HuangGradcurl3d(elementDOF, elements, elementFace, faces, elementEdge, edges);
	getFreenodesInfoHuangGradcurl3d(faces, edges, nodes, 1, elementDOF);

	getElementDOF_Lagrange3d(elementDOF+1, elements, elementFace, elementEdge, faces->row, edges->row, nodes->row, dop2);
	getFreenodesInfoLagrange3d(faces, edges, nodes, elementDOF+1);
		
	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_quadcurlHuang3d(&A, &b, &uh, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);

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
	// write_dvector4Matlab(&b, outputfileb); 
	
	check_symm(&A);
	check_diagpos(&A);
	check_diagdom(&A);

//	exit(1);/////////


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

	/* AMG solver */
	if (itsolver_type == 1) {
		printf("AMG preconditioned CG solver\n");
		classicAMG_PCG(&A, &b, &uh, &amgparam, print_level);
	}

	
	/* PCG+diag */
	else if (itsolver_type == 2) {
		printf("Diagonal preconditioned CG solver\n");
		diag_PCG(&A, &b, &uh, itsolver_maxit, itsolver_tol, print_level);
	}

	/* PCG+AMG */
	else if (itsolver_type == 3) {
		printf("AMG iterative solver\n");
		classicAMG(&A, &b, &uh, &amgparam);
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

	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_dvector(&uh, outputfile); */

	/** Step 5. Compute the error between numerical solution and exact solution */
	double errors[4];

	geterrorsQuadcurlHuang3d(errors, &uh, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);

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

	
	free_csr_matrix(&A);
	free_dvector(&b);
	free_dvector(&uh);
	free_elementDOF(elementDOF);
	free_elementDOF(elementDOF+1);
}

/**
* \fn void quadcurlHuang3d_msm(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief Huang element method for quad-curl equation in Maxwell-Stokes-Maxwell form in three dimensions
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
void quadcurlHuang3d_msm(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i,j;
	dCSRmat AA[4],As[4], A;
	dvector b[2], wh, phih[2], uh;
	ELEMENT_DOF elementDOF[2], elementDOFs[2];
		
	/** Step 0. Read input parameters */
	int print_level = Input->print_level;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;
//	int dop1 = Input->dop1;
	int dop1 = 2;
	int dop2 = 1;

	ASP_param aspparam; /* parameters for ASP */

	aspparam.print_level = Input->print_level;
	aspparam.max_iter = Input->itsolver_maxit;
	aspparam.tol = Input->itsolver_tol;
	aspparam.restart = Input->restart;

	aspparam.precond_type = Input->precond_type;
	aspparam.precond_scale = Input->precond_scale;
	aspparam.smoother = Input->smoother;
	aspparam.schwarz_type = Input->schwarz_type;
	aspparam.smooth_iter = Input->smooth_iter;

	aspparam.levelNum = Input->glevelNum;
	aspparam.mg_max_iter = Input->MG_maxit;
	aspparam.mg_tol = Input->MG_tol;
	aspparam.mg_smoother = Input->MG_smoother;
	aspparam.mg_smooth_iter = Input->MG_smooth_iter;

	aspparam.elements = elements;
	aspparam.elementFace = elementFace;
	aspparam.faces = faces;
	aspparam.elementEdge = elementEdge;
	aspparam.edges = edges;
	aspparam.nodes = nodes;
//	aspparam.edgesTran = edgesTran;/////////
//	aspparam.nodeCEdge = nodeCEdge;/////////

	aspparam.lambda = Input->lambda;
	aspparam.mu = 0.5;// Input->mu;
	aspparam.t = Input->t;
	aspparam.nu = Input->nu;
	aspparam.problem_num = Input->problem_num;
	if (aspparam.problem_num == 1)
		aspparam.mu = 0.5;
	else
		aspparam.mu = (1 - aspparam.nu) / 2;

	aspparam.variationalform_type = Input->variationalform_type;
	aspparam.stress_fem_type = Input->stress_fem_type;

	aspparam.max_levels = Input->AMG_levels;
	aspparam.AMG_coarsening_type = Input->AMG_coarsening_type;
	aspparam.AMG_interpolation_type = Input->AMG_interpolation_type;
	aspparam.AMG_coarse_dof = Input->AMG_coarse_dof;
	aspparam.AMG_strong_threshold = Input->AMG_strong_threshold;
	aspparam.AMG_truncation_threshold = Input->AMG_truncation_threshold;
	aspparam.AMG_max_row_sum = Input->AMG_max_row_sum;

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

		
	/************ Eqn1: Maxwell equation disctetized by Huang element ************/
	printf("Eqn1: Maxwell equation disctetized by Huang element\n");
	/** Step 1. generate degrees of freedom */
	getElementDOF_HuangGradcurl3d(elementDOF, elements, elementFace, faces, elementEdge, edges);
	getFreenodesInfoHuangGradcurl3d(faces, edges, nodes, 1, elementDOF);
	getElementDOF_Lagrange3d(elementDOF+1, elements, elementFace, elementEdge, faces->row, edges->row, nodes->row, dop2);
	getFreenodesInfoLagrange3d(faces, edges, nodes, elementDOF+1);
	/** Step 2. assemble stiffmatrix and right hand side term */
	assembleBiCurlHuang3d(&AA[0], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	assembleHuangGradLagrange3d(&AA[2], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	getTransposeOfSparse(&AA[2], &AA[1]);
	assembleMassmatrixLagrange3d(&AA[3], elements, elementDOF+1);
	assemble_maxwellHuang3d0(AA, &A, &b[0], &wh, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, quadcurl3d_f);
	/** Step 3. Check matrix properties */
	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A.row, A.col, A.nnz);
//	print_dcsr_matrix(&A);
//	print_dvector(b[0].row, &b[0]);
	/** Step 4. Solve the system */
	printf("Solve...\n");
	if (print_level>0) {
		printf("Maximal iteration number = %d\n", Input->itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", Input->itsolver_tol);
	}
	printf("Solving Huang by AMG preconditioned CG solver\n");
	classicAMG_PCG(&A, &b[0], &wh, &amgparam, print_level);
	free_csr_matrix(&A);
	free_dvector(&b[0]);

	/************ Eqn2: Stokes equation disctetized by nonconforming P1-P0 element ************/
	printf("Eqn2: Stokes equation disctetized by nonconforming P1-P0 element\n");
	/** Step 1. generate degrees of freedom */
	getElementDOF_NoncfmP13d(&elementDOFs[0], elementFace, faces->row);
	getFreenodesInfoNoncfmP1Vector3d(faces, &elementDOFs[0]);
	getElementDOFdg(&elementDOFs[1], elements->row, 1,0);
	getFreenodesInfo(&elementDOFs[1]);
	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_quadcurl_StokesNcP1P03d(As, b, phih, elements, elementFace, faces, elementEdge, edges, nodes, elementDOFs, &wh, elementDOF);
	free_dvector(&wh);
	/** Step 3. Check matrix properties */
	printf("As.row=%d, As.col=%d, As.nnz=%d\n", As[0].row, As[0].col, As[0].nnz);
	printf("Bs.row=%d, Bs.col=%d, Bs.nnz=%d\n", As[2].row, As[2].col, As[2].nnz);
	/** Step 4. Solve the system */
	aspparam.elementDOF = elementDOFs;
	printf("Solve...\n");
	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}
	if (itsolver_type == 1)
	{
		printf("\nASP Approximate Block Factorization preconditioned GMRES solver with auxiliary space method.\n\n");
		printf("Auxiliary space: P1 Lagrangian element.\n\n");
		AbfpAsP1StokesNcP1P0_GMRES(As, b, phih, &aspparam, print_level);
	//	Abfp2AsP1StokesNcP1P0_GMRES(As, b, phih, &aspparam, print_level);
	}
	else
	{
		printf("\nASP Block diagonal preconditioned MINRES solver with auxiliary space method\n");
		printf("Auxiliary space: P1 Lagrangian element.\n\n");
		DiagAsP1StokesNcP1P0_MINRES(As, b, phih, &aspparam, print_level);
	//	DiagAMGStokesNcP1P0_MINRES(As, b, phih, &aspparam, print_level); // sometimes better
	}

	for (i = 0; i<4; i++)
		free_csr_matrix(&As[i]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);

	/************ Eqn3: Maxwell equation disctetized by Huang element ************/
	printf("Eqn3: Maxwell equation disctetized by Huang element\n");
	/** Step 1. assemble the right hand side term */
	assemble_quadcurl_maxwellHuang3d(AA, &A, &b[0], &uh, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, &phih[0], &elementDOFs[0]);
	free_dvector(&phih[0]);
	free_dvector(&phih[1]);
	free_elementDOF(&elementDOFs[0]);
	free_elementDOF(&elementDOFs[1]);
	/** Step 2. Solve the system */
	printf("Solve...\n");
	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}
	printf("Solving Morley by AMG preconditioned CG solver\n");
	classicAMG_PCG(&A, &b[0], &uh, &amgparam, print_level);
	free_dvector(&b[0]);
	
	for(i=0;i<4;i++) free_csr_matrix(AA+i);
	
	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_dvector(&uh, outputfile); */

	/** Step 5. Compute the error between numerical solution and exact solution */
	double errors[4];

	geterrorsQuadcurlHuang3d(errors, &uh, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);

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

	
	free_csr_matrix(&A);
	free_dvector(&uh);
	free_elementDOF(elementDOF);
	free_elementDOF(elementDOF+1);
}

/**
* \fn void quadcurlHuangZhang3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
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
void quadcurlHuangZhang3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i,j;
	dCSRmat A;
	dvector b, uh;
	ELEMENT_DOF elementDOF[2];
		
	/** Step 0. Read input parameters */
	int print_level = Input->print_level;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;
//	int dop1 = Input->dop1;
	int dop1 = 2;
	int dop2 = 2;
		
	/** Step 1. generate degrees of freedom */
	getElementDOF_HuangZhang3d(elementDOF, elements, elementFace, faces, elementEdge, edges);
	getFreenodesInfoHuangZhang3d(faces, edges, nodes, 1, elementDOF);

	getElementDOF_Lagrange3d(elementDOF+1, elements, elementFace, elementEdge, faces->row, edges->row, nodes->row, dop2);
	getFreenodesInfoLagrange3d(faces, edges, nodes, elementDOF+1);

	/***************************Generate coefficient of basis functions**************************/
	ddenmat3 basisCoeffs;
	create_dden_matrix3(elements->row, elementDOF->col, elementDOF->col, &basisCoeffs);
	getHuangZhangBasisCoeffs(&basisCoeffs, elements, elementFace, faces, elementEdge, edges);
	/********************************************************************************************/

	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_quadcurlHuangZhang3d(&A, &b, &uh, &basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);

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
	// write_dvector4Matlab(&b, outputfileb); 
	
	check_symm(&A);
	check_diagpos(&A);
	check_diagdom(&A);

//	exit(1);/////////


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

		read_dvector4Matlab(&uh, outputfileuh);
	}

	/* AMG solver */
	if (itsolver_type == 1) {
		printf("AMG preconditioned CG solver\n");
		classicAMG_PCG(&A, &b, &uh, &amgparam, print_level);
	}

	
	/* PCG+diag */
	else if (itsolver_type == 2) {
		printf("Diagonal preconditioned CG solver\n");
		diag_PCG(&A, &b, &uh, itsolver_maxit, itsolver_tol, print_level);
	}

	/* PCG+AMG */
	else if (itsolver_type == 3) {
		printf("AMG iterative solver\n");
		classicAMG(&A, &b, &uh, &amgparam);
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

	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_dvector(&uh, outputfile); */

	/** Step 5. Compute the error between numerical solution and exact solution */
	double errors[4];

	geterrorsQuadcurlHuangZhang3d(errors, &uh, &basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);

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
* \fn void quadcurlDistribMixedFEM3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief Distributional mixed finite element method for quad-curl equation in three dimensions
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
void quadcurlDistribMixedFEM3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i,j;
	dCSRmat A, Asigma;
	dvector b, uh, sigmah, lambdah, ulambdah;
	ELEMENT_DOF elementDOF[5]; // sigma, u, lambda, phi
		
	/** Step 0. Read input parameters */
	int print_level = Input->print_level;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;
	int dop1 = Input->dop1;
	int dop2 = Input->dop2;
	int curlfemtype;
	if(dop1+1 == dop2)
		curlfemtype = 1;
	else // dop1+2 == dop2
		curlfemtype = 2;
		
	/** Step 1. generate degrees of freedom */
	// elementDOF: 0 sigma, 1 u, 2 lambda, 3 phi
	getElementDOFdg(elementDOF, elements->row, (dop1+1)*(dop1+2)*(dop1+3)/6, dop1);
	// getFreenodesInfo(elementDOF);
	if(curlfemtype == 1){
		getElementDOF_Nedelec1st3d(elementDOF+1, elements, elementFace, faces, elementEdge, edges, dop1+1);
		getFreenodesInfoNedelec1st3d(faces, edges, nodes, elementDOF+1);
	}
	else{ // curlfemtype == 2
		getElementDOF_Nedelec2nd3d(elementDOF+1, elements, elementFace, faces, elementEdge, edges, dop1+1);
		getFreenodesInfoNedelec2nd3d(faces, edges, nodes, elementDOF+1);
	}

	getElementDOFdg(elementDOF+2, faces->row, (dop1+1)*(dop1+2)/2, dop1);
	getFreenodesInfoFacesVector(faces, elementDOF+2, 2);
	getElementDOF_Lagrange3d(elementDOF+3, elements, elementFace, elementEdge, faces->row, edges->row, nodes->row, dop2);
	getFreenodesInfoLagrange3d(faces, edges, nodes, elementDOF+3);

	combineFreenodesInfo(&elementDOF[1].nfFlag, &elementDOF[2].nfFlag, elementDOF+4);

	printf(" sigmah: DG, dop = %d\n", elementDOF[0].dop);
	if(curlfemtype == 1)
		printf("     uh: the first kind Nedelec element, dop = %d\n", elementDOF[1].dop);
	else
		printf("     uh: the second kind Nedelec element, dop = %d\n", elementDOF[1].dop);
	printf("lambdah: DG, dop = %d\n", elementDOF[2].dop);
	printf("   phih: Lagrange element, dop = %d\n", elementDOF[3].dop);
	
	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_quadcurlDistribMixedFEM3d(&A, &Asigma, &b, &ulambdah, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, curlfemtype);

	free_ivector(&elementDOF[4].nfFlag);
	free_ivector(&elementDOF[4].freenodes);
	free_ivector(&elementDOF[4].nfreenodes);
	free_ivector(&elementDOF[4].index);

	
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

	/* AMG solver */
	if (itsolver_type == 1) {
		printf("AMG preconditioned CG solver\n");
		classicAMG_PCG(&A, &b, &ulambdah, &amgparam, print_level);
	}

	
	/* PCG+diag */
	else if (itsolver_type == 2) {
		printf("Diagonal preconditioned CG solver\n");
		diag_PCG(&A, &b, &ulambdah, itsolver_maxit, itsolver_tol, print_level);
	}

	/* PCG+AMG */
	else if (itsolver_type == 3) {
		printf("AMG iterative solver\n");
		classicAMG(&A, &b, &ulambdah, &amgparam);
	}

	/* CG */
	else if (itsolver_type == 4) {
		printf("Classical CG solver\n");
		standard_CG(&A, &b, &ulambdah, itsolver_maxit, itsolver_tol, print_level);
	}

	/* GMRES+AMG */
	else if (itsolver_type == 5) {
		printf("AMG preconditioned GMRES solver\n");
		classicAMG_GMRES(&A, &b, &ulambdah, &amgparam, print_level);
	}

	create_dvector(Asigma.row, &sigmah);
	sparse_mv0(-1.0, &Asigma, ulambdah.val, sigmah.val);
	free_csr_matrix(&Asigma);
	create_dvector(elementDOF[1].dof, &uh);
	create_dvector(elementDOF[2].dof*2, &lambdah);
	copy_array(uh.row, ulambdah.val, uh.val);
	copy_array(lambdah.row, ulambdah.val+uh.row, lambdah.val);
	free_dvector(&ulambdah);


	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_dvector(&uh, outputfile); */

	/** Step 5. Postprocessing */
	dvector uhstar;
	ELEMENT_DOF elementDOFuhstar;
	getElementDOFdg(&elementDOFuhstar, elements->row, (dop1 + 2)*(dop1 + 4)*(dop1 + 5) / 2, dop1 + 2);
	create_dvector(elementDOFuhstar.dof, &uhstar);
	postprocessQuadcurlDistribMixedFEM3d(&uhstar, &sigmah, &uh, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, &elementDOFuhstar, curlfemtype);	
	

	/** Step 6. Compute the error between numerical solution and exact solution */
	double errors[7];

	geterrorsQuadcurlDistribMixedFEM3d(errors, &sigmah, &uh, &lambdah, &uhstar, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, &elementDOFuhstar, curlfemtype);

	printf("\nA Priori Errors:\n");
	printf("L2 norm of sigma-sigma_h = %e\n", errors[0]);
	printf("L2 norm of u-u_h = %e\n", errors[1]);
	printf("L2 norm of curl(u-u_h) = %e\n", errors[2]);
	printf("L2 norm of gradcurl(u-u_h) = %e\n", errors[3]);
	printf("L2 norm of u-u_h* = %e\n", errors[4]);
	printf("L2 norm of curl(u-u_h*) = %e\n", errors[5]);
	printf("L2 norm of gradcurl(u-u_h*) = %e\n", errors[6]);
	
	/*********************************************************************************************/
	FILE *outputFile;
	outputFile = fopen("output/error.dat", "w");
	fprintf(outputFile, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", errors[0], errors[1], errors[2], errors[3], errors[4], errors[5], errors[6]);
	fclose(outputFile);
	/********************************************************************************************/

	
	free_csr_matrix(&A);
	free_dvector(&b);
	free_dvector(&sigmah);
	free_dvector(&uh);
	free_dvector(&lambdah);
	free_dvector(&uhstar);
	for(i=0;i<4;i++) free_elementDOF(elementDOF+i);
	free_elementDOF(&elementDOFuhstar);
}

/**
 * \fn void assemble_quadcurlHuang3d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
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
 * \param lambda Lame constant
 * \param mu Lame constant or Poisson ratio of plate
 * \return void
 */
void assemble_quadcurlHuang3d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Dx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	int i;
	dCSRmat AA[4], AB[4];
	dvector D, bb[2], wh[2];
	assembleBiGradcurlHuang3d(&AA[0], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	assembleHuangGradLagrange3d(&AA[2], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	getTransposeOfSparse(&AA[2], &AA[1]);
	assembleMassmatrixLagrange3d(&AA[3], elements, elementDOF+1);
	assembleRHSHuang3d(&bb[0], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, quadcurl3d_f);
	create_dvector(AA[1].col, bb+1);
	create_dvector(AA[1].row, wh);
	create_dvector(AA[1].col, wh+1);

	// Apply boundary condition
	updateFreenodes2bRHS(AA, bb, wh, elementDOF,elementDOF+1);
	updateFreenodes2bMatrix11(AA, AB, elementDOF,elementDOF+1);
	for(i=0;i<4;i++) free_csr_matrix(AA+i);
	getdiag(AB[3].row, &AB[3], &D);
	free_csr_matrix(&AB[3]);

	// Schur Complement : A = AB[0] + AB[1] Dinv AB[2]
	create_dvector(bb[0].row, b);
	copy_dvector(&bb[0], b);
	dotdiv_dvector(&D, &bb[1]);
	sparse_mv(1.0, &AB[1], bb[1].val, b->val);
	free_dvector(&bb[0]); free_dvector(&bb[1]);

	create_dvector(wh[0].row, uh);
	copy_dvector(&wh[0], uh);
	free_dvector(&wh[0]); free_dvector(&wh[1]);

	dDiagVectorInvMultiplydCSR(&D, &AB[2], &AA[2]);
	free_dvector(&D); free_csr_matrix(&AB[2]);
	sparseMultiplication(&AB[1], &AA[2], &AA[0]);
	free_csr_matrix(&AB[1]); free_csr_matrix(&AA[2]);
	sparseAddition(&AB[0], &AA[0], A);
	free_csr_matrix(&AB[0]); free_csr_matrix(&AA[0]);

/**********************************
	double EPS = 1e-300;
	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A.row, A.col, A.nnz);
	compress_dcsr(&A, ptr_A, EPS);
	free_csr_matrix(&A);
	printf("A.row=%d, A.col=%d, A.nnz=%d\n", ptr_A->row, ptr_A->col, ptr_A->nnz);
	*******************************************/
}

/**
 * \fn void assemble_quadcurl_StokesNcP1P03d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, dvector *wh, ELEMENT_DOF *elementDOFwh)
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
 * \param lambda Lame constant
 * \param mu Lame constant or Poisson ratio of plate
 * \return void
 */
void assemble_quadcurl_StokesNcP1P03d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, dvector *wh, ELEMENT_DOF *elementDOFwh)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Cx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	int i;
	dCSRmat AA[4];

	assembleBiGradVectorNcP13d(&AA[0], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	assembleDivNcP1P03d(&AA[2], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	getTransposeOfSparse(&AA[2], &AA[1]);
	create_csr_matrix(elementDOF[1].dof, elementDOF[1].dof, 0, &AA[3]);
	assembleRHSCurlHuangNcP13d4Stokes(b, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, wh, elementDOFwh);
    create_dvector(AA[1].row, uh);
	create_dvector(AA[1].col, uh+1);
	
	// Apply boundary condition
	updateFreenodes2bRHS(AA, b, uh, elementDOF, elementDOF+1);
	updateFreenodes2bMatrix11(AA, A, elementDOF, elementDOF+1);
	for(i=0;i<4;i++) free_csr_matrix(AA+i);	
}

/**
 * \fn void assemble_quadcurl_maxwellHuang3d(dCSRmat *AA, dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, dvector *uhs, ELEMENT_DOF *elementDOFs)
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
 * \param lambda Lame constant
 * \param mu Lame constant or Poisson ratio of plate
 * \return void
 */
void assemble_quadcurl_maxwellHuang3d(dCSRmat *AA, dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, dvector *uhs, ELEMENT_DOF *elementDOFs)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Dx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	int i;
	dCSRmat AB[4], tempA[2];
	dvector D, bb[2], wh[2];
	assembleRHSCurlHuangNcP13d4Maxwell(bb, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, uhs, elementDOFs);
	create_dvector(AA[1].col, bb+1);
	create_dvector(AA[1].row, wh);
	create_dvector(AA[1].col, wh+1);

	// Apply boundary condition
	updateFreenodes2bRHS(AA, bb, wh, elementDOF,elementDOF+1);
	updateFreenodes2bMatrix11(AA, AB, elementDOF,elementDOF+1);
	getdiag(AB[3].row, &AB[3], &D);
	free_csr_matrix(&AB[3]);

	// Schur Complement : A = AB[0] + AB[1] Dinv AB[2]
	create_dvector(bb[0].row, b);
	copy_dvector(&bb[0], b);
	dotdiv_dvector(&D, &bb[1]);
	sparse_mv(1.0, &AB[1], bb[1].val, b->val);
	free_dvector(&bb[0]); free_dvector(&bb[1]);

	create_dvector(wh[0].row, uh);
	copy_dvector(&wh[0], uh);
	free_dvector(&wh[0]); free_dvector(&wh[1]);

	dDiagVectorInvMultiplydCSR(&D, &AB[2], &tempA[0]);
	free_dvector(&D); free_csr_matrix(&AB[2]);
	sparseMultiplication(&AB[1], &tempA[0], &tempA[1]);
	free_csr_matrix(&AB[1]); free_csr_matrix(&tempA[0]);
	sparseAddition(&AB[0], &tempA[1], A);
	free_csr_matrix(&AB[0]); free_csr_matrix(&tempA[1]);
		
/**********************************
	double EPS = 1e-300;
	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A.row, A.col, A.nnz);
	compress_dcsr(&A, ptr_A, EPS);
	free_csr_matrix(&A);
	printf("A.row=%d, A.col=%d, A.nnz=%d\n", ptr_A->row, ptr_A->col, ptr_A->nnz);
	*******************************************/
}

/**
 * \fn void assemble_quadcurlHuang3d(dCSRmat *A, dvector *b, dvector *uh, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
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
 * \param lambda Lame constant
 * \param mu Lame constant or Poisson ratio of plate
 * \return void
 */
void assemble_quadcurlHuangZhang3d(dCSRmat *A, dvector *b, dvector *uh, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Dx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	int i;
	dCSRmat AA[4], AB[4];
	dvector D, bb[2], wh[2];
	assembleBiGradcurlHuangZhang3d(&AA[0], basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	assembleHuangZhangGradLagrange3d(&AA[2], basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	getTransposeOfSparse(&AA[2], &AA[1]);
	assembleMassmatrixLagrange3d(&AA[3], elements, elementDOF+1);
	assembleRHSHuangZhang3d(&bb[0], basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, quadcurl3d_f);
	create_dvector(AA[1].col, bb+1);
	create_dvector(AA[1].row, wh);
	create_dvector(AA[1].col, wh+1);

	// Apply boundary condition
	updateFreenodes2bRHS(AA, bb, wh, elementDOF,elementDOF+1);
	updateFreenodes2bMatrix11(AA, AB, elementDOF,elementDOF+1);
	for(i=0;i<4;i++) free_csr_matrix(AA+i);
	getdiag(AB[3].row, &AB[3], &D);
	free_csr_matrix(&AB[3]);

	// Schur Complement : A = AB[0] + AB[1] Dinv AB[2]
	create_dvector(bb[0].row, b);
	copy_dvector(&bb[0], b);
	dotdiv_dvector(&D, &bb[1]);
	sparse_mv(1.0, &AB[1], bb[1].val, b->val);
	free_dvector(&bb[0]); free_dvector(&bb[1]);

	create_dvector(wh[0].row, uh);
	copy_dvector(&wh[0], uh);
	free_dvector(&wh[0]); free_dvector(&wh[1]);

	dDiagVectorInvMultiplydCSR(&D, &AB[2], &AA[2]);
	free_dvector(&D); free_csr_matrix(&AB[2]);
	sparseMultiplication(&AB[1], &AA[2], &AA[0]);
	free_csr_matrix(&AB[1]); free_csr_matrix(&AA[2]);
	sparseAddition(&AB[0], &AA[0], A);
	free_csr_matrix(&AB[0]); free_csr_matrix(&AA[0]);

/**********************************
	double EPS = 1e-300;
	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A.row, A.col, A.nnz);
	compress_dcsr(&A, ptr_A, EPS);
	free_csr_matrix(&A);
	printf("A.row=%d, A.col=%d, A.nnz=%d\n", ptr_A->row, ptr_A->col, ptr_A->nnz);
	*******************************************/
}

/**
 * \fn void assemble_quadcurlDistribMixedFEM3d(dCSRmat *A, dCSRmat *Asigma, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, int curlfemtype)
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
 * \param curlfemtype type of Nedelec element
 * \return void
 */
void assemble_quadcurlDistribMixedFEM3d(dCSRmat *A, dCSRmat *Asigma, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, int curlfemtype)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Dx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	int i;
	dBDmat Ah, Ahinv; // hybridized A
	dCSRmat BCh[2]; // hybridized B, C
	dCSRmat AA[4], AB[4];
	dvector D, bb[2], wh[2];

	assembleMassmatrixTracelessDGPk3d(&Ah, elements, elementDOF);
	inverse_dBDmat(&Ah, &Ahinv);

/////////////////////////////////////////////
	// int j,k,k1,k2;
	// printf("Ah:\n");
	// for (k = 0; k<elements->row; k++){
	// 	for(i=0;i<8;i++){
	// 		printf("k = %d, i = %d \n", k,i);
	// 		for (k1 = 0; k1<elementDOF->col; k1++){
	// 			for (k2 = 0; k2<elementDOF->col; k2++)
	// 				printf("%e ",Ah.blk[k].val[k1][k2]);
	// 			printf("\n");
	// 		}
	// 	}
	// }
	// printf("Ahinv:\n");
	// for (k = 0; k<elements->row; k++){
	// 	for(i=0;i<8;i++){
	// 		printf("k = %d, i = %d \n", k,i);
	// 		for (k1 = 0; k1<elementDOF->col; k1++){
	// 			for (k2 = 0; k2<elementDOF->col; k2++)
	// 				printf("%e ",Ahinv.blk[k].val[k1][k2]);
	// 			printf("\n");
	// 		}
	// 	}
	// }
/////////////////////////////////////////////


	free_dbd_matrix(&Ah);	

	assembleDistribCurldivTracelessDGNedelecHybrid3d(&BCh[0], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, curlfemtype);
	getTransposeOfSparse(&BCh[0], &BCh[1]);
	

	dBDMultiplydCSR(&Ahinv, &BCh[1], Asigma);
	free_dbd_matrix(&Ahinv);
	free_csr_matrix(&BCh[1]);
	sparseMultiplication(&BCh[0], Asigma, &AA[0]);
	free_csr_matrix(&BCh[0]);

	if(curlfemtype == 1)
		assembleNedelec1stGradLagrange3d(&AA[2], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF+1, elementDOF+3);
	else // curlfemtype == 2
		assembleNedelec2ndGradLagrange3d(&AA[2], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF+1, elementDOF+3);
	AA[2].col += elementDOF[2].dof*2;
	getTransposeOfSparse(&AA[2], &AA[1]);
	assembleMassmatrixLagrange3d(&AA[3], elements, elementDOF+3);

	assembleRHSNedelec3d(&bb[0], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF+1, elementDOF[2].dof*2, curlfemtype, quadcurl3d_f);
	create_dvector(AA[1].col, bb+1);
	create_dvector(AA[1].row, wh);
	create_dvector(AA[1].col, wh+1);

	// Apply boundary condition
	updateFreenodes2bRHS(AA, bb, wh, elementDOF+4, elementDOF+3);
	updateFreenodes2bMatrix11(AA, AB, elementDOF+4, elementDOF+3);

	for(i=0;i<4;i++) free_csr_matrix(AA+i);
	getdiag(AB[3].row, &AB[3], &D);
	free_csr_matrix(&AB[3]);

	/* output A, b to a diskfile */
	// char *outputfileA1="output/A1.dat";
	// char *outputfileA2="output/A2.dat";
	// char *outputfileD="output/D.dat";
	// char *outputfileb="output/b.dat";
	// write_IJ_dCSRmat4Matlab(&AB[1], outputfileA1); 
	// write_IJ_dCSRmat4Matlab(&AB[2], outputfileA2); 
	// write_dvector4Matlab(&D, outputfileD); 
	// write_dvector4Matlab(&bb[0], outputfileb); 

	// Schur Complement : A = AB[0] + AB[1] Dinv AB[2]
	create_dvector(bb[0].row, b);
	copy_dvector(&bb[0], b);
	dotdiv_dvector(&D, &bb[1]);
	sparse_mv(1.0, &AB[1], bb[1].val, b->val);
	free_dvector(&bb[0]); free_dvector(&bb[1]);

	create_dvector(wh[0].row, uh);
	copy_dvector(&wh[0], uh);
	free_dvector(&wh[0]); free_dvector(&wh[1]);


	dDiagVectorInvMultiplydCSR(&D, &AB[2], &AA[2]);
	free_dvector(&D); free_csr_matrix(&AB[2]);
	sparseMultiplication(&AB[1], &AA[2], &AA[0]);
	free_csr_matrix(&AB[1]); free_csr_matrix(&AA[2]);
	sparseAddition(&AB[0], &AA[0], A);
	free_csr_matrix(&AB[0]); free_csr_matrix(&AA[0]);

/**********************************
	double EPS = 1e-300;
	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A.row, A.col, A.nnz);
	compress_dcsr(&A, ptr_A, EPS);
	free_csr_matrix(&A);
	printf("A.row=%d, A.col=%d, A.nnz=%d\n", ptr_A->row, ptr_A->col, ptr_A->nnz);
	*******************************************/
}

/**
 * \fn void assembleRHSCurlHuangNcP13d4Stokes(dvector *b, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, dvector *wh, ELEMENT_DOF *elementDOFwh)
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *BT pointer to stiffness matrix
 * \param *C pointer to stiffness matrix
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
 * \param lambda Lame constant
 * \param mu Lame constant or Poisson ratio of plate
 * \return void
 */
void assembleRHSCurlHuangNcP13d4Stokes(dvector *b, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, dvector *wh, ELEMENT_DOF *elementDOFwh)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Cx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[3], phi2[3];
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], **grd_lambda, **nv, *nvf[4], *etv[6], **vertices;
	double vol, s[4], h[6];
	int rowstart[3], row31[3];
	int count;
	short *forien, *eorien;
	int **fpermi;
	double *lambdaConst;

	int num_qp;
	double lambdas[100][4], lambdas2[100][3], weight[100];
	
	dvector lwh;
	create_dvector(elementDOFwh[0].col, &lwh);
	dvector lb;
	create_dvector(3*elementDOF[0].col, &lb);
	
	 /************************************************** right hand side b *****************************************************************/
	create_dvector(elementDOF[0].dof * 3, b);
	create_dvector(elementDOF[1].dof, b+1);
	num_qp = getNumQuadPoints_ShunnWilliams(3, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vol = elements->vol[k];
		grd_lambda = elements->gradLambda[k];
		vertices = elements->vertices[k];
		nv = elements->nvector[k];
		forien = elements->forien[k];
		eorien = elements->eorien[k];
		fpermi = elements->fpermi[k];
		lambdaConst = elements->lambdaConst[k];
		for (i = 0; i < elementFace->col; i++)
		{
			face = elementFace->val[k][i];
			nvf[i] = faces->nvector[face];
			s[i] = faces->area[face];
		}
		// for (i = 0; i<6; i++)
		// {
		// 	edge = elementEdge->val[k][i];
		// 	etv[i] = edges->tvector[edge];
		// 	h[i] = edges->length[edge];
		// }
		// end set parameters

		for (i = 0; i<elementDOFwh[0].col; i++)
		{
			j = elementDOFwh[0].val[k][i];
			lwh.val[i] = wh->val[j];
		}

		init_dvector(&lb, 0.0);
		for (i = 0; i<elementDOF[0].col; i++)
		{
			for (i1 = 0; i1<num_qp; i1++)
			{
				ncp13d_basis(lambdas[i1], i, &phi);
				init_array(3, phi2, 0.0);
				for (j = 0; j<elementDOFwh[0].col; j++)
				{
					huangQuadcurl3d_basisCurl(lambdas[i1], grd_lambda, nvf, eorien, fpermi, j, phi1);
					axpy_array(3, lwh.val[j], phi1, phi2);
				}
				lb.val[i] += vol*weight[i1] * phi2[0] * phi;
				lb.val[i + elementDOF[0].col] += vol*weight[i1] * phi2[1] * phi;
				lb.val[i + elementDOF[0].col*2] += vol*weight[i1] * phi2[2] * phi;
			} // i1
		} // k1 

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			b->val[i] += lb.val[k1];
			b->val[i + elementDOF[0].dof] += lb.val[k1+elementDOF[0].col];
			b->val[i + elementDOF[0].dof*2] += lb.val[k1+elementDOF[0].col*2];
		} // k1
	} // k
	free_dvector(&lwh);
	free_dvector(&lb);
}

/**
 * \fn void assembleRHSCurlHuangNcP13d4Maxwell(dvector *b, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, dvector *uhs, ELEMENT_DOF *elementDOFs)
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *BT pointer to stiffness matrix
 * \param *C pointer to stiffness matrix
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
 * \param lambda Lame constant
 * \param mu Lame constant or Poisson ratio of plate
 * \return void
 */
void assembleRHSCurlHuangNcP13d4Maxwell(dvector *b, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, dvector *uhs, ELEMENT_DOF *elementDOFs)
{
	int i,j,k,l;

	int nvertices=nodes->row;
	int nedges=edges->row;
	int nfaces=faces->row;
	int element, face, edge, node;
	
	double phi, phi1[3], phi2[3];
	int k1,k2,i1,j1,l1,l2,i2,ej;
	double x[3], **grd_lambda, **nv, *nvf[4], *etv[6], **vertices;
	double vol, s[4], h[6];
	int rowstart[3], row31[3];
	int count;
	short *forien, *eorien;
	int **fpermi;
	double *lambdaConst;

	int num_qp;
	double lambdas[100][4], lambdas2[100][3], weight[100];
		
	dvector luhs;
	create_dvector(elementDOFs[0].col*3, &luhs);
	
	dvector lb;
	create_dvector(elementDOF[0].col, &lb);
	/************************************************** right hand side b *****************************************************************/
	// dvector b;
	create_dvector(elementDOF[0].dof, b);
	num_qp = getNumQuadPoints_ShunnWilliams(3, 3); // the number of numerical intergation points
	init_ShunnWilliams3d(num_qp, lambdas, weight); // Shunn-Williams intergation initial
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		vol = elements->vol[k];
		grd_lambda = elements->gradLambda[k];
		vertices = elements->vertices[k];
		nv = elements->nvector[k];
		forien = elements->forien[k];
		eorien = elements->eorien[k];
		fpermi = elements->fpermi[k];
		lambdaConst = elements->lambdaConst[k];
		for (i = 0; i < elementFace->col; i++){
			face = elementFace->val[k][i];
			nvf[i] = faces->nvector[face];
			s[i] = faces->area[face];
		}
		// end set parameters

		for (i = 0; i<elementDOFs[0].col; i++){
			j = elementDOFs[0].val[k][i];
			luhs.val[i] = uhs[0].val[j];
			luhs.val[i+elementDOFs[0].col] = uhs[0].val[j+elementDOFs[0].dof];
			luhs.val[i+elementDOFs[0].col*2] = uhs[0].val[j+elementDOFs[0].dof*2];
		}

		init_dvector(&lb, 0.0);
		for (i = 0; i<elementDOF[0].col; i++)
		{
			for (i1 = 0; i1<num_qp; i1++)
			{
				huangQuadcurl3d_basisCurl(lambdas[i1], grd_lambda, nvf, eorien, fpermi, i, phi1);
				init_array(3, phi2, 0.0);
				for (j = 0; j<elementDOFs[0].col; j++)
				{
					ncp13d_basis(lambdas[i1], j, &phi);
					phi2[0] += phi * luhs.val[j];
					phi2[1] += phi * luhs.val[j + elementDOFs[0].col];
					phi2[2] += phi * luhs.val[j + elementDOFs[0].col*2];
				}
				lb.val[i] += vol*weight[i1] * dot_array(3, phi1, phi2);
			} // i1
		} // k1
	
		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			b->val[i] += lb.val[k1];
		} // k1 
	} // k
	free_dvector(&luhs);
	free_dvector(&lb);
}

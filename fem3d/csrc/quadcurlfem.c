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
	else
	{
		printf("Please set variationalform_type = 1 or 2!\n");
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
	dvector b, uh, _uh;
	ELEMENT_DOF elementDOF[2];
	iCSRmat elementdofTran[2];
		
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
	getTransposeOfelementDoF(elementDOF, elementdofTran, 0);

	getElementDOF_Lagrange3d(elementDOF+1, elements, elementFace, faces, elementEdge, edges, nodes->row, dop2);
	getFreenodesInfoLagrange3d(faces, edges, nodes, elementDOF+1);
	getTransposeOfelementDoF(elementDOF+1, elementdofTran+1, 0);
		
	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_quadcurlHuang3d(&A, &b, &uh, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, elementdofTran);
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
	dCSRmat As[4], A;
	dvector b[2], wh, phih[2], uh, _uh[2];
	ELEMENT_DOF elementDOF[2], elementDOFs[2];
	iCSRmat elementdofTran[2], elementdofTrans;
		
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
	/** Step 1. generate degrees of freedom */
	getElementDOF_HuangGradcurl3d(elementDOF, elements, elementFace, faces, elementEdge, edges);
	getFreenodesInfoHuangGradcurl3d(faces, edges, nodes, 1, elementDOF);
	getTransposeOfelementDoF(elementDOF, elementdofTran, 0);
	getElementDOF_Lagrange3d(elementDOF+1, elements, elementFace, faces, elementEdge, edges, nodes->row, dop2);
	getFreenodesInfoLagrange3d(faces, edges, nodes, elementDOF+1);
	getTransposeOfelementDoF(elementDOF+1, elementdofTran+1, 0);
	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_maxwellHuang3d(&A, &b[0], &wh, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, elementdofTran);
	for(i=0;i<2;i++)
		free_icsr_matrix(elementdofTran+i);
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
	create_dvector(b[0].row, &_uh[0]);
	init_dvector(&_uh[0], 0.0);
	printf("Solving Huang by AMG preconditioned CG solver\n");
	classicAMG_PCG(&A, &b[0], &_uh[0], &amgparam, print_level);
	for (i = 0; i < _uh[0].row; i++)
		wh.val[elementDOF[0].freenodes.val[i]] = _uh[0].val[i];

	free_dvector(&_uh[0]);
	free_dvector(&b[0]);
//	exit(1);////
	/************ Eqn2: Stokes equation disctetized by nonconforming P1-P0 element ************/
	/** Step 1. generate degrees of freedom */
	getElementDOF_NoncfmP13d(&elementDOFs[0], elementFace, faces->row);
	getFreenodesInfoNoncfmP1Vector3d(faces, &elementDOFs[0]);
	getElementDOF3d(&elementDOFs[1], elements->row, 0);
	getTransposeOfelementDoF(&elementDOFs[0], elementdofTran, 0);
	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_quadcurl_StokesNcP1P03d(As, b, phih, elements, elementFace, faces, elementEdge, edges, nodes, elementDOFs, elementdofTran, &wh, elementDOF);
	free_dvector(&wh);
	free_icsr_matrix(elementdofTran);
	/** Step 3. Check matrix properties */
	printf("As.row=%d, As.col=%d, As.nnz=%d\n", As[0].row, As[0].col, As[0].nnz);
	printf("Bs.row=%d, Bs.col=%d, Bs.nnz=%d\n", As[2].row, As[2].col, As[2].nnz);
	/** Step 4. Solve the system */
	aspparam.elementDOF = elementDOFs;
//	aspparam.elementdofTran = elementdofTran;
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	create_dvector(b[1].row, &_uh[1]);
	init_dvector(&_uh[0], 1);///////////////////////////////////

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}
	if (itsolver_type == 1)
	{
		printf("\nASP Approximate Block Factorization preconditioned GMRES solver with auxiliary space method.\n\n");
		printf("Auxiliary space: P1 Lagrangian element.\n\n");
		AbfpAsP1StokesNcP1P0_GMRES(As, b, _uh, &aspparam, print_level);
	//	Abfp2AsP1StokesNcP1P0_GMRES(As, b, _uh, &aspparam, print_level);
	}
	else
	{
		printf("\nASP Block diagonal preconditioned MINRES solver with auxiliary space method\n");
		printf("Auxiliary space: P1 Lagrangian element.\n\n");
		DiagAsP1StokesNcP1P0_MINRES(As, b, _uh, &aspparam, print_level);
	//	DiagAMGStokesNcP1P0_MINRES(As, b, _uh, &aspparam, print_level); // sometimes better
	}
	for (i = 0; i < _uh[0].row; i++)
		phih[0].val[elementDOFs[0].freenodes.val[i]] = _uh[0].val[i];
	/*	for (i = 0; i < _uh[1].row; i++)
	{
	if (problem_num == 1)
	uhs[1].val[i] = _uh[1].val[i];
	else
	uhs[1].val[i] = _uh[1].val[i] * (1 - nu);
	}*/
//	for (i = 0; i < _uh[1].row; i++)
//		uhs[1].val[i] = _uh[1].val[i];
	free_dvector(&_uh[0]);
	free_dvector(&_uh[1]);
	for (i = 0; i<3; i++)
		free_csr_matrix(&As[i]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);

	/************ Eqn3: Maxwell equation disctetized by Huang element ************/
	/** Step 1. assemble the right hand side term */
	create_dvector(elementDOF[0].dof, &uh);
	assembleRHSCurlHuangNcP13d4Maxwell(b, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, &phih[0], &elementDOFs[0]);
	free_dvector(&phih[0]);
	free_dvector(&phih[1]);
	free_elementDOF(&elementDOFs[0]);
	free_elementDOF(&elementDOFs[1]);
	/** Step 2. Solve the system */
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}
	printf("Solving Morley by AMG preconditioned CG solver\n");
	classicAMG_PCG(&A, &b[0], &_uh[0], &amgparam, print_level);
	for (i = 0; i < _uh[0].row; i++)
		uh.val[elementDOF[0].freenodes.val[i]] = _uh[0].val[i];
	free_dvector(&_uh[0]);
	free_dvector(&b[0]);
	


	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */

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
	dvector b, uh, _uh;
	ELEMENT_DOF elementDOF[2];
	iCSRmat elementdofTran[2];
		
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
	assemble_quadcurlHuangZhang3d(&A, &b, &uh, &basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, elementdofTran);
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
 * \fn void assemble_quadcurlHuang3d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran)
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
void assemble_quadcurlHuang3d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Dx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	dCSRmat AA, BB, CC, tempA;
	dCSRmat B, BT, C;
	dvector bb, D;

	assembleBiGradcurlHuang3d(&AA, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, elementdofTran);
	assembleHuangGradLagrange3d(&BB, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, elementdofTran);
	assembleMassmatrixLagrange3d(&CC, elements, elementDOF+1, elementdofTran+1);
    assembleRHSHuang3d(&bb, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, elementdofTran, quadcurl3d_f);
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

/**
 * \fn void assemble_quadcurl_StokesNcP1P03d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, dvector *wh, ELEMENT_DOF *elementDOFwh)
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
void assemble_quadcurl_StokesNcP1P03d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, dvector *wh, ELEMENT_DOF *elementDOFwh)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Cx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	dCSRmat AA, BB;
	dvector bb;

	assembleBiGradVectorNcP13d(&AA, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, elementdofTran);
	assembleDivNcP1P03d(&BB, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, elementdofTran);
	A[3].row = 0;
	A[3].col = 0;
	A[3].IA = NULL;
	A[3].JA = NULL;
	A[3].val = NULL;
	assembleRHSCurlHuangNcP13d4Stokes(&bb, &b[1], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, elementdofTran, wh, elementDOFwh);
    // initial solution
	create_dvector(bb.row, &uh[0]);
	create_dvector(b[1].row, &uh[1]);	
	
	// extract
	extractFreenodesVector2StokesDirichlet(&AA, &BB, &bb, b, &elementDOF[0], &uh[0]);
	free_dvector(&bb);
	extractFreenodesMatrix11(&AA, &A[0], &elementDOF[0], &elementDOF[0]);
	free_csr_matrix(&AA);
	extractFreenodesMatrix1c(&BB, &A[2], &elementDOF[0]);
	free_csr_matrix(&BB);
	getTransposeOfSparse(&A[2], &A[1]);
}

/**
 * \fn void assemble_quadcurlHuang3d(dCSRmat *A, dvector *b, dvector *uh, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran)
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
void assemble_quadcurlHuangZhang3d(dCSRmat *A, dvector *b, dvector *uh, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Dx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	dCSRmat AA, BB, CC, tempA;
	dCSRmat B, BT, C;
	dvector bb, D;

	assembleBiGradcurlHuangZhang3d(&AA, basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, elementdofTran);
	assembleHuangZhangGradLagrange3d(&BB, basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, elementdofTran);
	assembleMassmatrixLagrange3d(&CC, elements, elementDOF+1, elementdofTran+1);
    assembleRHSHuangZhang3d(&bb, basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, elementdofTran, quadcurl3d_f);
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

/**
 * \fn void assembleRHSCurlHuangNcP13d4Stokes(dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, dvector *wh, ELEMENT_DOF *elementDOFwh)
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
 * \param *elementdofTran pointer to transpose of elementDOF
 * \param lambda Lame constant
 * \param mu Lame constant or Poisson ratio of plate
 * \return void
 */
void assembleRHSCurlHuangNcP13d4Stokes(dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, dvector *wh, ELEMENT_DOF *elementDOFwh)
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
	
	int *index;
	int istart;
	
	 /************************************************** right hand side b *****************************************************************/
	create_dvector(elementDOF[0].dof * 3, b1);
	create_dvector(elementDOF[1].dof, b2);
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
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (i1 = 0; i1<num_qp; i1++)
			{
				ncp13d_basis(lambdas[i1], k1, &phi);
				init_array(3, phi2, 0.0);
				for (j = 0; j<elementDOFwh[0].col; j++)
				{
					huangQuadcurl3d_basisCurl(lambdas[i1], grd_lambda, nvf, eorien, fpermi, j, phi1);
					j1 = elementDOFwh[0].val[k][j];
					axpy_array(3, wh[0].val[j1], phi1, phi2);
				}
				b1->val[i] += vol*weight[i1] * phi2[0] * phi;
				b1->val[i + elementDOF[0].dof] += vol*weight[i1] * phi2[1] * phi;
				b1->val[i + elementDOF[0].dof*2] += vol*weight[i1] * phi2[2] * phi;
			} // i1
		} // k1 
	} // k
}

/**
 * \fn void assembleRHSCurlHuangNcP13d4Maxwell(dvector *rhs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, dvector *uhs, ELEMENT_DOF *elementDOFs)
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
 * \param *elementdofTran pointer to transpose of elementDOF
 * \param lambda Lame constant
 * \param mu Lame constant or Poisson ratio of plate
 * \return void
 */
void assembleRHSCurlHuangNcP13d4Maxwell(dvector *rhs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, dvector *uhs, ELEMENT_DOF *elementDOFs)
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
	
	int *index;
	int istart;
	
	/************************************************** right hand side b *****************************************************************/
	dvector b;
	create_dvector(elementDOF[0].dof, &b);
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
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (i1 = 0; i1<num_qp; i1++)
			{
				huangQuadcurl3d_basisCurl(lambdas[i1], grd_lambda, nvf, eorien, fpermi, k1, phi1);
				init_array(3, phi2, 0.0);
				for (j = 0; j<elementDOFs[0].col; j++)
				{
					ncp13d_basis(lambdas[i1], j, &phi);
					j1 = elementDOFs[0].val[k][j];
					phi2[0] += phi * uhs[0].val[j1];
					phi2[1] += phi * uhs[0].val[j1 + elementDOFs[0].dof];
					phi2[2] += phi * uhs[0].val[j1 + elementDOFs[0].dof*2];
				}
				b.val[i] += vol*weight[i1] * dot_array(3, phi1, phi2);
			} // i1
		} // k1 
	} // k

	ivector *freenodes = &elementDOF[0].freenodes;
	create_dvector(freenodes->row, rhs);
	for (i1 = 0; i1 < rhs->row; i1++)
	{
		i = freenodes->val[i1];
		rhs->val[i1] = b.val[i];
	}
	free_dvector(&b);
}

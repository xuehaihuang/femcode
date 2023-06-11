/*
 *  stokesfem3d.c
 *
 *  Created by Xuehai Huang on May 8, 2023.
 *  Copyright 2023 SUFE. All rights reserved.
 *
 */

/*! \file stokesfem3d.c
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
 * \fn void stokesfem3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
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
void stokesfem3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i;
	int variationalform_type = Input->variationalform_type;
	int fem_num = Input->fem_num;
	if (variationalform_type == 1) // primal formulation
	{
		printf("Discretize Stokes equation in primal formulation.\n");
		if(fem_num == 1)
		{
			printf("Discretized by Nonconforming P1-P0 element.\n");
			stokesNcP1P03d(elements, elementFace, faces, elementEdge, edges, nodes, Input);
		}
		
		if(fem_num == 2)
		{
			// printf("Discretized by Lagrange element.\n");
			// poissonLagrange3d(elements, elementFace, faces, elementEdge, edges, nodes, Input);
			printf("To be developed!\n");
			exit(0);
		}
	}
	else if (variationalform_type == 2) // decoupled mixed formulation
	{
//		printf("Discretize Poission equation in mixed formulation.\n");
//		poissonfem_mixed(uh, elements, elementEdge, edges, nodes, edgesTran, nodeCEdge, elementDOF, Input);
		printf("To be developed!\n");
		exit(0);
	}
	else
	{
		printf("Please set variationalform_type = 1 or 2!\n");
		exit(0);
	}
}

/**
* \fn void stokesNcP1P03d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief Second kind of Nedelec element method for Maxwell equation in three dimensions
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
void stokesNcP1P03d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i,j;
	dCSRmat A[4];
	dvector b[2], uh[2];
	ELEMENT_DOF elementDOF[2];
		
	/** Step 0. Read input parameters */
	int print_level = Input->print_level;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;
	int dop = Input->dop1;
		
	/** Step 1. generate degrees of freedom */
	getElementDOF_NoncfmP13d(&elementDOF[0], elementFace, faces->row);
	getFreenodesInfoNoncfmP1Vector3d(faces, &elementDOF[0]);

	getElementDOFdg(&elementDOF[1], elements->row, 1, 0);
	getFreenodesInfo(&elementDOF[1]);
			
	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_stokesNcP1P03d(A, b, uh, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);

//	print_dcsr_matrix(&A);///////////
//	print_dvector(0, &b);///////
	//print_darray(20, A.val);///
	
	/** Step 3. Check matrix properties */
	// check_symm(&A);
	// check_diagpos(&A);
	// check_diagdom(&A);
	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A[0].row, A[0].col, A[0].nnz);
	printf("B.row=%d, B.col=%d, B.nnz=%d\n", A[2].row, A[2].col, A[2].nnz);
	

	/** Step 4. Solve the system */

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", Input->itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", Input->itsolver_tol);
	}

	// create_dvector(b.row, &_uh);
	// init_dvector(&_uh, 0.0);

	// printf("A.row=%d, A.col=%d, A.nnz=%d\n", A.row, A.col, A.nnz);

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
	aspparam.elementDOF = elementDOF;
	
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

	if (itsolver_type == 1)
	{
		printf("\nASP Approximate Block Factorization preconditioned GMRES solver with auxiliary space method.\n\n");
		printf("Auxiliary space: P1 Lagrangian element.\n\n");
		AbfpAsP1StokesNcP1P0_GMRES(A, b, uh, &aspparam, print_level);
	//	Abfp2AsP1StokesNcP1P0_GMRES(A, b, uh, &aspparam, print_level);
	}
	else
	{
		printf("\nASP Block diagonal preconditioned MINRES solver with auxiliary space method\n");
		printf("Auxiliary space: P1 Lagrangian element.\n\n");
		DiagAsP1StokesNcP1P0_MINRES(A, b, uh, &aspparam, print_level);
	//	DiagAMGStokesNcP1P0_MINRES(A, b, uh, &aspparam, print_level); // sometimes better
	}
////////////////////////////////////////////////


	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */

	/** Step 5. Compute the error between numerical solution and exact solution */
	double errors[4];

	geterrorsStokesNcP1P03d(errors, uh, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);

	printf("\nA Priori Errors:\n");
	printf("L2 norm of u-u_h = %e\n", errors[0]);
	printf("H1 semi-norm of u-u_h = %e\n", errors[1]);
	printf("H1 norm of u-u_h = %e\n", errors[2]);
	printf("L2 norm of p-p_h = %e\n", errors[3]);
	
	/*********************************************************************************************/
	FILE *outputFile;
	outputFile = fopen("output/error.dat", "w");
	fprintf(outputFile, "%e\t%e\t%e\t%e\n", errors[0], errors[1], errors[2], errors[3]);
	fclose(outputFile);
	/********************************************************************************************/

	
	for(i=0;i<4;i++) free_csr_matrix(&A[i]);
	free_dvector(&b[0]); free_dvector(&b[1]);
	free_dvector(&uh[0]); free_dvector(&uh[1]);
	free_elementDOF(&elementDOF[0]); free_elementDOF(&elementDOF[1]);
}

/**
 * \fn void assemble_stokesNcP1P03d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
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
 * \return void
 */
void assemble_stokesNcP1P03d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1         = 0
	where x1: u_h, x2: p_h
	**/
	int i;
	dCSRmat AA[4];

	assembleBiGradVectorNcP13d(&AA[0], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	assembleDivNcP1P03d(&AA[2], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	getTransposeOfSparse(&AA[2], &AA[1]);
	create_csr_matrix(elementDOF[1].dof, elementDOF[1].dof, 0, &AA[3]);
	assembleRHSVectorNcP13d(b, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, stokes3d_f);
	assembleRHSdgPoly3d(b+1, elements, elementDOF+1, stokes3d_g);
	create_dvector(AA[1].row, uh);
	create_dvector(AA[1].col, uh+1);
	
	// Apply boundary condition
	updateFreenodes2bRHS(AA, b, uh, elementDOF, elementDOF+1);
	updateFreenodes2bMatrix11(AA, A, elementDOF, elementDOF+1);
	for(i=0;i<4;i++) free_csr_matrix(AA+i);	
}

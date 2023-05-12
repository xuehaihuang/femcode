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
 * \fn void maxwellfem(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
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
void maxwellfem(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i;
	int variationalform_type = Input->variationalform_type;
	int fem_num = Input->fem_num;
//	Input->itsolver_type = 3; //////////////////
	if (variationalform_type == 1) // primal formulation
	{
		printf("Discretize maxwell equation in primal formulation.\n");
		if(fem_num == 1)
		{
			printf("Discretized by the first kind of Nedelec element with k = %d.\n", Input->dop1);
			maxwellNedelec1st3d(elements, elementFace, faces, elementEdge, edges, nodes, Input);
		}
		
		if(fem_num == 2)
		{
			printf("Discretized by the second kind of Nedelec element.\n");
			maxwellNedelec2nd3d(elements, elementFace, faces, elementEdge, edges, nodes, Input);
		}

		if(fem_num == 3)
		{
			printf("Discretized by the Huang-Zhang element.\n");
			maxwellHuangZhang3d(elements, elementFace, faces, elementEdge, edges, nodes, Input);
		}

		if(fem_num == 0)
		{
			printf("Discretized by Huang element.\n");
			maxwellHuang3d(elements, elementFace, faces, elementEdge, edges, nodes, Input);
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
* \fn void maxwellNedelec1st3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief First kind of Nedelec element method for Maxwell equation in three dimensions
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
void maxwellNedelec1st3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
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
	int dop1 = Input->dop1;
	int dop2 = dop1;
		
	/** Step 1. generate degrees of freedom */
	getElementDOF_Nedelec1st3d(elementDOF, elements, elementFace, faces, elementEdge, edges, dop1);
	getFreenodesInfoNedelec1st3d(faces, edges, nodes, elementDOF);

	getElementDOF_Lagrange3d(elementDOF+1, elements, elementFace, faces, elementEdge, edges, nodes->row, dop2);
	getFreenodesInfoLagrange3d(faces, edges, nodes, elementDOF+1);
		
	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_maxwellNedelec1st3d(&A, &b, &uh, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);

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

	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */

	/** Step 5. Compute the error between numerical solution and exact solution */
	double errors[3];

	geterrorsMaxwellNedelec1st3d(errors, &uh, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);

	printf("\nA Priori Errors:\n");
	printf("L2 norm of u-u_h = %e\n", errors[0]);
	printf("L2 norm of curl(u-u_h) = %e\n", errors[1]);
	printf("Energy norm of u-u_h = %e\n", errors[2]);

	/*********************************************************************************************/
	FILE *outputFile;
	outputFile = fopen("output/error.dat", "w");
	fprintf(outputFile, "%e\t%e\t%e\n", errors[0], errors[1], errors[2]);
	fclose(outputFile);
	/********************************************************************************************/

	
	free_csr_matrix(&A);
	free_dvector(&b);
	free_dvector(&uh);
	free_elementDOF(elementDOF);
	free_elementDOF(elementDOF+1);
}

/**
* \fn void maxwellNedelec2nd3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
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
void maxwellNedelec2nd3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
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
	int dop1 = Input->dop1;
	int dop2 = dop1+1;
		
	/** Step 1. generate degrees of freedom */
	getElementDOF_Nedelec2nd3d(elementDOF, elements, elementFace, faces, elementEdge, edges, dop1);
	getFreenodesInfoNedelec2nd3d(faces, edges, nodes, elementDOF);

	getElementDOF_Lagrange3d(elementDOF+1, elements, elementFace, faces, elementEdge, edges, nodes->row, dop2);
	getFreenodesInfoLagrange3d(faces, edges, nodes, elementDOF+1);
		
	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_maxwellNedelec2nd3d(&A, &b, &uh, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);

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

	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */

	/** Step 5. Compute the error between numerical solution and exact solution */
	double errors[3];

	geterrorsMaxwellNedelec2nd3d(errors, &uh, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);

	printf("\nA Priori Errors:\n");
	printf("L2 norm of u-u_h = %e\n", errors[0]);
	printf("L2 norm of curl(u-u_h) = %e\n", errors[1]);
	printf("Energy norm of u-u_h = %e\n", errors[2]);

	/*********************************************************************************************/
	FILE *outputFile;
	outputFile = fopen("output/error.dat", "w");
	fprintf(outputFile, "%e\t%e\t%e\n", errors[0], errors[1], errors[2]);
	fclose(outputFile);
	/********************************************************************************************/

	
	free_csr_matrix(&A);
	free_dvector(&b);
	free_dvector(&uh);
	free_elementDOF(elementDOF);
	free_elementDOF(elementDOF+1);
}

/**
* \fn void maxwellHuangZhang3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief Huang-Zhang element method for Maxwell equation in three dimensions
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
void maxwellHuangZhang3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
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
	int dop1 = 2;
	int dop2 = 2;
		
	/** Step 1. generate degrees of freedom */
	getElementDOF_HuangZhang3d(elementDOF, elements, elementFace, faces, elementEdge, edges);
	getFreenodesInfoHuangZhang3d(faces, edges, nodes, 2, elementDOF);
	
	getElementDOF_Lagrange3d(elementDOF+1, elements, elementFace, faces, elementEdge, edges, nodes->row, dop2);
	getFreenodesInfoLagrange3d(faces, edges, nodes, elementDOF+1);

	/***************************Generate coefficient of basis functions**************************/
	ddenmat3 basisCoeffs;
	create_dden_matrix3(elements->row, elementDOF->col, elementDOF->col, &basisCoeffs);
	getHuangZhangBasisCoeffs(&basisCoeffs, elements, elementFace, faces, elementEdge, edges);
	/********************************************************************************************/

	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_maxwellHuangZhang3d(&A, &b, &uh, &basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);

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

	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */

	/** Step 5. Compute the error between numerical solution and exact solution */
	double errors[3];

	geterrorsMaxwellHuangZhang3d(errors, &uh, &basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);

	printf("\nA Priori Errors:\n");
	printf("L2 norm of u-u_h = %e\n", errors[0]);
	printf("L2 norm of curl(u-u_h) = %e\n", errors[1]);
	printf("Energy norm of u-u_h = %e\n", errors[2]);

	/*********************************************************************************************/
	FILE *outputFile;
	outputFile = fopen("output/error.dat", "w");
	fprintf(outputFile, "%e\t%e\t%e\n", errors[0], errors[1], errors[2]);
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
* \fn void maxwellHuang3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief Huang element method for Maxwell equation in three dimensions
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
void maxwellHuang3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
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
	int dop1 = 2;
	int dop2 = 1;
		
	/** Step 1. generate degrees of freedom */
	getElementDOF_HuangGradcurl3d(elementDOF, elements, elementFace, faces, elementEdge, edges);
	getFreenodesInfoHuangGradcurl3d(faces, edges, nodes, 2, elementDOF);

	getElementDOF_Lagrange3d(elementDOF+1, elements, elementFace, faces, elementEdge, edges, nodes->row, dop2);
	getFreenodesInfoLagrange3d(faces, edges, nodes, elementDOF+1);
		
	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_maxwellHuang3d(&A, &b, &uh, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, maxwell3d_f);

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

	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */

	/** Step 5. Compute the error between numerical solution and exact solution */
	double errors[3];

	geterrorsMaxwellHuang3d(errors, &uh, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);

	printf("\nA Priori Errors:\n");
	printf("L2 norm of u-u_h = %e\n", errors[0]);
	printf("L2 norm of curl(u-u_h) = %e\n", errors[1]);
	printf("Energy norm of u-u_h = %e\n", errors[2]);

	/*********************************************************************************************/
	FILE *outputFile;
	outputFile = fopen("output/error.dat", "w");
	fprintf(outputFile, "%e\t%e\t%e\n", errors[0], errors[1], errors[2]);
	fclose(outputFile);
	/********************************************************************************************/

	
	free_csr_matrix(&A);
	free_dvector(&b);
	free_dvector(&uh);
	free_elementDOF(elementDOF);
	free_elementDOF(elementDOF+1);
}

/**
 * \fn void assemble_maxwellNedelec1st3d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
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
void assemble_maxwellNedelec1st3d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Dx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	int i;
	dCSRmat AA[4], AB[4];
	dvector D, bb[2], wh[2];
	assembleBiCurlNedelec1st3d(&AA[0], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	assembleNedelec1stGradLagrange3d(&AA[2], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	getTransposeOfSparse(&AA[2], &AA[1]);
	assembleMassmatrixLagrange3d(&AA[3], elements, elementDOF+1);
	assembleRHSNedelec1st3d(&bb[0], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, maxwell3d_f);
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
 * \fn void assemble_maxwellNedelec2nd3d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
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
void assemble_maxwellNedelec2nd3d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Dx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	int i;
	dCSRmat AA[4], AB[4];
	dvector D, bb[2], wh[2];
	assembleBiCurlNedelec2nd3d(&AA[0], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	assembleNedelec2ndGradLagrange3d(&AA[2], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	getTransposeOfSparse(&AA[2], &AA[1]);
	assembleMassmatrixLagrange3d(&AA[3], elements, elementDOF+1);
	assembleRHSNedelec2nd3d(&bb[0], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, maxwell3d_f);
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
 * \fn void assemble_maxwellHuangZhang3d(dCSRmat *A, dvector *b, dvector *uh, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
 * \brief assemble stiffness matrix *C and righ hand side *ptr_b
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side
 * \param *uh pointer to initial value
 * \param *basisCoeffs pointer to coefficients of basis
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
void assemble_maxwellHuangZhang3d(dCSRmat *A, dvector *b, dvector *uh, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Dx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	int i;
	dCSRmat AA[4], AB[4];
	dvector D, bb[2], wh[2];
	assembleBiCurlHuangZhang3d(&AA[0], basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	assembleHuangZhangGradLagrange3d(&AA[2], basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	getTransposeOfSparse(&AA[2], &AA[1]);
	assembleMassmatrixLagrange3d(&AA[3], elements, elementDOF+1);
	assembleRHSHuangZhang3d(&bb[0], basisCoeffs, elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, maxwell3d_f);
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
 * \fn void assemble_maxwellHuang3d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double (*f)(double *))
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
void assemble_maxwellHuang3d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, void (*f)(double *, double *))
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Dx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	int i;
	dCSRmat AA[4], AB[4], tempA[2];
	dvector D, bb[2], wh[2];
	assembleBiCurlHuang3d(&AA[0], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	assembleHuangGradLagrange3d(&AA[2], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	getTransposeOfSparse(&AA[2], &AA[1]);
	assembleMassmatrixLagrange3d(&AA[3], elements, elementDOF+1);
	assembleRHSHuang3d(&bb[0], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, f);
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
 * \fn void assemble_maxwellHuang3d0(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double (*f)(double *))
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
void assemble_maxwellHuang3d0(dCSRmat *AA, dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, void (*f)(double *, double *))
{
	/**	
	Ax1 + B^Tx2 = b
	Bx1 - Dx2   = 0
	where x1: u_h, x2: lambda_h
	**/
	int i;
	dCSRmat AB[4], tempA[2];
	dvector D, bb[2], wh[2];
	assembleBiCurlHuang3d(&AA[0], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	assembleHuangGradLagrange3d(&AA[2], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF);
	getTransposeOfSparse(&AA[2], &AA[1]);
	assembleMassmatrixLagrange3d(&AA[3], elements, elementDOF+1);
	assembleRHSHuang3d(&bb[0], elements, elementFace, faces, elementEdge, edges, nodes, elementDOF, f);
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
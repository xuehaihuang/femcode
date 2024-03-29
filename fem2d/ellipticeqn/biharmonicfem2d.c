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
 * \fn void biharmonicfem2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
 * \brief finite element methods for biharmonic equation in two dimensions
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
void biharmonicfem2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i;
	int variationalform_type = Input->variationalform_type;
	int fem_num = Input->fem_num;
	int glevelNum = Input->glevelNum;
	short extrap = Input->extrap; // for extrapolation

	if (variationalform_type == 1) // primal formulation
	{
		printf("Discretize biharmonic equation in primal formulation.\n");
		if(fem_num == 1)
		{
			printf("Discretized by Morley element.\n");
			// biharmonicMorley2d(elements, elementEdge, edges, nodes, Input);
			biharmonicMorley2d(&elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], Input);
		}
		
		if(fem_num == 2)
		{
			printf("Discretized by C0 IPDG with Lagrange element.\n");
			if(extrap == 1){
				if(glevelNum>1){
					printf("Extrapolation is used.\n");
					int cl = glevelNum - 2;
					biharmonicC0ipdgExtrap2d(elements+cl, elementEdge+cl, edges+cl, nodes+cl, Input);
				}
				else{
					printf("Extrapolation: Level of meshes must be greater than 1.\n");
				}
			}
			else{
				// biharmonicC0ipdg2d(elements, elementEdge, edges, nodes, Input);
				biharmonicC0ipdg2d(&elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], Input);
			}
			
		}
	}
	else if (variationalform_type == 2) // decoupled mixed formulation
	{
		printf("Discretize biharmonic equation in decoupled Poisson-Stokes-Poisson formulation.\n");
		if(fem_num == 1)
		{
			printf("Discretized by Morley : Nonconforming P1-P0 : Morley element.\n");
			biharmonicPSP_MorleyNcfmP1P0Morley2d(&elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], Input);
		}
	}
	else
	{
		printf("Please set variationalform_type = 1 or 2!\n");
		exit(0);
	}
}

/**
* \fn void biharmonicMorley2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief Morley element method for Biharmonic equation in two dimensions
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
void biharmonicMorley2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i,j;
	dvector uh;
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
	// int dop = Input->dop1;
		
	/** Step 1. generate degrees of freedom */
	getElementDOF_Morley2d(&elementDOF, elements, elementEdge, edges, nodes->row);
	getFreenodesInfoMorley2d(edges, nodes, &elementDOF);
	// getTransposeOfelementDoF(&elementDOF, &elementdofTran, 0);

	/** Step 2. assemble stiffmatrix and right hand side term */
	solve_biharmonicMorley2d(&uh, &elementDOF, elements, elementEdge, edges, nodes, Input);

	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */

	/** Step 3. Compute the error between numerical solution and exact solution */
	double errors[4];

	geterrorsBiharmonicMorley2d(errors, &uh, elements, elementEdge, edges, nodes, &elementDOF);

	printf("\nA Priori Errors:\n");
	printf("L2 norm of u-u_h = %e\n", errors[0]);
	printf("H1 semi-norm of u-u_h = %e\n", errors[1]);
	printf("H2 semi-norm of u-u_h = %e\n", errors[2]);
	printf("H2 norm of u-u_h = %e\n", errors[3]);
	
	/*********************************************************************************************/
	FILE *outputFile;
	outputFile = fopen("output/error.dat", "w");
	fprintf(outputFile, "%e\t%e\t%e\t%e\n", errors[0], errors[1], errors[2], errors[3]);
	fclose(outputFile);
	/********************************************************************************************/

	
	// free_csr_matrix(&A);
	// free_dvector(&b);
	free_dvector(&uh);
	free_elementDOF(&elementDOF);
}

/**
* \fn void solve_biharmonicMorley2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief Morley element method for Biharmonic equation in two dimensions
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
void solve_biharmonicMorley2d(dvector *uh, ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i,j;
	dCSRmat A;
	dvector b;
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
		
	/** Step 1. generate degrees of freedom */
	// getElementDOF_Morley2d(&elementDOF, elements, elementEdge, edges, nodes->row);
	// getFreenodesInfoMorley2d(edges, nodes, &elementDOF);
	// getTransposeOfelementDoF(elementDOF, &elementdofTran, 0);

	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_biharmonicMorley2d(&A, &b, uh, elements, elementEdge, edges, nodes, elementDOF);
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

	create_dvector(b.row, uh);
	init_dvector(uh, 0.0);

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

	/* PCG+diag */
	if (itsolver_type == 1) {
		printf("Diagonal preconditioned CG solver\n");
		diag_PCG(&A, &b, uh, itsolver_maxit, itsolver_tol, print_level);
	}
	
	/* PCG+AMG */
	else if (itsolver_type == 2) {
		printf("AMG iterative solver\n");
		classicAMG_PCG(&A, &b, uh, &amgparam, print_level);
	}

	/* AMG solver */
	else if (itsolver_type == 3) {
		printf("AMG preconditioned CG solver\n");
		classicAMG(&A, &b, uh, &amgparam);
	}


	/* CG */
	else if (itsolver_type == 4) {
		printf("Classical CG solver\n");
		standard_CG(&A, &b, uh, itsolver_maxit, itsolver_tol, print_level);
	}

	/* GMRES+AMG */
	else if (itsolver_type == 5) {
		printf("AMG preconditioned GMRES solver\n");
		classicAMG_GMRES(&A, &b, uh, &amgparam, print_level);
	}

	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */
	
	free_csr_matrix(&A);
	free_dvector(&b);
}

/**
 * \fn void assemble_biharmonicMorley2d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
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
void assemble_biharmonicMorley2d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	/**	
	Ax = b
	**/
	dCSRmat AA;

	assembleBiHessMorley2d(&AA, elements, elementEdge, edges, nodes, elementDOF);
	assembleRHSMorley2d(b, elements, elementEdge, edges, elementDOF, biharmonic2d_f, NULL);
    // initial solution
	create_dvector(b->row, uh);
	// Apply boundary condition
	updateFreenodesRHS(&AA, b, uh, elementDOF, elementDOF, 1);
	updateFreenodesMatrix11(&AA, A, elementDOF, elementDOF, 1);

	// extractFreenodesVector(&AA, &bb, b, elementDOF, uh);
	// // free_dvector(&bb);
	// extractFreenodesMatrix11(&AA, A, elementDOF, elementDOF);

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
* \fn void biharmonicC0ipdgExtrap2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief Extrapolation of C0 IPDG method for Biharmonic equation in two dimensions
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
void biharmonicC0ipdgExtrap2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i,j,l;
	dvector uh[2];
	ELEMENT_DOF elementDOF[2];
	iCSRmat elementdofTran[2];
		
	/** Step 0. Read input parameters */
	int print_level = Input->print_level;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;
	double lambda = Input->lambda;
	double mu = Input->mu;
	double parapenalty = Input->parapenalty;
	int glevelNum = Input->glevelNum;
	int CglevelNum = Input->CglevelNum;
	int FglevelNum = Input->FglevelNum;
	int dop = Input->dop1;

	double errors[9];

	FILE *outputFile;
	
	printf("k = %d,  parapenalty = %f\n", dop, parapenalty);

	for(l=0;l<2;l++){
		/** Step 1. generate degrees of freedom */
		getElementDOF_Lagrange2d(elementDOF+l, elements+l, elementEdge+l, edges+l, nodes[l].row, dop);
		getFreenodesInfoLagrange2d(edges+l, nodes+l, elementDOF+l);
		getTransposeOfelementDoF(elementDOF+l, elementdofTran+l, 0);

		/** Step 2. assemble and solve */
		solve_biharmonicC0ipdg2d(uh+l, elementDOF+l, elements+l, elementEdge+l, edges+l, nodes+l, elementdofTran+l, Input);

		/** Step 3. Compute the error between numerical solution and exact solution */
		
		geterrorsBiharmonicC0ipdg2d(errors, uh+l, elements+l, elementEdge+l, edges+l, nodes+l, elementDOF+l);

		printf("\nA Priori Errors:\n");
		printf("L2 norm of u-u_h = %e\n", errors[0]);
		printf("H1 semi-norm of u-u_h = %e\n", errors[1]);
		printf("H2 semi-norm of u-u_h = %e\n", errors[2]);
		printf("H2 norm of u-u_h = %e\n", errors[3]);
		printf("jump of u_h = %e\n", errors[4]);
		printf("Energy norm of u-u_h = %e\n", errors[5]);
	
		/*********************************************************************************************/
		if(l==0)
			outputFile = fopen("output/error.dat", "w");
		else
			outputFile = fopen("output/error.dat", "a");
		fprintf(outputFile, "%e\t%e\t%e\t%e\t%e\t%e\n", errors[0], errors[1], errors[2], errors[3], errors[4], errors[5]);
		fclose(outputFile);
		/********************************************************************************************/
	}

	double c0 =0, c1 = 1;
	if(dop == 2){
		c0 = -1.0/3.0;
		c1 = 4.0/3.0;
		// c0 = -1.0;
		// c1 = 2.0;
	}
	else if(dop == 3){
		c0 = -1.0/15.0;
		c1 = 16.0/15.0;
		// c0 = -1.0/3.0;
		// c1 = 4.0/3.0;
	}

	// double maxuh=0, val, *x;
	// for(i=0;i<nodes[0].row;i++){
	// 	x = nodes[0].val[i];
	// 	if(fabs(x[0]-0.5)<0.25+1e-10 && fabs(x[1]-0.5)<0.25+1e-10){
	// 		val = c0*uh[0].val[i] + c1*uh[1].val[i];
	// 		val -= biharmonic2d_u(x, NULL);
	// 		val = fabs(val);
	// 		if(val>maxuh) maxuh=val;
	// 	}
	// }
	// printf("maxuh = %e\n", maxuh);

	geterrorsBiharmonicC0ipdgExtrap2d(errors, uh, elements, elementEdge, edges, nodes, elementDOF, elementdofTran);

	printf("\nExtrapolation Errors:\n");
	printf("max norm of u-u_h = %e\n", errors[0]);
	printf("max norm of grad(u-u_h) = %e\n", errors[1]);
	printf("max norm of hess(u-u_h) = %e\n", errors[2]);
	printf("\nNon-extrapolation Errors:\n");
	printf("max norm of u-u_h = %e\n", errors[3]);
	printf("max norm of grad(u-u_h) = %e\n", errors[4]);
	printf("max norm of hess(u-u_h) = %e\n", errors[5]);
	printf("max norm of u-u_h = %e\n", errors[6]);
	printf("max norm of grad(u-u_h) = %e\n", errors[7]);
	printf("max norm of hess(u-u_h) = %e\n", errors[8]);

	/*********************************************************************************************/
	outputFile = fopen("output/error.dat", "a");	
	fprintf(outputFile, "%e\t%e\t%e\n", errors[0], errors[1], errors[2]);
	fprintf(outputFile, "%e\t%e\t%e\n", errors[3], errors[4], errors[5]);
	fprintf(outputFile, "%e\t%e\t%e\n", errors[6], errors[7], errors[8]);
	// fprintf(outputFile, "%e\n", &maxuh);
	fclose(outputFile);
	/********************************************************************************************/


	
	for(l=0;l<2;l++){
		free_dvector(uh+l);
		free_elementDOF(elementDOF+l);
		free_icsr_matrix(elementdofTran+l);
	}	
}

/**
* \fn void biharmonicC0ipdg2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief C0 IPDG method for Biharmonic equation in two dimensions
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
void biharmonicC0ipdg2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i,j;
	dvector uh;
	ELEMENT_DOF elementDOF;
	iCSRmat elementdofTran;
		
	/** Step 0. Read input parameters */
	int print_level = Input->print_level;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;
	double lambda = Input->lambda;
	double mu = Input->mu;
	double parapenalty = Input->parapenalty;
	int glevelNum = Input->glevelNum;
	int CglevelNum = Input->CglevelNum;
	int FglevelNum = Input->FglevelNum;
	int dop = Input->dop1;
	
	printf("k = %d,  parapenalty = %f\n", dop, parapenalty);

	/** Step 1. generate degrees of freedom */
	getElementDOF_Lagrange2d(&elementDOF, elements, elementEdge, edges, nodes->row, dop);
	getFreenodesInfoLagrange2d(edges, nodes, &elementDOF);
	getTransposeOfelementDoF(&elementDOF, &elementdofTran, 0);

	/** Step 2. assemble and solve */
	solve_biharmonicC0ipdg2d(&uh, &elementDOF, elements, elementEdge, edges, nodes, &elementdofTran, Input);
	


	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */

	/** Step 3. Compute the error between numerical solution and exact solution */
	double errors[6];

	geterrorsBiharmonicC0ipdg2d(errors, &uh, elements, elementEdge, edges, nodes, &elementDOF);

	printf("\nA Priori Errors:\n");
	printf("L2 norm of u-u_h = %e\n", errors[0]);
	printf("H1 semi-norm of u-u_h = %e\n", errors[1]);
	printf("H2 semi-norm of u-u_h = %e\n", errors[2]);
	printf("H2 norm of u-u_h = %e\n", errors[3]);
	printf("jump of u_h = %e\n", errors[4]);
	printf("Energy norm of u-u_h = %e\n", errors[5]);
	
	/*********************************************************************************************/
	FILE *outputFile;
	outputFile = fopen("output/error.dat", "w");
	fprintf(outputFile, "%e\t%e\t%e\t%e\t%e\t%e\n", errors[0], errors[1], errors[2], errors[3], errors[4], errors[5]);
	fclose(outputFile);
	/********************************************************************************************/

	
	// free_csr_matrix(&A);
	// free_dvector(&b);
	free_dvector(&uh);
	free_elementDOF(&elementDOF);
	free_icsr_matrix(&elementdofTran);
}

/**
* \fn void solve_biharmonicC0ipdg2d(dvector *uh, ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *elementdofTran, Input_data *Input)
* \brief C0 IPDG method for Biharmonic equation in two dimensions
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
void solve_biharmonicC0ipdg2d(dvector *uh, ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *elementdofTran, Input_data *Input)
{
	int i,j;
	dCSRmat A;
	dvector b;
	// ELEMENT_DOF elementDOF;
	// iCSRmat elementdofTran;
		
	/** Step 0. Read input parameters */
	int print_level = Input->print_level;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;
	double lambda = Input->lambda;
	double mu = Input->mu;
	double parapenalty = Input->parapenalty;
	int glevelNum = Input->glevelNum;
	int CglevelNum = Input->CglevelNum;
	int FglevelNum = Input->FglevelNum;
	int dop = Input->dop1;
	
	// printf("k = %d,  parapenalty = %f\n", dop, parapenalty);

	/** Step 1. generate degrees of freedom */
	// getElementDOF_Lagrange2d(&elementDOF, elements, elementEdge, edges, nodes->row, dop);
	// getFreenodesInfoLagrange2d(edges, nodes, &elementDOF);
	// getTransposeOfelementDoF(elementDOF, &elementdofTran, 0);
	
	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_biharmonicC0ipdg2d(&A, &b, uh, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, parapenalty);
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


	create_dvector(b.row, uh);
	init_dvector(uh, 0.0);


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

		read_dvector4Matlab(uh, outputfileuh);
	}

	/* PCG+diag */
	if (itsolver_type == 1) {
		printf("Diagonal preconditioned CG solver\n");
		diag_PCG(&A, &b, uh, itsolver_maxit, itsolver_tol, print_level);
	}

	/* PCG+AMG */
	else if (itsolver_type == 2) {
		printf("AMG preconditioned CG solver\n");
		classicAMG_PCG(&A, &b, uh, &amgparam, print_level);
	}

	/* AMG solver */
	else if (itsolver_type == 3) {
		printf("AMG iterative solver\n");
		classicAMG(&A, &b, uh, &amgparam);
	}

	
	/* CG */
	else if (itsolver_type == 4) {
		printf("Classical CG solver\n");
		standard_CG(&A, &b, uh, itsolver_maxit, itsolver_tol, print_level);
	}

	/* GMRES+AMG */
	else if (itsolver_type == 5) {
		printf("AMG preconditioned GMRES solver\n");
		classicAMG_GMRES(&A, &b, uh, &amgparam, print_level);
	}

	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */
	
	printf("\n A.row=%d, A.col=%d, A.nnz=%d\n", A.row, A.col, A.nnz);

	free_csr_matrix(&A);
	free_dvector(&b);
}

/**
 * \fn void assemble_biharmonicC0ipdg2d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double parapenalty)
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
void assemble_biharmonicC0ipdg2d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double parapenalty)
{
	/**	
	Ax = b
	**/
	dCSRmat AA;

	assembleBiHessC0ipdg2d(&AA, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, parapenalty);
	assembleRHSLagrange2d(b, elements, elementDOF, biharmonic2d_f, NULL);
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
* \fn void biharmonicPSP_MorleyNcfmP1P0Morley2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief Morley : Nonconforming P1-P0 : Morley element element method for the decoupled Poisson-Stokes-Poisson formulation of biharmonic equation in two dimensions
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
void biharmonicPSP_MorleyNcfmP1P0Morley2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
{
	int i,j;
	dvector uh;
	ELEMENT_DOF elementDOF[4];
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
	getElementDOF_CrouzeixRaviart2d(elementDOF+2, elements, elementEdge, edges);
	getFreenodesInfoCrouzeixRaviart2d(edges, elementDOF+2);
	repeatElementDoF(elementDOF+2, elementDOF+1, 2);
	getElementDOFdg(elementDOF+3, elements->row, 1, 0);
	getFreenodesInfoDG(elementDOF+3);

	/** Step 2. assemble stiffmatrix and right hand side term */
	solve_biharmonicPSP_MorleyNcfmP1P0Morley2d(&uh, elementDOF, elements, elementEdge, edges, nodes, Input);

	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */

	/** Step 3. Compute the error between numerical solution and exact solution */
	double errors[4];

	geterrorsBiharmonicMorley2d(errors, &uh, elements, elementEdge, edges, nodes, elementDOF);
	// geterrorsBiharmonicPSP_MorleyNcfmP1P0Morley2d(errors, &uh, elements, elementEdge, edges, nodes, elementDOF);

	printf("\nA Priori Errors:\n");
	printf("L2 norm of u-u_h = %e\n", errors[0]);
	printf("H1 semi-norm of u-u_h = %e\n", errors[1]);
	printf("H2 semi-norm of u-u_h = %e\n", errors[2]);
	printf("H2 norm of u-u_h = %e\n", errors[3]);
	
	/*********************************************************************************************/
	FILE *outputFile;
	outputFile = fopen("output/error.dat", "w");
	fprintf(outputFile, "%e\t%e\t%e\t%e\n", errors[0], errors[1], errors[2], errors[3]);
	fclose(outputFile);
	/********************************************************************************************/

	
	// free_csr_matrix(&A);
	// free_dvector(&b);
	free_dvector(&uh);
	for(i=0;i<4;i++) free_elementDOF(elementDOF+i);
}

/**
* \fn void solve_biharmonicPSP_MorleyNcfmP1P0Morley2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
* \brief Morley element method for Biharmonic equation in two dimensions
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
void solve_biharmonicPSP_MorleyNcfmP1P0Morley2d(dvector *uh, ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
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
	/************ Eqn1: Poisson equation disctetized by Morley element ************/
	assembleBiGradMorley2d(&AAm, elements, elementEdge, edges, nodes, elementDOF);
	assembleRHSMorley2d(b, elements, elementEdge, edges, elementDOF, biharmonic2d_f, NULL);
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
	
	/************ Eqn2: Stokes equation disctetized by nonconforming P1-P0 element ************/
	assembleBiGradCrouzeixRaviartDiagBlockRepeat2d(&AA[0], elements, elementEdge, edges, nodes, elementDOF+2, 0.5);
	assembleDivCrouzeixRaviartL2poly2d(&AA[1], elements, elementDOF+2, elementDOF+3);
	getTransposeOfSparse(&AA[1], &AA[2]);
	create_csr_matrix(AA[1].col, AA[1].col, 0, &AA[3]);
	assembleRHScurlMorleyCrouzeixRaviart2d(b, elements, elementEdge, edges, elementDOF, elementDOF+2, wh);
	free_dvector(wh);
	create_dvector(AA[1].col, b+1);
	create_dvector(AA[1].row, wh);
	create_dvector(AA[1].col, wh+1);
	// Apply boundary condition
	updateFreenodes2bRHS(AA, b, wh, elementDOF+1,elementDOF+3);
	updateFreenodes2bMatrix11(AA, A, elementDOF+1,elementDOF+3);
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
		AbfpAsP1StokesNcP1P0_GMRES(A, b, wh, &aspparam, print_level);
	//	Abfp2AsP1StokesNcP1P0_GMRES(A, b, wh, &aspparam, print_level);
	}
	else{
		printf("\nASP Block diagonal preconditioned MINRES solver with auxiliary space method\n");
		printf("Auxiliary space: P1 Lagrangian element.\n\n");
		DiagAsP1StokesNcP1P0_MINRES(A, b, wh, &aspparam, print_level);
	//	DiagAMGStokesNcP1P0_MINRES(A, b, wh, &aspparam, print_level); // sometimes better
	}

	for(i=0;i<4;i++) free_csr_matrix(A+i);
	free_dvector(&b[0]); free_dvector(&b[1]);

	/************ Eqn3: Poisson equation disctetized by Morley element ************/
	assembleRHSCrouzeixRaviartcurlMorley2d(b, elements, elementEdge, edges, elementDOF+2, elementDOF, wh);
	free_dvector(&wh[0]); free_dvector(&wh[1]);
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

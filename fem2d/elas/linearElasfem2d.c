/*
 *  linearElasfem2d.c
 *
 *  Created by Xuehai Huang on Jun 06, 2022.
 *  Copyright 2022 SUFE. All rights reserved.
 *
 */

/*! \file linearElasfem2d.c
 *  \brief FEM for linear elasticity
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "header.h"
#include "matvec.h"
#include "checkmat.h"

 /**
 * \fn void linearElas2dfem(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, Input_data *Input)
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
void linearElas2dfem(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, Input_data *Input)
{
	int i;
	int variationalform_type = Input->variationalform_type;
	int fem_num = Input->fem_num;
	if (variationalform_type == 1) // primal formulation
	{
		printf("Discretize Linear Elasticity in primal formulation.\n");
		if(fem_num == 1)
		{
			// printf("Discretized by Lagrange element.\n");
			// poissonLagrange3d(elements, elementFace, faces, elementEdge, edges, nodes, Input);	
		}
		printf("To be developed!\n");
		exit(0);
	}
	else if (variationalform_type == 2) // Stress-Displacement formulation
	{
		printf("Discretize Linear Elasticity in Stress-Displacement formulation.\n");
		if(fem_num == 1)
		{
			printf("Discretized by Hu-Zhang element.\n");
			linearElasHuZhang2d_mfem(elements, elementEdge, edges, nodes, edgesTran, nodeCEdge, Input);	
		}
		else{
			printf("To be developed!\n");
			exit(0);
		}
		
	}
	else if (variationalform_type == 3) // Lame system
	{
//		printf("Discretize Lame system in Stress-Displacement formulation.\n");
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
* \fn void linearElasHuZhang2d_mfem(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input)
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
void linearElasHuZhang2d_mfem(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, Input_data *Input)
{
	int i,j;
	dCSRmat A[4];
	dvector b[2], uh[2];
	ELEMENT_DOF elementDOF[3];
	iCSRmat elementdofTran;
		
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
	int dop2 = dop1-1;
	if(dop2 < 0) dop2 = 0;

	double paras[2];
	paras[0] = lambda;
	paras[1] = mu;

	printf("k = %d,  lambda = %f,  mu = %f\n", dop1, lambda, mu);

	/** Step 1. generate degrees of freedom */
	getElementDOF_HuZhang(&elementDOF[0], elements, elementEdge, edges, nodes->row, dop1);
	getElementDOF(&elementDOF[1], elements->row, dop2);
	getElementDOF(&elementDOF[2], elements->row, dop2+2);

	getTransposeOfelementDoF(elementDOF, &elementdofTran, 0);

		
	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_linearElasHuZhang2d(A, b, elements, elementEdge, edges, nodes, elementDOF, &elementdofTran, lambda, mu);

	// print_dcsr_matrix(&A[0]);///////////
	// print_dcsr_matrix(&A[1]);///////////
	// print_dvector(0, &b[0]);///////
	// print_dvector(0, &b[1]);///////
	// print_darray(20, A[0].val);///
	
	/** Step 3. Check matrix properties */
	for (i = 0; i<3; i++)
	{
		check_symm(&A[i]);
		check_diagpos(&A[i]);
		check_diagdom(&A[i]);
	}

	/** Step 4. Solve the system */

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", Input->itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", Input->itsolver_tol);
	}

	// create_dvector(b.row, &_uh);
	// init_dvector(&_uh, 0.0);
	create_dvector(b[0].row, &uh[0]);
	create_dvector(b[1].row, &uh[1]);
	init_dvector(&uh[0], 0.0);
	init_dvector(&uh[1], 0.0);

	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A[0].row, A[0].col, A[0].nnz);
	printf("B.row=%d, B.col=%d, B.nnz=%d\n", A[1].row, A[1].col, A[1].nnz);
	if (A[3].row>0)
		printf("C.row=%d, C.col=%d, C.nnz=%d\n", A[3].row, A[3].col, A[3].nnz);
	// printf("A.row=%d, A.col=%d, A.nnz=%d\n", A.row, A.col, A.nnz);

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
	aspparam.edgesTran = edgesTran;
	aspparam.nodeCEdge = nodeCEdge;
	aspparam.elementDOF = elementDOF;
	aspparam.elementdofTran = &elementdofTran;

	aspparam.lambda = Input->lambda;
	aspparam.mu = Input->mu;

	aspparam.max_levels = Input->AMG_levels;
	aspparam.AMG_coarsening_type = Input->AMG_coarsening_type;
	aspparam.AMG_interpolation_type = Input->AMG_interpolation_type;
	aspparam.AMG_coarse_dof = Input->AMG_coarse_dof;
	aspparam.AMG_strong_threshold = Input->AMG_strong_threshold;
	aspparam.AMG_truncation_threshold = Input->AMG_truncation_threshold;
	aspparam.AMG_max_row_sum = Input->AMG_max_row_sum;

	/* GMRES + Block triangular preconditioner */
	if (itsolver_type == 1)
	{
		printf("Block triangular preconditioned GMRES solver with auxiliary space method\n");
		TriAsP1ElasDG_GMRES(A, b, uh, &aspparam, print_level);
	}
	/* MINRES + Block diagonal preconditioner */
	else
	{
		//	aspparam.mu = 0.5;
		printf("Block diagonal preconditioned MINRES solver with auxiliary space method\n");
		DiagAsP1ElasDG_MINRES(A, b, uh, &aspparam, print_level);
	}

	/** Postprocessing */
	dvector Qhu, uhstar;
	create_dvector(uh[1].row, &Qhu);
	printf("Postprocessing...\n");
	projL2PkVec2d(&Qhu, elements, nodes, &elementDOF[1], linearElas2d_u, paras);
	postprocess2newDisplacement(&uhstar, &uh[0], &uh[1], elements, nodes, elementDOF, lambda, mu);
	printf("Postprocessing finish!\n");

	/* output solution to a diskfile */
	/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */

	/** Step 5. Compute the error between numerical solution and exact solution */
	double errors[10], posterrors[7];

	// geterrorsPoissonLagrange3d(errors, &uh, elements, elementFace, faces, elementEdge, edges, nodes, &elementDOF);
	geterrorslinearElasHuZhang2d(errors, &uh[0], &uh[1], &Qhu, &uhstar, elements, elementEdge, edges, nodes, elementDOF, lambda, mu);
	getposteriorierrorslinearElasHuZhang2d(posterrors, &uh[0], &uh[1], &uhstar, elements, elementEdge, edges, nodes, elementDOF, lambda, mu);

	printf("\nA Priori Errors:\n\n");
	printf("L2 norm of sigma-sigma_h = %e\n", errors[9]);
	printf("A norm of sigma-sigma_h = %e\n", errors[0]);
	printf("L2 norm of divergence of sigma-sigma_h = %e\n", errors[1]);
	printf("Jump error of u_h = %e\n", errors[2]);
	printf("L2 norm of u-u_h = %e\n", errors[3]);
	printf("Energy norm = %e\n", errors[4]);
	printf("H1 seminorm of Qhu-u_h = %e\n", errors[5]);
	printf("L2 norm of Qhu-u_h = %e\n", errors[6]);
	printf("H1 seminorm of u-uhstar = %e\n", errors[7]);
	printf("L2 norm of u-uhstar = %e\n", errors[8]);

	printf("\nA Posteriori Errors:\n\n");
	printf("rotrot(A sigma_h) = %e\n", posterrors[0]);
	printf("L2 norm of Asigmah - varepsilon_h(uhstar) = %e\n", posterrors[1]);
	printf("Oscillation: f-Qhf = %e\n", posterrors[2]);
	printf("Jump error of M_{tt}(A sigma_h) = %e\n", posterrors[3]);
	printf("Jump error of rot(Asigmah) t - dt(M_{nt}(A sigma_h)) = %e\n", posterrors[4]);
	printf("Eta = %e\n", posterrors[5]);
	printf("Eta + Oscillation = %e\n", posterrors[6]);
	
	/*********************************************************************************************/
	FILE *outputFile;
	outputFile = fopen("output/error.dat", "w");
	fprintf(outputFile, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", errors[9], errors[0], errors[1], errors[2], errors[3], errors[4], errors[5], errors[6], errors[7], errors[8]);
	fprintf(outputFile, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", posterrors[0], posterrors[1], posterrors[2], posterrors[3], posterrors[4], posterrors[5], posterrors[6]);
	fclose(outputFile);
	/********************************************************************************************/

	
	for (i = 0; i<3; i++)
		free_csr_matrix(&A[i]);
	if(A[3].row>0)
		free_csr_matrix(&A[3]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);
	free_dvector(&uh[0]);
	free_dvector(&uh[1]);
	free_dvector(&Qhu);
	free_dvector(&uhstar);
	free_elementDOF(elementDOF);
	free_elementDOF(elementDOF + 1);
	free_elementDOF(elementDOF + 2);

	free(elementdofTran.IA);
	free(elementdofTran.JA);
}

/**
* \fn void assemble_linearElasHuZhang2d(dCSRmat *ptr_A, dvector *ptr_b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu)
* \brief assemble stiffness matrix *C and righ hand side *ptr_b
* \param *ptr_A pointer to stiffness matrix
* \param *ptr_b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assemble_linearElasHuZhang2d(dCSRmat *ptr_A, dvector *ptr_b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu)
{
	/**	A[0]=A, A[1]=B^T, A{2]=-B
	Ax1 + B^Tx2 = 0
	-Bx1 +       = b
	where x1: sigma_h, x2: u_h
	**/
	double paras[2];
	paras[0] = lambda;
	paras[1] = mu;

	assembleweightedMassatrixHuZhang2d(&ptr_A[0], elements, elementDOF, elementdofTran, lambda, mu);
	assembleDivHuZhangL2poly2d(&ptr_A[1], elements, elementDOF, elementdofTran);
	assembleRHSHuZhang2d(&ptr_b[1], elements, nodes, elementDOF+1, elementdofTran+1, linearElas2d_f, paras);

	if (elementDOF[0].dop > 2)
	{
		ptr_A[3].row = 0;
		ptr_A[3].col = 0;
		ptr_A[3].IA = NULL;
		ptr_A[3].JA = NULL;
		ptr_A[3].val = NULL;
	}
	else
		assembleJumpL2poly2d(&ptr_A[3], elements, elementEdge, edges, nodes, elementDOF+1, -1.0);

	getTransposeOfSparse(&ptr_A[1], &ptr_A[2]);
	ax_array(ptr_b[1].row, -1.0, ptr_b[1].val);
	create_dvector(ptr_A[0].row, &ptr_b[0]);
}
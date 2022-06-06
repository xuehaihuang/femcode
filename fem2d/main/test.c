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
	dCSRmat A[4];
	dvector b[2], uh[2];
	ELEMENT *elements;
	idenmat *elementEdge;
	EDGE *edges;
	dennode *nodes;
	iCSRmat *edgesTran;
	ivector nodeCEdge;
	ELEMENT_DOF elementDOF[3];
	iCSRmat elementdofTran;



	/** Step 0. Read input parameters */
	char *inputfile = "ini/input.dat";
	Input_data Input;
	read_Input_data(inputfile, &Input);

	int print_level = Input.print_level;
	int problem_num = Input.problem_num;
	int itsolver_type = Input.itsolver_type;
	int itsolver_maxit = Input.itsolver_maxit;
	double itsolver_tol = Input.itsolver_tol;

	double lambda = Input.lambda;
	double mu = Input.mu;
	int glevelNum = Input.glevelNum;
	int CglevelNum = Input.CglevelNum;
	int FglevelNum = Input.FglevelNum;
	int dop1 = Input.dop1;

	int dop2=dop1-1;

	elements = (ELEMENT*)calloc(glevelNum, sizeof(ELEMENT));
	elementEdge = (idenmat*)calloc(glevelNum, sizeof(idenmat));
	edges = (EDGE*)calloc(glevelNum, sizeof(EDGE));
	nodes = (dennode*)calloc(glevelNum, sizeof(dennode));
	edgesTran = (iCSRmat*)calloc(glevelNum, sizeof(iCSRmat));


	/** Step 1. generate mesh */
	if (getmesh(elements, elementEdge, edges, nodes, edgesTran, &nodeCEdge, glevelNum) == 0)
	{
		printf("It's fail to generate mesh.\n");
		return 0;
	}

	int i;

	getElementEdgeGeoInfo(&elements[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1]);
//	for (i = 0; i < glevelNum; i++)
//		getElementEdgeGeoInfo(&elements[i], &edges[i], &nodes[i]);

//	getElementDOF(&elementDOFipdg0, elements[glevelNum - 1].row, dop1);

	getElementDOF_HuZhang(&elementDOF[0], &elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], nodes[glevelNum - 1].row, dop1);
	getElementDOF(&elementDOF[1], elements[glevelNum - 1].row, dop2);
	getElementDOF(&elementDOF[2], elements[glevelNum - 1].row, dop2+2);

	getTransposeOfelementDoF(elementDOF, &elementdofTran, 0);



	printf("h=%f\n", edges[glevelNum - 1].length[1]);
	printf("k=%d,  lambda=%f,  mu=%f\n", dop1, lambda, mu);


	/***************************************output edge******************************************************/
/*	printf("Edge: %d\n", edges.row);
	for(i=0;i<edges.row;i++)
		printf("%d %d          %d %d\n",edges.val[i][0]+1,edges.val[i][1]+1,edges.val[i][2]+1,edges.val[i][3]+1);
	printf("Node:%d %d\n", nodes.row,2);
	for(i=0;i<nodes.row;i++)
		printf("%f %f\n",nodes.val[i][0],nodes.val[i][1]);
	printf("Elements: %d %d\n",elements.row,3);
	for(i=0;i<elements.row;i++)
		printf("%d %d %d\n",elements.val[i][0]+1,elements.val[i][1]+1,elements.val[i][2]+1);	*/
		/********************************************************************************************/



		/***************************************output mesh******************************************************/
/*	FILE *meshFile;
	meshFile = fopen("output/mesh.dat", "w");
	fprintf(meshFile, "%d %d\n", nodes[glevelNum - 1].row, 2);
	for (i = 0; i < nodes[glevelNum - 1].row; i++)
		fprintf(meshFile, "%f %f\n", nodes[glevelNum - 1].val[i][0], nodes[glevelNum - 1].val[i][1]);
	fprintf(meshFile, "%d %d\n", elements[glevelNum - 1].row, 3);
	for (i = 0; i < elements[glevelNum - 1].row; i++)
		fprintf(meshFile, "%d %d %d\n", elements[glevelNum - 1].val[i][0] + 1, elements[glevelNum - 1].val[i][1] + 1, elements[glevelNum - 1].val[i][2] + 1);
	fclose(meshFile);*/
	/********************************************************************************************/



	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble(A, b, &elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], elementDOF, &elementdofTran, lambda, mu);

	/***************************write matrix to file**********************************************
	write_IJ_matrix(&A[0], "output/A11.dat");
	write_IJ_matrix(&A[1], "output/A12.dat");
	*********************************************************************************************/


	// initial solution
	create_dvector(b[0].row, &uh[0]);
	create_dvector(b[1].row, &uh[1]);

//	init_dvector(&uh[0], 1.0);
//	init_dvector(&uh[1], 1.0);
//	init_dvector(&b[0], 1.0);
//	init_dvector(&b[1], 1.0);


	// initial solution
//	create_dvector(elementDOF[glevelNum - 1].dof * 2, &uh);
//	create_dvector(elementDOFipdg0.dof * 2, &uh);

//	assembleStiffmatrixElasIPDG0(&A, &b, &elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], &elementDOFipdg0);

/*	int j;
	printf("A:\n");
	for (i = 0; i < Adg.row; i++)
	{
		for (j = Adg.IA[i]; j < Adg.IA[i + 1];j++)
			printf("%d, %f; ", Adg.JA[j], Adg.val[j]);
		printf("\n");
	}
	printf("b:\n");
	for (i = 0; i < bdg.row; i++)
		printf("%f; ", bdg.val[i]);
	printf("\n");
*/

	
	/** Step 3. Check matrix properties */
	for (i = 0; i<3; i++)
	{
		check_symm(&A[i]);
		check_diagpos(&A[i]);
		check_diagdom(&A[i]);
	}

//	return 0;///////////////////////////////////////

	/** Step 4. Solve the system */ 

	AMG_param amgparam; /* parameters for AMG */

	amgparam.print_level = Input.print_level;
	amgparam.restart = Input.restart;
	amgparam.max_iter = Input.itsolver_maxit;
	amgparam.tol = Input.itsolver_tol;
	amgparam.AMG_max_iter = Input.MG_maxit;
	amgparam.AMG_tol = Input.MG_tol;

	amgparam.max_levels = Input.AMG_levels;
	amgparam.smoother = Input.MG_smoother;
	amgparam.presmooth_iter = Input.MG_smooth_iter;
	amgparam.postsmooth_iter = Input.MG_smooth_iter;

	amgparam.coarsening_type = Input.AMG_coarsening_type;
	amgparam.interpolation_type = Input.AMG_interpolation_type;
	amgparam.coarse_dof = Input.AMG_coarse_dof;
	amgparam.strong_threshold = Input.AMG_strong_threshold;
	amgparam.truncation_threshold = Input.AMG_truncation_threshold;
	amgparam.max_row_sum = Input.AMG_max_row_sum;

//	cg_pcg(A, &b[1], &uh[0], &uh[1], itsolver_maxit, itsolver_tol, &amgparam, itsolver_type, print_level);


/*	if (print_level>0) {
		printf("Maximal iteration number = %d\n", Input.itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", Input.itsolver_tol);
	}

	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A.row, A.col, A.nnz);
	*/
	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A[0].row, A[0].col, A[0].nnz);
	printf("B.row=%d, B.col=%d, B.nnz=%d\n", A[1].row, A[1].col, A[1].nnz);
	if (A[3].row>0)
		printf("C.row=%d, C.col=%d, C.nnz=%d\n", A[3].row, A[3].col, A[3].nnz);

	ASP_param aspparam; /* parameters for AMG */

	aspparam.print_level = Input.print_level;
	aspparam.max_iter = Input.itsolver_maxit;
	aspparam.tol = Input.itsolver_tol;
	aspparam.restart = Input.restart;

	aspparam.precond_type = Input.precond_type;
	aspparam.precond_scale = Input.precond_scale;
	aspparam.smoother = Input.smoother;
	aspparam.smooth_iter = Input.smooth_iter;
	
	aspparam.levelNum = glevelNum;
	aspparam.mg_max_iter = Input.MG_maxit;
	aspparam.mg_tol = Input.MG_tol;
	aspparam.mg_smoother = Input.MG_smoother;
	aspparam.mg_smooth_iter = Input.MG_smooth_iter;

	aspparam.elements = &elements[glevelNum - 1];
	aspparam.elementEdge = &elementEdge[glevelNum - 1];
	aspparam.edges = &edges[glevelNum - 1];
	aspparam.nodes = &nodes[glevelNum - 1];
	aspparam.edgesTran = &edgesTran[glevelNum - 1];
	aspparam.nodeCEdge = &nodeCEdge;
	aspparam.elementDOF = elementDOF;
	aspparam.elementdofTran = &elementdofTran;

	aspparam.lambda = Input.lambda;
	aspparam.mu = Input.mu;

	aspparam.max_levels = Input.AMG_levels;
	aspparam.AMG_coarsening_type = Input.AMG_coarsening_type;
	aspparam.AMG_interpolation_type = Input.AMG_interpolation_type;
	aspparam.AMG_coarse_dof = Input.AMG_coarse_dof;
	aspparam.AMG_strong_threshold = Input.AMG_strong_threshold;
	aspparam.AMG_truncation_threshold = Input.AMG_truncation_threshold;
	aspparam.AMG_max_row_sum = Input.AMG_max_row_sum;


/*********************************************************
	dCSRmat Adg;
	assembleStiffmatrixElasIPDG0(&Adg, &elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], &elementDOF[1]);
	aspparam.elementDOF = &elementDOF[1]; 
	printf("Adg.row=%d, Adg.col=%d, Adg.nnz=%d\n", Adg.row, Adg.col, Adg.nnz);
	printf("Auxiliary space preconditioned CG solver\n");
	asP1ElasDG_PCG(&Adg, &b[1], &uh[1], &aspparam, print_level);
	return 0;
	/*************************************************************/
	if (itsolver_type == 1)
	{
//		aspparam.mu = 0.5;
		printf("Block diagonal preconditioned MINRES solver with auxiliary space method\n");
		DiagAsP1ElasDG_MINRES(A, b, uh, &aspparam, print_level);
	}
	else
	{
		printf("Block triangular preconditioned GMRES solver with auxiliary space method\n");
		TriAsP1ElasDG_GMRES(A, b, uh, &aspparam, print_level);
	}

	
	/* output solution to a diskfile */ 
/*	char *outputfile="output/sol.out";
	write_IJ_dvector(&uh, outputfile); */

	/** Postprocessing */
	dvector Qhu, uhstar;
	create_dvector(uh[1].row, &Qhu);
	printf("Postprocessing...\n");
	projPiecewiseLagrangeDisplacement(&Qhu, &elements[glevelNum - 1], &nodes[glevelNum - 1], &elementDOF[1], lambda, mu);
	postprocess2newDisplacement(&uhstar, &uh[0], &uh[1], &elements[glevelNum - 1], &nodes[glevelNum - 1], elementDOF, lambda, mu);
	printf("Postprocessing finish!\n");
	/** Step 5. Compute the error between numerical solution and exact solution */
	double errors[10], posterrors[7];

	geterrors(errors, &uh[0], &uh[1], &Qhu, &uhstar, &elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], elementDOF, lambda, mu);
	getposteriorierrors(posterrors, &uh[0], &uh[1], &uhstar, &elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], elementDOF, lambda, mu);
	

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

	//	for(i=0;i<elements.row*6;i++)
	//		printf("%e,",uh.val[i]);

	/*********************************************************************************************/
	FILE *outputFile;
	outputFile = fopen("output/error.dat", "w");
	fprintf(outputFile, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", errors[9], errors[0], errors[1], errors[2], errors[3], errors[4], errors[5], errors[6], errors[7], errors[8]);
	fprintf(outputFile, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", posterrors[0], posterrors[1], posterrors[2], posterrors[3], posterrors[4], posterrors[5], posterrors[6]);
	fclose(outputFile);
	/********************************************************************************************/

	for (i = 0; i < glevelNum; i++)
	{
		free_ELEMENT(&elements[i]);
		free_iden_matrix(&elementEdge[i]);
		free_EDGE(&edges[i]);
		free_dennode(&nodes[i]);
		free_icsr_matrix(&edgesTran[i]);
	}
	free_ivector(&nodeCEdge);
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

	free(elements);
	free(elementEdge);
	free(edges);
	free(nodes);
	free(edgesTran);
	
	return 1;
}

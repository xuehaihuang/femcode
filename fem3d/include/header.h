/*
 *  header.h
 *  MGV
 *
 *  Created by Xuehai Huang on 03/27/2009.
 *  Modified by Chensong Zhang on 03/29/2009.
 *
 *  Copyright 2008 PSU. All rights reserved.
 *
 */

/*! \file header.h
 *  \brief main header file for the MSC package
 */

/** 
 * \brief Definition of max, min, abs
 */
#define max(a,b) (((a)>(b))?(a):(b)) /**< bigger one in a and b */
#define min(a,b) (((a)<(b))?(a):(b)) /**< smaller one in a and b */
#define abs(a) (((a)>0.0)?(a):-(a)) /**< absolute value of a */

#define PI 3.1415926535897932384626
#define pi 3.1415926535897932384626
#define pii 3.1415926535897932384626

#define DBL_DIG         15                      /* # of decimal digits of precision */
#define DBL_EPSILON     2.2204460492503131e-016 /* smallest such that 1.0+DBL_EPSILON != 1.0 */
#define DBL_MANT_DIG    53                      /* # of bits in mantissa */
#define DBL_MAX         1.7976931348623158e+308 /* max value */
#define DBL_MAX_10_EXP  308                     /* max decimal exponent */
#define DBL_MAX_EXP     1024                    /* max binary exponent */
#define DBL_MIN         2.2250738585072014e-308 /* min positive value */
#define DBL_MIN_10_EXP  (-307)                  /* min decimal exponent */
#define DBL_MIN_EXP     (-1021)                 /* min binary exponent */
#define _DBL_RADIX      2                       /* exponent radix */
#define _DBL_ROUNDS     1                       /* addition rounding: near */

/** 
 * \brief Smoother type
 */
#define JACOBI 1  /**< Jacobi smoother */
#define GS     2  /**< Gauss-Seidel smoother */
#define SGS    3  /**< symm Gauss-Seidel smoother */
#define MSWZ   4  /**< multiplicative Schwarz smoother */
#define ASWZ   5  /**< additive Schwarz smoother */
#define SMSWZ   6  /**< symm multiplicative Schwarz smoother */

/**
 * \brief input parameters 
 *
 * Input parameters, reading from disk file
 */
typedef struct {
	
	// problem parameters	
	int problem_num;	/**< problem number to solve */
	
	int domain_num; /**< domain number */
					
	// parameters for iterative solvers
	int itsolver_type; /**< type of iterative solvers */
	double itsolver_tol; /**< tolerance for iterative linear solver */
	int itsolver_maxit; /**< maximal number of iterations for iterative solvers */
	int restart; /**< restart number used in GMRES */

	int precond_type; /**< type of precondition in ASP: 1 additive; 2 multiplicative  */
	double precond_scale[2]; /**< scale for the preconditioned variables in precondition */
	int smoother; /**< type of smoother in ASP */
	int schwarz_type; /**< type of Schwarz smoother */
	int smooth_iter; /**< number of smoothing in ASP */

	// parameters for MG
	double MG_tol; /**< tolerance for MG if used as preconditioner */
	int MG_maxit; /**< max number of iterations for MG if used as preconditioner */
	int MG_smoother; /**< type of smoother in MG */
	int MG_smooth_iter; /**< number of smoothing in MG */
//	int MG_presmooth_iter; /**< number of presmoothing */
//	int MG_postsmooth_iter; /**< number of postsmoothing */
	int AMG_levels; /**< maximal number of levels */
	int AMG_coarsening_type; /**< coarsening type */
	int AMG_interpolation_type; /**< interpolation type */
	int AMG_coarse_dof;	/**< minimal coarsest level dof */
	double AMG_strong_threshold; /**< strong threshold for coarsening */
	double AMG_truncation_threshold; /**< truncation factor for interpolation */
	double AMG_max_row_sum; /**< maximal row sum */

	// parameters for AFEM
	double AFEM_tol; /**< tolerance for AFEM */
	int AFEM_maxit; /**< max number of iterations for AFEM */
	double AFEM_mark_threshold; /**< mark threshold for AFEM */
	
	// output flags
	int print_level; /**< print level */	
	
	// temp
	double nu;
	double lambda;
	double mu;
	double t;
	double paraeps;
	double alpha1;
	double alpha2;
	double alpha3;
	double beta1;
	double beta2;
	double beta3;
	int glevelNum;
	int CglevelNum;
	int FglevelNum;
	int cDOP;
	int rDOP;
	int dop1; /**< degree of polynomial for stress tensor */
	int dop2; /**< degree of polynomial for displacement */
	int dop3; /**< degree of polynomial for trace of normal derivative of displacement */
	int dop4; /**< degree of polynomial for trace of displacement */

	int variationalform_type; /**< type of variational form */
	int stress_fem_type; /**< type of fem space for stress */
	int fem_num; /**< number of fem */
	short nitsche; /**< Use the Nitsche’s technique or not. 1: Yes, 0: No */
	
} Input_data;

/** 
 * \brief AMG_param: parameters for AMG solver.
 *
 * This is needed for the AMG solver.
 *
 * What coarsening type available? 
 * What are the parameters in detail?
 */
typedef struct {

	int print_level;								/**< print level for AMG */
	int max_levels;								/**< setup param, max number of levels */
	int coarse_dof;								/**< minimal coarsest level dof */
	int max_iter;									/**< solve params, max number of iterations */
	double tol;										/**< solve params, tolerance for solver */
	int AMG_max_iter;							/**< AMG params, max number of iterations */
	double AMG_tol;								/**< AMG params, tolerance for solver */
	int restart;                                /**< restart number used in GMRES */

	int smoother;									/**< smoother type */
	int presmooth_iter;						/**< number of presmoothers */
	int postsmooth_iter;						/**< number of postsmoothers */
	
	int coarsening_type;						/**< coarsening type */
	int interpolation_type;			  /**< interpolation type */	
	double strong_threshold;				/**< strong connection threshold for coarsening */
	double max_row_sum;						/**< maximal row sum parameter */
	double truncation_threshold;		/**< truncation threshold */

} AMG_param;

/** 
 * \brief Dense matrix of double type.
 *
 * A dense double matrix
 */ 
typedef struct ddenmat{
	/** number of rows */
	int row;	  
	/** number of columns */
	int col;	
	/** actual matrix entries */
	double **val;
}ddenmat;

/** 
 * \brief Dense matrix of int type.
 *
 * A dense int matrix
 */ 
typedef struct idenmat{
	/** number of rows */
	int row;	  
	/** number of columns */
	int col;	
	/** actual matrix entries */
	int **val;
}idenmat;

/** 
 * \brief ELEMENT struct.
 *
 * ELEMENT struct
 */ 
typedef struct ELEMENT{
	/** number of rows */
	int row;	  
	/** number of columns */
	int col;	
	/** the first 3 columns store the indexes of vertices,
	    the last 3 columns store the indexes of edge's midpoints */
	int **val;
	/** parent element of current element, which will be the root if its parent = -1*/
	int *parent;
	/** type of current element, which will be used for triangulation in 3d*/
	int *type;
	/** coordinate of vertices of the element */
	double ***vertices;
	/** volume of the element*/
	double *vol;
	/** barycenter of the element*/
	double **barycenter;
	/** barycenter of four faces of the element*/
	double ***bcFace;
	/** value of the barycentric coordinates lambda at the origin of coordinate (0,0,0) */
	double **lambdaConst;
	/** gradient of the barycentric coordinates lambda */
	double ***gradLambda;
	/** unit normal vector of four faces*/
	double ***nvector;
	/** permutation of four faces*/
	int ***fperm;
	/** inverse of the permutation fperm*/
	int ***fpermi;
	/** orientation of four faces*/
	short **forien;
	/** permutation of six edges*/
	int ***eperm;
	/** orientation of six edges*/
	short **eorien;
}ELEMENT;

/** 
 * \brief FACE struct.
 *
 * FACE struct
 */ 
typedef struct FACE{
	/** number of rows */
	int row;	  
	/** number of columns */
	int col;	
	/** the first three columns store the three vertice, the fourth and fifth columns store the affiliated elements, the fifth column stores -1 if the face is on boundary */
	int **val;
	/** unit normal vector of face*/
	double **nvector;
	/** unit tangential vector of face*/
	double **t1vector;
	double **t2vector;
	/** tangential vector of face*/
	double **t01;
	double **t02;
	double **t12;
	/** area of face*/
	double *area;
	/** barycenter of face*/
	double **barycenter;
	/** boundary type of face
	0 : non - boundary, i.e., an interior edge or face.
	1 : first type, i.e., a Dirichlet boundary edge or face.
	2 : second type, i.e., a Neumann boundary edge or face.
	3 : third type, i.e., a Robin boundary edge or face.
	**/
	int *bdFlag;
}FACE;

/** 
 * \brief EDGE struct.
 *
 * EDGE struct
 */ 
typedef struct EDGE{
	/** number of rows */
	int row;	  
	/** number of columns */
	int col;	
	/** the first two columns store the two vertice */
	int **val;
	/** unit normal vector of edge*/
	double **n1vector;
	double **n2vector;
	/** unit tangential vector of edge*/
	double **tvector;
	/** length of edge*/
	double *length;
	/** boundary type of edge
	0 : non - boundary, i.e., an interior edge or face.
	1 : first type, i.e., a Dirichlet boundary edge or face.
	2 : second type, i.e., a Neumann boundary edge or face.
	3 : third type, i.e., a Robin boundary edge or face.
	**/
	int *bdFlag;
}EDGE;

/** 
 * \brief Sparse matrix of double type in CSR format.
 *
 * CSR Format (IA,JA,A)
 *
 * The starting index of A is 0, other data stuctures also take this convention.  
 */
typedef struct dCSRmat{
	//! row number of matrix A, m
	int row;   
	//! column of matrix A, n
	int col;   
	//! number of nonzero entries
	int nnz;
	//! integer array of row pointers, the size is m+1
	int *IA;   
	//! integer array of column indexes, the size is nnz
	int *JA;    
	//! nonzero entries of A
	double *val;
}dCSRmat;

/** 
 * \brief Vector with n entries of double type.
 */
typedef struct dvector{
  //! number of rows
	int row;
  //! actual vector entries
	double *val;
}dvector;

/** 
 * \brief Vector with n entries of int type.
 */
typedef struct ivector{
  //! number of rows
	int row;
  //! actual vector entries
	int *val;
}ivector;

/** 
 * \brief Sparse matrix of int type in CSR format.
 *
 * CSR Format (IA,JA,A)
 *
 * The starting index of A is 0, other data stuctures also take this convention.  
 */
typedef struct iCSRmat{
	//! row number of matrix A, m
	int row;   
	//! column of matrix A, n
	int col;   
	//! number of nonzero entries
	int nnz;
	//! integer array of row pointers, the size is m+1
	int *IA;   
	//! integer array of column indexes, the size is nnz
	int *JA;    
	//! nonzero entries of A
	int *val;
}iCSRmat;

/** 
 * \brief Block diagonal matrix of double type.
 *
 * A dense double block diagonal matrix
 */ 
typedef struct dBDmat{
	/** number of rows */
	int row;	  
	/** number of columns */
	int col;	
	/** number of blocks */
	int nb;	
	/** blocks */
	ddenmat *blk;
}dBDmat;

/**
* \brief Overlapped block diagonal matrix of double type.
*
* A dense overlapped double block diagonal matrix
*/
typedef struct dOBDmat {
	/** number of rows */
	int row;
	/** number of columns */
	int col;
	/** number of blocks */
	int nb;
	/** blocks */
	ddenmat *blk;
	/** indices to local rows */
	ivector *rindices;
	/** indices to local columns */
	ivector *cindices;
}dOBDmat;

/** 
 * \brief Sparse matrix of double type in IJ format.
 *
 * Coordinate Format (I,J,A)
 *
 */
typedef struct dIJmat{
	//! row number of matrix A, m
	int row;   
	//! column of matrix A, n
	int col;   
	//! number of nonzero entries
	int nnz;
	//! integer array of row indices, the size is nnz
	int *I;   
	//! integer array of column indices, the size is nnz
	int *J;    
	//! nonzero entries of A
	double *val;
}dIJmat;

/** 
 * \brief Dennode type.
 */ 
typedef struct dennode{
	/** number of rows */
	int row;	  
	/** number of columns */
	int col;	
	/** actual matrix entries */
	double **val;
	/** boundary type of vetex
	0 : non - boundary, i.e., an interior vertex.
	1 : first type, i.e., a Dirichlet boundary vertex.
	2 : second type, i.e., a Neumann boundary vertex.
	3 : third type, i.e., a Robin boundary vertex.
	12: a Dirichlet-Neumann boundary vertex.
	22: a Neumann-Neumann boundary vertex.
	**/
	int *bdFlag;
}dennode;

/** 
 * \brief ELEMENT_DOF type.
 */ 
typedef struct ELEMENT_DOF{
	/** number of degree of polynomial */
	int dop;	  
	/** number of global degree of freedom */
	int dof;	  
	/** number of rows */
	int row;	  
	/** number of columns */
	int col;	
	/** actual matrix entries */
	int **val;
	/** non-free flag
	1 : non-free node
	0 : free node
	**/
	ivector nfFlag;
	/** free nodes */
	ivector freenodes;
	/** non-free nodes */
	ivector nfreenodes;
	/** index of free and non-free nodes */
	ivector index;
}ELEMENT_DOF;

/**
* \brief AS_param: parameters for auxiliary space preconditioner.
*
* This is needed for the AS preconditioner.
*
* What coarsening type available?
* What are the parameters in detail?
*/
typedef struct {

	int problem_num;
	int print_level;								/**< print level for ASP */
	int max_iter;									/**< solve params, max number of iterations */
	double tol;										/**< solve params, tolerance for solver */
	int restart;                                /**< restart number used in GMRES */

	int mass_precond_type; /**< type of preconditioning mass matrix: 1 diagonal; 2 full  */
	int precond_type; /**< type of precondition in ASP: 1 additive; 2 multiplicative  */
	double *precond_scale; /**< scale for the preconditioned variables in precondition  */
	int smoother;									/**< smoother type */
	int schwarz_type;                       /**< type of Schwarz smoother */
	int smooth_iter;						/**< number of smoothers */

	int levelNum;								/**< setup param, max number of levels */
	int mg_max_iter;							/**< MG params, max number of iterations */
	double mg_tol;								/**< MG params, tolerance for solver */
	int mg_smoother;									/**< smoother type */
	int mg_smooth_iter;						/**< number of smoothers */

											/************* Geometric information *************/
	ELEMENT *elements;
	idenmat *elementFace;
	FACE *faces;
	idenmat *elementEdge;
	EDGE *edges;
	dennode *nodes;
	iCSRmat *edgesTran;
	ivector *nodeCEdge;
	ELEMENT_DOF *elementDOF;
	iCSRmat *elementdofTran;

	/************* Material information *************/
	double lambda; /**< Lame constant */
	double mu; 
	double t;
	double nu;

	int variationalform_type; /**< type of variational form */
	int stress_fem_type; /**< type of fem space for stress */

	int max_levels;								/**< setup param, max number of levels */
	int AMG_levels; /**< maximal number of levels */
	int AMG_coarsening_type; /**< coarsening type */
	int AMG_interpolation_type; /**< interpolation type */
	int AMG_coarse_dof;	/**< minimal coarsest level dof */
	double AMG_strong_threshold; /**< strong threshold for coarsening */
	double AMG_truncation_threshold; /**< truncation factor for interpolation */
	double AMG_max_row_sum; /**< maximal row sum */
} ASP_param;

/**
* \brief AFEM_param: parameters for adaptive finite element method.
*
* This is needed for adaptive finite element method.
*
*/
typedef struct {

	AMG_param *amgparam;                             /**< AMG params */
	ASP_param *aspparam;                             /**< AMG params */
	double tol;										/**< tolerance for AFEM */
	int max_iter;							        /**< max number of iterations for AFEM */
	int solver_type;   					        /**< type of solver */
	double mark_threshold;			        	/**< marking threshold for AFEM */
} AFEM_param;

/** 
 * \brief Dense three-dimensional matrix of double type.
 *
 * A dense three-dimensional double matrix
 */ 
typedef struct ddenmat3{
	/** number of pages */
	int pag;	  
	/** number of rows */
	int row;	  
	/** number of columns */
	int col;	
	/** actual matrix entries */
	double ***val;
}ddenmat3;


/** 
 * \brief Links
 */
struct linked_list
{
  //! data
	int                        data;
	//! next element
	struct linked_list *next_elt;
	//! previous element
	struct linked_list *prev_elt;
	//! starting position
	int                        head;
	//! ending position
	int                        tail;
};
typedef struct linked_list ListElement;

/** 
 * \brief List of links
 */
typedef ListElement  *LinkList;  



/*--- Function declaration ---*/

/* input.c */
int read_Input_data(char *filenm, Input_data *Input);

/* coarsening.c */
void coarsening(dCSRmat *A, ivector *vertices, dCSRmat *P, AMG_param *param);
void generateS(dCSRmat *A, iCSRmat *S, AMG_param *param);
void generateSRS(dCSRmat *A, iCSRmat *S, double epsilon_str, int coarsening_type);
void getCoarseStiffMatrix(dCSRmat *PT, dCSRmat *A, dCSRmat *P, dCSRmat *Ac);

/* classicAMG.c */
void classicAMG(dCSRmat *A, dvector *b, dvector *x, AMG_param *param);

/* classicAMG_setup.c */
void classicAMG_setup(dCSRmat *A, dCSRmat *P, dCSRmat *PT, int *levels, AMG_param *param);

/* classicAMG_solve.c */
void classicAMG_solve(dCSRmat *A, dvector *b, dvector *x, dCSRmat *P, dCSRmat *PT, int levelNum, AMG_param *param);

/* interpolation.c */
void interpVecP1toNcP1_3d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFcr, FACE *faces);
void interpolation(dCSRmat *A, ivector *vertices, dCSRmat *P, AMG_param *param);
void interpolationRS(dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param);
int getiteval(dCSRmat *A, dCSRmat *it);

/* smoother.c */
void jacobi(dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, int L);
void gs(dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, int L);
void mulschwarz(dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, dOBDmat *B, int L);
void getSchwarzblocks_vertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, int nvertices, ELEMENT_DOF *elementDOF);
void getSchwarzblocks_edge(dOBDmat *B, dCSRmat *A, EDGE *edges, ELEMENT_DOF *elementDOF);
void getSchwarzblocks_edgevertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, EDGE *edges, int nvertices, ELEMENT_DOF *elementDOF);
void getSchwarzblocks_element(dOBDmat *B, dCSRmat *A, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF);
void getSchwarzblocks_elementvertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, int nvertices, ELEMENT_DOF *elementDOF);
void getSchwarzblocksVec2_vertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, int nvertices, ELEMENT_DOF *elementDOF);
void getSchwarzblocksVec2_edge(dOBDmat *B, dCSRmat *A, EDGE *edges, ELEMENT_DOF *elementDOF);
void getSchwarzblocksVec2_edgevertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, EDGE *edges, int nvertices, ELEMENT_DOF *elementDOF);
void getSchwarzblocksVec2_element(dOBDmat *B, dCSRmat *A, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF);
void getSchwarzblocksVec2_elementvertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, int nvertices, ELEMENT_DOF *elementDOF);

/* multigrid.c */
void multigrid(dCSRmat *A, dvector *b, dvector *x, dCSRmat *R, dCSRmat *P, int level, int levelNum, int smoother, int m1, int m2, int m0);

/* itsolver.c */
int standard_CG(dCSRmat *A, dvector *b, dvector *x, int MaxIt, double tol, int print_level);
int diag_PCG(dCSRmat *A, dvector *b, dvector *x, int MaxIt, double tol, int print_level);
int classicAMG_PCG(dCSRmat *A, dvector *b, dvector *x, AMG_param *param, int print_level);
int classicAMG_GMRES(dCSRmat *A, dvector *b, dvector *x, AMG_param *param, int print_level);
int DiagAsP1StokesNcP1P0_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int AbfpAsP1StokesNcP1P0_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int cg_pcg(dCSRmat *A, dvector *b, dvector *sigma, dvector *u, int MaxIt, double tol, AMG_param *param, int itsolver_type, int print_level);

/* linklist.c */
void dispose_elt ( LinkList element_ptr );
void remove_point ( LinkList *LoL_head_ptr , LinkList *LoL_tail_ptr , int measure , int index , int *lists , int *where );
LinkList create_elt ( int Item );
void enter_on_lists ( LinkList *LoL_head_ptr , LinkList *LoL_tail_ptr , int measure , int index , int *lists , int *where );

/* sparse.c */
int get_block(dCSRmat *A, int m, int n, int *rows, int *cols, double *Aloc, int *mask);
int get_block_dden(dCSRmat *A, int m, int n, int *rows, int *cols, ddenmat *Aloc, int *mask);
iCSRmat getCSRfromIntT(idenmat *T, int nNode);
iCSRmat getTransposeOfA(iCSRmat *A);
void getTransposeOfSparse(dCSRmat *A,dCSRmat *AT);
void getTransposeOfiden(idenmat *A, iCSRmat *AT, int cols, int nNode);
void getTransposeOfELEMENT(ELEMENT *A, iCSRmat *AT, int cols, int nNode);
void getTransposeOfelementDoF(ELEMENT_DOF *A, iCSRmat *AT, int cols);
int sparseAddition(dCSRmat *A, dCSRmat *B, dCSRmat *C);
int sparseSubtraction(dCSRmat *A, dCSRmat *B, dCSRmat *C);
int dCSRPlusdDiagVector(dCSRmat *A, dvector *B, dCSRmat *C);
int dCSRPlusdBD(dCSRmat *A, dBDmat *B, dCSRmat *C);
int dBDMinusdCSR(dBDmat *A, dCSRmat *B, dCSRmat *C);
void sparseMultiplication(dCSRmat *A, dCSRmat *B, dCSRmat *C);
void dBDMultiplydCSR(dBDmat *A, dCSRmat *B, dCSRmat *C);
void dBD2MultiplydCSR(dBDmat *A, int srow, int scol, dCSRmat *B, dCSRmat *C);
void dCSRMultiplydBD(dCSRmat *A, dBDmat *B, dCSRmat *C);
int dCSRMultiplydDiagVector(dCSRmat *A, dvector *D, dCSRmat *C);
int dCSRMultiplydDiagVectorInv(dCSRmat *A, dvector *D, dCSRmat *C);
int dDiagVectorMultiplydCSR(dvector *D, dCSRmat *A, dCSRmat *C);
int dDiagVectorInvMultiplydCSR(dvector *D, dCSRmat *A, dCSRmat *C);
void sparseTripleMultiplication(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B);
void sparseTripleMultiplication1(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B);
void sparseTripleMultiplication2(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B);
void sparseTripleMultiplication3(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B);
int dIJtoCSR(dIJmat *A, dCSRmat *B);
int dCSRtoIJ(dCSRmat *A, dIJmat *B);

/* poisson.c */
double poisson3d_f(double *x);
double poisson3d_u(double *x);
void poisson3d_gradu(double *x, double *val);

/* poissonfem3d.c */
void poissonfem3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void poissonLagrange3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void assemble_poissonLagrange3d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);

/* poissonerror3d.c */
void geterrorsPoissonLagrange3d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);

/* maxwell.c */
void maxwell3d_f(double *x, double *val);
void maxwell3d_u(double *x, double *val);
void maxwell3d_curlu(double *x, double *val);
void maxwell3d_gradcurlu(double *x, double *val);
void maxwell3d_gradu(double *x, double *val);
void maxwell3d_gradu1(double *x, double *val);
void maxwell3d_gradu2(double *x, double *val);
void maxwell3d_gradu3(double *x, double *val);
double maxwell3d_u1_x(double x, double y, double z);
double maxwell3d_u1_y(double x, double y, double z);
double maxwell3d_u1_z(double x, double y, double z);
double maxwell3d_u2_x(double x, double y, double z);
double maxwell3d_u2_y(double x, double y, double z);
double maxwell3d_u2_z(double x, double y, double z);
double maxwell3d_u3_x(double x, double y, double z);
double maxwell3d_u3_y(double x, double y, double z);
double maxwell3d_u3_z(double x, double y, double z);

/* maxwellfem.c */
void maxwellfem(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void maxwellNedelec1st3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void maxwellNedelec2nd3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void maxwellHuangZhang3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void maxwellHuang3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void assemble_maxwellNedelec1st3d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assemble_maxwellNedelec2nd3d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assemble_maxwellHuangZhang3d(dCSRmat *A, dvector *b, dvector *uh, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assemble_maxwellHuang3d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);

/* maxwellerror.c */
void geterrorsMaxwellNedelec1st3d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void geterrorsMaxwellNedelec2nd3d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void geterrorsMaxwellCHHcurlHermite3d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void geterrorsMaxwellHuangZhang3d(double *errors, dvector *uh, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void geterrorsMaxwellHuang3d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);

/* quadcurl.c */
void quadcurl3d_f(double *x, double *val);
void quadcurl3d_u(double *x, double *val);
void quadcurl3d_curlu(double *x, double *val);
void quadcurl3d_gradcurlu(double *x, double *val);

/* quadcurlfem.c */
void quadcurlfem(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void quadcurlHuang3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void quadcurlHuang3d_msm(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void quadcurlHuangZhang3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void assemble_quadcurlHuang3d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assemble_quadcurl_StokesNcP1P03d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, dvector *wh, ELEMENT_DOF *elementDOFwh);
void assemble_quadcurlHuangZhang3d(dCSRmat *A, dvector *b, dvector *uh, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assembleRHSCurlHuangNcP13d4Stokes(dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, dvector *wh, ELEMENT_DOF *elementDOFwh);
void assembleRHSCurlHuangNcP13d4Maxwell(dvector *rhs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, dvector *uhs, ELEMENT_DOF *elementDOFs);

/* quadcurlerror.c */
void geterrorsQuadcurlHuang3d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void geterrorsQuadcurlHuangZhang3d(double *errors, dvector *uh, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);

/* quadcurlperturbfem.c */
void quadcurlperturbfem(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void quadcurlperturbHuangZhang3d(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void assemble_quadcurlperturbHuangZhang3d(dCSRmat *A, dvector *b, dvector *uh, double paraeps, short nitsche, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);

/* quadcurlperturberror.c */
void geterrorsQuadcurlperturbHuangZhang3d(double *errors, double paraeps, short nitsche, dvector *uh, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);

/* assemble.c */
int getmesh(int domain_num, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, int levelNum);
void getElementDOF1d(ELEMENT_DOF *elementDOF, int ne, int dop);
void getElementDOF1d_Continue(ELEMENT_DOF *edgeDOF, ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop);
void getElementDOF(ELEMENT_DOF *elementDOF, int nt, int dop);
void getElementDOF_Lagrange(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop);
void getElementDOF3d(ELEMENT_DOF *elementDOF, int nt, int dop);
void getElementDOF_Lagrange3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, int nvertices, int dop);
void getElementDOF_NoncfmP13d(ELEMENT_DOF *elementDOF, idenmat *elementFace, int nf);
void getElementDOF_Morley3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges);
void getElementDOF_HuZhang(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop);
void getElementDOF_HuZhang3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, int nvertices, int dop);
void getElementDOF_Nedelec1st3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, int dop);
void getElementDOF_Nedelec2nd3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, int dop);
void getElementDOF_CHHcurlHermite3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, int nvertices, int dop);
void getElementDOF_HuangGradcurl3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges);
void getElementDOF_HuangZhang3d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges);
void assembleMassmatrixLagrange3d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assembleMassmatrixP03d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF);
void assembleBiGradLagrange3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assembleRHSLagrange3d(dvector *b, ELEMENT *elements, ELEMENT_DOF *elementDOF, double (*f)(double *));
void assembleBiGradVectorNcP13d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assembleDivNcP1P03d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assembleBiCurlNedelec1st3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assembleNedelec1stGradLagrange3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assembleRHSNedelec1st3d(dvector *b, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, void (*f)(double *, double *));
void assembleBiCurlNedelec2nd3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assembleNedelec2ndGradLagrange3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assembleRHSNedelec2nd3d(dvector *b, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, void (*f)(double *, double *));
void assembleBiCurlHuangZhang3d(dCSRmat *A, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assembleBiGradcurlHuangZhang3d(dCSRmat *A, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assembleBiGradcurlperturbHuangZhang3d(dCSRmat *A, double paraeps, short nitsche, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assembleHuangZhangGradLagrange3d(dCSRmat *A, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assembleRHSHuangZhang3d(dvector *b, ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, void (*f)(double *, double *));
void assembleBiCurlHuang3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assembleBiGradcurlHuang3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assembleHuangGradLagrange3d(dCSRmat *A, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assembleRHSHuang3d(dvector *b, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, void (*f)(double *, double *));
void jumpOperatorMorley3d(double *prt_lambdas, int face, ELEMENT *elements, idenmat *elementFace, FACE *faces, ELEMENT_DOF *elementDOF, int node, double *jump);
void averageOperatorNormalDerivativeMorley3d(double *prt_lambdas, int face, ELEMENT *elements, idenmat *elementFace, FACE *faces, ELEMENT_DOF *elementDOF, int node, double *ave);
void jumpOperatorVectorTensor3d(double *prt_lambdas, int face, ELEMENT *elements, idenmat *elementFace, FACE *faces, ELEMENT_DOF *elementDOF, int node, double *jump);
void jumpOperatorVector3d(double *prt_lambdas, int face, ELEMENT *elements, idenmat *elementFace, FACE *faces, ELEMENT_DOF *elementDOF, int node, double *jump);
void jumpOperatorVector3dOld(double *prt_lambdas, int face, ELEMENT *elements, idenmat *elementFace, FACE *faces, ELEMENT_DOF *elementDOF, int node, double *jump);
void getElementFaceEdgeGeoInfo(ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes);
void getFreenodesInfoLagrange3d(FACE *faces, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void getFreenodesInfoMorley3d(FACE *faces, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void getFreenodesInfoNoncfmP1Vector3d(FACE *faces, ELEMENT_DOF *elementDOF);
void getFreenodesInfoNedelec1st3d(FACE *faces, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void getFreenodesInfoNedelec2nd3d(FACE *faces, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void getFreenodesInfoCHHcurlHermite3d(FACE *faces, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void getFreenodesInfoHuangGradcurl3d(FACE *faces, EDGE *edges, dennode *nodes, int type, ELEMENT_DOF *elementDOF);
void getFreenodesInfoHuangZhang3d(FACE *faces, EDGE *edges, dennode *nodes, int type, ELEMENT_DOF *elementDOF);
void getFaceEdgeNTtensor(double **fnv, double **ft1v, double **ft2v, double **etv, double **en1v, double **en2v, double (*fnn)[6], double (*fnt1)[6], double (*fnt2)[6], double (*ft1t1)[6], double (*ft2t2)[6], double (*ft1t2)[6], double (*ett)[6], double (*etn1)[6], double (*etn2)[6], double (*en1n1)[6], double (*en2n2)[6], double (*en1n2)[6]);
void getBoundaryInfoNormalTrace(EDGE *edges, ELEMENT_DOF *elementDOF, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index);
void getBoundaryInfoEdge(EDGE *edges, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index);


/* assembleInfo.c */
void getFaceInfo(ELEMENT *elements, iCSRmat *elementsTran, FACE *faces, dennode *nodes);
void getFaceBdflag(FACE *faces, dennode *nodes);
void getEdgeInfo(idenmat *elements, iCSRmat *elementsTran, idenmat *edges, iCSRmat *edgesTran);
int getCoarseInfo(int domain_num, dennode *nodes, ELEMENT *elements, FACE *faces, EDGE *edges, iCSRmat *elementsTran, idenmat *elementEdge);
void uniformrefine(dennode *Cnodes, ELEMENT *Celements, FACE *Cfaces, EDGE *Cedges, iCSRmat *CelementsTran, dennode *Fnodes, ELEMENT *Felements, FACE *Ffaces, EDGE *Fedges, iCSRmat *FelementsTran, idenmat *FelementEdge);
void extractNondirichletMatrix(dCSRmat *A, dCSRmat *A1, ivector *dirichlet, ivector *nondirichlet);
void extractNondirichletVector(dvector *b, dvector *b1, ivector *dirichlet, ivector *nondirichlet);
void extractFreenodesVector2StokesDirichlet(dCSRmat *A, dCSRmat *B, dvector *b, dvector *b1, ELEMENT_DOF *elementDOF, dvector *uh);
void extractFreenodesVector(dCSRmat *A, dvector *b, dvector *b1, ELEMENT_DOF *elementDOF, dvector *uh);
void extractFreenodesMatrix11(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc);
void extractFreenodesMatrix1r(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF);
void extractFreenodesMatrix1c(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF);
void extractFreenodesMatrix1cBlock(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF);
void extractFreenodesMatrix11cBlock(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc);
void extractFreenodesMatrix11cBlock3d(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc);
int getEdgeDOFsTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row31, int *row32);
int getEdgeDOFsScalarTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int rowstart, int row31, int row32);
int getEdgeDOFsVectorTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row31, int *row32);
int getEdgeDOFsVector(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row21);
int getFaceDOFsVector(dCSRmat *A, int count, int element, int face, idenmat *elementFace, ELEMENT_DOF *elementDOF, int *rowstart, int *row31);
int getFaceDOFsVectorOld(dCSRmat *A, int count, int element, int face, idenmat *elementFace, ELEMENT_DOF *elementDOF, int *rowstart, int *row31);
int getFaceDOFsVectorlocalFace(dCSRmat *A, int count, int element, int lface, ELEMENT_DOF *elementDOF, int *rowstart, int *row31);
int getFaceDOFsVectorlocalFaceOld(dCSRmat *A, int count, int element, int lface, ELEMENT_DOF *elementDOF, int *rowstart, int *row31);
int getFaceDOFs(int *A, int count, int element, int face, idenmat *elementFace, ELEMENT_DOF *elementDOF);
void patchtogether22(dCSRmat *A11, dCSRmat *A12, dCSRmat *A21, dCSRmat *A22, dCSRmat *A);
void patchtogether33(dCSRmat *A11, dCSRmat *A12, dCSRmat *A13, dCSRmat *A21, dCSRmat *A22, dCSRmat *A23, dCSRmat *A31, dCSRmat *A32, dCSRmat *A33, dCSRmat *A);
void getPenaltyParameters(ddenmat *etas, idenmat *elementEdge, EDGE *edges, double *alpha, double *beta);
void getEdgeInfo3(ELEMENT *elements, iCSRmat *elementsTran, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes);
void getEdgeBdflag3(EDGE *edges, ELEMENT *elements, FACE *faces, idenmat *elementEdge, dennode *nodes);
void getPermutation(int *a, int *b, int *perm, int n);
int vv2edge3d(int v1, int v2);
void edge2vv3d(int edge, int *v);
void face2vertices3d(int face, int *v);

/* basiscoeff.c */
void getHuangZhangBasisCoeffs(ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges);

/* basicData.c */
void morley_basis(double lambda1, double lambda2, double lambda3, double s, double elen[3], double eta[3], double xi[3], double sij[3], double orient[3], int index, double *phi);
void morley_basis1(double lambda1, double lambda2, double lambda3, double s, double elen[3], double eta[3], double xi[3], double sij[3], double orient[3], int index, double phi[2]);
void morley_basis2(double lambda1, double lambda2, double lambda3, double s, double elen[3], double eta[3], double xi[3], double sij[3], double orient[3], int index, double phi[3]);
void lagrange1D_basis(double lambda, int index, int dop, double *phi);
void lagrange1D_basis1(double lambda, int index, int dop, double h, double *phi);
void lagrange_basis(double *lambda, int index, int dop, double *phi);
void lagrange_basis1(double *lambda, double s, double eta[3], double xi[3], int index, int dop, double phi[2]);
void lagrange_basis2(double *lambda, double s, double eta[3], double xi[3], int index, int dop, double phi[3]);
void lagrange3d_basis(double *lambda, int index, int dop, double *phi);
void lagrange3d_basis1(double *lambda, double **grd_lambda, int index, int dop, double phi[3]);
void ncp13d_basis(double *lambda, int index, double *phi);
void ncp13d_basis1(double **grd_lambda, int index, double *phi);
void morley3d_basis(double *lambda, double v, double **grd_lambda, double **nv, double **nvf, int index, double *phi);
void morley3d_basis1(double *lambda, double v, double **grd_lambda, double **nv, double **nvf, int index, double phi[3]);
void morley3d_basis2(double *lambda, double v, double **grd_lambda, double **nv, double **nvf, int index, double phi[6]);
void rt_basis(double x, double y, double (*T)[2], double s, double elen[3], double eta[3], double xi[3], double orient[3], int index, int dop, double phi[2]);
void rt_basis1(double x, double y, double (*T)[2], double s, double elen[3], double eta[3], double xi[3], double orient[3], int index, int dop, double phi[4]);
void arnoldwinther_basis(double *lambda, double *x, double *y, ddenmat3 *basisCoeffs, int element, int index, double *phi);
void arnoldwinther_basisDIV(double *lambda, ddenmat3 *basisCoeffs, int element, double s, double eta[3], double xi[3], int index, double *phi);
void arnoldwinther_basisDIV2(ddenmat3 *basisCoeffs, int element, double s, double eta[3], double xi[3], int index, double *phi);
void huzhang_basis(double *lambda, double **nv, double **tv, int index, int dop, double phi[3]);
void huzhang_basisDIV(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, int index, int dop, double phi[2]);
void huzhang3d_basis(double *lambda, double (*fnn)[6], double (*fnt1)[6], double (*fnt2)[6], double (*ft1t1)[6], double (*ft2t2)[6], double (*ft1t2)[6], double (*ett)[6], double (*etn1)[6], double (*etn2)[6], double (*en1n1)[6], double (*en2n2)[6], double (*en1n2)[6], int index, int dop, double phi[6]);
void huzhang3d_basisDIV(double *lambda, double v, double **grd_lambda, double (*fnn)[6], double (*fnt1)[6], double (*fnt2)[6], double (*ft1t1)[6], double (*ft2t2)[6], double (*ft1t2)[6], double (*ett)[6], double (*etn1)[6], double (*etn2)[6], double (*en1n1)[6], double (*en2n2)[6], double (*en1n2)[6], int index, int dop, double phi[3]);
void nedelec1st3d_basis(double *lambda, double **grd_lambda, short *eorien, int **fpermi, int index, int dop, double phi[3]);
void nedelec1st3d_basisCurl(double *lambda, double **grd_lambda, short *eorien, int **fpermi, int index, int dop, double phi[3]);
void nedelec2nd3d_basis(double *lambda, double **grd_lambda, int **eperm, int index, int dop, double phi[3]);
void nedelec2nd3d_basisCurl(double *lambda, double **grd_lambda, int **eperm, int index, int dop, double phi[3]);
void chhcurlHermite3d_basis(double *lambda, double **grd_lambda, double **etv, int index, int dop, double phi[3]);
void chhcurlHermite3d_basisCurl(double *lambda, double **grd_lambda, double **etv, int index, int dop, double phi[3]);
void huangQuadcurl3d_facebasis(double *x, double *xK, double *lambda, double **grd_lambda, double **vertices, double *nvf[4], int l, int li, double phi[3]);
void huangQuadcurl3d_facebasisCurl(double *lambda, double **grd_lambda, double *nvf[4], int l, int li, double phi[3]);
void huangQuadcurl3d_facebasisGradCurl(double **grd_lambda, double *nvf[4], int l, int li, double phi[9]);
void huangQuadcurl3d_basis(double *x, double *xK, double *lambda, double **grd_lambda, double **vertices, double *nvf[4], short *eorien, int **fpermi, int index, double phi[3]);
void huangQuadcurl3d_basisCurl(double *lambda, double **grd_lambda, double *nvf[4], short *eorien, int **fpermi, int index, double phi[3]);
void huangQuadcurl3d_basisGradCurl(double **grd_lambda, double *nvf[4], short *eorien, int **fpermi, int index, double phi[9]);
void huangzhang03d_basis(double *lambda, double **grd_lambda, int index, double phi[3]);
void huangzhang03d_basisCurl(double *lambda, double **grd_lambda, int index, double phi[3]);
void huangzhang03d_basisGradCurl(double *lambda, double **grd_lambda, int index, double phi[9]);
double volume(double (*tet)[3]);

/* quadrature.c */
int getNumQuadPoints(int dop, int dim);
int getNumQuadPoints_ShunnWilliams(int dop, int dim);
void init_quadrature(int num_qp, int ncoor, double (*gauss)[3]);
void init_Gauss(int num_qp, int ncoor, double (*gauss)[3]);
void init_Gauss1D(int num_qp, int ncoor, double (*gauss)[2]);
void init_NewtonCotes1D(int num_qp, int ncoor, double (*newtoncotes)[2]);
void init_ShunnWilliams2d(int num_qp, double (*lambdas)[3], double *weight);
void init_ShunnWilliams3d(int num_qp, double (*lambdas)[4], double *weight);

/* lu.c */
int LU_Decomp(double *A, int pivot[], int n);
int LU_Solve(double *A, double b[], int pivot[], double x[], int n);

/* output.c */
int writeElementsNodes(ELEMENT *elements, dennode *nodes, char *fname1, char *fname2);
int write_IJ_dCSRmat(dCSRmat *A, char *fname);
void write_IJ_dCSRmat4Matlab(dCSRmat *A, char *fname);
void write_dvector4Matlab(dvector *vec, char *fname);
void read_dvector4Matlab(dvector *vec, char *fname);

/* xuludmil.for (Xiangtan energy minimization code) */
//void get_p_xuludmil_(int *ia,int *ja,double *a, int *n, int *nc,int *ip,int *jp,double *pn,int *ipt,int *jpt,int *mf);

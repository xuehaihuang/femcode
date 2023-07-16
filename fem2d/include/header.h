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
	double parapenalty;
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
	short extrap; /**< type of extrapolation. 0: No extrapolation*/

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
	/** coordinate of vertices of the element */
	double ***vertices;
	/** volume of the element*/
	double *vol;
	/** barycenter of the element*/
	double **barycenter;
	double **eta; 
	double **xi;
	/** gradient of the barycentric coordinates lambda */
	double ***gradLambda;
	/** unit normal vector of three edges*/
	double ***nvector;
	/** unit tangential vector of three edges*/
	double ***tvector;
	/** tangential vector of three edges*/
	double ***tij;
	/** length of three edges*/
	double **edgeslength;
	/** height to three edges*/
	double **height;
	/** permutation of six edges*/
	int ***eperm;
	/** orientation of six edges*/
	short **eorien;
}ELEMENT;

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
	/** the first two columns store the two vertice, the third and fourth columns store the affiliated elements，
        the fourth column stores -1 if the edge is on boundary */
	int **val;
	/** sizes of the matrix entries*/
	double *h;
	/** x1-x2*/
	// double *xi;
	/** y1-y2*/
	// double *eta;
	/** unit normal vector of edge*/
	double **nvector;
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

	int print_level;								/**< print level for ASP */
	int max_iter;									/**< solve params, max number of iterations */
	double tol;										/**< solve params, tolerance for solver */
	int restart;                                /**< restart number used in GMRES */

	int precond_type; /**< type of precondition in ASP: 1 additive; 2 multiplicative  */
	double *precond_scale; /**< scale for the preconditioned variables in precondition  */
	int smoother;									/**< smoother type */
	int smooth_iter;						/**< number of smoothers */

	int levelNum;								/**< setup param, max number of levels */
	int mg_max_iter;							/**< MG params, max number of iterations */
	double mg_tol;								/**< MG params, tolerance for solver */
	int mg_smoother;									/**< smoother type */
	int mg_smooth_iter;						/**< number of smoothers */

											/************* Geometric information *************/
	ELEMENT *elements;
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
void interpP1toDG2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFdg);
void interpP1toP2_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFp2, EDGE *edges);
void interpVecP1toDG2d(dCSRmat *P, int n, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFdg);
void interpVecP1toNcP1_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFcr, EDGE *edges);
void interpVecP1toMINI_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFmini);
void interpStensorP1toMINI_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFmini);
void interpVecP1toCR_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFcr, EDGE *edges);
void interpVecP2toCR_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp2, ELEMENT_DOF *elementDOFcr);
void interpolation(dCSRmat *A, ivector *vertices, dCSRmat *P, AMG_param *param);
void interpolationRS(dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param);
int getiteval(dCSRmat *A, dCSRmat *it);

/* smoother.c */
void jacobi(dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, int L);
void gs(dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, int L);

/* multigrid.c */
void multigrid(dCSRmat *A, dvector *b, dvector *x, dCSRmat *R, dCSRmat *P, int level, int levelNum, int smoother, int m1, int m2, int m0);
double mgvVectorP1_solve(dCSRmat *A, dvector *b, dvector *x, EDGE *edges, iCSRmat *edgesTran, ivector *nodeCEdge, ivector *isInNode, ivector *nondirichlet, ivector *index, int levelNum, int smoother, int m);
double mgvP1_solve(dCSRmat *A, dvector *b, dvector *x, EDGE *edges, iCSRmat *edgesTran, ivector *nodeCEdge, ivector *isInNode, ivector *nondirichlet, ivector *index, int levelNum, int smoother, int m);
void interpolationPvector2d(EDGE *Cedges, ivector *CisInNode, ivector *Cindex, ivector *Fnondirichlet, ivector *nodeCEdge, int cnn, dvector *e, dvector *Pe);
void interpolationPvector2d2023(EDGE *Cedges, ivector *nodeCEdge, int cnn, dvector *e, dvector *Pe);
void interpolationP2d(EDGE *Cedges, ivector *CisInNode, ivector *Cindex, ivector *Fnondirichlet, ivector *nodeCEdge, int cnn, dvector *e, dvector *Pe);
void restrictionPTvector2d(iCSRmat *FedgesTran, ivector *Findex, ivector *Cnondirichlet, dvector *r, dvector *PTr);
void restrictionPTvector2d2023(iCSRmat *FedgesTran, dvector *r, dvector *PTr);
void restrictionPT2d(iCSRmat *FedgesTran, ivector *Findex, ivector *Cnondirichlet, dvector *r, dvector *PTr);
double mgvVectorP1_solveOld(dCSRmat *A, dvector *b, dCSRmat *A11, dvector *b1, dvector *x, EDGE *edges, iCSRmat *edgesTran, ivector *nodeCEdge, ivector *nondirichlet, ivector *index, int levelNum, int m);
void interpolationPvector2dOld(EDGE *Cedges, ivector *nodeCEdge, dvector *e, dvector *Pe);
void restrictionPTvector2dOld(iCSRmat *edgesTran, dvector *r, dvector *PTr);

/* itsolver.c */
int standard_CG(dCSRmat *A, dvector *b, dvector *x, int MaxIt, double tol, int print_level);
int diag_PCG(dCSRmat *A, dvector *b, dvector *x, int MaxIt, double tol, int print_level);
int classicAMG_PCG(dCSRmat *A, dvector *b, dvector *x, AMG_param *param, int print_level);
int DiagAsP1PoissonDG_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int TriAsP1PoissonDG_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int DiagAsP1StokesNcP1P0_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int AbfpAsP1StokesNcP1P0_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int DiagAsP1symStokesMINI_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int AbfpAsP1symStokesMINI_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int asP1ElasDG_PCG(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int DiagAsP1ElasDG_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level, void(*assembleMassmatrix)(dCSRmat *, ELEMENT *, ELEMENT_DOF *, double, double));
int TriAsP1ElasDG_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level, void(*assembleMassmatrix)(dCSRmat *, ELEMENT *, ELEMENT_DOF *, double, double));
int classicAMG_GMRES(dCSRmat *A, dvector *b, dvector *x, AMG_param *param, int print_level);
int cg_pcg(dCSRmat *A, dvector *b, dvector *sigma, dvector *u, int MaxIt, double tol, AMG_param *param, int itsolver_type, int print_level);

/* asElasDG.c */
int mgvVectorP1AsElasDG(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
double mgvVectorP1As_solve(dCSRmat *A, dvector *b, dvector *x, void *data);

/* linklist.c */
void dispose_elt ( LinkList element_ptr );
void remove_point ( LinkList *LoL_head_ptr , LinkList *LoL_tail_ptr , int measure , int index , int *lists , int *where );
LinkList create_elt ( int Item );
void enter_on_lists ( LinkList *LoL_head_ptr , LinkList *LoL_tail_ptr , int measure , int index , int *lists , int *where );

/* sparse.c */
int get_block(dCSRmat *A, int m, int n, int *rows, int *cols, double *Aloc, int *mask);
iCSRmat getCSRfromIntT(idenmat *T, int nNode);
iCSRmat getTransposeOfA(iCSRmat *A);
void getTransposeOfSparse(dCSRmat *A,dCSRmat *AT);
void getTransposeOfiden(idenmat *A, iCSRmat *AT, int cols, int nNode);
void getTransposeOfELEMENT(ELEMENT *A, iCSRmat *AT, int cols, int nNode);
void getTransposeOfelementDoF(ELEMENT_DOF *A, iCSRmat *AT, int cols);
int sparseAddition(dCSRmat *A, dCSRmat *B, dCSRmat *C);
int sparseSubtraction(dCSRmat *A, dCSRmat *B, dCSRmat *C);
int dCSRPlusdDiagVector(dCSRmat *A, dvector *B, dCSRmat *C);
int dCSRMultiplydDiagVector(dCSRmat *A, dvector *D, dCSRmat *C);
int dDiagVectorMultiplydCSR(dvector *D, dCSRmat *A, dCSRmat *C);
int dCSRPlusdBD(dCSRmat *A, dBDmat *B, dCSRmat *C);
int dBDMinusdCSR(dBDmat *A, dCSRmat *B, dCSRmat *C);
void sparseMultiplication(dCSRmat *A, dCSRmat *B, dCSRmat *C);
void dBDMultiplydCSR(dBDmat *A, dCSRmat *B, dCSRmat *C);
void dBD2MultiplydCSR(dBDmat *A, int srow, int scol, dCSRmat *B, dCSRmat *C);
void dCSRMultiplydBD(dCSRmat *A, dBDmat *B, dCSRmat *C);
void sparseTripleMultiplication(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B);
void sparseTripleMultiplication1(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B);
void sparseTripleMultiplication2(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B);
void sparseTripleMultiplication3(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B);
int sorteddIJtoCSR(dIJmat *A, dCSRmat *B);
int dIJtoCSR(dCSRmat *A, int *ia, int *ja, double *val, int N, int row, int col);
int dIJtoCSReps(dCSRmat *A, int *ia, int *ja, double *val, int N, int row, int col, double eps);
int dCSRtoIJ(dCSRmat *A, dIJmat *B);

/* poisson2d.c */
void poisson2d_f(double *x, double *val, double *paras);
void poisson2d_u(double *x, double *val, double *paras);
void poisson2d_gradu(double *x, double *val, double *paras);

/* poissonfem2d.c */
void poissonfem2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void poissonLagrange2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void assemble_poissonLagrange2d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void poissonCrouzeixRaviart2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void assemble_poissonCrouzeixRaviart2d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void poissonRaviartThomas2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void assemble_poissonRaviartThomas2d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);


/* poissonerror2d.c */
void geterrorsPoissonLagrange2d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void geterrorsPoissonCrouzeixRaviart2d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void geterrorsPoissonRaviartThomas2d(double *errors, dvector *uh, dvector *Qhu, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);

/* biharmonic2d.c */
void biharmonic2d_f(double *x, double *val, double *paras);
void biharmonic2d_u(double *x, double *val, double *paras);
void biharmonic2d_gradu(double *x, double *val, double *paras);
void biharmonic2d_hessu(double *x, double *val, double *paras);

/* biharmonicfem2d.c */
void biharmonicfem2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void biharmonicMorley2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void solve_biharmonicMorley2d(dvector *uh, ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void assemble_biharmonicMorley2d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void biharmonicC0ipdgExtrap2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void biharmonicC0ipdg2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void solve_biharmonicC0ipdg2d(dvector *uh, ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *elementdofTran, Input_data *Input);
void assemble_biharmonicC0ipdg2d(dCSRmat *A, dvector *b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double parapenalty);
void biharmonicPSP_MorleyNcfmP1P0Morley2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void solve_biharmonicPSP_MorleyNcfmP1P0Morley2d(dvector *uh, ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);

/* biharmonicerror2d.c */
void geterrorsBiharmonicMorley2d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void geterrorsBiharmonicC0ipdg2d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void geterrorsBiharmonicC0ipdgExtrap2d(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);

/* triharmonic2d.c */
void triharmonic2d_f(double *x, double *val, double *paras);
void triharmonic2d_u(double *x, double *val, double *paras);
void triharmonic2d_gradu(double *x, double *val, double *paras);
void triharmonic2d_hessu(double *x, double *val, double *paras);
void triharmonic2d_grad3u(double *x, double *val, double *paras);

/* triharmonicfem2d.c */
void triharmonicfem2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void triharmonicBSB_MorleyAnHuangMorley2d(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);
void solve_triharmonicBSB_MorleyAnHuangMorley2d(dvector *uh, dvector *sigmah, ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, Input_data *Input);

/* triharmonicerror2d.c */
void geterrorsTriharmonicMorleyMINIMorley2d(double *errors, dvector *uh, dvector *sigmah, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1);

/* linearElas2d.c */
void linearElas2d_f(double *x, double *val, double *paras);
void linearElas2d_u(double *x, double *val, double *paras);
void linearElas2d_gradu(double *x, double *val, double *paras);
void linearElas2d_sigma(double *x, double *val, double *paras);

/* linearElasfem2d.c */
void linearElas2dfem(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, Input_data *Input);
void linearElasHuZhang2d_mfem(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, Input_data *Input);
void assemble_linearElasHuZhang2d(dCSRmat *ptr_A, dvector *ptr_b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu);
void linearElasHuangZhou2d_mfem(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, Input_data *Input);
void assemble_linearElasHuangZhou2d(dCSRmat *ptr_A, dvector *ptr_b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu);

/* linearElaserror.c */
void geterrorslinearElasHuZhang2d(double *errors, dvector *sigmah, dvector *uh, dvector *Qhu, dvector *uhstar, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu);
void getposteriorierrorslinearElasHuZhang2d(double *errors, dvector *sigmah, dvector *uh, dvector *uhstar, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu);
void geterrorslinearElasHuangZhou2d(double *errors, dvector *sigmah, dvector *uh, dvector *Qhu, dvector *uhstar, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu);


/* assemble.c */
int getmesh(int domain_num, ELEMENT *ptr_elements, idenmat *ptr_elementEdge, EDGE *ptr_edges, dennode *ptr_nodes, iCSRmat *edgesTran, ivector *nodeCEdge, int levelNum);
void getElementDOFdg(ELEMENT_DOF *elementDOF, int nt, int ldof, int dop);
void getElementDOF1D_Continue(ELEMENT_DOF *edgeDOF, ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop);
void getElementDOF_Lagrange2d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop);
void getElementDOF_CrouzeixRaviart2d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges);
void getElementDOF_RaviartThomas2d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int dop);
void getElementDOF_BrezziDouglasMarini2d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int dop);
void getElementDOF_MINI2d(ELEMENT_DOF *elementDOF, ELEMENT *elements, int nvertices);
void getElementDOF_Morley2d(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices);
void getElementDOF_HuZhang(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop);
void getElementDOF_HuangZhou(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices);
void getFreenodesInfoLagrange2d(EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void getFreenodesInfoCrouzeixRaviart2d(EDGE *edges, ELEMENT_DOF *elementDOF);
void getFreenodesInfoMINI2d(dennode *nodes, ELEMENT_DOF *elementDOF);
void getFreenodesInfoMorley2d(EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void getFreenodesInfoDG(ELEMENT_DOF *elementDOF);
void assembleMassmatrixLagrange2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double mu);
void assembleMassmatrixLagrangeDiagBlockRepeat2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double mu);
void assembleBiGradLagrange2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double mu);
void assembleRHSLagrange2d(dvector *b, ELEMENT *elements, ELEMENT_DOF *elementDOF, void (*f)(double *, double *, double *), double *paras);
void assembleBiGradCrouzeixRaviart2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double mu);
void assembleRHSCrouzeixRaviart2d(dvector *b, ELEMENT *elements, ELEMENT_DOF *elementDOF, void (*f)(double *, double *, double *), double *paras);
void assembleMassmatrixRaviartThomas2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF);
void assembleDivRaviartThomasL2poly2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1);
void assembleRHSdgPoly2d(dvector *b, ELEMENT *elements, ELEMENT_DOF *elementDOF, void (*f)(double *, double *, double *), double *paras);
void assembleBiGradMorley2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void assembleBiHessMorley2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void assembleRHSMorley2d(dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, void (*f)(double *, double *, double *), double *paras);
void assembleBiHessC0ipdg2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double parapenalty);
void assembleBiGradCrouzeixRaviartDiagBlockRepeat2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double mu);
void assembleDivCrouzeixRaviartL2poly2d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1);
void assembleRHScurlMorleyCrouzeixRaviart2d(dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1, dvector *uh);
void assembleRHSCrouzeixRaviartcurlMorley2d(dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1, dvector *uh);
void assembleMassmatrixP02d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF);
void assembleBiGradMINIsymtensorDiagBlockRepeat2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double mu);
void assembleRotSMINILagrange2d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1);
void assembleRHShessMorleySMINI2d(dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1, dvector *uh);
void assembleRHSSMINIhessMorley2d(dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1, dvector *uh);
void assembleweightedMassmatrixHuZhang2d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF, double lambda, double mu);
void assembleDivHuZhangL2poly2d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF);
void assembleRHSdgPolyVector2d(dvector *b, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, void (*f)(double *, double *, double *), double *paras);
void assembleJumpL2poly2d(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double parapenalty);
void assembleweightedMassmatrixHuangZhou2d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF, double lambda, double mu);
void assembleDivHuangZhouL2poly2d(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF);
// void assembleStiffmatrixHuZhangIP2d(dCSRmat *A, dCSRmat *BT, dCSRmat *C, dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu);
void assembleStiffmatrixElasLagrange(dCSRmat *A, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu);
void assembleStiffmatrixVecLagrange(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu);
void sumNormalDerivativeEta(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, ddenmat *etas, double *sum);
void jumpOperatorVectorTensor(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump);
void jumpOperatorVector(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump);
void jumpOperator(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump);
void getinfo4jump(double lambda1, double lambda2, int *edgeNode, double elen, int element, ELEMENT *elements, int dop, int index, double *jump);
void TangentDerivative4Edge(double lambda1, double lambda2, int edge, int element, ELEMENT *elements, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *val);
void averageNormalDerivative2Lagrange2d(double *lambdas1d, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *avg);
void jumpNormalDerivativeLagrange2d(double *lambdas1d, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump);
void jumpOperatorTensor(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump);
void jumpOperatorTensorNormal(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump);
void jumpOperatorTensorTangent(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump);
void jumpOperatorTensorDIVPt(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump);
void jumpOperatorTensorDIV(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump);
void jumpOperatorATensorTangent2(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double lambda, double mu, double *jump);
void jumpOperatorRotATensorTangentPt(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double lambda, double mu, double *jump);
void getElementEdgeGeoInfo(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes);
void getBoundaryInfo(EDGE *edges, dennode *nodes, int dof, int dop, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index);
void getBoundaryInfoVector2d(EDGE *edges, dennode *nodes, int dof, int dop, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index);
void getBoundaryInfoNormalTrace(EDGE *edges, ELEMENT_DOF *elementDOF, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index);
void getBoundaryInfoEdge(EDGE *edges, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index);


/* assembleInfo.c */
void getEdgeInfo(ELEMENT *elements, iCSRmat *elementsTran, EDGE *edges, iCSRmat *edgesTran, dennode *nodes);
void getEdgeBdflag(EDGE *edges, dennode *nodes);
int getCoarseInfo(int domain_num, dennode *nodes, ELEMENT *elements, EDGE *edges, iCSRmat *elementsTran, iCSRmat *edgesTran, ivector *nodeCEdge);
void uniformrefine(dennode *Cnodes, ELEMENT *Celements, EDGE *Cedges, iCSRmat *CelementsTran, dennode *Fnodes, ELEMENT *Felements, EDGE *Fedges, iCSRmat *FelementsTran, iCSRmat *FedgesTran, ivector *nodeCEdge);
void repeatElementDoF(ELEMENT_DOF *elementDOFs, ELEMENT_DOF *elementDOFv, int n);
void extractNondirichletMatrixVector(dCSRmat *A, dvector *b, dCSRmat *A11, dvector *b1, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index, dvector *uh);
void extractNondirichletMatrix11(dCSRmat *A, dCSRmat *A11, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index);
void extractNondirichletMatrix1r(dCSRmat *A, dCSRmat *A1, ivector *dirichlet, ivector *nondirichlet);
void extractNondirichletMatrix1c(dCSRmat *A, dCSRmat *A1, ivector *isInNode, ivector *dirichlet, ivector *index);
void extractNondirichletMatrix1cBlock(dCSRmat *A, dCSRmat *A1, ivector *isInNode, ivector *dirichlet, ivector *index);
void extractNondirichletVector(dCSRmat *A, dvector *b, dvector *b1, ivector *dirichlet, ivector *nondirichlet, dvector *uh);
void extractFreenodesVector2StokesDirichlet(dCSRmat *A, dCSRmat *B, dvector *b, dvector *b1, ELEMENT_DOF *elementDOF, dvector *uh);
void extractFreenodesVector(dCSRmat *A, dvector *b, dvector *b1, ELEMENT_DOF *elementDOF, dvector *uh);
void updateFreenodesRHS(dCSRmat *A, dvector *b, dvector *x, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc, int flag);
void updateFreenodes2bRHS(dCSRmat *A, dvector *b, dvector *x, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1);
void extractFreenodesMatrix11(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc);
void updateFreenodesMatrix11(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc, int flag);
void updateFreenodes2bMatrix11(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOF0, ELEMENT_DOF *elementDOF1);
void extractFreenodesMatrix1r(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF);
void extractFreenodesMatrix1c(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF);
void extractFreenodesMatrix1cBlock(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF);
void extractFreenodesMatrix11cBlock(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc);
void extractFreenodesMatrix11cBlock3d(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc);
int getEdgeDOFsTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row31, int *row32);
int getEdgeDOFsScalarTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int rowstart, int row31, int row32);
int getEdgeDOFsVectorTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row31, int *row32);
int getEdgeDOFsVector(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart);
void patchtogether22(dCSRmat *A11, dCSRmat *A12, dCSRmat *A21, dCSRmat *A22, dCSRmat *A);
void patchtogether33(dCSRmat *A11, dCSRmat *A12, dCSRmat *A13, dCSRmat *A21, dCSRmat *A22, dCSRmat *A23, dCSRmat *A31, dCSRmat *A32, dCSRmat *A33, dCSRmat *A);
void getPenaltyParameters(ddenmat *etas, idenmat *elementEdge, EDGE *edges, double *alpha, double *beta);

/* basicData.c */
void cartToPol2d(double x, double y, double *r, double *theta);
void cartToBary2d(double *x, double *lambda, double **tri);
void baryToCart2d(double *lambda, double *x, double **tri);
void morley_basis(double *lambda, double **gradLambda, double *nve[3], int index, double *phi);
void morley_basis1(double *lambda, double **gradLambda, double *nve[3], int index, double phi[2]);
void morley_basis2(double **gradLambda, double *nve[3], int index, double phi[3]);
void lagrange1D_basis(double lambda, int index, int dop, double *phi);
void lagrange1D_basis1(double lambda, int index, int dop, double h, double *phi);
void lagrange_basis(double *lambda, int index, int dop, double *phi);
void lagrange_basis1(double *lambda, double **gradLambda, int index, int dop, double phi[2]);
void lagrange_basis2(double *lambda, double **gradLambda, int index, int dop, double phi[3]);
void bernstein2d_basis(double *lambda, int index, int dop, double *phi);
void bernstein2d_basis1(double *lambda, double **gradLambda, int index, int dop, double phi[2]);
void cr_basis(int n, double *lambda, int index, double *phi);
void cr_basis1(int n, double **gradLambda, int index, double *phi);
void rt_basis(double *lambda, double *height, double **tij, short *eorien, int index, int dop, double phi[2]);
void rt_basisDIV(double *lambda, double **grd_lambda, double *height, double **tij, short *eorien, int index, int dop, double *phi);
void rt_basis1(double *lambda, double **grd_lambda, double *height, double **tij, short *eorien, int index, int dop, double phi[4]);
void bdm_basis(double *lambda, double s, double eta[3], double xi[3], double **nv, double **nve, int index, int dop, double phi[2]);
void bdm_basis1(double *lambda, double s, double eta[3], double xi[3], double **nv, double **nve, int index, int dop, double phi[4]);
void arnoldwinther_basis(double *lambda, double *x, double *y, ddenmat3 *basisCoeffs, int element, int index, double *phi);
void arnoldwinther_basisDIV(double *lambda, ddenmat3 *basisCoeffs, int element, double s, double eta[3], double xi[3], int index, double *phi);
void arnoldwinther_basisDIV2(ddenmat3 *basisCoeffs, int element, double s, double eta[3], double xi[3], int index, double *phi);
void huzhang_basis(double *lambda, double **nv, double **tv, int index, int dop, double phi[3]);
void huzhang_basis1(double *lambda, double **gradLambda, double **nv, double **tv, int index, int dop, double(*phi)[2]);
void huzhang_basisDIV(double *lambda, double **gradLambda, double **nv, double **tv, int index, int dop, double phi[2]);
void huzhang_basisROT(double *lambda, double **gradLambda, double **nv, double **tv, int index, int dop, double phi[2]);
void huzhang_basisCurlTrace(double *lambda, double **gradLambda, double **nv, double **tv, int index, int dop, double phi[2]);
void huzhang_basisROTROT(double *lambda, double **gradLambda, double **nv, double **tv, int index, int dop, double *phi);
void huzhang_basisLaplaceTrace(double *lambda, double **gradLambda, double **nv, double **tv, int index, int dop, double *phi);
void divS_huangzhou_prebasis(double *lambda, double **tij, int i, double phi[3]);
void divS_huangzhou_prebasisDIV(double *lambda, double **tij, int i, double phi[2]);
void divS_huangzhou_basis(double *lambda, double s, double **nv, double **tv, double **tij, int index, double phi[3]);
void divS_huangzhou_basisDIV(double *lambda, double **gradLambda, double s, double **nv, double **tv, double **tij, int index, double phi[2]);
void mini_basis(double *lambda, int index, double *phi);
void mini_basis1(double *lambda, double **gradLambda, int index, double phi[2]);
double area(double **tri);
double area0(double *v0, double *v1, double *v2);

/* basiscoeff.c */
void generateBasisCoeffs(ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes);
void generateBasisCoeffsEdgeProj(dBDmat *B, dBDmat *bC, EDGE *edges, int dopl, int dopk);

/* quadrature.c */
int getNumQuadPoints(int dop, int dim);
int getNumQuadPoints_ShunnWilliams(int dop, int dim);
void init_quadrature(int num_qp, int ncoor, double (*gauss)[3]);
void init_Gauss(int num_qp, int ncoor, double (*gauss)[3]);
void init_Gauss2d(int num_qp, double (*gauss)[3], double *weight);
void init_Gauss1D(int num_qp, int ncoor, double (*gauss)[2]);
void init_Gauss1d(int num_qp, double (*gauss)[2], double *weight);
void init_NewtonCotes1D(int num_qp, int ncoor, double (*newtoncotes)[2]);
void init_ShunnWilliams2d(int num_qp, double (*lambdas)[3], double *weight);
void init_ShunnWilliams3d(int num_qp, double (*lambdas)[4], double *weight);

/* lu.c */
int LU_Decomp(double *A, int pivot[], int n);
int LU_Solve(double *A, double b[], int pivot[], double x[], int n);

/* post_processing.c */
void projL2Pk2d(dvector *Qhu, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, void (*u)(double *, double *, double *), double *paras);
void projL2PkVec2d(dvector *Qhu, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, void (*u)(double *, double *, double *), double *paras);

/* linearElaspostprocessing.c */
void postprocess2newDisplacementHuZhang2d(dvector *uhstar, dvector *sigmah, dvector *uh, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu);



/* output.c */
int getElementsNodes4Matlab(ELEMENT *elements, dennode *nodes, dvector *uh);
int writeElementsNodes(ELEMENT *elements, dennode *nodes, char *fname1, char *fname2);
int write_IJ_dCSRmat(dCSRmat *A, char *fname);
void write_IJ_dCSRmat4Matlab(dCSRmat *A, char *fname);
void write_dvector4Matlab(dvector *vec, char *fname);
void write_ivector4Matlab(ivector *vec, char *fname);
void read_dvector4Matlab(dvector *vec, char *fname);
void read_dvector4Matlab2b(dvector *vec, char *fname);


/* xuludmil.for (Xiangtan energy minimization code) */
//void get_p_xuludmil_(int *ia,int *ja,double *a, int *n, int *nc,int *ip,int *jp,double *pn,int *ipt,int *jpt,int *mf);

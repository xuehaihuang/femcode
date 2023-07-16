/*
 *  triharmonic.c
 *
 *  Created by Xuehai Huang on 30 Jan, 2023.
 *  Copyright 2023 SUFE. All rights reserved.
 *
 */

/*! \file triharmonic.c
 *  \brief Triharmonic equation
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "header.h"
#include "matvec.h"
#include "checkmat.h"

/**
 * \fn double triiharmonic2d_f(double *x, double *paras)
 * \brief load f, i.e. right hand side when the true solution u is sin(pi*x)*sin(pi*y)
 *		  -\Delta^3 u = f
 * \param *x the cooridates of the point in three dimensions
 * \return function value
 */
void triharmonic2d_f(double *x, double *val, double *paras)
{
	double t0, t2, t3, t4, t5, t6, t7;
	t2 = pi*x[0];
  	t3 = pi*x[1];
  	t4 = sin(t2);
  	t5 = sin(t3);
  	t6 = t4*t4;
  	t7 = t5*t5;
  	t0 = (pi*pi*pi*pi*pi*pi)*t4*t5*(t6*1.51E+2+t7*1.51E+2-t6*t7*2.43E+2-9.0E+1)*-2.4E+1;
	*val = t0;
}

/**
 * \fn double triharmonic2d_u(double *x, double *paras)
 * \brief true solution u
 * \param *x the cooridates of the point in three dimensions
 * \return function value 
 */
void triharmonic2d_u(double *x, double *val, double *paras)
{
	*val = pow(sin(pi*x[0]),3.0)*pow(sin(pi*x[1]),3.0);
}

/**
 * \fn void triharmonic2d_gradu(double *x, double *val, double *paras)
 * \brief the gradient of true solution u
 * \param *x the cooridates of the point in three dimensions
 * \param *val the gradient of u
 * \return the gradient of true solution u 
 */
void triharmonic2d_gradu(double *x, double *val, double *paras)
{
	double t0, t2, t3, t4, t5;
	t2 = pi*x[0];
	t3 = pi*x[1];
  	t4 = sin(t2);
  	t5 = sin(t3);
  	val[0] = pi*(t4*t4)*(t5*t5*t5)*cos(t2)*3.0;
  	val[1] = pi*(t4*t4*t4)*(t5*t5)*cos(t3)*3.0;
}

/**
 * \fn void triharmonic2d_hessu(double *x, double *val, double *paras)
 * \brief the hessian of true solution u
 * \param *x the cooridates of the point in three dimensions
 * \param *val the hessian of u
 * \return u_xx, u_yy, u_xy
 */
void triharmonic2d_hessu(double *x, double *val, double *paras)
{
	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
	t2 = pi*x[0];
	t3 = pi*x[1];
  	t4 = pi*pi;
  	t5 = cos(t2);
  	t6 = cos(t3);
  	t7 = sin(t2);
  	t8 = sin(t3);
  	t9 = t7*t7*t7;
  	t10 = t8*t8*t8;
  	t11 = t4*t9*t10*3.0;
  	t12 = -t11;
  	val[0] = t12+t4*(t5*t5)*t7*t10*6.0;
  	val[1] = t12+t4*(t6*t6)*t8*t9*6.0;
  	val[2] = t4*t5*t6*(t7*t7)*(t8*t8)*9.0;
}

/**
 * \fn void triharmonic2d_grad3u(double *x, double *val, double *paras)
 * \brief the third-order grad of true solution u
 * \param *x the cooridates of the point in three dimensions
 * \param *val the third-order grad of u
 * \return u_xxx, u_yyy, u_xxy, u_xyy 
 */
void triharmonic2d_grad3u(double *x, double *val, double *paras)
{
	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
	t2 = pi*x[0];
  	t3 = pi*x[1];
  	t4 = pi*pi*pi;
  	t5 = cos(t2);
  	t6 = cos(t3);
  	t7 = sin(t2);
  	t8 = sin(t3);
  	t9 = t7*t7;
  	t10 = t7*t7*t7;
  	t11 = t8*t8;
  	t12 = t8*t8*t8;
  	val[0] = t4*(t5*t5*t5)*t12*6.0-t4*t5*t9*t12*2.1E+1;
  	val[1] = t4*(t6*t6*t6)*t10*6.0-t4*t6*t10*t11*2.1E+1;
  	val[2] = t4*t6*t10*t11*-9.0+t4*(t5*t5)*t6*t7*t11*1.8E+1;
  	val[3] = t4*t5*t9*t12*-9.0+t4*t5*(t6*t6)*t8*t9*1.8E+1;
}

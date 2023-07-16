/*
 *  poisson.c
 *
 *  Created by Xuehai Huang on Apr 30, 2022.
 *  Copyright 2022 SUFE. All rights reserved.
 *
 */

/*! \file biharmonic.c
 *  \brief Biharmonic equation
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "header.h"
#include "matvec.h"
#include "checkmat.h"

/**
 * \fn double biharmonic2d_f(double *x, double *paras)
 * \brief load f, i.e. right hand side when the true solution u is sin(pi*x)*sin(pi*y)
 *		  -\Delta u = f
 * \param *x the cooridates of the point in three dimensions
 * \return function value
 */
void biharmonic2d_f(double *x, double *val, double *paras)
{
	double v = 8*pi*pi*pi*pi*(cos(2*pi*x[0])*cos(2*pi*x[1]) - cos(2*pi*x[0])*sin(pi*x[1])*sin(pi*x[1]) - cos(2*pi*x[1])*sin(pi*x[0])*sin(pi*x[0]));
	*val = v;
}

/**
 * \fn double biharmonic2d_u(double *x, double *paras)
 * \brief true solution u
 * \param *x the cooridates of the point in three dimensions
 * \return function value 
 */
void biharmonic2d_u(double *x, double *val, double *paras)
{
	*val = sin(pi*x[0])*sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[1]);
}

/**
 * \fn void biharmonic2d_gradu(double *x, double *val, double *paras)
 * \brief the gradient of true solution u
 * \param *x the cooridates of the point in three dimensions
 * \param *val the gradient of u
 * \return the gradient of true solution u 
 */
void biharmonic2d_gradu(double *x, double *val, double *paras)
{
	val[0] = pi*sin(2*pi*x[0])*sin(pi*x[1])*sin(pi*x[1]);
	val[1] = pi*sin(2*pi*x[1])*sin(pi*x[0])*sin(pi*x[0]);
}

/**
 * \fn void biharmonic2d_hessu(double *x, double *val, double *paras)
 * \brief the hessian of true solution u
 * \param *x the cooridates of the point in three dimensions
 * \param *val the hessian of u
 * \return the gradient of true solution u 
 */
void biharmonic2d_hessu(double *x, double *val, double *paras)
{
	val[0] = 2*pi*pi*cos(2*pi*x[0])*sin(pi*x[1])*sin(pi*x[1]);
	val[1] = 2*pi*pi*cos(2*pi*x[1])*sin(pi*x[0])*sin(pi*x[0]);
	val[2] = pi*pi*sin(2*pi*x[0])*sin(2*pi*x[1]);
}

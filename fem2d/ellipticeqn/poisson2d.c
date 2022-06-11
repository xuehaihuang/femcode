/*
 *  poisson.c
 *
 *  Created by Xuehai Huang on Apr 30, 2022.
 *  Copyright 2022 SUFE. All rights reserved.
 *
 */

/*! \file poisson.c
 *  \brief Poisson equation
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "header.h"
#include "matvec.h"
#include "checkmat.h"

/**
 * \fn double poisson2d_f(double *x, double *paras)
 * \brief load f, i.e. right hand side when the true solution u is sin(pi*x)*sin(pi*y)
 *		  -\Delta u = f
 * \param *x the cooridates of the point in three dimensions
 * \return function value
 */
double poisson2d_f(double *x, double *paras)
{
	double uxx=-pi*pi*sin(pi*x[0])*sin(pi*x[1]);
	double uyy=-pi*pi*sin(pi*x[0])*sin(pi*x[1]);
	return -uxx-uyy;
}

/**
 * \fn double poisson2d_u(double *x, double *paras)
 * \brief true solution u
 * \param *x the cooridates of the point in three dimensions
 * \return function value 
 */
double poisson2d_u(double *x, double *paras)
{
	return sin(pi*x[0])*sin(pi*x[1]);
}

/**
 * \fn void poisson2d_gradu(double *x, double *val, double *paras)
 * \brief the gradient of true solution u
 * \param *x the cooridates of the point in three dimensions
 * \param gradu the gradient of u
 * \return the gradient of true solution u 
 */
void poisson2d_gradu(double *x, double *val, double *paras)
{
	val[0] = pi*cos(pi*x[0])*sin(pi*x[1]);
	val[1] = pi*sin(pi*x[0])*cos(pi*x[1]);
}

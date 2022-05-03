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
 * \fn double poisson3d_f(double *x)
 * \brief load f, i.e. right hand side when the true solution u is (0.25-x*x)*(0.25-y*y)*(1-z*z)
 *		  -\Delta u = f
 * \param *x the cooridates of the point in three dimensions
 * \return function value
 */
double poisson3d_f(double *x)
{
	double uxx=-pi*pi*sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[2]);
	double uyy=-pi*pi*sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[2]);
	double uzz=-pi*pi*sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[2]);
	return -uxx-uyy-uzz;
}

/**
 * \fn double poisson3d_u(double *x)
 * \brief true solution u
 * \param *x the cooridates of the point in three dimensions
 * \return function value 
 */
double poisson3d_u(double *x)
{
	return sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[2]);
}

/**
 * \fn void poisson3d_gradu(double *x, double *val)
 * \brief the gradient of true solution u
 * \param *x the cooridates of the point in three dimensions
 * \param gradu the gradient of u
 * \return the gradient of true solution u 
 */
void poisson3d_gradu(double *x, double *val)
{
	val[0] = pi*cos(pi*x[0])*sin(pi*x[1])*sin(pi*x[2]);
	val[1] = pi*sin(pi*x[0])*cos(pi*x[1])*sin(pi*x[2]);
	val[2] = pi*sin(pi*x[0])*sin(pi*x[1])*cos(pi*x[2]);
}

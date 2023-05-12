/*
 *  stokes3d.c
 *
 *  Created by Xuehai Huang on MAy 8, 2023.
 *  Copyright 2023 SUFE. All rights reserved.
 *
 */

/*! \file stokes3d.c
 *  \brief Stokes equation
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "header.h"
#include "matvec.h"
#include "checkmat.h"

/**
 * \fn void stokes3d_f(double *x, double *val)
 * \brief load f, i.e. right hand side when the true solution u is (0.25-x*x)*(0.25-y*y)*(1-z*z)
 *		  -\Delta u - \grad p = f
 * \param *x the coordinates of the point in three dimensions
 * \return function value
 */
void stokes3d_f(double *x, double *val)
{
	double u1xx=2*pi*pi*cos(2*pi*x[0])*sin(2*pi*x[1])*sin(pi*x[2])*sin(pi*x[2]);
	double u1yy=-4*pi*pi*sin(pi*x[0])*sin(pi*x[0])*sin(2*pi*x[1])*sin(pi*x[2])*sin(pi*x[2]);
	double u1zz=2*pi*pi*sin(pi*x[0])*sin(pi*x[0])*sin(2*pi*x[1])*cos(2*pi*x[2]);

	double u2xx=4*pi*pi*sin(2*pi*x[0])*sin(pi*x[1])*sin(pi*x[1])*sin(pi*x[2])*sin(pi*x[2]);
	double u2yy=-2*pi*pi*sin(2*pi*x[0])*cos(2*pi*x[1])*sin(pi*x[2])*sin(pi*x[2]);
	double u2zz=-2*pi*pi*sin(2*pi*x[0])*sin(pi*x[1])*sin(pi*x[1])*cos(2*pi*x[2]);

	val[0] = -u1xx-u1yy-u1zz;
	val[1] = -u2xx-u2yy-u2zz;
	val[2] = 0;
	// double uxx=-pi*pi*sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[2]);
	// double uyy=-pi*pi*sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[2]);
	// double uzz=-pi*pi*sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[2]);

	// val[0] = -uxx-uyy-uzz;
	// val[1] = 2*val[0];
	// val[2] = 3*val[0];
}

/**
 * \fn void stokes3d_g(double *x)
 * \brief load g = div p
 * \param *x the coordinates of the point in three dimensions
 * \return function value
 */
double stokes3d_g(double *x)
{
	return 0;	
}

/**
 * \fn void stokes3d_u(double *x, double *val)
 * \brief true solution u
 * \param *x the coordinates of the point in three dimensions
 * \return function value 
 */
void stokes3d_u(double *x, double *val)
{
	val[0] = sin(pi*x[0])*sin(pi*x[0])*sin(2*pi*x[1])*sin(pi*x[2])*sin(pi*x[2]);
	val[1] = -sin(2*pi*x[0])*sin(pi*x[1])*sin(pi*x[1])*sin(pi*x[2])*sin(pi*x[2]);
	val[2] = 0;
	// val[0] = sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[2]);
	// val[1] = 2*sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[2]);
	// val[2] = 3*sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[2]);
}

/**
 * \fn double stokes3d_p(double *x)
 * \brief true solution p
 * \param *x the coordinates of the point in three dimensions
 * \return function value
 */
double stokes3d_p(double *x)
{
	return 0;
}

/**
 * \fn void stokes3d_gradu(double *x, double *val)
 * \brief the gradient of true solution u
 * \param *x the cooridates of the point in three dimensions
 * \param gradu the gradient of u
 * \return the gradient of true solution u 
 */
void stokes3d_gradu(double *x, double *val)
{
	val[0] = pi*sin(2*pi*x[0])*sin(2*pi*x[1])*sin(pi*x[2])*sin(pi*x[2]);
	val[1] = 2*pi*sin(pi*x[0])*sin(pi*x[0])*cos(2*pi*x[1])*sin(pi*x[2])*sin(pi*x[2]);
	val[2] = pi*sin(pi*x[0])*sin(pi*x[0])*sin(2*pi*x[1])*sin(2*pi*x[2]);

	val[3] = -2*pi*cos(2*pi*x[0])*sin(pi*x[1])*sin(pi*x[1])*sin(pi*x[2])*sin(pi*x[2]);
	val[4] = -pi*sin(2*pi*x[0])*sin(2*pi*x[1])*sin(pi*x[2])*sin(pi*x[2]);
	val[5] = -pi*sin(2*pi*x[0])*sin(pi*x[1])*sin(pi*x[1])*sin(2*pi*x[2]);

	val[6] = 0;
	val[7] = 0;
	val[8] = 0;
	// val[0] = pi*cos(pi*x[0])*sin(pi*x[1])*sin(pi*x[2]);
	// val[1] = pi*sin(pi*x[0])*cos(pi*x[1])*sin(pi*x[2]);
	// val[2] = pi*sin(pi*x[0])*sin(pi*x[1])*cos(pi*x[2]);

	// val[3] = 2*pi*cos(pi*x[0])*sin(pi*x[1])*sin(pi*x[2]);
	// val[4] = 2*pi*sin(pi*x[0])*cos(pi*x[1])*sin(pi*x[2]);
	// val[5] = 2*pi*sin(pi*x[0])*sin(pi*x[1])*cos(pi*x[2]);

	// val[6] = 3*pi*cos(pi*x[0])*sin(pi*x[1])*sin(pi*x[2]);
	// val[7] = 3*pi*sin(pi*x[0])*cos(pi*x[1])*sin(pi*x[2]);
	// val[8] = 3*pi*sin(pi*x[0])*sin(pi*x[1])*cos(pi*x[2]);
}

/*
 *  quadcurl.c
 *
 *  Created by Xuehai Huang on May 3, 2022.
 *  Copyright 2022 SUFE. All rights reserved.
 *
 */

/*! \file quadcurl.c
 *  \brief Quad curl equation
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "header.h"
#include "matvec.h"
#include "checkmat.h"

double quadcurl3d_phi1(double x, double y, double z);
double quadcurl3d_phi2(double x, double y, double z);
double quadcurl3d_phi3(double x, double y, double z);
double quadcurl3d_phi1_x(double x, double y, double z);
double quadcurl3d_phi1_y(double x, double y, double z);
double quadcurl3d_phi1_z(double x, double y, double z);
double quadcurl3d_phi2_x(double x, double y, double z);
double quadcurl3d_phi2_y(double x, double y, double z);
double quadcurl3d_phi2_z(double x, double y, double z);
double quadcurl3d_phi3_x(double x, double y, double z);
double quadcurl3d_phi3_y(double x, double y, double z);
double quadcurl3d_phi3_z(double x, double y, double z);

/**
* \fn void linearElas2d_f(double *x, double *val, double *paras)
* \brief f = (f1, f2)
* \param x pointer to the coordinate value of the point
* \param paras pointer to (lambda, mu)
* \return f = (f1, f2)
*/
void linearElas2d_f(double *x, double *val, double *paras)
{
	double lambda = paras[0], mu = paras[1];
	val[0] = -mu*pi*pi*pi * sin(2 * pi*x[1])*(2 * cos(2 * pi*x[0]) - 1);
    val[1] = mu*pi*pi*pi * sin(2 * pi*x[0])*(2 * cos(2 * pi*x[1]) - 1);
}

/**
* \fn void linearElas2d_u(double *x, double *val, double *paras)
* \brief u = (u1, u2)
* \param x pointer to the coordinate value of the point
* \param paras pointer to (lambda, mu)
* \return u = (u1, u2)
*/
void linearElas2d_u(double *x, double *val, double *paras)
{
	double lambda = paras[0], mu = paras[1];
	val[0] = pi*cos(pi*x[1])*sin(pi*x[0])*sin(pi*x[0]) * sin(pi*x[1]);
    val[1] = -pi*cos(pi*x[0])*sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[1]);
}

/**
* \fn void linearElas2d_gradu(double *x, double *val, double *paras)
* \brief (u1x, u1y, u2x, u2y)
* \param x pointer to the coordinate value of the point
* \param paras pointer to (lambda, mu)
* \return gradient of u
*/
void linearElas2d_gradu(double *x, double *val, double *paras)
{
	double lambda = paras[0], mu = paras[1];
	val[0] = 2 * pi * pi * cos(pi*x[0])*cos(pi*x[1])*sin(pi*x[0])*sin(pi*x[1]);
    val[1] = pi * pi * cos(2 * pi*x[1])*sin(pi*x[0])*sin(pi*x[0]);
	val[2] = -pi * pi * cos(2 * pi*x[0])*sin(pi*x[1])*sin(pi*x[1]);
	val[3] = -2 * pi * pi * cos(pi*x[0])*cos(pi*x[1])*sin(pi*x[0])*sin(pi*x[1]);
}

/**
* \fn void linearElas2d_sigma(double *x, double *val, double *paras)
* \brief sigma = (sigma11, sigma22, sigma12)
* \param x pointer to the coordinate value of the point
* \param paras pointer to (lambda, mu)
* \return sigma
*/
void linearElas2d_sigma(double *x, double *val, double *paras)
{
	double lambda = paras[0], mu = paras[1];
	val[0] = 4 * mu*pi*pi * cos(pi*x[0])*cos(pi*x[1])*sin(pi*x[0])*sin(pi*x[1]);
    val[1] = -4 * mu*pi*pi * cos(pi*x[0])*cos(pi*x[1])*sin(pi*x[0])*sin(pi*x[1]);
	val[2] = -(mu*pi*pi * (cos(2 * pi*x[0]) - cos(2 * pi*x[1]))) / 2;
}
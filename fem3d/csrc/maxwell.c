/*
 *  maxwell.c
 *
 *  Created by Xuehai Huang on Apr 30, 2022.
 *  Copyright 2022 SUFE. All rights reserved.
 *
 */

/*! \file maxwell.c
 *  \brief Maxwell equation
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "header.h"
#include "matvec.h"
#include "checkmat.h"

/**
 * \fn void maxwell3d_f(double *x, double *val)
 * \brief the vector load f
 * \param *x the cooridates of the point in three dimensions
 * \return function value
 */
void maxwell3d_f(double *x, double *val)
{
	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
	t2 = pii*x[0];
	t3 = pii*x[1];
  	t4 = pii*x[2];
  	t5 = cos(t2);
  	t6 = cos(t4);
  	t7 = t5*t5;
  	t8 = t6*t6;
  	val[0] = (pii*pii*pii)*cos(t3)*sin(t3)*(t7*5.0+t8*5.0-t7*t8*6.0-4.0)*-4.0;

	t2 = pii*x[0];
  	t3 = pii*x[1];
  	t4 = pii*x[2];
  	t5 = cos(t3);
  	t6 = cos(t4);
  	t7 = t5*t5;
  	t8 = t6*t6;
  	val[1] = (pii*pii*pii)*cos(t2)*sin(t2)*(t7*5.0+t8*5.0-t7*t8*6.0-4.0)*4.0;

	val[2] = 0.0;
}

/**
 * \fn void maxwell3d_u(double *x, double *val)
 * \brief the true solution u
 * \param *x the cooridates of the point in three dimensions
 * \return function value
 */
void maxwell3d_u(double *x, double *val)
{
	double t0, t2;
	t2 = pii*x[1];
	val[0] = pii*pow(sin(pii*x[0]),2.0)*pow(sin(pii*x[2]),2.0)*cos(t2)*sin(t2)*2.0;
	t2 = pii*x[0];
  	val[1] = pii*pow(sin(pii*x[1]),2.0)*pow(sin(pii*x[2]),2.0)*cos(t2)*sin(t2)*-2.0;
	val[2] = 0.0;
}

/**
 * \fn void maxwell3d_curlu(double *x, double *val)
 * \brief curl of the true solution u
 * \param *x the cooridates of the point in three dimensions
 * \return function value
 */
void maxwell3d_curlu(double *x, double *val)
{
	val[0] = maxwell3d_u3_y(x[0], x[1], x[2]) - maxwell3d_u2_z(x[0], x[1], x[2]);
	val[1] = maxwell3d_u1_z(x[0], x[1], x[2]) - maxwell3d_u3_x(x[0], x[1], x[2]);
	val[2] = maxwell3d_u2_x(x[0], x[1], x[2]) - maxwell3d_u1_y(x[0], x[1], x[2]);
}

/**
 * \fn void maxwell3d_gradu(double *x, double *val)
 * \brief gradient of the true solution u
 * \param *x the cooridates of the point in three dimensions
 * \return function value
 */
void maxwell3d_gradu(double *x, double *val)
{
	maxwell3d_gradu1(x, val);
	maxwell3d_gradu2(x, val+3);
	maxwell3d_gradu3(x, val+6);
}

/**
 * \fn void maxwell3d_gradu1(double *x, double *val)
 * \brief gradient of the true solution u1
 * \param *x the cooridates of the point in three dimensions
 * \return function value
 */
void maxwell3d_gradu1(double *x, double *val)
{
	val[0] = maxwell3d_u1_x(x[0], x[1], x[2]);
	val[1] = maxwell3d_u1_y(x[0], x[1], x[2]);
	val[2] = maxwell3d_u1_z(x[0], x[1], x[2]);
}

/**
 * \fn void maxwell3d_gradu2(double *x, double *val)
 * \brief gradient of the true solution u2
 * \param *x the cooridates of the point in three dimensions
 * \return function value
 */
void maxwell3d_gradu2(double *x, double *val)
{
	val[0] = maxwell3d_u2_x(x[0], x[1], x[2]);
	val[1] = maxwell3d_u2_y(x[0], x[1], x[2]);
	val[2] = maxwell3d_u2_z(x[0], x[1], x[2]);
}

/**
 * \fn void maxwell3d_gradu3(double *x, double *val)
 * \brief gradient of the true solution u3
 * \param *x the cooridates of the point in three dimensions
 * \return function value
 */
void maxwell3d_gradu3(double *x, double *val)
{
	val[0] = maxwell3d_u3_x(x[0], x[1], x[2]);
	val[1] = maxwell3d_u3_y(x[0], x[1], x[2]);
	val[2] = maxwell3d_u3_z(x[0], x[1], x[2]);
}

/**
* \fn double maxwell3d_u1_x(double x, double y, double z)
* \brief the x-directional partial derivative of true solution u1
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param z the z-axis value of the point
* \return x-directional partial derivative of true solution u1
*/
double maxwell3d_u1_x(double x, double y, double z)
{
	double t0, t2, t3;
	t2 = pii*x;
  	t3 = pii*y;
  	t0 = (pii*pii)*pow(sin(pii*z),2.0)*cos(t2)*cos(t3)*sin(t2)*sin(t3)*4.0;

	return t0;
}

/**
* \fn double maxwell3d_u1_y(double x, double y, double z)
* \brief the y-directional partial derivative of true solution u1
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param z the z-axis value of the point
* \return y-directional partial derivative of true solution u1
*/
double maxwell3d_u1_y(double x, double y, double z)
{
	double t0, t2, t3;
	t0 = (pii*pii)*cos(pii*y*2.0)*pow(sin(pii*x),2.0)*pow(sin(pii*z),2.0)*2.0;

	return t0;
}

/**
* \fn double maxwell3d_u1_z(double x, double y, double z)
* \brief the z-directional partial derivative of true solution u1
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param z the z-axis value of the point
* \return z-directional partial derivative of true solution u1
*/
double maxwell3d_u1_z(double x, double y, double z)
{
	double t0, t2, t3;
	t2 = pii*y;
  	t3 = pii*z;
  	t0 = (pii*pii)*pow(sin(pii*x),2.0)*cos(t2)*cos(t3)*sin(t2)*sin(t3)*4.0;

	return t0;
}

/**
* \fn double maxwell3d_u2_x(double x, double y, double z)
* \brief the x-directional partial derivative of true solution u2
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param z the z-axis value of the point
* \return x-directional partial derivative of true solution u2
*/
double maxwell3d_u2_x(double x, double y, double z)
{
	double t0, t2, t3;
	t0 = (pii*pii)*cos(pii*x*2.0)*pow(sin(pii*y),2.0)*pow(sin(pii*z),2.0)*-2.0;

	return t0;
}

/**
* \fn double maxwell3d_u2_y(double x, double y, double z)
* \brief the y-directional partial derivative of true solution u2
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param z the z-axis value of the point
* \return y-directional partial derivative of true solution u2
*/
double maxwell3d_u2_y(double x, double y, double z)
{
	double t0, t2, t3;
	t2 = pii*x;
  	t3 = pii*y;
  	t0 = (pii*pii)*pow(sin(pii*z),2.0)*cos(t2)*cos(t3)*sin(t2)*sin(t3)*-4.0;

	return t0;
}

/**
* \fn double maxwell3d_u2_z(double x, double y, double z)
* \brief the z-directional partial derivative of true solution u2
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param z the z-axis value of the point
* \return z-directional partial derivative of true solution u2
*/
double maxwell3d_u2_z(double x, double y, double z)
{
	double t0, t2, t3;
	t2 = pii*x;
  	t3 = pii*z;
  	t0 = (pii*pii)*pow(sin(pii*y),2.0)*cos(t2)*cos(t3)*sin(t2)*sin(t3)*-4.0;

	return t0;
}

/**
* \fn double maxwell3d_u3_x(double x, double y, double z)
* \brief the x-directional partial derivative of true solution u3
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param z the z-axis value of the point
* \return x-directional partial derivative of true solution u3
*/
double maxwell3d_u3_x(double x, double y, double z)
{
	return 0.0;
}

/**
* \fn double maxwell3d_u3_y(double x, double y, double z)
* \brief the y-directional partial derivative of true solution u3
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param z the z-axis value of the point
* \return y-directional partial derivative of true solution u3
*/
double maxwell3d_u3_y(double x, double y, double z)
{
	return 0.0;
}

/**
* \fn double maxwell3d_u3_z(double x, double y, double z)
* \brief the z-directional partial derivative of true solution u3
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param z the z-axis value of the point
* \return z-directional partial derivative of true solution u3
*/
double maxwell3d_u3_z(double x, double y, double z)
{
	return 0.0;
}

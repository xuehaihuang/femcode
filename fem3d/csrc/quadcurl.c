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
 * \fn void quadcurl3d_f(double *x, double *val)
 * \brief the vector load f
 * \param *x the cooridates of the point in three dimensions
 * \return function value
 */
void quadcurl3d_f(double *x, double *val)
{
	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
	t2 = pii*x[0];
	t3 = pii*x[1];
  	t4 = pii*x[2];
  	t5 = cos(t2);
  	t6 = cos(t3);
  	t7 = cos(t4);
  	t8 = t5*t5;
  	t9 = t6*t6;
  	t10 = t7*t7;
  	val[0] = (pii*pii*pii*pii*pii)*t6*sin(t2)*sin(t4)*(t8*3.85E+2+t9*2.49E+2+t10*3.85E+2-t8*t9*4.53E+2-t8*t10*6.37E+2-t9*t10*4.53E+2+t8*t9*t10*7.29E+2-2.05E+2)*-3.0;
	  
	// t2 = pii*x;
  	// t3 = pii*y;
  	// t4 = pii*z;
  	// t5 = cos(t2);
 	// t6 = cos(t3);
  	// t7 = cos(t4);
  	// t8 = t5*t5;
  	// t9 = t6*t6;
  	// t10 = t7*t7;
  	val[1] = (pii*pii*pii*pii*pii)*t5*sin(t3)*sin(t4)*(t8*2.49E+2+t9*3.85E+2+t10*3.85E+2-t8*t9*4.53E+2-t8*t10*4.53E+2-t9*t10*6.37E+2+t8*t9*t10*7.29E+2-2.05E+2)*3.0;
  	
	val[2] = 0.0;
}

/**
 * \fn void quadcurl3d_u(double *x, double *val)
 * \brief the true solution u
 * \param *x the cooridates of the point in three dimensions
 * \return function value
 */
void quadcurl3d_u(double *x, double *val)
{
	double t0, t2;
	t2 = pii*x[1];
	val[0] = pii*pow(sin(pii*x[0]),3.0)*pow(sin(pii*x[2]),3.0)*cos(t2)*pow(sin(t2),2.0)*3.0;
	t2 = pii*x[0];
  	val[1] = pii*pow(sin(pii*x[1]),3.0)*pow(sin(pii*x[2]),3.0)*cos(t2)*pow(sin(t2),2.0)*-3.0;
	val[2] = 0.0;
}

/**
 * \fn void quadcurl3d_curlu(double *x, double *val)
 * \brief curl of the true solution u
 * \param *x the cooridates of the point in three dimensions
 * \return function value
 */
void quadcurl3d_curlu(double *x, double *val)
{
	double t0, t2, t3, t4, t5, t6, t7;
	t2 = pii*x[0];
  	t3 = pii*x[2];
  	val[0] = (pii*pii)*pow(sin(pii*x[1]),3.0)*cos(t2)*cos(t3)*pow(sin(t2),2.0)*pow(sin(t3),2.0)*9.0;

	t2 = pii*x[1];
  	t3 = pii*x[2];
  	val[1] = (pii*pii)*pow(sin(pii*x[0]),3.0)*cos(t2)*cos(t3)*pow(sin(t2),2.0)*pow(sin(t3),2.0)*9.0;

	t2 = pii*x[0];
  	t3 = pii*x[1];
  	t4 = sin(t2);
  	t5 = sin(t3);
  	t6 = t4*t4;
  	t7 = t5*t5;
  	val[2] = (pii*pii)*t4*t5*pow(sin(pii*x[2]),3.0)*(t6+t7-t6*t7*3.0)*-6.0;
}

/**
 * \fn void quadcurl3d_gradcurlu(double *x, double *val)
 * \brief GradCurl of the true solution u
 * \param *x the cooridates of the point in three dimensions
 * \return function value
 */
void quadcurl3d_gradcurlu(double *x, double *val)
{
	val[0] = quadcurl3d_phi1_x(x[0], x[1], x[2]);
	val[1] = quadcurl3d_phi1_y(x[0], x[1], x[2]);
	val[2] = quadcurl3d_phi1_z(x[0], x[1], x[2]);
	val[3] = quadcurl3d_phi2_x(x[0], x[1], x[2]);
	val[4] = quadcurl3d_phi2_y(x[0], x[1], x[2]);
	val[5] = quadcurl3d_phi2_z(x[0], x[1], x[2]);
	val[6] = quadcurl3d_phi3_x(x[0], x[1], x[2]);
	val[7] = quadcurl3d_phi3_y(x[0], x[1], x[2]);
	val[8] = quadcurl3d_phi3_z(x[0], x[1], x[2]);
}

/**
* \fn double quadcurl3d_phi1_x(double x, double y, double z)
* \brief the x-directional partial derivative of true solution phi1
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param z the z-axis value of the point
* \return x-directional partial derivative of true solution curlu1
*/
double quadcurl3d_phi1_x(double x, double y, double z)
{
	double t0, t2;
	t2 = pii*z;
  	t0 = (pii*pii*pii)*sin(pii*x)*pow(sin(pii*y),3.0)*cos(t2)*pow(sin(t2),2.0)*(cos(pii*x*2.0)*(3.0/2.0)+1.0/2.0)*9.0;

	return t0;
}

/**
* \fn double quadcurl3d_phi1_y(double x, double y, double z)
* \brief the y-directional partial derivative of true solution phi1
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param z the z-axis value of the point
* \return y-directional partial derivative of true solution curlu1
*/
double quadcurl3d_phi1_y(double x, double y, double z)
{
	double t0, t2, t3, t4;
	t2 = pii*x;
  	t3 = pii*y;
  	t4 = pii*z;
  	t0 = (pii*pii*pii)*cos(t2)*cos(t3)*cos(t4)*pow(sin(t2),2.0)*pow(sin(t3),2.0)*pow(sin(t4),2.0)*2.7E+1;

	return t0;
}

/**
* \fn double quadcurl3d_phi1_z(double x, double y, double z)
* \brief the z-directional partial derivative of true solution phi1
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param z the z-axis value of the point
* \return z-directional partial derivative of true solution curlu1
*/
double quadcurl3d_phi1_z(double x, double y, double z)
{
	double t0, t2;
	t2 = pii*x;
  	t0 = (pii*pii*pii)*pow(sin(pii*y),3.0)*sin(pii*z)*cos(t2)*pow(sin(t2),2.0)*(cos(pii*z*2.0)*(3.0/2.0)+1.0/2.0)*9.0;

	return t0;
}

/**
* \fn double quadcurl3d_phi2_x(double x, double y, double z)
* \brief the x-directional partial derivative of true solution phi2
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param z the z-axis value of the point
* \return x-directional partial derivative of true solution curlu2
*/
double quadcurl3d_phi2_x(double x, double y, double z)
{
	double t0, t2, t3, t4;
	t2 = pii*x;
  	t3 = pii*y;
  	t4 = pii*z;
  	t0 = (pii*pii*pii)*cos(t2)*cos(t3)*cos(t4)*pow(sin(t2),2.0)*pow(sin(t3),2.0)*pow(sin(t4),2.0)*2.7E+1;

	return t0;
}

/**
* \fn double quadcurl3d_phi2_y(double x, double y, double z)
* \brief the y-directional partial derivative of true solution phi2
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param z the z-axis value of the point
* \return y-directional partial derivative of true solution curlu2
*/
double quadcurl3d_phi2_y(double x, double y, double z)
{
	double t0, t2;
	t2 = pii*z;
  	t0 = (pii*pii*pii)*pow(sin(pii*x),3.0)*sin(pii*y)*cos(t2)*pow(sin(t2),2.0)*(cos(pii*y*2.0)*(3.0/2.0)+1.0/2.0)*9.0;

	return t0;
}

/**
* \fn double quadcurl3d_phi2_z(double x, double y, double z)
* \brief the z-directional partial derivative of true solution phi2
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param z the z-axis value of the point
* \return z-directional partial derivative of true solution curlu2
*/
double quadcurl3d_phi2_z(double x, double y, double z)
{
	double t0, t2;
	t2 = pii*y;
  	t0 = (pii*pii*pii)*pow(sin(pii*x),3.0)*sin(pii*z)*cos(t2)*pow(sin(t2),2.0)*(cos(pii*z*2.0)*(3.0/2.0)+1.0/2.0)*9.0;

	return t0;
}

/**
* \fn double quadcurl3d_phi3_x(double x, double y, double z)
* \brief the x-directional partial derivative of true solution phi3
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param z the z-axis value of the point
* \return x-directional partial derivative of true solution curlu3
*/
double quadcurl3d_phi3_x(double x, double y, double z)
{
	double t0, t2, t3, t4, t5, t6, t7;
	t2 = pii*x;
  	t3 = pii*y;
  	t4 = sin(t2);
  	t5 = sin(t3);
  	t6 = t4*t4;
  	t7 = t5*t5;
  	t0 = (pii*pii*pii)*t5*pow(sin(pii*z),3.0)*cos(t2)*(t6*3.0+t7-t6*t7*9.0)*-6.0;

	return t0;
}

/**
* \fn double quadcurl3d_phi3_y(double x, double y, double z)
* \brief the y-directional partial derivative of true solution phi3
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param z the z-axis value of the point
* \return y-directional partial derivative of true solution curlu3
*/
double quadcurl3d_phi3_y(double x, double y, double z)
{
	double t0, t2, t3, t4, t5, t6, t7;
	t2 = pii*x;
  	t3 = pii*y;
  	t4 = sin(t2);
  	t5 = sin(t3);
  	t6 = t4*t4;
  	t7 = t5*t5;
  	t0 = (pii*pii*pii)*t4*pow(sin(pii*z),3.0)*cos(t3)*(t6+t7*3.0-t6*t7*9.0)*-6.0;

	return t0;
}

/**
* \fn double quadcurl3d_phi3_z(double x, double y, double z)
* \brief the z-directional partial derivative of true solution phi3
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param z the z-axis value of the point
* \return z-directional partial derivative of true solution curlu3
*/
double quadcurl3d_phi3_z(double x, double y, double z)
{
	double t0, t2, t3, t4, t5, t6, t7, t8;
	t2 = pii*x;
  	t3 = pii*y;
  	t4 = pii*z;
  	t5 = sin(t2);
  	t6 = sin(t3);
  	t7 = t5*t5;
  	t8 = t6*t6;
  	t0 = (pii*pii*pii)*t5*t6*cos(t4)*pow(sin(t4),2.0)*(t7+t8-t7*t8*3.0)*-1.8E+1;

	return t0;
}
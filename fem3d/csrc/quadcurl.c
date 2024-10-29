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

/**
 * \fn void quadcurl3d_f(double *x, double *val)
 * \brief the vector load f
 * \param *x the cooridates of the point in three dimensions
 * \return function value
 */
void quadcurl3d_f(double *x, double *val)
{

  double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31;
	
  t2 = pii*x[0];
  t3 = pii*x[1];
  t4 = pii*x[2];
  t5 = pii*pii*pii*pii*pii;
  t6 = cos(t2);
  t7 = cos(t3);
  t8 = cos(t4);
  t9 = sin(t2);
  t10 = sin(t3);
  t11 = sin(t4);
  t12 = t6*t6;
  t13 = t6*t6*t6;
  t14 = t7*t7;
  t15 = t7*t7*t7;
  t16 = t8*t8;
  t17 = t9*t9;
  t18 = t10*t10;
  t19 = t11*t11;
  t20 = t17*2.4E+1;
  t21 = t18*2.4E+1;
  t22 = t19*7.2E+1;
  t23 = t17*t18*9.2E+1;
  t24 = t17*t19*2.76E+2;
  t25 = t18*t19*2.76E+2;
  t29 = t17*t18*t19*7.29E+2;
  t26 = -t23;
  t27 = -t24;
  t28 = -t25;
  t30 = t20+t21+t22+t26+t27+t28+t29;
  t31 = t5*t8*t9*t10*t30*3.0;
  val[0] = -t31;
  val[1] = t31;
  val[2] = t5*t11*(t6*t10*2.05E+2-t7*t9*2.05E+2-t10*t13*2.49E+2+t9*t15*2.49E+2+t7*t9*t12*3.85E+2-t6*t10*t14*3.85E+2-t6*t10*t16*3.85E+2+t7*t9*t16*3.85E+2-t9*t12*t15*4.53E+2+t10*t13*t14*4.53E+2+t10*t13*t16*4.53E+2-t9*t15*t16*4.53E+2-t7*t9*t12*t16*6.37E+2+t6*t10*t14*t16*6.37E+2+t9*t12*t15*t16*7.29E+2-t10*t13*t14*t16*7.29E+2)*3.0;

  // t2 = pii*x[0];
  // t3 = pii*x[1];
  // t4 = pii*x[2];
  // t5 = pii*pii*pii*pii*pii;
  // t6 = cos(t2);
  // t7 = cos(t3);
  // t8 = cos(t4);
  // t9 = sin(t4);
  // t10 = t6*t6;
  // t11 = t7*t7;
  // t12 = t8*t8;
  // t13 = t12*3.85E+2;
  // t14 = t10*t11*4.53E+2;
  // t16 = t10*t11*t12*7.29E+2;
  // t15 = -t14;
  // val[0] = t5*t7*t9*sin(t2)*(t10*3.85E+2+t11*2.49E+2+t13+t15+t16-t10*t12*6.37E+2-t11*t12*4.53E+2-2.05E+2)*-3.0;
  // val[1] = t5*t6*t9*sin(t3)*(t10*2.49E+2+t11*3.85E+2+t13+t15+t16-t10*t12*4.53E+2-t11*t12*6.37E+2-2.05E+2)*3.0;
  // val[2] = 0;
 
}

/**
 * \fn void quadcurl3d_u(double *x, double *val)
 * \brief the true solution u
 * \param *x the cooridates of the point in three dimensions
 * \return function value
 */
void quadcurl3d_u(double *x, double *val)
{
	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13;

  t2 = pii*x[0];
  t3 = pii*x[1];
  t4 = pii*x[2];
  t5 = cos(t4);
  t6 = sin(t2);
  t7 = sin(t3);
  t8 = sin(t4);
  t9 = t6*t6*t6;
  t10 = t7*t7*t7;
  t11 = t8*t8;
  t12 = t8*t8*t8;
  t13 = pii*t5*t9*t10*t11*3.0;
  val[0] = -t13;
  val[1] = t13;
  val[2] = pii*(t6*t6)*t10*t12*cos(t2)*3.0-pii*(t7*t7)*t9*t12*cos(t3)*3.0;

  // t2 = pii*x[0];
  // t3 = pii*x[1];
  // t4 = pii*x[2];
  // t5 = sin(t2);
  // t6 = sin(t3);
  // t7 = sin(t4);
  // t8 = t7*t7*t7;
  // val[0] = pii*(t5*t5*t5)*(t6*t6)*t8*cos(t3)*3.0;
  // val[1] = pii*(t5*t5)*(t6*t6*t6)*t8*cos(t2)*-3.0;
  // val[2] = 0;

}

/**
 * \fn void quadcurl3d_curlu(double *x, double *val)
 * \brief curl of the true solution u
 * \param *x the cooridates of the point in three dimensions
 * \return function value
 */
void quadcurl3d_curlu(double *x, double *val)
{
  double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16;
	
  t2 = pii*x[0];
  t3 = pii*x[1];
  t4 = pii*x[2];
  t5 = pii*pii;
  t6 = cos(t2);
  t7 = cos(t3);
  t8 = cos(t4);
  t9 = sin(t2);
  t10 = sin(t3);
  t11 = sin(t4);
  t12 = t8*t8;
  t13 = t9*t9;
  t14 = t10*t10;
  t15 = t11*t11;
  val[0] = t5*t10*t11*t13*((t7*t7)*t9*t15*2.0+t9*t12*t14*2.0-t9*t14*t15*2.0-t6*t7*t10*t15*3.0)*-3.0;
  val[1] = t5*t9*t11*t14*((t6*t6)*t10*t15*2.0+t10*t12*t13*2.0-t10*t13*t15*2.0-t6*t7*t9*t15*3.0)*-3.0;
  val[2] = t5*t8*t13*t14*t15*sin(pii*(x[0]+x[1]))*9.0;

  // t2 = pii*x[0];
  // t3 = pii*x[1];
  // t4 = pii*x[2];
  // t5 = pii*pii;
  // t6 = cos(t4);
  // t7 = sin(t2);
  // t8 = sin(t3);
  // t9 = sin(t4);
  // t10 = t7*t7;
  // t11 = t8*t8;
  // t12 = t9*t9;
  // val[0] = t5*t6*(t8*t8*t8)*t10*t12*cos(t2)*9.0;
  // val[1] = t5*t6*(t7*t7*t7)*t11*t12*cos(t3)*9.0;
  // val[2] = t5*t7*t8*(t9*t9*t9)*(t10+t11-t10*t11*3.0)*-6.0;


}

/**
 * \fn void quadcurl3d_gradcurlu(double *x, double *val)
 * \brief GradCurl of the true solution u
 * \param *x the cooridates of the point in three dimensions
 * \return function value
 */
void quadcurl3d_gradcurlu(double *x, double *val)
{
   double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24;
	
  t2 = pii*x[0];
  t3 = pii*x[1];
  t4 = pii*x[2];
  t5 = pii*pii*pii;
  t6 = cos(t2);
  t7 = cos(t3);
  t8 = cos(t4);
  t9 = sin(t2);
  t10 = sin(t3);
  t11 = sin(t4);
  t12 = t6*t6;
  t13 = t7*t7;
  t14 = t8*t8;
  t15 = t9*t9;
  t16 = t10*t10;
  t17 = t11*t11;
  t18 = t6*t9*t13*t17*2.0;
  t19 = t7*t10*t12*t17*2.0;
  val[0] = t5*t9*t10*t11*(t18-t19+t6*t9*t14*t16*2.0-t6*t9*t16*t17*2.0+t7*t10*t15*t17)*-9.0;
  val[1] = t5*t11*t15*(t6*(t10*t10*t10)*t17*3.0+(t7*t7*t7)*t9*t17*2.0-t6*t10*t13*t17*6.0+t7*t9*t14*t16*6.0-t7*t9*t16*t17*1.0E+1)*-3.0;
  val[2] = t5*t8*t10*t15*(t9*t13*t17*6.0+t9*t14*t16*2.0-t9*t16*t17*1.0E+1-t6*t7*t10*t17*9.0)*-3.0;
  val[3] = t5*t11*t16*(t7*(t9*t9*t9)*t17*3.0+(t6*t6*t6)*t10*t17*2.0+t6*t10*t14*t15*6.0-t7*t9*t12*t17*6.0-t6*t10*t15*t17*1.0E+1)*-3.0;
  val[4] = t5*t9*t10*t11*(-t18+t19+t7*t10*t14*t15*2.0+t6*t9*t16*t17-t7*t10*t15*t17*2.0)*-9.0;
  val[5] = t5*t8*t9*t16*(t10*t12*t17*6.0+t10*t14*t15*2.0-t10*t15*t17*1.0E+1-t6*t7*t9*t17*9.0)*-3.0;
  val[6] = t5*t8*t9*t16*t17*(t10/2.0+sin(pii*(x[0]*2.0+x[1]))*(3.0/2.0))*9.0;
  val[7] = t5*t8*t10*t15*t17*(t9/2.0+sin(pii*(x[0]+x[1]*2.0))*(3.0/2.0))*9.0;
  val[8] = t5*t11*t15*t16*sin(pii*(x[0]+x[1]))*(cos(t4*2.0)*(3.0/2.0)+1.0/2.0)*9.0;

  // t2 = pii*x[0];
  // t3 = pii*x[1];
  // t4 = pii*x[2];
  // t5 = pii*pii*pii;
  // t6 = t4*2.0;
  // t7 = cos(t2);
  // t8 = cos(t3);
  // t9 = cos(t4);
  // t10 = sin(t2);
  // t11 = sin(t3);
  // t12 = sin(t4);
  // t13 = cos(t6);
  // t14 = t10*t10;
  // t15 = t10*t10*t10;
  // t16 = t11*t11;
  // t17 = t11*t11*t11;
  // t18 = t12*t12;
  // t19 = t12*t12*t12;
  // t20 = t13*(3.0/2.0);
  // t22 = t14*t16*9.0;
  // t24 = t5*t7*t8*t9*t14*t16*t18*2.7E+1;
  // t21 = t20+1.0/2.0;
  // t23 = -t22;
  // val[0] = t5*t9*t10*t17*t18*(cos(t2*2.0)*(3.0/2.0)+1.0/2.0)*9.0;
  // val[1] = t24;
  // val[2] = t5*t7*t12*t14*t17*t21*9.0;
  // val[3] = t24;
  // val[4] = t5*t9*t11*t15*t18*(cos(t3*2.0)*(3.0/2.0)+1.0/2.0)*9.0;
  // val[5] = t5*t8*t12*t15*t16*t21*9.0;
  // val[6] = t5*t7*t11*t19*(t14*3.0+t16+t23)*-6.0;
  // val[7] = t5*t8*t10*t19*(t14+t16*3.0+t23)*-6.0;
  // val[8] = t5*t9*t10*t11*t18*(t14+t16-t14*t16*3.0)*-1.8E+1;


}
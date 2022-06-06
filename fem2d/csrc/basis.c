/*
 *  basis.c
 *
 *  Created by Xuehai Huang on 3/29/09.
 *  Copyright 2009 PSU. All rights reserved.
 *
 */

/*! \file basis.c
 *  \brief Basis functions
 */
 
#include <math.h>
#include "header.h"

 /**
 * \fn double f(double x, double y)
 * \brief load f, i.e. right hand side when the true solution u is x*x*(x-1)*(x-1)*y*y*(y-1)*(y-1)
 *		  \Delta^2 u = f
 * \param x the x-axis value of the point
 * \param y the y-axis value of the point
 * \return function value
 */
double f(double x, double y)
{
	double lambda = 0.3;
	double nu = 0.35;

	return 0;
}

double f1(double x, double y, double lambda, double mu)
{
	return -mu*pi*pi*pi * sin(2 * pi*y)*(2 * cos(2 * pi*x) - 1);
}

double f2(double x, double y, double lambda, double mu)
{
	return mu*pi*pi*pi * sin(2 * pi*x)*(2 * cos(2 * pi*y) - 1);
}

/**
* \fn double u1(double x, double y, double lambda, double mu)
* \brief true solution u1
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param nu Lame constant or Poisson ratio of plate
* \return function value
*/
double u1(double x, double y, double lambda, double mu)
{
	return pi*cos(pi*y)*sin(pi*x)*sin(pi*x) * sin(pi*y);
}

/**
* \fn double u2(double x, double y, double lambda, double mu)
* \brief true solution u2
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return function value
*/
double u2(double x, double y, double lambda, double mu)
{
	return -pi*cos(pi*x)*sin(pi*x)*sin(pi*y)*sin(pi*y);
}

/**
* \fn double u1_x(double x, double y, double lambda, double mu)
* \brief the x-directional partial derivative of true solution u1
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return x-directional partial derivative of true solution u
*/
double u1_x(double x, double y, double lambda, double mu)
{
	return 2 * pi * pi * cos(pi*x)*cos(pi*y)*sin(pi*x)*sin(pi*y);
}

/**
* \fn double u1_y(double x, double y, double lambda, double mu)
* \brief the y-directional partial derivative of true solution u
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return y-directional partial derivative of true solution u
*/
double u1_y(double x, double y, double lambda, double mu)
{
	return pi * pi * cos(2 * pi*y)*sin(pi*x)*sin(pi*x);
}

/**
* \fn double u2_x(double x, double y, double lambda, double mu)
* \brief the x-directional partial derivative of true solution u1
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return x-directional partial derivative of true solution u
*/
double u2_x(double x, double y, double lambda, double mu)
{
	return -pi * pi * cos(2 * pi*x)*sin(pi*y)*sin(pi*y);
}

/**
* \fn double u2_y(double x, double y, double lambda, double mu)
* \brief the y-directional partial derivative of true solution u
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return y-directional partial derivative of true solution u
*/
double u2_y(double x, double y, double lambda, double mu)
{
	return -2 * pi * pi * cos(pi*x)*cos(pi*y)*sin(pi*x)*sin(pi*y);
}

/**
* \fn double sigma11(double x, double y, double lambda, double mu)
* \brief sigma11
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return sigma11
*/
double sigma11(double x, double y, double lambda, double mu)
{
	return 4 * mu*pi*pi * cos(pi*x)*cos(pi*y)*sin(pi*x)*sin(pi*y);
}

/**
* \fn double sigma22(double x, double y, double lambda, double mu)
* \brief sigma22
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return sigma22
*/
double sigma22(double x, double y, double lambda, double mu)
{
	return -4 * mu*pi*pi * cos(pi*x)*cos(pi*y)*sin(pi*x)*sin(pi*y);
}

/**
* \fn double sigma12(double x, double y, double lambda, double mu)
* \brief sigma12
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return sigma12
*/
double sigma12(double x, double y, double lambda, double mu)
{
	return -(mu*pi*pi * (cos(2 * pi*x) - cos(2 * pi*y))) / 2;
}

/**
* \fn double sigma21(double x, double y, double lambda, double mu)
* \brief sigma21
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return sigma21
*/
double sigma21(double x, double y, double lambda, double mu)
{
	return -(mu*pi*pi * (cos(2 * pi*x) - cos(2 * pi*y))) / 2;
}

/** 
 * \fn void morley_basis(double lambda1, double lambda2, double lambda3, double s, double elen[3], double eta[3], double xi[3], double sij[3], double orient[3], int index, double *phi)
 * \brief basis function of Morley element
 * \param lambda1 the first area coordiante
 * \param lambda2 the second area coordiante
 * \param lambda3 the third area coordiante
 * \param s the area of the triangule
 * \param elen[3] length of three edge
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param sij[3] some auxiliary parameter
 * \param orient[3] some auxiliary parameter
 * \param index the indicator of the basis function
 * \param *phi basis function
 * \return void
 */
void morley_basis(double lambda1, double lambda2, double lambda3, double s, double elen[3], double eta[3], double xi[3], double sij[3], double orient[3], int index, double *phi)
{
	if(index==0)
	{
		*phi=lambda1*lambda1 + (sij[2]/(elen[1]*elen[1])+sij[1]/(elen[2]*elen[2]))*lambda2*lambda3 + sij[1]/(elen[2]*elen[2])*lambda3*lambda1 + sij[2]/(elen[1]*elen[1])*lambda1*lambda2;
	}
	else if(index==1)
	{
		*phi=lambda2*lambda2 + sij[0]/(elen[2]*elen[2])*lambda2*lambda3 + (sij[0]/(elen[2]*elen[2])+sij[2]/(elen[0]*elen[0]))*lambda3*lambda1 + sij[2]/(elen[0]*elen[0])*lambda1*lambda2;
	}
	else if(index==2)
	{
		*phi=lambda3*lambda3 + sij[0]/(elen[1]*elen[1])*lambda2*lambda3 + sij[1]/(elen[0]*elen[0])*lambda3*lambda1 + (sij[1]/(elen[0]*elen[0])+sij[0]/(elen[1]*elen[1]))*lambda1*lambda2;
	}
	else if(index==3)
	{
		*phi=2*s/elen[0]*lambda1*(lambda1-1)*orient[0];
	}
	else if(index==4)
	{
		*phi=2*s/elen[1]*lambda2*(lambda2-1)*orient[1];
	}
	else if(index==5)
	{
		*phi=2*s/elen[2]*lambda3*(lambda3-1)*orient[2];
	}
	else
	{
		*phi=0;
	}
}

/** 
 * \fn void morley_basis1(double lambda1, double lambda2, double lambda3, double s, double elen[3], double eta[3], double xi[3], double sij[3], double orient[3], int index, double phi[2])
 * \brief the first order derivative of Morley element basis function: (\partial_{x}phi, \partial_{y}phi)
 * \param lambda1 the first area coordiante
 * \param lambda2 the second area coordiante
 * \param lambda3 the third area coordiante
 * \param s the area of the triangule
 * \param elen[3] length of three edge
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param sij[3] some auxiliary parameter
 * \param orient[3] some auxiliary parameter
 * \param index the indicator of the basis function
 * \param phi[2] the first order derivative of Morley element basis function: (\partial_{x}phi, \partial_{y}phi)
 * \return void
 */
void morley_basis1(double lambda1, double lambda2, double lambda3, double s, double elen[3], double eta[3], double xi[3], double sij[3], double orient[3], int index, double phi[2])
{
	if(index==0)
	{
		phi[0]=(2*lambda1*eta[0]+sij[2]/(elen[1]*elen[1])*eta[1]*(1-2*lambda2)+sij[1]/(elen[2]*elen[2])*eta[2]*(1-2*lambda3))/(2*s);
		phi[1]=-(2*lambda1*xi[0]+sij[2]/(elen[1]*elen[1])*xi[1]*(1-2*lambda2)+sij[1]/(elen[2]*elen[2])*xi[2]*(1-2*lambda3))/(2*s);
	}
	else if(index==1)
	{
		phi[0]=(2*lambda2*eta[1]+sij[0]/(elen[2]*elen[2])*eta[2]*(1-2*lambda3)+sij[2]/(elen[0]*elen[0])*eta[0]*(1-2*lambda1))/(2*s);
		phi[1]=-(2*lambda2*xi[1]+sij[0]/(elen[2]*elen[2])*xi[2]*(1-2*lambda3)+sij[2]/(elen[0]*elen[0])*xi[0]*(1-2*lambda1))/(2*s);
	}
	else if(index==2)
	{
		phi[0]=(2*lambda3*eta[2]+sij[1]/(elen[0]*elen[0])*eta[0]*(1-2*lambda1)+sij[0]/(elen[1]*elen[1])*eta[1]*(1-2*lambda2))/(2*s);
		phi[1]=-(2*lambda3*xi[2]+sij[1]/(elen[0]*elen[0])*xi[0]*(1-2*lambda1)+sij[0]/(elen[1]*elen[1])*xi[1]*(1-2*lambda2))/(2*s);
	}
	else if(index==3)
	{
		phi[0]=eta[0]/elen[0]*(2*lambda1-1)*orient[0];
		phi[1]=-xi[0]/elen[0]*(2*lambda1-1)*orient[0];
	}
	else if(index==4)
	{
		phi[0]=eta[1]/elen[1]*(2*lambda2-1)*orient[1];
		phi[1]=-xi[1]/elen[1]*(2*lambda2-1)*orient[1];
	}
	else if(index==5)
	{
		phi[0]=eta[2]/elen[2]*(2*lambda3-1)*orient[2];
		phi[1]=-xi[2]/elen[2]*(2*lambda3-1)*orient[2];
	}
	else
	{
		phi[0]=0;
		phi[1]=0;
	}
}

/** 
 * \fn void morley_basis2(double lambda1, double lambda2, double lambda3, double s, double elen[3], double eta[3], double xi[3], double sij[3], double orient[3], int index, double phi[3])
 * \brief the second order derivative of Morley element basis function: (\partial_{xx}phi, \partial_{yy}phi, \partial_{xy}phi)
 * \param lambda1 the first area coordiante
 * \param lambda2 the second area coordiante
 * \param lambda3 the third area coordiante
 * \param s the area of the triangule
 * \param elen[3] length of three edge
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param sij[3] some auxiliary parameter
 * \param orient[3] some auxiliary parameter
 * \param index the indicator of the basis function
 * \param phi[3] the second order derivative of Morley element basis function: (\partial_{xx}phi, \partial_{yy}phi, \partial_{xy}phi)
 * \return void
 */
void morley_basis2(double lambda1, double lambda2, double lambda3, double s, double elen[3], double eta[3], double xi[3], double sij[3], double orient[3], int index, double phi[3])
{
	if(index==0)
	{
		phi[0]=(eta[0]*eta[0]-sij[2]/(elen[1]*elen[1])*eta[1]*eta[1]-sij[1]/(elen[2]*elen[2])*eta[2]*eta[2])/(2*s*s);
		phi[1]=(xi[0]*xi[0]-sij[2]/(elen[1]*elen[1])*xi[1]*xi[1]-sij[1]/(elen[2]*elen[2])*xi[2]*xi[2])/(2*s*s);
		phi[2]=-(xi[0]*eta[0]-sij[2]/(elen[1]*elen[1])*xi[1]*eta[1]-sij[1]/(elen[2]*elen[2])*xi[2]*eta[2])/(2*s*s);
	}
	else if(index==1)
	{
		phi[0]=(eta[1]*eta[1]-sij[0]/(elen[2]*elen[2])*eta[2]*eta[2]-sij[2]/(elen[0]*elen[0])*eta[0]*eta[0])/(2*s*s);
		phi[1]=(xi[1]*xi[1]-sij[0]/(elen[2]*elen[2])*xi[2]*xi[2]-sij[2]/(elen[0]*elen[0])*xi[0]*xi[0])/(2*s*s);
		phi[2]=-(xi[1]*eta[1]-sij[0]/(elen[2]*elen[2])*xi[2]*eta[2]-sij[2]/(elen[0]*elen[0])*xi[0]*eta[0])/(2*s*s);
	}
	else if(index==2)
	{
		phi[0]=(eta[2]*eta[2]-sij[1]/(elen[0]*elen[0])*eta[0]*eta[0]-sij[0]/(elen[1]*elen[1])*eta[1]*eta[1])/(2*s*s);
		phi[1]=(xi[2]*xi[2]-sij[1]/(elen[0]*elen[0])*xi[0]*xi[0]-sij[0]/(elen[1]*elen[1])*xi[1]*xi[1])/(2*s*s);
		phi[2]=-(xi[2]*eta[2]-sij[1]/(elen[0]*elen[0])*xi[0]*eta[0]-sij[0]/(elen[1]*elen[1])*xi[1]*eta[1])/(2*s*s);
	}
	else if(index==3)
	{
		phi[0]=eta[0]*eta[0]/(s*elen[0])*orient[0];
		phi[1]=xi[0]*xi[0]/(s*elen[0])*orient[0];
		phi[2]=-xi[0]*eta[0]/(s*elen[0])*orient[0];
	}
	else if(index==4)
	{
		phi[0]=eta[1]*eta[1]/(s*elen[1])*orient[1];
		phi[1]=xi[1]*xi[1]/(s*elen[1])*orient[1];
		phi[2]=-xi[1]*eta[1]/(s*elen[1])*orient[1];
	}
	else if(index==5)
	{
		phi[0]=eta[2]*eta[2]/(s*elen[2])*orient[2];
		phi[1]=xi[2]*xi[2]/(s*elen[2])*orient[2];
		phi[2]=-xi[2]*eta[2]/(s*elen[2])*orient[2];
	}
	else
	{
		phi[0]=0;
		phi[1]=0;
		phi[2]=0;
	}
}

/** 
 * \fn void lagrange1D_basis(double lambda, int index, int dop, double *phi)
 * \brief basis function of Lagrange element in 1d space
 * \param lambda  area coordiante
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param *phi basis function
 * \return void
 */
void lagrange1D_basis(double lambda, int index, int dop, double *phi)
{
	int in, ie, ii;
	int dofs = dop+1; // degrees of freedom
	if(index>= dofs || index<0)
	{
		*phi=0;
		return;
	}

	if(dop==0)
	{
		*phi = 1;
	}

	else if(dop==1)
	{
		switch(index)
		{
		case 0: *phi = lambda; break;
		case 1: *phi = 1-lambda; break;
		default: *phi = 0;
		}
	} // dop=1
	
	else if(dop==2)
	{
		switch(index)
		{
		case 0: *phi = lambda*(2*lambda-1); break;
		case 1: *phi = (1-lambda)*(1-2*lambda); break;
		case 2: *phi = 4*lambda*(1-lambda); break;
		default: *phi = 0;
		}
	} // dop=2
	
	else if(dop==3)
	{
		switch(index)
		{
		case 0: *phi = lambda*(3*lambda-1)*(3*lambda-2)/2.0; break;
		case 1: *phi = -(3*lambda-1)*(3*lambda-2)*(lambda-1)/2.0; break;
		case 2: *phi = -9*lambda*(3*lambda-1)*(lambda-1)/2.0; break;
		case 3: *phi = 9*lambda*(3*lambda-2)*(lambda-1)/2.0; break;
		default: *phi = 0;
		}
	} // dop=3

	else if(dop==4)
	{
		switch(index)
		{
		case 0: *phi = lambda*(4*lambda-1)*(2*lambda-1)*(4*lambda-3)/3.0; break;
		case 1: *phi = (4*lambda-1)*(2*lambda-1)*(4*lambda-3)*(lambda-1)/3.0; break;
		case 2: *phi = -16*lambda*(4*lambda-1)*(2*lambda-1)*(lambda-1)/3.0; break;
		case 3: *phi = 4*lambda*(4*lambda-1)*(4*lambda-3)*(lambda-1); break;
		case 4: *phi = -16*lambda*(2*lambda-1)*(4*lambda-3)*(lambda-1)/3.0; break;
		default: *phi = 0;
		}
	} // dop=4

	else if(dop==5)
	{
		switch(index)
		{
		case 0: *phi = lambda*(5*lambda-1)*(5*lambda-2)*(5*lambda-3)*(5*lambda-4)/24.0; break;
		case 1: *phi = -(5*lambda-1)*(5*lambda-2)*(5*lambda-3)*(5*lambda-4)*(lambda-1)/24.0; break;
		case 2: *phi = -25*lambda*(5*lambda-1)*(5*lambda-2)*(5*lambda-3)*(lambda-1)/24.0; break;
		case 3: *phi = 25*lambda*(5*lambda-1)*(5*lambda-2)*(5*lambda-4)*(lambda-1)/12.0; break;
		case 4: *phi = -25*lambda*(5*lambda-1)*(5*lambda-3)*(5*lambda-4)*(lambda-1)/12.0; break;
		case 5: *phi = 25*lambda*(5*lambda-2)*(5*lambda-3)*(5*lambda-4)*(lambda-1)/24.0; break;
		default: *phi = 0;
		}
	} // dop=5
}

/** 
 * \fn void lagrange1D_basis1(double lambda, int index, int dop, double h, double *phi)
 * \brief the first order derivative of basis function of Lagrange element in 1d space
 * \param lambda  area coordiante
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param h length of edge
 * \param *phi basis function
 * \return void
 */
void lagrange1D_basis1(double lambda, int index, int dop, double h, double *phi)
{
	int in, ie, ii;
	int dofs = dop+1; // degrees of freedom
	if(index>= dofs || index<0)
	{
		*phi=0;
		return;
	}

	double slp=-1.0/h; // slope

	if(dop==0)
	{
		*phi = 0;
	}

	else if(dop==1)
	{
		switch(index)
		{
		case 0: *phi = slp; break;
		case 1: *phi = -slp; break;
		default: *phi = 0;
		}
	} // dop=1
	
	else if(dop==2)
	{
		switch(index)
		{
		case 0: *phi = slp*(2*lambda-1) + lambda*2*slp; break;
		case 1: *phi = -slp*(1-2*lambda) - (1-lambda)*2*slp; break;
		case 2: *phi = 4*slp*(1-lambda) - 4*lambda*slp; break;
		default: *phi = 0;
		}
	} // dop=2
	
	else if(dop==3)
	{
		switch(index)
		{
		case 0: *phi = slp*(3*lambda-1)*(3*lambda-2)/2.0 + lambda*3*slp*(3*lambda-2)/2.0 + lambda*(3*lambda-1)*3*slp/2.0; break;
		case 1: *phi = -3*slp*(3*lambda-2)*(lambda-1)/2.0 - (3*lambda-1)*3*slp*(lambda-1)/2.0 - (3*lambda-1)*(3*lambda-2)*slp/2.0; break;
		case 2: *phi = -9*slp*(3*lambda-1)*(lambda-1)/2.0 - 9*lambda*3*slp*(lambda-1)/2.0 - 9*lambda*(3*lambda-1)*slp/2.0; break;
		case 3: *phi = 9*slp*(3*lambda-2)*(lambda-1)/2.0 + 9*lambda*3*slp*(lambda-1)/2.0 + 9*lambda*(3*lambda-2)*slp/2.0; break;
		default: *phi = 0;
		}
	} // dop=3

	else if(dop==4)
	{
		switch(index)
		{
		case 0: *phi = slp*(4*lambda-1)*(2*lambda-1)*(4*lambda-3)/3.0 + lambda*4*slp*(2*lambda-1)*(4*lambda-3)/3.0 + lambda*(4*lambda-1)*2*slp*(4*lambda-3)/3.0 + lambda*(4*lambda-1)*(2*lambda-1)*4*slp/3.0; break;
		case 1: *phi = 4*slp*(2*lambda-1)*(4*lambda-3)*(lambda-1)/3.0 + (4*lambda-1)*2*slp*(4*lambda-3)*(lambda-1)/3.0 + (4*lambda-1)*(2*lambda-1)*4*slp*(lambda-1)/3.0 + (4*lambda-1)*(2*lambda-1)*(4*lambda-3)*slp/3.0; break;
		case 2: *phi = -16*slp*(4*lambda-1)*(2*lambda-1)*(lambda-1)/3.0 - 16*lambda*4*slp*(2*lambda-1)*(lambda-1)/3.0 - 16*lambda*(4*lambda-1)*2*slp*(lambda-1)/3.0 - 16*lambda*(4*lambda-1)*(2*lambda-1)*slp/3.0; break;
		case 3: *phi = 4*slp*(4*lambda-1)*(4*lambda-3)*(lambda-1) + 4*lambda*4*slp*(4*lambda-3)*(lambda-1) + 4*lambda*(4*lambda-1)*4*slp*(lambda-1) + 4*lambda*(4*lambda-1)*(4*lambda-3)*slp; break;
		case 4: *phi = -16*slp*(2*lambda-1)*(4*lambda-3)*(lambda-1)/3.0 - 16*lambda*2*slp*(4*lambda-3)*(lambda-1)/3.0 - 16*lambda*(2*lambda-1)*4*slp*(lambda-1)/3.0 - 16*lambda*(2*lambda-1)*(4*lambda-3)*slp/3.0; break;
		default: *phi = 0;
		}
	} // dop=4

	else if(dop==5)
	{
		switch(index)
		{
		case 0: *phi = slp*(5*lambda-1)*(5*lambda-2)*(5*lambda-3)*(5*lambda-4)/24.0 + lambda*5*slp*(5*lambda-2)*(5*lambda-3)*(5*lambda-4)/24.0 + lambda*(5*lambda-1)*5*slp*(5*lambda-3)*(5*lambda-4)/24.0 + lambda*(5*lambda-1)*(5*lambda-2)*5*slp*(5*lambda-4)/24.0 + lambda*(5*lambda-1)*(5*lambda-2)*(5*lambda-3)*5*slp/24.0; break;
		case 1: *phi = -5*slp*(5*lambda-2)*(5*lambda-3)*(5*lambda-4)*(lambda-1)/24.0 - (5*lambda-1)*5*slp*(5*lambda-3)*(5*lambda-4)*(lambda-1)/24.0 - (5*lambda-1)*(5*lambda-2)*5*slp*(5*lambda-4)*(lambda-1)/24.0 - (5*lambda-1)*(5*lambda-2)*(5*lambda-3)*5*slp*(lambda-1)/24.0 - (5*lambda-1)*(5*lambda-2)*(5*lambda-3)*(5*lambda-4)*slp/24.0; break;
		case 2: *phi = -25*slp*(5*lambda-1)*(5*lambda-2)*(5*lambda-3)*(lambda-1)/24.0 - 25*lambda*5*slp*(5*lambda-2)*(5*lambda-3)*(lambda-1)/24.0 - 25*lambda*(5*lambda-1)*5*slp*(5*lambda-3)*(lambda-1)/24.0 - 25*lambda*(5*lambda-1)*(5*lambda-2)*5*slp*(lambda-1)/24.0 - 25*lambda*(5*lambda-1)*(5*lambda-2)*(5*lambda-3)*slp/24.0; break;
		case 3: *phi = 25*slp*(5*lambda-1)*(5*lambda-2)*(5*lambda-4)*(lambda-1)/12.0 + 25*lambda*5*slp*(5*lambda-2)*(5*lambda-4)*(lambda-1)/12.0 + 25*lambda*(5*lambda-1)*5*slp*(5*lambda-4)*(lambda-1)/12.0 + 25*lambda*(5*lambda-1)*(5*lambda-2)*5*slp*(lambda-1)/12.0 + 25*lambda*(5*lambda-1)*(5*lambda-2)*(5*lambda-4)*slp/12.0; break;
		case 4: *phi = -25*slp*(5*lambda-1)*(5*lambda-3)*(5*lambda-4)*(lambda-1)/12.0 - 25*lambda*5*slp*(5*lambda-3)*(5*lambda-4)*(lambda-1)/12.0 - 25*lambda*(5*lambda-1)*5*slp*(5*lambda-4)*(lambda-1)/12.0 - 25*lambda*(5*lambda-1)*(5*lambda-3)*5*slp*(lambda-1)/12.0 - 25*lambda*(5*lambda-1)*(5*lambda-3)*(5*lambda-4)*slp/12.0; break;
		case 5: *phi = 25*slp*(5*lambda-2)*(5*lambda-3)*(5*lambda-4)*(lambda-1)/24.0 + 25*lambda*5*slp*(5*lambda-3)*(5*lambda-4)*(lambda-1)/24.0 + 25*lambda*(5*lambda-2)*5*slp*(5*lambda-4)*(lambda-1)/24.0 + 25*lambda*(5*lambda-2)*(5*lambda-3)*5*slp*(lambda-1)/24.0 + 25*lambda*(5*lambda-2)*(5*lambda-3)*(5*lambda-4)*slp/24.0; break;
		default: *phi = 0;
		}
	} // dop=5
}

/** 
 * \fn void lagrange_basis(double *lambda, int index, int dop, double *phi)
 * \brief basis function of Lagrange element
 * \param *lambda pointer to the area coordiante
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param *phi basis function
 * \return void
 */
void lagrange_basis(double *lambda, int index, int dop, double *phi)
{
	int in, ie, ii;
	int dofs = (dop+1)*(dop+2)/2; // degrees of freedom
	if(index>= dofs || index<0)
	{
		*phi=0;
		return;
	}

	if(dop==0)
	{
		*phi = 1;
	}

	else if(dop==1)
	{
		*phi = lambda[index];
	} // dop=1
	
	else if(dop==2)
	{
		if(index<3)
		{
			*phi = lambda[index]*(2*lambda[index]-1);
		}
		else
		{
			*phi = 4.0*lambda[(index+1)%3]*lambda[(index+2)%3];
		}
	} // dop=2
	
	else if(dop==3)
	{
		if(index<3)
		{
			*phi = lambda[index]*(3*lambda[index]-1)*(3*lambda[index]-2)/2.0;
		}
		else if(index < 3*dop) 
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1);
			*phi = lambda[(ie+1)%3]*lambda[(ie+2)%3]*(3*lambda[(ie+1+ii)%3]-1)*9.0/2.0;
		}
		else
		{
			*phi = 27.0*lambda[0]*lambda[1]*lambda[2];
		}
	} // dop=3

	else if(dop==4)
	{
		if(index<3)
		{
			*phi = lambda[index]*(4*lambda[index]-1)*(4*lambda[index]-2)*(4*lambda[index]-3)/6.0;
		}
		else if(index < 3*dop)
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1);
			switch(ii)
			{
			case 0:
				*phi = lambda[(ie+1)%3]*lambda[(ie+2)%3]*(4*lambda[(ie+1)%3]-1)*(4*lambda[(ie+1)%3]-2)*8.0/3.0;
				break; 
			case 1:
				*phi = lambda[(ie+1)%3]*lambda[(ie+2)%3]*(4*lambda[(ie+1)%3]-1)*(4*lambda[(ie+2)%3]-1)*4.0;
				break;
			case 2:
				*phi = lambda[(ie+1)%3]*lambda[(ie+2)%3]*(4*lambda[(ie+2)%3]-1)*(4*lambda[(ie+2)%3]-2)*8.0/3.0;
				break;
			default:
				*phi = 0;
			}
		}
		else
		{
			in = index-3*dop;
			*phi = 32.0*lambda[0]*lambda[1]*lambda[2]*(4*lambda[in]-1);
		}
	} // dop=4

	else if(dop==5)
	{
		if(index<3)
		{
			*phi = lambda[index]*(5*lambda[index]-1)*(5*lambda[index]-2)*(5*lambda[index]-3)*(5*lambda[index]-4)/24.0;
		}
		else if(index < 3*dop)
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1);
			switch(ii)
			{
			case 0:
				*phi = lambda[(ie+1)%3]*lambda[(ie+2)%3]*(5*lambda[(ie+1)%3]-1)*(5*lambda[(ie+1)%3]-2)*(5*lambda[(ie+1)%3]-3)*25.0/24.0;
				break; 
			case 1:
				*phi = lambda[(ie+1)%3]*lambda[(ie+2)%3]*(5*lambda[(ie+1)%3]-1)*(5*lambda[(ie+2)%3]-1)*(5*lambda[(ie+1)%3]-2)*25.0/12.0;
				break;
			case 2:
				*phi = lambda[(ie+1)%3]*lambda[(ie+2)%3]*(5*lambda[(ie+1)%3]-1)*(5*lambda[(ie+2)%3]-1)*(5*lambda[(ie+2)%3]-2)*25.0/12.0;
				break;
			case 3:
				*phi = lambda[(ie+1)%3]*lambda[(ie+2)%3]*(5*lambda[(ie+2)%3]-1)*(5*lambda[(ie+2)%3]-2)*(5*lambda[(ie+2)%3]-3)*25.0/24.0;
				break;
			default:
				*phi = 0;
			}
		}
		else if(index < 3*dop+3)
		{
			in = index-3*dop;
			*phi = lambda[0]*lambda[1]*lambda[2]*(5*lambda[in]-1)*(5*lambda[in]-2)*125.0/6.0;
		}
		else
		{
			in = index-3*dop-3;
			*phi = lambda[0]*lambda[1]*lambda[2]*(5*lambda[(in+1)%3]-1)*(5*lambda[(in+2)%3]-1)*125.0/4.0;
		}
	} // dop=5
}

/** 
 * \fn void lagrange_basis1(double *lambda, double s, double eta[3], double xi[3], int index, int dop, double phi[2])
 * \brief the first order derivative of Lagrange element basis function: (\partial_{x}phi, \partial_{y}phi)
 * \param *lambda pointer to the area coordiante
 * \param s the area of the triangule
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param phi[2] the first order derivative of Morley element basis function: (\partial_{x}phi, \partial_{y}phi)
 * \return void
 */
void lagrange_basis1(double *lambda, double s, double eta[3], double xi[3], int index, int dop, double phi[2])
{
	int in, ie, ii, i1, i2, i3;
	int dofs = (dop+1)*(dop+2)/2; // degrees of freedom
	if(index>= dofs || index<0)
	{
		phi[0]=0;
		phi[1]=0;
		return;
	}

	if(dop==0)
	{
		phi[0]=0;
		phi[1]=0;
	} // dop=0

	else if(dop==1)
	{
		phi[0]=eta[index]/(2.0*s);
		phi[1]=-xi[index]/(2.0*s);
	} // dop=1

	else if(dop==2)
	{
		if(index<3)
		{
			phi[0]=(4*lambda[index]-1)*eta[index]/(2.0*s);
			phi[1]=-(4*lambda[index]-1)*xi[index]/(2.0*s);
		}
		else
		{
			i1=(index+1)%3;
			i2=(index+2)%3;
			phi[0]=4*(lambda[i1]*eta[i2]+lambda[i2]*eta[i1])/(2.0*s);
			phi[1]=-4*(lambda[i1]*xi[i2]+lambda[i2]*xi[i1])/(2.0*s);
		}
	} // dop=2

	else if(dop==3)
	{
		if(index<3)
		{
			phi[0]=((3*lambda[index]-1)*(3*lambda[index]-2) + 3*lambda[index]*(3*lambda[index]-2) + 3*lambda[index]*(3*lambda[index]-1))*eta[index]/(4.0*s);
			phi[1]=-((3*lambda[index]-1)*(3*lambda[index]-2) + 3*lambda[index]*(3*lambda[index]-2) + 3*lambda[index]*(3*lambda[index]-1))*xi[index]/(4.0*s);
		}
		else if(index < 3*dop) 
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1);
			i1 = (ie+1)%3;
			i2 = (ie+2)%3;
			i3 = (ie+1+ii)%3;
			phi[0]=(lambda[i2]*(3*lambda[i3]-1)*eta[i1] + lambda[i1]*(3*lambda[i3]-1)*eta[i2] + 3*lambda[i1]*lambda[i2]*eta[i3])*9.0/(4.0*s);
			phi[1]=-(lambda[i2]*(3*lambda[i3]-1)*xi[i1] + lambda[i1]*(3*lambda[i3]-1)*xi[i2] + 3*lambda[i1]*lambda[i2]*xi[i3])*9.0/(4.0*s);
		}
		else
		{
			phi[0]=27*(lambda[1]*lambda[2]*eta[0] + lambda[2]*lambda[0]*eta[1] + lambda[0]*lambda[1]*eta[2])/(2.0*s);
			phi[1]=-27*(lambda[1]*lambda[2]*xi[0] + lambda[2]*lambda[0]*xi[1] + lambda[0]*lambda[1]*xi[2])/(2.0*s);
		}
	} // dop=3

	else if(dop==4)
	{
		if(index<3)
		{
			phi[0]=((4*lambda[index]-1)*(4*lambda[index]-2)*(4*lambda[index]-3) + 4*lambda[index]*(4*lambda[index]-2)*(4*lambda[index]-3) + 4*lambda[index]*(4*lambda[index]-1)*(4*lambda[index]-3) + 4*lambda[index]*(4*lambda[index]-1)*(4*lambda[index]-2))*eta[index]/(12.0*s);
			phi[1]=-((4*lambda[index]-1)*(4*lambda[index]-2)*(4*lambda[index]-3) + 4*lambda[index]*(4*lambda[index]-2)*(4*lambda[index]-3) + 4*lambda[index]*(4*lambda[index]-1)*(4*lambda[index]-3) + 4*lambda[index]*(4*lambda[index]-1)*(4*lambda[index]-2))*xi[index]/(12.0*s);
		}
		else if(index < 3*dop)
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1);
			i1 = (ie+1)%3;
			i2 = (ie+2)%3;
			switch(ii)
			{
			case 0:
				phi[0]=(eta[i1]*lambda[i2]*(4*lambda[i1]-1)*(4*lambda[i1]-2) + lambda[i1]*eta[i2]*(4*lambda[i1]-1)*(4*lambda[i1]-2) + lambda[i1]*lambda[i2]*4*eta[i1]*(4*lambda[i1]-2) + lambda[i1]*lambda[i2]*(4*lambda[i1]-1)*4*eta[i1])*4.0/(3.0*s);
				phi[1]=-(xi[i1]*lambda[i2]*(4*lambda[i1]-1)*(4*lambda[i1]-2) + lambda[i1]*xi[i2]*(4*lambda[i1]-1)*(4*lambda[i1]-2) + lambda[i1]*lambda[i2]*4*xi[i1]*(4*lambda[i1]-2) + lambda[i1]*lambda[i2]*(4*lambda[i1]-1)*4*xi[i1])*4.0/(3.0*s);
				break; 
			case 1:
				phi[0]=(eta[i1]*lambda[i2]*(4*lambda[i1]-1)*(4*lambda[i2]-1) + lambda[i1]*eta[i2]*(4*lambda[i1]-1)*(4*lambda[i2]-1) + lambda[i1]*lambda[i2]*4*eta[i1]*(4*lambda[i2]-1) + lambda[i1]*lambda[i2]*(4*lambda[i1]-1)*4*eta[i2])*2.0/s;
				phi[1]=-(xi[i1]*lambda[i2]*(4*lambda[i1]-1)*(4*lambda[i2]-1) + lambda[i1]*xi[i2]*(4*lambda[i1]-1)*(4*lambda[i2]-1) + lambda[i1]*lambda[i2]*4*xi[i1]*(4*lambda[i2]-1) + lambda[i1]*lambda[i2]*(4*lambda[i1]-1)*4*xi[i2])*2.0/s;
				break;
			case 2:
				phi[0]=(eta[i1]*lambda[i2]*(4*lambda[i2]-1)*(4*lambda[i2]-2) + lambda[i1]*eta[i2]*(4*lambda[i2]-1)*(4*lambda[i2]-2) + lambda[i1]*lambda[i2]*4*eta[i2]*(4*lambda[i2]-2) + lambda[i1]*lambda[i2]*(4*lambda[i2]-1)*4*eta[i2])*4.0/(3.0*s);
				phi[1]=-(xi[i1]*lambda[i2]*(4*lambda[i2]-1)*(4*lambda[i2]-2) + lambda[i1]*xi[i2]*(4*lambda[i2]-1)*(4*lambda[i2]-2) + lambda[i1]*lambda[i2]*4*xi[i2]*(4*lambda[i2]-2) + lambda[i1]*lambda[i2]*(4*lambda[i2]-1)*4*xi[i2])*4.0/(3.0*s);
				break;
			default:
				phi[0] = 0;
				phi[1] = 0;
			}
		}
		else
		{
			in = index-3*dop;
			phi[0]=(eta[0]*lambda[1]*lambda[2]*(4*lambda[in]-1) + lambda[0]*eta[1]*lambda[2]*(4*lambda[in]-1) + lambda[0]*lambda[1]*eta[2]*(4*lambda[in]-1) + lambda[0]*lambda[1]*lambda[2]*4*eta[in])*16.0/s;
			phi[1]=-(xi[0]*lambda[1]*lambda[2]*(4*lambda[in]-1) + lambda[0]*xi[1]*lambda[2]*(4*lambda[in]-1) + lambda[0]*lambda[1]*xi[2]*(4*lambda[in]-1) + lambda[0]*lambda[1]*lambda[2]*4*xi[in])*16.0/s;
		}
	} // dop=4

	else if(dop==5)
	{
		if(index<3)
		{
			phi[0]=(eta[index]*(5*lambda[index]-1)*(5*lambda[index]-2)*(5*lambda[index]-3)*(5*lambda[index]-4) + lambda[index]*5*eta[index]*(5*lambda[index]-2)*(5*lambda[index]-3)*(5*lambda[index]-4) + lambda[index]*(5*lambda[index]-1)*5*eta[index]*(5*lambda[index]-3)*(5*lambda[index]-4) + lambda[index]*(5*lambda[index]-1)*(5*lambda[index]-2)*5*eta[index]*(5*lambda[index]-4) + lambda[index]*(5*lambda[index]-1)*(5*lambda[index]-2)*(5*lambda[index]-3)*5*eta[index])/(48.0*s);
			phi[1]=-(xi[index]*(5*lambda[index]-1)*(5*lambda[index]-2)*(5*lambda[index]-3)*(5*lambda[index]-4) + lambda[index]*5*xi[index]*(5*lambda[index]-2)*(5*lambda[index]-3)*(5*lambda[index]-4) + lambda[index]*(5*lambda[index]-1)*5*xi[index]*(5*lambda[index]-3)*(5*lambda[index]-4) + lambda[index]*(5*lambda[index]-1)*(5*lambda[index]-2)*5*xi[index]*(5*lambda[index]-4) + lambda[index]*(5*lambda[index]-1)*(5*lambda[index]-2)*(5*lambda[index]-3)*5*xi[index])/(48.0*s);
		}
		else if(index < 3*dop)
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1);
			i1 = (ie+1)%3;
			i2 = (ie+2)%3;
			switch(ii)
			{
			case 0:
				phi[0]=(eta[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i1]-2)*(5*lambda[i1]-3) + lambda[i1]*eta[i2]*(5*lambda[i1]-1)*(5*lambda[i1]-2)*(5*lambda[i1]-3) + lambda[i1]*lambda[i2]*5*eta[i1]*(5*lambda[i1]-2)*(5*lambda[i1]-3) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*5*eta[i1]*(5*lambda[i1]-3) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i1]-2)*5*eta[i1])*25.0/(48.0*s);
				phi[1]=-(xi[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i1]-2)*(5*lambda[i1]-3) + lambda[i1]*xi[i2]*(5*lambda[i1]-1)*(5*lambda[i1]-2)*(5*lambda[i1]-3) + lambda[i1]*lambda[i2]*5*xi[i1]*(5*lambda[i1]-2)*(5*lambda[i1]-3) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*5*xi[i1]*(5*lambda[i1]-3) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i1]-2)*5*xi[i1])*25.0/(48.0*s);
				break; 
			case 1:
				phi[0]=(eta[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*(5*lambda[i1]-2) + lambda[i1]*eta[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*(5*lambda[i1]-2) + lambda[i1]*lambda[i2]*5*eta[i1]*(5*lambda[i2]-1)*(5*lambda[i1]-2) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*5*eta[i2]*(5*lambda[i1]-2) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*5*eta[i1])*25.0/(24.0*s);
				phi[1]=-(xi[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*(5*lambda[i1]-2) + lambda[i1]*xi[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*(5*lambda[i1]-2) + lambda[i1]*lambda[i2]*5*xi[i1]*(5*lambda[i2]-1)*(5*lambda[i1]-2) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*5*xi[i2]*(5*lambda[i1]-2) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*5*xi[i1])*25.0/(24.0*s);
				break;
			case 2:
				phi[0]=(eta[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*(5*lambda[i2]-2) + lambda[i1]*eta[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*(5*lambda[i2]-2) + lambda[i1]*lambda[i2]*5*eta[i1]*(5*lambda[i2]-1)*(5*lambda[i2]-2) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*5*eta[i2]*(5*lambda[i2]-2) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*5*eta[i2])*25.0/(24.0*s);
				phi[1]=-(xi[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*(5*lambda[i2]-2) + lambda[i1]*xi[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*(5*lambda[i2]-2) + lambda[i1]*lambda[i2]*5*xi[i1]*(5*lambda[i2]-1)*(5*lambda[i2]-2) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*5*xi[i2]*(5*lambda[i2]-2) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*5*xi[i2])*25.0/(24.0*s);
				break;
			case 3:
				phi[0]=(eta[i1]*lambda[i2]*(5*lambda[i2]-1)*(5*lambda[i2]-2)*(5*lambda[i2]-3) + lambda[i1]*eta[i2]*(5*lambda[i2]-1)*(5*lambda[i2]-2)*(5*lambda[i2]-3) + lambda[i1]*lambda[i2]*5*eta[i2]*(5*lambda[i2]-2)*(5*lambda[i2]-3) + lambda[i1]*lambda[i2]*(5*lambda[i2]-1)*5*eta[i2]*(5*lambda[i2]-3) + lambda[i1]*lambda[i2]*(5*lambda[i2]-1)*(5*lambda[i2]-2)*5*eta[i2])*25.0/(48.0*s);
				phi[1]=-(xi[i1]*lambda[i2]*(5*lambda[i2]-1)*(5*lambda[i2]-2)*(5*lambda[i2]-3) + lambda[i1]*xi[i2]*(5*lambda[i2]-1)*(5*lambda[i2]-2)*(5*lambda[i2]-3) + lambda[i1]*lambda[i2]*5*xi[i2]*(5*lambda[i2]-2)*(5*lambda[i2]-3) + lambda[i1]*lambda[i2]*(5*lambda[i2]-1)*5*xi[i2]*(5*lambda[i2]-3) + lambda[i1]*lambda[i2]*(5*lambda[i2]-1)*(5*lambda[i2]-2)*5*xi[i2])*25.0/(48.0*s);
				break;
			default:
				phi[0] = 0;
				phi[1] = 0;
			}
		}
		else if(index < 3*dop+3)
		{
			in = index-3*dop;
			phi[0]=(eta[0]*lambda[1]*lambda[2]*(5*lambda[in]-1)*(5*lambda[in]-2) + lambda[0]*eta[1]*lambda[2]*(5*lambda[in]-1)*(5*lambda[in]-2) + lambda[0]*lambda[1]*eta[2]*(5*lambda[in]-1)*(5*lambda[in]-2) + lambda[0]*lambda[1]*lambda[2]*5*eta[in]*(5*lambda[in]-2) + lambda[0]*lambda[1]*lambda[2]*(5*lambda[in]-1)*5*eta[in])*125.0/(12.0*s);
			phi[1]=-(xi[0]*lambda[1]*lambda[2]*(5*lambda[in]-1)*(5*lambda[in]-2) + lambda[0]*xi[1]*lambda[2]*(5*lambda[in]-1)*(5*lambda[in]-2) + lambda[0]*lambda[1]*xi[2]*(5*lambda[in]-1)*(5*lambda[in]-2) + lambda[0]*lambda[1]*lambda[2]*5*xi[in]*(5*lambda[in]-2) + lambda[0]*lambda[1]*lambda[2]*(5*lambda[in]-1)*5*xi[in])*125.0/(12.0*s);
		}
		else
		{
			in = index-3*dop-3;
			i1 = (in+1)%3;
			i2 = (in+2)%3;
			phi[0]=(eta[0]*lambda[1]*lambda[2]*(5*lambda[i1]-1)*(5*lambda[i2]-1) + lambda[0]*eta[1]*lambda[2]*(5*lambda[i1]-1)*(5*lambda[i2]-1) + lambda[0]*lambda[1]*eta[2]*(5*lambda[i1]-1)*(5*lambda[i2]-1) + lambda[0]*lambda[1]*lambda[2]*5*eta[i1]*(5*lambda[i2]-1) + lambda[0]*lambda[1]*lambda[2]*(5*lambda[i1]-1)*5*eta[i2])*125.0/(8.0*s);
			phi[1]=-(xi[0]*lambda[1]*lambda[2]*(5*lambda[i1]-1)*(5*lambda[i2]-1) + lambda[0]*xi[1]*lambda[2]*(5*lambda[i1]-1)*(5*lambda[i2]-1) + lambda[0]*lambda[1]*xi[2]*(5*lambda[i1]-1)*(5*lambda[i2]-1) + lambda[0]*lambda[1]*lambda[2]*5*xi[i1]*(5*lambda[i2]-1) + lambda[0]*lambda[1]*lambda[2]*(5*lambda[i1]-1)*5*xi[i2])*125.0/(8.0*s);
		}
	} // dop=5
}

/** 
 * \fn void lagrange_basis2(double *lambda, double s, double eta[3], double xi[3], int index, int dop, double phi[3])
 * \brief the second order derivative of Lagrange element basis function: (\partial_{xx}phi, \partial_{yy}phi, \partial_{xy}phi)
 * \param *lambda pointer to the area coordiante
 * \param s the area of the triangule
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param phi[2] the first order derivative of Morley element basis function: (\partial_{xx}phi, \partial_{yy}phi, \partial_{xy}phi)
 * \return void
 */
void lagrange_basis2(double *lambda, double s, double eta[3], double xi[3], int index, int dop, double phi[3])
{
	int in, ie, ii, i1, i2, i3;
	int dofs = (dop+1)*(dop+2)/2; // degrees of freedom
	if(index>= dofs || index<0)
	{
		phi[0]=0;
		phi[1]=0;
		phi[2]=0;
		return;
	}

	if(dop==0)
	{
		phi[0]=0;
		phi[1]=0;
		phi[2]=0;
	}

	if(dop==1)
	{
		phi[0]=0;
		phi[1]=0;
		phi[2]=0;
	}

	else if(dop==2)
	{
		if(index<3)
		{
			phi[0]=eta[index]*eta[index]/(s*s);
			phi[1]=xi[index]*xi[index]/(s*s);
			phi[2]=-eta[index]*xi[index]/(s*s);
		}
		else
		{
			i1=(index+1)%3;
			i2=(index+2)%3;
			phi[0]=2*eta[i1]*eta[i2]/(s*s);
			phi[1]=2*xi[i1]*xi[i2]/(s*s);
			phi[2]=-(eta[i1]*xi[i2]+eta[i2]*xi[i1])/(s*s);
		}
	} // dop=2

	else if(dop==3)
	{
		if(index<3)
		{
			phi[0]=9*(3*lambda[index]-1)*eta[index]*eta[index]/(4*s*s);
			phi[1]=9*(3*lambda[index]-1)*xi[index]*xi[index]/(4*s*s);
			phi[2]=-9*(3*lambda[index]-1)*eta[index]*xi[index]/(4*s*s);
		}
		else if(index < 3*dop) 
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1);
			i1 = (ie+1)%3;
			i2 = (ie+2)%3;
			i3 = (ie+1+ii)%3;
			phi[0]=9*(3*lambda[i1]*eta[i2]*eta[i3]+3*lambda[i2]*eta[i3]*eta[i1]+(3*lambda[i3]-1)*eta[i1]*eta[i2])/(4*s*s);
			phi[1]=9*(3*lambda[i1]*xi[i2]*xi[i3]+3*lambda[i2]*xi[i3]*xi[i1]+(3*lambda[i3]-1)*xi[i1]*xi[i2])/(4*s*s);
			phi[2]=-9*(3*lambda[i1]*(eta[i2]*xi[i3]+eta[i3]*xi[i2]) + 3*lambda[i2]*(eta[i3]*xi[i1]+eta[i1]*xi[i3]) + (3*lambda[i3]-1)*(eta[i1]*xi[i2]+eta[i2]*xi[i1]))/(8*s*s);
		}
		else
		{
			phi[0]=27*(lambda[0]*eta[1]*eta[2]+lambda[1]*eta[2]*eta[0]+lambda[2]*eta[0]*eta[1])/(2*s*s);
			phi[1]=27*(lambda[0]*xi[1]*xi[2]+lambda[1]*xi[2]*xi[0]+lambda[2]*xi[0]*xi[1])/(2*s*s);
			phi[2]=-27*(lambda[0]*(eta[1]*xi[2]+eta[2]*xi[1])+lambda[1]*(eta[2]*xi[0]+eta[0]*xi[2])+lambda[2]*(eta[0]*xi[1]+eta[1]*xi[0]))/(4*s*s);
		}
	} // dop=3

	else if(dop==4)
	{
		if(index<3)
		{
			phi[0]=(128*lambda[index]*lambda[index]-96*lambda[index]+44.0/3.0)*eta[index]*eta[index]/(4*s*s);
			phi[1]=(128*lambda[index]*lambda[index]-96*lambda[index]+44.0/3.0)*xi[index]*xi[index]/(4*s*s);
			phi[2]=-(128*lambda[index]*lambda[index]-96*lambda[index]+44.0/3.0)*eta[index]*xi[index]/(4*s*s);
		}
		else if(index < 3*dop)
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1);
			i1 = (ie+1)%3;
			i2 = (ie+2)%3;
			switch(ii)
			{
			case 0:
				phi[0]=(32.0/3.0*(24*lambda[i1]*lambda[i1]-12*lambda[i1]+1)*eta[i1]*eta[i2] + 64*lambda[i2]*(4*lambda[i1]-1)*eta[i1]*eta[i1])/(4*s*s);
				phi[1]=(32.0/3.0*(24*lambda[i1]*lambda[i1]-12*lambda[i1]+1)*xi[i1]*xi[i2] + 64*lambda[i2]*(4*lambda[i1]-1)*xi[i1]*xi[i1])/(4*s*s);
				phi[2]=-(16.0/3.0*(24*lambda[i1]*lambda[i1]-12*lambda[i1]+1)*(eta[i1]*xi[i2]+eta[i2]*xi[i1]) + 64*lambda[i2]*(4*lambda[i1]-1)*eta[i1]*xi[i1])/(4*s*s);
				break; 
			case 1:
				phi[0]=(32*(4*lambda[i2]*lambda[i2]-lambda[i2])*eta[i1]*eta[i1] + 8*(8*lambda[i1]-1)*(8*lambda[i2]-1)*eta[i1]*eta[i2] + 32*(4*lambda[i1]*lambda[i1]-lambda[i1])*eta[i2]*eta[i2])/(4*s*s);
				phi[1]=(32*(4*lambda[i2]*lambda[i2]-lambda[i2])*xi[i1]*xi[i1] + 8*(8*lambda[i1]-1)*(8*lambda[i2]-1)*xi[i1]*xi[i2] + 32*(4*lambda[i1]*lambda[i1]-lambda[i1])*xi[i2]*xi[i2])/(4*s*s);
				phi[2]=-(32*(4*lambda[i2]*lambda[i2]-lambda[i2])*eta[i1]*xi[i1] + 4*(8*lambda[i1]-1)*(8*lambda[i2]-1)*(eta[i1]*xi[i2]+eta[i2]*xi[i1]) + 32*(4*lambda[i1]*lambda[i1]-lambda[i1])*eta[i2]*xi[i2])/(4*s*s);
				break;
			case 2:
				phi[0]=(32.0/3.0*(24*lambda[i2]*lambda[i2]-12*lambda[i2]+1)*eta[i1]*eta[i2] + 64*lambda[i1]*(4*lambda[i2]-1)*eta[i2]*eta[i2])/(4*s*s);
				phi[1]=(32.0/3.0*(24*lambda[i2]*lambda[i2]-12*lambda[i2]+1)*xi[i1]*xi[i2] + 64*lambda[i1]*(4*lambda[i2]-1)*xi[i2]*xi[i2])/(4*s*s);
				phi[2]=-(16.0/3.0*(24*lambda[i2]*lambda[i2]-12*lambda[i2]+1)*(eta[i1]*xi[i2]+eta[i2]*xi[i1]) + 64*lambda[i1]*(4*lambda[i2]-1)*eta[i2]*xi[i2])/(4*s*s);
				break;
			default:
				phi[0] = 0;
				phi[1] = 0;
				phi[2] = 0;
			}
		}
		else
		{
			in = index-3*dop;
			phi[0]=(eta[1]*eta[2]*lambda[0]*(4*lambda[in]-1) + eta[2]*eta[0]*lambda[1]*(4*lambda[in]-1) + eta[0]*eta[1]*lambda[2]*(4*lambda[in]-1) + 4*eta[in]*eta[0]*lambda[1]*lambda[2] + 4*eta[in]*eta[1]*lambda[2]*lambda[0] + 4*eta[in]*eta[2]*lambda[0]*lambda[1])*16.0/(s*s);
			phi[1]=(xi[1]*xi[2]*lambda[0]*(4*lambda[in]-1) + xi[2]*xi[0]*lambda[1]*(4*lambda[in]-1) + xi[0]*xi[1]*lambda[2]*(4*lambda[in]-1) + 4*xi[in]*xi[0]*lambda[1]*lambda[2] + 4*xi[in]*xi[1]*lambda[2]*lambda[0] + 4*xi[in]*xi[2]*lambda[0]*lambda[1])*16.0/(s*s);
			phi[2]=-((eta[1]*xi[2]+eta[2]*xi[1])*lambda[0]*(4*lambda[in]-1) + (eta[2]*xi[0]+eta[0]*xi[2])*lambda[1]*(4*lambda[in]-1) + (eta[0]*xi[1]+eta[1]*xi[0])*lambda[2]*(4*lambda[in]-1) + 4*(eta[in]*xi[0]+eta[0]*xi[in])*lambda[1]*lambda[2] + 4*(eta[in]*xi[1]+eta[1]*xi[in])*lambda[2]*lambda[0] + 4*(eta[in]*xi[2]+eta[2]*xi[in])*lambda[0]*lambda[1])*8.0/(s*s);
		}
	} // dop=4

	else if(dop==5)
	{
		if(index<3)
		{
			phi[0]=(3125.0/6.0*lambda[index]*lambda[index]*lambda[index]-625*lambda[index]*lambda[index]+875.0/4.0*lambda[index]-125.0/6.0)*eta[index]*eta[index]/(4*s*s);
			phi[1]=(3125.0/6.0*lambda[index]*lambda[index]*lambda[index]-625*lambda[index]*lambda[index]+875.0/4.0*lambda[index]-125.0/6.0)*xi[index]*xi[index]/(4*s*s);
			phi[2]=-(3125.0/6.0*lambda[index]*lambda[index]*lambda[index]-625*lambda[index]*lambda[index]+875.0/4.0*lambda[index]-125.0/6.0)*eta[index]*xi[index]/(4*s*s);
		}
		else if(index < 3*dop)
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1);
			i1 = (ie+1)%3;
			i2 = (ie+2)%3;
			switch(ii)
			{
			case 0:
				phi[0]= 25.0/12.0*((500*lambda[i1]*lambda[i1]*lambda[i1]-450*lambda[i1]*lambda[i1]+110*lambda[i1]-6)*eta[i1]*eta[i2] + (750*lambda[i1]*lambda[i1]-450*lambda[i1]+55)*lambda[i2]*eta[i1]*eta[i1])/(4*s*s);
				phi[1]= 25.0/12.0*((500*lambda[i1]*lambda[i1]*lambda[i1]-450*lambda[i1]*lambda[i1]+110*lambda[i1]-6)*xi[i1]*xi[i2] + (750*lambda[i1]*lambda[i1]-450*lambda[i1]+55)*lambda[i2]*xi[i1]*xi[i1])/(4*s*s);
				phi[2]= -25.0/12.0*((250*lambda[i1]*lambda[i1]*lambda[i1]-225*lambda[i1]*lambda[i1]+55*lambda[i1]-3)*(eta[i1]*xi[i2]+eta[i2]*xi[i1]) + (750*lambda[i1]*lambda[i1]-450*lambda[i1]+55)*lambda[i2]*eta[i1]*xi[i1])/(4*s*s);
				break; 
			case 1:
				phi[0]= 25.0/12.0*((150*lambda[i1]-30)*(5*lambda[i2]*lambda[i2]-lambda[i2])*eta[i1]*eta[i1] + (75*lambda[i1]*lambda[i1]-30*lambda[i1]+2)*(10*lambda[i2]-1)*2*eta[i1]*eta[i2] + (25*lambda[i1]*lambda[i1]*lambda[i1]-15*lambda[i1]*lambda[i1]+2*lambda[i1])*10*eta[i2]*eta[i2])/(4*s*s);
				phi[1]= 25.0/12.0*((150*lambda[i1]-30)*(5*lambda[i2]*lambda[i2]-lambda[i2])*xi[i1]*xi[i1] + (75*lambda[i1]*lambda[i1]-30*lambda[i1]+2)*(10*lambda[i2]-1)*2*xi[i1]*xi[i2] + (25*lambda[i1]*lambda[i1]*lambda[i1]-15*lambda[i1]*lambda[i1]+2*lambda[i1])*10*xi[i2]*xi[i2])/(4*s*s);
				phi[2]= -25.0/12.0*((150*lambda[i1]-30)*(5*lambda[i2]*lambda[i2]-lambda[i2])*eta[i1]*xi[i1] + (75*lambda[i1]*lambda[i1]-30*lambda[i1]+2)*(10*lambda[i2]-1)*(eta[i1]*xi[i2]+eta[i2]*xi[i1]) + (25*lambda[i1]*lambda[i1]*lambda[i1]-15*lambda[i1]*lambda[i1]+2*lambda[i1])*10*eta[i2]*xi[i2])/(4*s*s);
				break;
			case 2:
				phi[0]= 25.0/12.0*((150*lambda[i2]-30)*(5*lambda[i1]*lambda[i1]-lambda[i1])*eta[i2]*eta[i2] + (75*lambda[i2]*lambda[i2]-30*lambda[i2]+2)*(10*lambda[i1]-1)*2*eta[i1]*eta[i2] + (25*lambda[i2]*lambda[i2]*lambda[i2]-15*lambda[i2]*lambda[i2]+2*lambda[i2])*10*eta[i1]*eta[i1])/(4*s*s);
				phi[1]= 25.0/12.0*((150*lambda[i2]-30)*(5*lambda[i1]*lambda[i1]-lambda[i1])*xi[i2]*xi[i2] + (75*lambda[i2]*lambda[i2]-30*lambda[i2]+2)*(10*lambda[i1]-1)*2*xi[i1]*xi[i2] + (25*lambda[i2]*lambda[i2]*lambda[i2]-15*lambda[i2]*lambda[i2]+2*lambda[i2])*10*xi[i1]*xi[i1])/(4*s*s);
				phi[2]= -25.0/12.0*((150*lambda[i2]-30)*(5*lambda[i1]*lambda[i1]-lambda[i1])*eta[i2]*xi[i2] + (75*lambda[i2]*lambda[i2]-30*lambda[i2]+2)*(10*lambda[i1]-1)*(eta[i1]*xi[i2]+eta[i2]*xi[i1]) + (25*lambda[i2]*lambda[i2]*lambda[i2]-15*lambda[i2]*lambda[i2]+2*lambda[i2])*10*eta[i1]*xi[i1])/(4*s*s);
				break;
			case 3:
				phi[0]= 25.0/12.0*((500*lambda[i2]*lambda[i2]*lambda[i2]-450*lambda[i2]*lambda[i2]+110*lambda[i2]-6)*eta[i1]*eta[i2] + (750*lambda[i2]*lambda[i2]-450*lambda[i2]+55)*lambda[i1]*eta[i2]*eta[i2])/(4*s*s);
				phi[1]= 25.0/12.0*((500*lambda[i2]*lambda[i2]*lambda[i2]-450*lambda[i2]*lambda[i2]+110*lambda[i2]-6)*xi[i1]*xi[i2] + (750*lambda[i2]*lambda[i2]-450*lambda[i2]+55)*lambda[i1]*xi[i2]*xi[i2])/(4*s*s);
				phi[2]= -25.0/12.0*((250*lambda[i2]*lambda[i2]*lambda[i2]-225*lambda[i2]*lambda[i2]+55*lambda[i2]-3)*(eta[i1]*xi[i2]+eta[i2]*xi[i1]) + (750*lambda[i2]*lambda[i2]-450*lambda[i2]+55)*lambda[i1]*eta[i2]*xi[i2])/(4*s*s);
				break;
			default:
				phi[0] = 0;
				phi[1] = 0;
				phi[2] = 0;
			}
		}
		else if(index < 3*dop+3)
		{
			in = index-3*dop;
			phi[0]= 125.0/6.0*(2*(eta[0]*eta[1]*lambda[2]+eta[1]*eta[2]*lambda[0]+eta[2]*eta[0]*lambda[1])*(25*lambda[in]*lambda[in]-15*lambda[in]+2) + 2*(eta[0]*lambda[1]*lambda[2]+eta[1]*lambda[2]*lambda[0]+eta[2]*lambda[0]*lambda[1])*(50*lambda[in]-15)*eta[in] + lambda[0]*lambda[1]*lambda[2]*50*eta[in]*eta[in])/(4*s*s);
			phi[1]= 125.0/6.0*(2*(xi[0]*xi[1]*lambda[2]+xi[1]*xi[2]*lambda[0]+xi[2]*xi[0]*lambda[1])*(25*lambda[in]*lambda[in]-15*lambda[in]+2) + 2*(xi[0]*lambda[1]*lambda[2]+xi[1]*lambda[2]*lambda[0]+xi[2]*lambda[0]*lambda[1])*(50*lambda[in]-15)*xi[in] + lambda[0]*lambda[1]*lambda[2]*50*xi[in]*xi[in])/(4*s*s);
			phi[2]= -125.0/6.0*(((eta[0]*xi[1]+eta[1]*xi[0])*lambda[2]+(eta[1]*xi[2]+eta[2]*xi[1])*lambda[0]+(eta[2]*xi[0]+eta[0]*xi[2])*lambda[1])*(25*lambda[in]*lambda[in]-15*lambda[in]+2) + (xi[0]*lambda[1]*lambda[2]+xi[1]*lambda[2]*lambda[0]+xi[2]*lambda[0]*lambda[1])*(50*lambda[in]-15)*eta[in] + (eta[0]*lambda[1]*lambda[2]+eta[1]*lambda[2]*lambda[0]+eta[2]*lambda[0]*lambda[1])*(50*lambda[in]-15)*xi[in] + lambda[0]*lambda[1]*lambda[2]*50*eta[in]*xi[in])/(4*s*s);
		}
		else
		{
			in = index-3*dop-3;
			i1 = (in+1)%3;
			i2 = (in+2)%3;
			phi[0]= 125.0/4.0*(2*(eta[0]*eta[1]*lambda[2]+eta[1]*eta[2]*lambda[0]+eta[2]*eta[0]*lambda[1])*(5*lambda[i1]-1)*(5*lambda[i2]-1) + 2*(eta[0]*lambda[1]*lambda[2]+eta[1]*lambda[2]*lambda[0]+eta[2]*lambda[0]*lambda[1])*((25*lambda[i2]-5)*eta[i1]+(25*lambda[i1]-5)*eta[i2]) + lambda[0]*lambda[1]*lambda[2]*50*eta[i1]*eta[i2])/(4*s*s);
			phi[1]= 125.0/4.0*(2*(xi[0]*xi[1]*lambda[2]+xi[1]*xi[2]*lambda[0]+xi[2]*xi[0]*lambda[1])*(5*lambda[i1]-1)*(5*lambda[i2]-1) + 2*(xi[0]*lambda[1]*lambda[2]+xi[1]*lambda[2]*lambda[0]+xi[2]*lambda[0]*lambda[1])*((25*lambda[i2]-5)*xi[i1]+(25*lambda[i1]-5)*xi[i2]) + lambda[0]*lambda[1]*lambda[2]*50*xi[i1]*xi[i2])/(4*s*s);
			phi[2]= -125.0/4.0*(((eta[0]*xi[1]+eta[1]*xi[0])*lambda[2]+(eta[1]*xi[2]+eta[2]*xi[1])*lambda[0]+(eta[2]*xi[0]+eta[0]*xi[2])*lambda[1])*(5*lambda[i1]-1)*(5*lambda[i2]-1) + (xi[0]*lambda[1]*lambda[2]+xi[1]*lambda[2]*lambda[0]+xi[2]*lambda[0]*lambda[1])*((25*lambda[i2]-5)*eta[i1]+(25*lambda[i1]-5)*eta[i2]) + (eta[0]*lambda[1]*lambda[2]+eta[1]*lambda[2]*lambda[0]+eta[2]*lambda[0]*lambda[1])*((25*lambda[i2]-5)*xi[i1]+(25*lambda[i1]-5)*xi[i2]) + lambda[0]*lambda[1]*lambda[2]*25*(eta[i1]*xi[i2]+eta[i2]*xi[i1]))/(4*s*s);
		}
	} // dop=5
}


/** 
 * \fn void rt_basis(double x, double y, double (*T)[2], double s, double elen[3], double eta[3], double xi[3], double orient[3], int index, int dop, double phi[2])
 * \brief basis function of Raviart-Thomas element: (phi1, phi2)
 * \param x the x-axis coordiante
 * \param y the y-axis coordiante
 * \param (*T)[2] point the coordiantes of all vertices of current element
 * \param s the area of the triangule
 * \param elen[3] length of three edge
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param orient[3] some auxiliary parameter
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param phi[2] basis function of Raviart-Thomas element: (phi1, phi2)
 * \return void
 */
void rt_basis(double x, double y, double (*T)[2], double s, double elen[3], double eta[3], double xi[3], double orient[3], int index, int dop, double phi[2])
{
	double c1, c2, d1, d2, d3, a1, a2, a3, b1, b2, b3, aa1, aa2, aa3, bb1, bb2, bb3;
	double x1 = T[0][0];
	double x2 = T[1][0];
	double x3 = T[2][0];
	double y1 = T[0][1];
	double y2 = T[1][1];
	double y3 = T[2][1];
	if(dop==1)
	{
		if(index==0)
		{
			c1 = -elen[0]*eta[0]/(s*s);
			c2 = elen[0]*xi[0]/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1;
			a2 = -d2*x2-elen[0]*xi[2]/(2*s);
			a3 = -d3*x3+elen[0]*xi[1]/(2*s);
			b1 = -d1*y1;
			b2 = -d2*y2-elen[0]*eta[2]/(2*s);
			b3 = -d3*y3+elen[0]*eta[1]/(2*s);
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
			aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
			bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
			c1 *= orient[0];
			c2 *= orient[0];
			aa1 *= orient[0];
			aa2 *= orient[0];
			aa3 *= orient[0];
			bb1 *= orient[0];
			bb2 *= orient[0];
			bb3 *= orient[0];
		}
		else if(index==1)
		{
			c1 = -2*elen[0]*(eta[1]-eta[2])/(s*s);
			c2 = 2*elen[0]*(xi[1]-xi[2])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1;
			a2 = -d2*x2+3*elen[0]*xi[2]/s;
			a3 = -d3*x3+3*elen[0]*xi[1]/s;
			b1 = -d1*y1;
			b2 = -d2*y2+3*elen[0]*eta[2]/s;
			b3 = -d3*y3+3*elen[0]*eta[1]/s;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
			aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
			bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else if(index==2)
		{
			c1 = -elen[1]*eta[1]/(s*s);
			c2 = elen[1]*xi[1]/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1+elen[1]*xi[2]/(2*s);
			a2 = -d2*x2;
			a3 = -d3*x3-elen[1]*xi[0]/(2*s);
			b1 = -d1*y1+elen[1]*eta[2]/(2*s);
			b2 = -d2*y2;
			b3 = -d3*y3-elen[1]*eta[0]/(2*s);
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
			aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
			bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
			c1 *= orient[1];
			c2 *= orient[1];
			aa1 *= orient[1];
			aa2 *= orient[1];
			aa3 *= orient[1];
			bb1 *= orient[1];
			bb2 *= orient[1];
			bb3 *= orient[1];
		}
		else if(index==3)
		{
			c1 = -2*elen[1]*(eta[2]-eta[0])/(s*s);
			c2 = 2*elen[1]*(xi[2]-xi[0])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1+3*elen[1]*xi[2]/s;
			a2 = -d2*x2;
			a3 = -d3*x3+3*elen[1]*xi[0]/s;
			b1 = -d1*y1+3*elen[1]*eta[2]/s;
			b2 = -d2*y2;
			b3 = -d3*y3+3*elen[1]*eta[0]/s;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
			aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
			bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else if(index==4)
		{
			c1 = -elen[2]*eta[2]/(s*s);
			c2 = elen[2]*xi[2]/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1-elen[2]*xi[1]/(2*s);
			a2 = -d2*x2+elen[2]*xi[0]/(2*s);
			a3 = -d3*x3;
			b1 = -d1*y1-elen[2]*eta[1]/(2*s);
			b2 = -d2*y2+elen[2]*eta[0]/(2*s);
			b3 = -d3*y3;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
			aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
			bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
			c1 *= orient[2];
			c2 *= orient[2];
			aa1 *= orient[2];
			aa2 *= orient[2];
			aa3 *= orient[2];
			bb1 *= orient[2];
			bb2 *= orient[2];
			bb3 *= orient[2];
		}
		else if(index==5)
		{
			c1 = -2*elen[2]*(eta[0]-eta[1])/(s*s);
			c2 = 2*elen[2]*(xi[0]-xi[1])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1+3*elen[2]*xi[1]/s;
			a2 = -d2*x2+3*elen[2]*xi[0]/s;
			a3 = -d3*x3;
			b1 = -d1*y1+3*elen[2]*eta[1]/s;
			b2 = -d2*y2+3*elen[2]*eta[0]/s;
			b3 = -d3*y3;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
			aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
			bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else if(index==6)
		{
			c1 = -(eta[0]*eta[0]+eta[1]*eta[1]+eta[2]*eta[2])/(s*s);
			c2 = (xi[0]*eta[0]+xi[1]*eta[1]+xi[2]*eta[2])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1;
			a2 = -d2*x2;
			a3 = -d3*x3;
			b1 = -d1*y1;
			b2 = -d2*y2;
			b3 = -d3*y3;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
			aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
			bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else if(index==7)
		{
			c1 = (xi[0]*eta[0]+xi[1]*eta[1]+xi[2]*eta[2])/(s*s);
			c2 = -(xi[0]*xi[0]+xi[1]*xi[1]+xi[2]*xi[2])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1;
			a2 = -d2*x2;
			a3 = -d3*x3;
			b1 = -d1*y1;
			b2 = -d2*y2;
			b3 = -d3*y3;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
			aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
			bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else
		{
			c1 = 0;
			c2 = 0;
			d1 = 0;
			d2 = 0;
			d3 = 0;
			a1 = 0;
			a2 = 0;
			a3 = 0;
			b1 = 0;
			b2 = 0;
			b3 = 0;
			aa1 = 0;
			aa2 = 0;
			aa3 = 0;
			bb1 = 0;
			bb2 = 0;
			bb3 = 0;
		}
		phi[0] = aa1*x+aa2*y+aa3+x*(c1*x+c2*y);
		phi[1] = bb1*x+bb2*y+bb3+y*(c1*x+c2*y);
	} // dop=1
	else
	{
		phi[0]=0;
		phi[1]=0;
	}
}

/** 
 * \fn void rt_basis1(double x, double y, double (*T)[2], double s, double elen[3], double eta[3], double xi[3], double orient[3], int index, int dop, double phi[4])
 * \brief the first order derivative of Raviart-Thomas element basis function: (\partial_{x}phi1, \partial_{y}phi1, \partial_{x}phi2, \partial_{y}phi2)
 * \param x the x-axis coordiante
 * \param y the y-axis coordiante
 * \param (*T)[2] point the coordiantes of all vertices of current element
 * \param s the area of the triangule
 * \param elen[3] length of three edge
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param orient[3] some auxiliary parameter
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param phi[4] the first order derivative of Raviart-Thomas element basis function: (\partial_{x}phi1, \partial_{y}phi1, \partial_{x}phi2, \partial_{y}phi2)
 * \return void
 */
void rt_basis1(double x, double y, double (*T)[2], double s, double elen[3], double eta[3], double xi[3], double orient[3], int index, int dop, double phi[4])
{
	double c1, c2, d1, d2, d3, a1, a2, a3, b1, b2, b3, aa1, aa2, aa3, bb1, bb2, bb3;
	double x1 = T[0][0];
	double x2 = T[1][0];
	double x3 = T[2][0];
	double y1 = T[0][1];
	double y2 = T[1][1];
	double y3 = T[2][1];
	if(dop==1)
	{
		if(index==0)
		{
			c1 = -elen[0]*eta[0]/(s*s);
			c2 = elen[0]*xi[0]/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1;
			a2 = -d2*x2-elen[0]*xi[2]/(2*s);
			a3 = -d3*x3+elen[0]*xi[1]/(2*s);
			b1 = -d1*y1;
			b2 = -d2*y2-elen[0]*eta[2]/(2*s);
			b3 = -d3*y3+elen[0]*eta[1]/(2*s);
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
	//		aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
	//		bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
			c1 *= orient[0];
			c2 *= orient[0];
			aa1 *= orient[0];
			aa2 *= orient[0];
	//		aa3 *= orient[0];
			bb1 *= orient[0];
			bb2 *= orient[0];
	//		bb3 *= orient[0];
		}
		else if(index==1)
		{
			c1 = -2*elen[0]*(eta[1]-eta[2])/(s*s);
			c2 = 2*elen[0]*(xi[1]-xi[2])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1;
			a2 = -d2*x2+3*elen[0]*xi[2]/s;
			a3 = -d3*x3+3*elen[0]*xi[1]/s;
			b1 = -d1*y1;
			b2 = -d2*y2+3*elen[0]*eta[2]/s;
			b3 = -d3*y3+3*elen[0]*eta[1]/s;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
	//		aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
	//		bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else if(index==2)
		{
			c1 = -elen[1]*eta[1]/(s*s);
			c2 = elen[1]*xi[1]/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1+elen[1]*xi[2]/(2*s);
			a2 = -d2*x2;
			a3 = -d3*x3-elen[1]*xi[0]/(2*s);
			b1 = -d1*y1+elen[1]*eta[2]/(2*s);
			b2 = -d2*y2;
			b3 = -d3*y3-elen[1]*eta[0]/(2*s);
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
	//		aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
	//		bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
			c1 *= orient[1];
			c2 *= orient[1];
			aa1 *= orient[1];
			aa2 *= orient[1];
	//		aa3 *= orient[1];
			bb1 *= orient[1];
			bb2 *= orient[1];
	//		bb3 *= orient[1];
		}
		else if(index==3)
		{
			c1 = -2*elen[1]*(eta[2]-eta[0])/(s*s);
			c2 = 2*elen[1]*(xi[2]-xi[0])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1+3*elen[1]*xi[2]/s;
			a2 = -d2*x2;
			a3 = -d3*x3+3*elen[1]*xi[0]/s;
			b1 = -d1*y1+3*elen[1]*eta[2]/s;
			b2 = -d2*y2;
			b3 = -d3*y3+3*elen[1]*eta[0]/s;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
	//		aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
	//		bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else if(index==4)
		{
			c1 = -elen[2]*eta[2]/(s*s);
			c2 = elen[2]*xi[2]/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1-elen[2]*xi[1]/(2*s);
			a2 = -d2*x2+elen[2]*xi[0]/(2*s);
			a3 = -d3*x3;
			b1 = -d1*y1-elen[2]*eta[1]/(2*s);
			b2 = -d2*y2+elen[2]*eta[0]/(2*s);
			b3 = -d3*y3;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
	//		aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
	//		bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
			c1 *= orient[2];
			c2 *= orient[2];
			aa1 *= orient[2];
			aa2 *= orient[2];
	//		aa3 *= orient[2];
			bb1 *= orient[2];
			bb2 *= orient[2];
	//		bb3 *= orient[2];
		}
		else if(index==5)
		{
			c1 = -2*elen[2]*(eta[0]-eta[1])/(s*s);
			c2 = 2*elen[2]*(xi[0]-xi[1])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1+3*elen[2]*xi[1]/s;
			a2 = -d2*x2+3*elen[2]*xi[0]/s;
			a3 = -d3*x3;
			b1 = -d1*y1+3*elen[2]*eta[1]/s;
			b2 = -d2*y2+3*elen[2]*eta[0]/s;
			b3 = -d3*y3;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
	//		aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
	//		bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else if(index==6)
		{
			c1 = -(eta[0]*eta[0]+eta[1]*eta[1]+eta[2]*eta[2])/(s*s);
			c2 = (xi[0]*eta[0]+xi[1]*eta[1]+xi[2]*eta[2])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1;
			a2 = -d2*x2;
			a3 = -d3*x3;
			b1 = -d1*y1;
			b2 = -d2*y2;
			b3 = -d3*y3;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
	//		aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
	//		bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else if(index==7)
		{
			c1 = (xi[0]*eta[0]+xi[1]*eta[1]+xi[2]*eta[2])/(s*s);
			c2 = -(xi[0]*xi[0]+xi[1]*xi[1]+xi[2]*xi[2])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1;
			a2 = -d2*x2;
			a3 = -d3*x3;
			b1 = -d1*y1;
			b2 = -d2*y2;
			b3 = -d3*y3;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
	//		aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
	//		bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else
		{
			c1 = 0;
			c2 = 0;
			d1 = 0;
			d2 = 0;
			d3 = 0;
			a1 = 0;
			a2 = 0;
			a3 = 0;
			b1 = 0;
			b2 = 0;
			b3 = 0;
			aa1 = 0;
			aa2 = 0;
	//		aa3 = 0;
			bb1 = 0;
			bb2 = 0;
	//		bb3 = 0;
		}
		phi[0] = aa1+2*c1*x+c2*y;
		phi[1] = aa2+c2*x;
		phi[2] = bb1+c1*y;
		phi[3] = bb2+c1*x+2*c2*y;
	} // dop=1
	else
	{
		phi[0]=0;
		phi[1]=0;
		phi[2]=0;
		phi[3]=0;
	}
}

/** 
 * \fn void arnoldwinther_basis(double *lambda, double *x, double *y, ddenmat3 *basisCoeffs, int element, int index, double *phi)
 * \brief basis function of Arnold-Winther element
 * \param *lambda pointer to the area coordiante
 * \param *x pointer to the horizontal ordinates of three vertices
 * \param *y pointer to the longitudinal ordinates of three vertices
 * \param *basisCoeffs pointer to coefficients of basis functions
 * \param element the current element of triangulation
 * \param index the indicator of the basis function
 * \param *phi basis function
 * \return void
 */
void arnoldwinther_basis(double *lambda, double *x, double *y, ddenmat3 *basisCoeffs, int element, int index, double *phi)
{
	double xx=lambda[0]*x[0]+lambda[1]*x[1]+lambda[2]*x[2];
	double yy=lambda[0]*y[0]+lambda[1]*y[1]+lambda[2]*y[2];
	double *coeffs=basisCoeffs->val[element][index];
	
	phi[0] = coeffs[0]*lambda[0] + coeffs[3]*lambda[1] + coeffs[6]*lambda[2] + coeffs[9]*lambda[1]*lambda[2] + coeffs[12]*lambda[2]*lambda[0] + coeffs[15]*lambda[0]*lambda[1];
	phi[1] = coeffs[1]*lambda[0] + coeffs[4]*lambda[1] + coeffs[7]*lambda[2] + coeffs[10]*lambda[1]*lambda[2] + coeffs[13]*lambda[2]*lambda[0] + coeffs[16]*lambda[0]*lambda[1];
	phi[2] = coeffs[2]*lambda[0] + coeffs[5]*lambda[1] + coeffs[8]*lambda[2] + coeffs[11]*lambda[1]*lambda[2] + coeffs[14]*lambda[2]*lambda[0] + coeffs[17]*lambda[0]*lambda[1];

	phi[0] += coeffs[18]*pow(yy,3) - coeffs[21]*3*xx*yy*yy - coeffs[22]*pow(xx,3)/3.0 - coeffs[23]*xx*xx*yy;
	phi[1] += coeffs[19]*pow(xx,3) - coeffs[20]*3*xx*xx*yy - coeffs[22]*xx*yy*yy - coeffs[23]*pow(yy,3)/3.0;
	phi[2] += coeffs[20]*pow(xx,3) + coeffs[21]*pow(yy,3) + coeffs[22]*xx*xx*yy + coeffs[23]*xx*yy*yy;
}

/** 
 * \fn void arnoldwinther_basisDIV(double *lambda, ddenmat3 *basisCoeffs, int element, double s, double eta[3], double xi[3], int index, double *phi)
 * \brief divergence of basis function of Arnold-Winther element
 * \param *lambda pointer to the area coordiante
 * \param *basisCoeffs pointer to coefficients of basis functions
 * \param element the current element of triangulation
 * \param s the area of the triangule
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param index the indicator of the basis function
 * \param *phi divergence of basis function
 * \return void
 */
void arnoldwinther_basisDIV(double *lambda, ddenmat3 *basisCoeffs, int element, double s, double eta[3], double xi[3], int index, double *phi)
{
	double *coeffs=basisCoeffs->val[element][index];
	double lambdaGrad[3][2];
	int i;
	for(i=0;i<3;i++)
	{
		lambdaGrad[i][0]=eta[i]/(2.0*s);
		lambdaGrad[i][1]=-xi[i]/(2.0*s);
	}

	phi[0] = coeffs[0]*lambdaGrad[0][0] + coeffs[3]*lambdaGrad[1][0] + coeffs[6]*lambdaGrad[2][0];
	phi[0] += coeffs[9]*(lambda[1]*lambdaGrad[2][0]+lambda[2]*lambdaGrad[1][0]) + coeffs[12]*(lambda[2]*lambdaGrad[0][0]+lambda[0]*lambdaGrad[2][0]) + coeffs[15]*(lambda[0]*lambdaGrad[1][0]+lambda[1]*lambdaGrad[0][0]);
	phi[0] += coeffs[2]*lambdaGrad[0][1] + coeffs[5]*lambdaGrad[1][1] + coeffs[8]*lambdaGrad[2][1];
	phi[0] += coeffs[11]*(lambda[1]*lambdaGrad[2][1]+lambda[2]*lambdaGrad[1][1]) + coeffs[14]*(lambda[2]*lambdaGrad[0][1]+lambda[0]*lambdaGrad[2][1]) + coeffs[17]*(lambda[0]*lambdaGrad[1][1]+lambda[1]*lambdaGrad[0][1]);

	phi[1] = coeffs[1]*lambdaGrad[0][1] + coeffs[4]*lambdaGrad[1][1] + coeffs[7]*lambdaGrad[2][1];
	phi[1] += coeffs[10]*(lambda[1]*lambdaGrad[2][1]+lambda[2]*lambdaGrad[1][1]) + coeffs[13]*(lambda[2]*lambdaGrad[0][1]+lambda[0]*lambdaGrad[2][1]) + coeffs[16]*(lambda[0]*lambdaGrad[1][1]+lambda[1]*lambdaGrad[0][1]);
	phi[1] += coeffs[2]*lambdaGrad[0][0] + coeffs[5]*lambdaGrad[1][0] + coeffs[8]*lambdaGrad[2][0];
	phi[1] += coeffs[11]*(lambda[1]*lambdaGrad[2][0]+lambda[2]*lambdaGrad[1][0]) + coeffs[14]*(lambda[2]*lambdaGrad[0][0]+lambda[0]*lambdaGrad[2][0]) + coeffs[17]*(lambda[0]*lambdaGrad[1][0]+lambda[1]*lambdaGrad[0][0]);
}

/** 
 * \fn void arnoldwinther_basisDIV2(ddenmat3 *basisCoeffs, int element, double s, double eta[3], double xi[3], int index, double *phi)
 * \brief divergence of divergence of basis function of Arnold-Winther element
 * \param *basisCoeffs pointer to coefficients of basis functions
 * \param element the current element of triangulation
 * \param s the area of the triangule
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param index the indicator of the basis function
 * \param *phi divergence of divergence of basis function
 * \return void
 */
void arnoldwinther_basisDIV2(ddenmat3 *basisCoeffs, int element, double s, double eta[3], double xi[3], int index, double *phi)
{
	double *coeffs=basisCoeffs->val[element][index];
	double lambdaGrad[3][2];
	int i;
	for(i=0;i<3;i++)
	{
		lambdaGrad[i][0]=eta[i]/(2.0*s);
		lambdaGrad[i][1]=-xi[i]/(2.0*s);
	}

	*phi = coeffs[9]*2*lambdaGrad[1][0]*lambdaGrad[2][0] + coeffs[12]*2*lambdaGrad[2][0]*lambdaGrad[0][0] + coeffs[15]*2*lambdaGrad[0][0]*lambdaGrad[1][0];
	*phi += coeffs[10]*2*lambdaGrad[1][1]*lambdaGrad[2][1] + coeffs[13]*2*lambdaGrad[2][1]*lambdaGrad[0][1] + coeffs[16]*2*lambdaGrad[0][1]*lambdaGrad[1][1];
	*phi += coeffs[11]*2*(lambdaGrad[1][0]*lambdaGrad[2][1]+lambdaGrad[1][1]*lambdaGrad[2][0]);
	*phi += coeffs[14]*2*(lambdaGrad[2][0]*lambdaGrad[0][1]+lambdaGrad[2][1]*lambdaGrad[0][0]);
	*phi += coeffs[17]*2*(lambdaGrad[0][0]*lambdaGrad[1][1]+lambdaGrad[0][1]*lambdaGrad[1][0]);
}

/** 
 * \fn void huzhang_basis(double *lambda, double **nv, double **tv, int index, int dop, double phi[3])
 * \brief basis function of Hu-Zhang element
 * \param *lambda pointer to the area coordiante
 * \param **nv the unit normal vectors of the three edges
 * \param **tv the unit tangential vectors of the three edges
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param *phi basis function
 * \return void
 */
void huzhang_basis(double *lambda, double **nv, double **tv, int index, int dop, double phi[3])
{
	int dofs = (dop+1)*(dop+2)/2; // degrees of freedom

	phi[0]=0;
	phi[1]=0;
	phi[2]=0;
	if(index>= dofs*3 || index<0)
		return;

	double val;
	double nn[3][3], nt[3][3], tt[3][3];
	int i;
	for(i=0;i<3;i++)
	{
		nn[i][0]=nv[i][0]*nv[i][0]; nn[i][1]=nv[i][1]*nv[i][1]; nn[i][2]=nv[i][0]*nv[i][1];
		nt[i][0]=nv[i][0]*tv[i][0]*sqrt(2); nt[i][1]=nv[i][1]*tv[i][1]*sqrt(2); nt[i][2]=(nv[i][0]*tv[i][1]+nv[i][1]*tv[i][0])/sqrt(2);
		tt[i][0]=tv[i][0]*tv[i][0]; tt[i][1]=tv[i][1]*tv[i][1]; tt[i][2]=tv[i][0]*tv[i][1];
	}

	if(dop==0)
	{
		phi[index] = 1;
	}
	
	else // dop>=1
	{
		if(index<9)
		{
			lagrange_basis(lambda, index%3, dop, &val);
			phi[index/3] = val;
		}
		else if(index<9+(dop-1)*3)
		{
			i=index-9;
			lagrange_basis(lambda, 3+i, dop, &val);
			phi[0]=val*nn[i/(dop-1)][0];
			phi[1]=val*nn[i/(dop-1)][1];
			phi[2]=val*nn[i/(dop-1)][2];
		}
		else if(index<9+(dop-1)*6)
		{
			i=index-9-(dop-1)*3;
			lagrange_basis(lambda, 3+i, dop, &val);
			phi[0]=val*nt[i/(dop-1)][0];
			phi[1]=val*nt[i/(dop-1)][1];
			phi[2]=val*nt[i/(dop-1)][2];
		}
		else if(index<9+(dop-1)*9)
		{
			i=index-9-(dop-1)*6;
			lagrange_basis(lambda, 3+i, dop, &val);
			phi[0]=val*tt[i/(dop-1)][0];
			phi[1]=val*tt[i/(dop-1)][1];
			phi[2]=val*tt[i/(dop-1)][2];
		}
		else
		{
			i=index-dop*9;
			lagrange_basis(lambda, dop*3+i%(dofs-dop*3), dop, &val);
			phi[i/(dofs-dop*3)]=val;
		}
	} 
	
}

/**
* \fn void huzhang_basis1(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, int index, int dop, double(*phi)[2])
* \brief the first order derivative of basis function of Hu-Zhang element
* \param *lambda pointer to the area coordiante
* \param s the area of the triangule
* \param eta[3] some auxiliary parameter
* \param xi[3] some auxiliary parameter
* \param **nv the unit normal vectors of the three edges
* \param **tv the unit tangential vectors of the three edges
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param double(*phi)[2] the first order derivative of basis function of Hu-Zhang element
* \return void
*/
void huzhang_basis1(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, int index, int dop, double(*phi)[2])
{
	int dofs = (dop + 1)*(dop + 2) / 2; // degrees of freedom

	phi[0][0] = 0; phi[0][1] = 0;
	phi[1][0] = 0; phi[1][1] = 0;
	phi[2][0] = 0; phi[2][1] = 0;
	if (index >= dofs * 3 || index<0)
		return;

	double val[2];
	double nn[3][3], nt[3][3], tt[3][3];
	int i;
	for (i = 0; i<3; i++)
	{
		nn[i][0] = nv[i][0] * nv[i][0]; nn[i][1] = nv[i][1] * nv[i][1]; nn[i][2] = nv[i][0] * nv[i][1];
		nt[i][0] = nv[i][0] * tv[i][0] * sqrt(2); nt[i][1] = nv[i][1] * tv[i][1] * sqrt(2); nt[i][2] = (nv[i][0] * tv[i][1] + nv[i][1] * tv[i][0]) / sqrt(2);
		tt[i][0] = tv[i][0] * tv[i][0]; tt[i][1] = tv[i][1] * tv[i][1]; tt[i][2] = tv[i][0] * tv[i][1];
	}

	if (dop == 0)
	{
		return;
	}

	else // dop>=1
	{
		if (index<9)
		{
			lagrange_basis1(lambda, s, eta, xi, index % 3, dop, val);
			phi[index / 3][0] = val[0]; 
			phi[index / 3][1] = val[1];
		}
		else if (index<9 + (dop - 1) * 3)
		{
			i = index - 9;
			lagrange_basis1(lambda, s, eta, xi, 3 + i, dop, val);
			phi[0][0] = val[0] * nn[i / (dop - 1)][0];
			phi[0][1] = val[1] * nn[i / (dop - 1)][0];
			phi[1][0] = val[0] * nn[i / (dop - 1)][1];
			phi[1][1] = val[1] * nn[i / (dop - 1)][1];
			phi[2][0] = val[0] * nn[i / (dop - 1)][2];
			phi[2][1] = val[1] * nn[i / (dop - 1)][2];
		}
		else if (index<9 + (dop - 1) * 6)
		{
			i = index - 9 - (dop - 1) * 3;
			lagrange_basis1(lambda, s, eta, xi, 3 + i, dop, val);
			phi[0][0] = val[0] * nt[i / (dop - 1)][0];
			phi[0][1] = val[1] * nt[i / (dop - 1)][0];
			phi[1][0] = val[0] * nt[i / (dop - 1)][1];
			phi[1][1] = val[1] * nt[i / (dop - 1)][1];
			phi[2][0] = val[0] * nt[i / (dop - 1)][2];
			phi[2][1] = val[1] * nt[i / (dop - 1)][2];
		}
		else if (index<9 + (dop - 1) * 9)
		{
			i = index - 9 - (dop - 1) * 6;
			lagrange_basis1(lambda, s, eta, xi, 3 + i, dop, val);
			phi[0][0] = val[0] * tt[i / (dop - 1)][0];
			phi[0][1] = val[1] * tt[i / (dop - 1)][0];
			phi[1][0] = val[0] * tt[i / (dop - 1)][1];
			phi[1][1] = val[1] * tt[i / (dop - 1)][1];
			phi[2][0] = val[0] * tt[i / (dop - 1)][2];
			phi[2][1] = val[1] * tt[i / (dop - 1)][2];
		}
		else
		{
			i = index - dop * 9;
			lagrange_basis1(lambda, s, eta, xi, dop * 3 + i % (dofs - dop * 3), dop, val);
			phi[i / (dofs - dop * 3)][0] = val[0];
			phi[i / (dofs - dop * 3)][1] = val[1];
		}
	}

}

/** 
 * \fn void huzhang_basisDIV(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, int index, int dop, double phi[2])
 * \brief divergence of basis function of Hu-Zhang element
 * \param *lambda pointer to the area coordiante
 * \param s the area of the triangule
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param **nv the unit normal vectors of the three edges
 * \param **tv the unit tangential vectors of the three edges
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param *phi divergence of basis function
 * \return void
 */
void huzhang_basisDIV(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, int index, int dop, double phi[2])
{
	int dofs = (dop+1)*(dop+2)/2; // degrees of freedom

	phi[0]=0;
	phi[1]=0;
	if(index>= dofs*3 || index<0)
		return;

	double val[2];
	double nn[3][3], nt[3][3], tt[3][3];
	int i;
	for(i=0;i<3;i++)
	{
		nn[i][0]=nv[i][0]*nv[i][0]; nn[i][1]=nv[i][1]*nv[i][1]; nn[i][2]=nv[i][0]*nv[i][1];
		nt[i][0]=nv[i][0]*tv[i][0]*sqrt(2); nt[i][1]=nv[i][1]*tv[i][1]*sqrt(2); nt[i][2]=(nv[i][0]*tv[i][1]+nv[i][1]*tv[i][0])/sqrt(2);
		tt[i][0]=tv[i][0]*tv[i][0]; tt[i][1]=tv[i][1]*tv[i][1]; tt[i][2]=tv[i][0]*tv[i][1];
	}

	if(dop==0)
	{
		return;
	}
	
	else // dop>=1
	{
		if(index<3)
		{
			lagrange_basis1(lambda, s, eta, xi, index, dop, val);
			phi[0] = val[0];
			phi[1] = 0;
		}
		else if(index<6)
		{
			lagrange_basis1(lambda, s, eta, xi, index%3, dop, val);
			phi[0] = 0;
			phi[1] = val[1];
		}
		else if(index<9)
		{
			lagrange_basis1(lambda, s, eta, xi, index%3, dop, val);
			phi[0] = val[1];
			phi[1] = val[0];
		}
		else if(index<9+(dop-1)*3)
		{
			i=index-9;
			lagrange_basis1(lambda, s, eta, xi, 3+i, dop, val);
			phi[0]=val[0]*nn[i/(dop-1)][0]+val[1]*nn[i/(dop-1)][2];
			phi[1]=val[0]*nn[i/(dop-1)][2]+val[1]*nn[i/(dop-1)][1];
		}
		else if(index<9+(dop-1)*6)
		{
			i=index-9-(dop-1)*3;
			lagrange_basis1(lambda, s, eta, xi, 3+i, dop, val);
			phi[0]=val[0]*nt[i/(dop-1)][0]+val[1]*nt[i/(dop-1)][2];
			phi[1]=val[0]*nt[i/(dop-1)][2]+val[1]*nt[i/(dop-1)][1];
		}
		else if(index<9+(dop-1)*9)
		{
			i=index-9-(dop-1)*6;
			lagrange_basis1(lambda, s, eta, xi, 3+i, dop, val);
			phi[0]=val[0]*tt[i/(dop-1)][0]+val[1]*tt[i/(dop-1)][2];
			phi[1]=val[0]*tt[i/(dop-1)][2]+val[1]*tt[i/(dop-1)][1];
		}
		else if(index<dop*9+(dofs-dop*3))
		{
			i=index-dop*9;
			lagrange_basis1(lambda, s, eta, xi, dop*3+i, dop, val);
			phi[0] = val[0];
			phi[1] = 0;
		}
		else if(index<dop*9+(dofs-dop*3)*2)
		{
			i=index-dop*9-(dofs-dop*3);
			lagrange_basis1(lambda, s, eta, xi, dop*3+i, dop, val);
			phi[0] = 0;
			phi[1] = val[1];
		}
		else
		{
			i=index-dop*9-(dofs-dop*3)*2;
			lagrange_basis1(lambda, s, eta, xi, dop*3+i, dop, val);
			phi[0] = val[1];
			phi[1] = val[0];
		}
	} 
}

/**
* \fn void huzhang_basisROT(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, int index, int dop, double phi[2])
* \brief rotation of basis function of Hu-Zhang element
* \param *lambda pointer to the area coordiante
* \param s the area of the triangule
* \param eta[3] some auxiliary parameter
* \param xi[3] some auxiliary parameter
* \param **nv the unit normal vectors of the three edges
* \param **tv the unit tangential vectors of the three edges
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi divergence of basis function
* \return void
*/
void huzhang_basisROT(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, int index, int dop, double phi[2])
{
	int dofs = (dop + 1)*(dop + 2) / 2; // degrees of freedom

	phi[0] = 0;
	phi[1] = 0;
	if (index >= dofs * 3 || index<0)
		return;

	double val[2];
	double nn[3][3], nt[3][3], tt[3][3];
	int i;
	for (i = 0; i<3; i++)
	{
		nn[i][0] = nv[i][0] * nv[i][0]; nn[i][1] = nv[i][1] * nv[i][1]; nn[i][2] = nv[i][0] * nv[i][1];
		nt[i][0] = nv[i][0] * tv[i][0] * sqrt(2); nt[i][1] = nv[i][1] * tv[i][1] * sqrt(2); nt[i][2] = (nv[i][0] * tv[i][1] + nv[i][1] * tv[i][0]) / sqrt(2);
		tt[i][0] = tv[i][0] * tv[i][0]; tt[i][1] = tv[i][1] * tv[i][1]; tt[i][2] = tv[i][0] * tv[i][1];
	}

	if (dop == 0)
	{
		return;
	}

	else // dop>=1
	{
		if (index<3)
		{
			lagrange_basis1(lambda, s, eta, xi, index, dop, val);
			phi[0] = -val[1];
			phi[1] = 0;
		}
		else if (index<6)
		{
			lagrange_basis1(lambda, s, eta, xi, index % 3, dop, val);
			phi[0] = 0;
			phi[1] = val[0];
		}
		else if (index<9)
		{
			lagrange_basis1(lambda, s, eta, xi, index % 3, dop, val);
			phi[0] = val[0];
			phi[1] = -val[1];
		}
		else if (index<9 + (dop - 1) * 3)
		{
			i = index - 9;
			lagrange_basis1(lambda, s, eta, xi, 3 + i, dop, val);
			phi[0] = -val[1] * nn[i / (dop - 1)][0] + val[0] * nn[i / (dop - 1)][2];
			phi[1] = -val[1] * nn[i / (dop - 1)][2] + val[0] * nn[i / (dop - 1)][1];
		}
		else if (index<9 + (dop - 1) * 6)
		{
			i = index - 9 - (dop - 1) * 3;
			lagrange_basis1(lambda, s, eta, xi, 3 + i, dop, val);
			phi[0] = -val[1] * nt[i / (dop - 1)][0] + val[0] * nt[i / (dop - 1)][2];
			phi[1] = -val[1] * nt[i / (dop - 1)][2] + val[0] * nt[i / (dop - 1)][1];
		}
		else if (index<9 + (dop - 1) * 9)
		{
			i = index - 9 - (dop - 1) * 6;
			lagrange_basis1(lambda, s, eta, xi, 3 + i, dop, val);
			phi[0] = -val[1] * tt[i / (dop - 1)][0] + val[0] * tt[i / (dop - 1)][2];
			phi[1] = -val[1] * tt[i / (dop - 1)][2] + val[0] * tt[i / (dop - 1)][1];
		}
		else if (index<dop * 9 + (dofs - dop * 3))
		{
			i = index - dop * 9;
			lagrange_basis1(lambda, s, eta, xi, dop * 3 + i, dop, val);
			phi[0] = -val[1];
			phi[1] = 0;
		}
		else if (index<dop * 9 + (dofs - dop * 3) * 2)
		{
			i = index - dop * 9 - (dofs - dop * 3);
			lagrange_basis1(lambda, s, eta, xi, dop * 3 + i, dop, val);
			phi[0] = 0;
			phi[1] = val[0];
		}
		else
		{
			i = index - dop * 9 - (dofs - dop * 3) * 2;
			lagrange_basis1(lambda, s, eta, xi, dop * 3 + i, dop, val);
			phi[0] = val[0];
			phi[1] = -val[1];
		}
	}
}

/**
* \fn void huzhang_basisCurlTrace(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, int index, int dop, double phi[2])
* \brief curl of the trace of basis function of Hu-Zhang element
* \param *lambda pointer to the area coordiante
* \param s the area of the triangule
* \param eta[3] some auxiliary parameter
* \param xi[3] some auxiliary parameter
* \param **nv the unit normal vectors of the three edges
* \param **tv the unit tangential vectors of the three edges
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi divergence of basis function
* \return void
*/
void huzhang_basisCurlTrace(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, int index, int dop, double phi[2])
{
	int dofs = (dop + 1)*(dop + 2) / 2; // degrees of freedom

	phi[0] = 0;
	phi[1] = 0;
	if (index >= dofs * 3 || index<0)
		return;

	double val[2];
	double nn[3][3], nt[3][3], tt[3][3];
	int i;
	for (i = 0; i<3; i++)
	{
		nn[i][0] = nv[i][0] * nv[i][0]; nn[i][1] = nv[i][1] * nv[i][1]; nn[i][2] = nv[i][0] * nv[i][1];
		nt[i][0] = nv[i][0] * tv[i][0] * sqrt(2); nt[i][1] = nv[i][1] * tv[i][1] * sqrt(2); nt[i][2] = (nv[i][0] * tv[i][1] + nv[i][1] * tv[i][0]) / sqrt(2);
		tt[i][0] = tv[i][0] * tv[i][0]; tt[i][1] = tv[i][1] * tv[i][1]; tt[i][2] = tv[i][0] * tv[i][1];
	}

	if (dop == 0)
	{
		return;
	}

	else // dop>=1
	{
		if (index<3)
		{
			lagrange_basis1(lambda, s, eta, xi, index, dop, val);
			phi[0] = -val[1];
			phi[1] = val[0];
		}
		else if (index<6)
		{
			lagrange_basis1(lambda, s, eta, xi, index % 3, dop, val);
			phi[0] = -val[1];
			phi[1] = val[0];
		}
		else if (index<9)
		{
			lagrange_basis1(lambda, s, eta, xi, index % 3, dop, val);
			phi[0] = 0;
			phi[1] = 0;
		}
		else if (index<9 + (dop - 1) * 3)
		{
			i = index - 9;
			lagrange_basis1(lambda, s, eta, xi, 3 + i, dop, val);
			phi[0] = -val[1] * (nn[i / (dop - 1)][0] + nn[i / (dop - 1)][1]);
			phi[1] = val[0] * (nn[i / (dop - 1)][0] + nn[i / (dop - 1)][1]);
		}
		else if (index<9 + (dop - 1) * 6)
		{
			i = index - 9 - (dop - 1) * 3;
			lagrange_basis1(lambda, s, eta, xi, 3 + i, dop, val);
			phi[0] = -val[1] * (nt[i / (dop - 1)][0] + nt[i / (dop - 1)][1]);
			phi[1] = val[0] * (nt[i / (dop - 1)][0] + nt[i / (dop - 1)][1]);
		}
		else if (index<9 + (dop - 1) * 9)
		{
			i = index - 9 - (dop - 1) * 6;
			lagrange_basis1(lambda, s, eta, xi, 3 + i, dop, val);
			phi[0] = -val[1] * (tt[i / (dop - 1)][0] + tt[i / (dop - 1)][1]);
			phi[1] = val[0] * (tt[i / (dop - 1)][0] + tt[i / (dop - 1)][1]);
		}
		else if (index<dop * 9 + (dofs - dop * 3))
		{
			i = index - dop * 9;
			lagrange_basis1(lambda, s, eta, xi, dop * 3 + i, dop, val);
			phi[0] = -val[1];
			phi[1] = val[0];
		}
		else if (index<dop * 9 + (dofs - dop * 3) * 2)
		{
			i = index - dop * 9 - (dofs - dop * 3);
			lagrange_basis1(lambda, s, eta, xi, dop * 3 + i, dop, val);
			phi[0] = -val[1];
			phi[1] = val[0];
		}
		else
		{
			i = index - dop * 9 - (dofs - dop * 3) * 2;
			lagrange_basis1(lambda, s, eta, xi, dop * 3 + i, dop, val);
			phi[0] = 0;
			phi[1] = 0;
		}
	}
}

/**
* \fn void huzhang_basisROTROT(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, int index, int dop, double *phi)
* \brief rotrot of basis function of Hu-Zhang element
* \param *lambda pointer to the area coordiante
* \param s the area of the triangule
* \param eta[3] some auxiliary parameter
* \param xi[3] some auxiliary parameter
* \param **nv the unit normal vectors of the three edges
* \param **tv the unit tangential vectors of the three edges
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi divergence of basis function
* \return void
*/
void huzhang_basisROTROT(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, int index, int dop, double *phi)
{
	int dofs = (dop + 1)*(dop + 2) / 2; // degrees of freedom

	*phi = 0;
	if (index >= dofs * 3 || index<0)
		return;

	double val[3];
	double nn[3][3], nt[3][3], tt[3][3];
	int i;
	for (i = 0; i<3; i++)
	{
		nn[i][0] = nv[i][0] * nv[i][0]; nn[i][1] = nv[i][1] * nv[i][1]; nn[i][2] = nv[i][0] * nv[i][1];
		nt[i][0] = nv[i][0] * tv[i][0] * sqrt(2); nt[i][1] = nv[i][1] * tv[i][1] * sqrt(2); nt[i][2] = (nv[i][0] * tv[i][1] + nv[i][1] * tv[i][0]) / sqrt(2);
		tt[i][0] = tv[i][0] * tv[i][0]; tt[i][1] = tv[i][1] * tv[i][1]; tt[i][2] = tv[i][0] * tv[i][1];
	}

	if (dop < 2)
	{
		return;
	}

	else // dop>=2
	{
		if (index<3)
		{
			lagrange_basis2(lambda, s, eta, xi, index, dop, val);
			*phi = val[1];
		}
		else if (index<6)
		{
			lagrange_basis2(lambda, s, eta, xi, index % 3, dop, val);
			*phi = val[0];
		}
		else if (index<9)
		{
			lagrange_basis2(lambda, s, eta, xi, index % 3, dop, val);
			*phi = -2 * val[2];
		}
		else if (index<9 + (dop - 1) * 3)
		{
			i = index - 9;
			lagrange_basis2(lambda, s, eta, xi, 3 + i, dop, val);
			*phi = val[1] * nn[i / (dop - 1)][0] + val[0] * nn[i / (dop - 1)][1] - 2 * val[2] * nn[i / (dop - 1)][2];
		}
		else if (index<9 + (dop - 1) * 6)
		{
			i = index - 9 - (dop - 1) * 3;
			lagrange_basis2(lambda, s, eta, xi, 3 + i, dop, val);
			*phi = val[1] * nt[i / (dop - 1)][0] + val[0] * nt[i / (dop - 1)][1] - 2 * val[2] * nt[i / (dop - 1)][2];
		}
		else if (index<9 + (dop - 1) * 9)
		{
			i = index - 9 - (dop - 1) * 6;
			lagrange_basis2(lambda, s, eta, xi, 3 + i, dop, val);
			*phi = val[1] * tt[i / (dop - 1)][0] + val[0] * tt[i / (dop - 1)][1] - 2 * val[2] * tt[i / (dop - 1)][2];
		}
		else if (index<dop * 9 + (dofs - dop * 3))
		{
			i = index - dop * 9;
			lagrange_basis2(lambda, s, eta, xi, dop * 3 + i, dop, val);
			*phi = val[1];
		}
		else if (index<dop * 9 + (dofs - dop * 3) * 2)
		{
			i = index - dop * 9 - (dofs - dop * 3);
			lagrange_basis2(lambda, s, eta, xi, dop * 3 + i, dop, val);
			*phi = val[0];
		}
		else
		{
			i = index - dop * 9 - (dofs - dop * 3) * 2;
			lagrange_basis2(lambda, s, eta, xi, dop * 3 + i, dop, val);
			*phi = -2 * val[2];
		}
	}
}

/**
* \fn void huzhang_basisLaplaceTrace(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, int index, int dop, double *phi)
* \brief Laplace of the trace of basis function of Hu-Zhang element
* \param *lambda pointer to the area coordiante
* \param s the area of the triangule
* \param eta[3] some auxiliary parameter
* \param xi[3] some auxiliary parameter
* \param **nv the unit normal vectors of the three edges
* \param **tv the unit tangential vectors of the three edges
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi divergence of basis function
* \return void
*/
void huzhang_basisLaplaceTrace(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, int index, int dop, double *phi)
{
	int dofs = (dop + 1)*(dop + 2) / 2; // degrees of freedom

	*phi = 0;
	if (index >= dofs * 3 || index<0)
		return;

	double val[3];
	double nn[3][3], nt[3][3], tt[3][3];
	int i;
	for (i = 0; i<3; i++)
	{
		nn[i][0] = nv[i][0] * nv[i][0]; nn[i][1] = nv[i][1] * nv[i][1]; nn[i][2] = nv[i][0] * nv[i][1];
		nt[i][0] = nv[i][0] * tv[i][0] * sqrt(2); nt[i][1] = nv[i][1] * tv[i][1] * sqrt(2); nt[i][2] = (nv[i][0] * tv[i][1] + nv[i][1] * tv[i][0]) / sqrt(2);
		tt[i][0] = tv[i][0] * tv[i][0]; tt[i][1] = tv[i][1] * tv[i][1]; tt[i][2] = tv[i][0] * tv[i][1];
	}

	if (dop < 2)
	{
		return;
	}

	else // dop>=2
	{
		if (index<3)
		{
			lagrange_basis2(lambda, s, eta, xi, index, dop, val);
			*phi = val[0] + val[1];
		}
		else if (index<6)
		{
			lagrange_basis2(lambda, s, eta, xi, index % 3, dop, val);
			*phi = val[0] + val[1];
		}
		else if (index<9)
		{
			lagrange_basis2(lambda, s, eta, xi, index % 3, dop, val);
			*phi = 0;
		}
		else if (index<9 + (dop - 1) * 3)
		{
			i = index - 9;
			lagrange_basis2(lambda, s, eta, xi, 3 + i, dop, val);
			*phi = (val[0] + val[1]) * (nn[i / (dop - 1)][0] + nn[i / (dop - 1)][1]);
		}
		else if (index<9 + (dop - 1) * 6)
		{
			i = index - 9 - (dop - 1) * 3;
			lagrange_basis2(lambda, s, eta, xi, 3 + i, dop, val);
			*phi = (val[0] + val[1]) * (nt[i / (dop - 1)][0] + nt[i / (dop - 1)][1]);
		}
		else if (index<9 + (dop - 1) * 9)
		{
			i = index - 9 - (dop - 1) * 6;
			lagrange_basis2(lambda, s, eta, xi, 3 + i, dop, val);
			*phi = (val[0] + val[1]) * (tt[i / (dop - 1)][0] + tt[i / (dop - 1)][1]);
		}
		else if (index<dop * 9 + (dofs - dop * 3))
		{
			i = index - dop * 9;
			lagrange_basis2(lambda, s, eta, xi, dop * 3 + i, dop, val);
			*phi = val[0] + val[1];
		}
		else if (index<dop * 9 + (dofs - dop * 3) * 2)
		{
			i = index - dop * 9 - (dofs - dop * 3);
			lagrange_basis2(lambda, s, eta, xi, dop * 3 + i, dop, val);
			*phi = val[0] + val[1];
		}
		else
		{
			i = index - dop * 9 - (dofs - dop * 3) * 2;
			lagrange_basis2(lambda, s, eta, xi, dop * 3 + i, dop, val);
			*phi = 0;
		}
	}
}

/**
 * \fn double area(double x1,double x2,double x3,double y1,double y2,double y3)
 * \brief get area for triangle p1(x1,y1),p2(x2,y2),p3(x3,y3)
 * \param x1 the x-axis value of the point p1
 * \param x2 the x-axis value of the point p2
 * \param x3 the x-axis value of the point p3
 * \param y1 the y-axis value of the point p1
 * \param y2 the y-axis value of the point p2
 * \param y3 the y-axis value of the point p3
 * \return area of the trianle
 */
double area(double x1,double x2,double x3,double y1,double y2,double y3)
{
	return ((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))/2;
}

/** 
 * \fn void localb(double (*nodes)[2],double *b)
 * \brief get local right-hand side b from triangle nodes
 * \param (*nodes)[2] the vertice of the triangule
 * \param *b local right-hand side
 * \return void
 */
void localb(double (*nodes)[2],double *b)
{
	int i;
	double x,y,a;
	double s=area(nodes[0][0],nodes[1][0],nodes[2][0],nodes[0][1],nodes[1][1],nodes[2][1]);
	
	int num_qp=49; // the number of numerical intergation points
	double gauss[num_qp][3];
	init_Gauss(num_qp, 2, gauss); // gauss intergation initial	
	
	for(i=0;i<3;i++)
		b[i]=0;
	
	for(i=0;i<num_qp;i++)
	{
		x=nodes[0][0]*gauss[i][0]+nodes[1][0]*gauss[i][1]+nodes[2][0]*(1-gauss[i][0]-gauss[i][1]);
		y=nodes[0][1]*gauss[i][0]+nodes[1][1]*gauss[i][1]+nodes[2][1]*(1-gauss[i][0]-gauss[i][1]);
		a=f(x,y);
		
		b[0]+=2*s*gauss[i][2]*a*gauss[i][0];
		b[1]+=2*s*gauss[i][2]*a*gauss[i][1];
		b[2]+=2*s*gauss[i][2]*a*(1-gauss[i][0]-gauss[i][1]);
		
	}
}

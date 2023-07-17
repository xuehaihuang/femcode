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
#include <stdio.h>
#include "header.h"
#include "matvec.h"


/**
* \fn void cartToPol2d(x,y,r,theta)
* \brief Transform Cartesian coordinate to polar  coordinate
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param r the distance from the origin to a point in the x-y plane
* \param theta a counterclockwise angular displacement in radians from the positive x-axis
* \return void
*/
void cartToPol2d(double x, double y, double *r, double *theta)
{
	(*r) = sqrt(x*x + y*y);
	(*theta) = atan2(y, x);
	if ((*theta)<0)
		(*theta) += 2 * PI;
}

/**
* \fn void cartToBary2d(double *x, double *lambda, double **tri)
* \brief Transform Cartesian coordinate to barycentric coordinate
* \param x the Cartesian coordinate
* \param lambda the barycentric coordinate
* \param tri the Cartesian coordinate of three vertices of triangle
* \return void
*/
void cartToBary2d(double *x, double *lambda, double **tri)
{
	double s = area(tri); 

	lambda[0] = area0(x, tri[1], tri[2]) / s;
	lambda[1] = area0(tri[0], x, tri[2]) / s;
	lambda[2] = area0(tri[0], tri[1], x) / s;	
}

/**
* \fn void baryToCart2d(double *lambda, double *x, double **tri)
* \brief Transform barycentric to Cartesian coordinate coordinate
* \param lambda the barycentric coordinate
* \param x the Cartesian coordinate
* \param tri the Cartesian coordinate of three vertices of triangle
* \return void
*/
void baryToCart2d(double *lambda, double *x, double **tri)
{
	x[0] = lambda[0]*tri[0][0] + lambda[1]*tri[1][0] + lambda[2]*tri[2][0];
	x[1] = lambda[0]*tri[0][1] + lambda[1]*tri[1][1] + lambda[2]*tri[2][1];
}


/** 
 * \fn void morley_basis(double *lambda, double **gradLambda, double *nve[3], int index, double *phi)
 * \brief basis function of Morley element
 * \param lambda the barycentric coordinate
 * \param gradLambda the gradient of  barycentric coordinate
 * \param nve[3] the unit normal vectors of three edge
 * \param index the indicator of the basis function
 * \param *phi basis function
 * \return void
 */
void morley_basis(double *lambda, double **gradLambda, double *nve[3], int index, double *phi)
{
	int i,j,k;
	double c1, c2;
	if(index<0 || index>5) *phi=0.0;

	if(index<3){
		i = index;
		j = (i+1)%3;
		k = (i+2)%3;
		c1 = dot_array(2, gradLambda[i], gradLambda[j])/ lpnormp_array(2, gradLambda[j], 2);
		c2 = dot_array(2, gradLambda[i], gradLambda[k])/ lpnormp_array(2, gradLambda[k], 2);
		*phi = lambda[i]*lambda[i] + c1*lambda[j]*(lambda[j]-1) + c2*lambda[k]*(lambda[k]-1);
	}
	else{
		i = index-3;
		c1 = -1.0 / dot_array(2, nve[i], gradLambda[i]);
		*phi = c1*lambda[i]*(lambda[i]-1);
	}
}

/** 
 * \fn void morley_basis1(double *lambda, double **gradLambda, double *nve[3], int index, double phi[2])
 * \brief the gradient of basis function of Morley element
 * \param lambda the barycentric coordinate
 * \param gradLambda the gradient of  barycentric coordinate
 * \param nve[3] the unit normal vectors of three edge
 * \param index the indicator of the basis function
 * \param *phi basis function
 * \return void
 */
void morley_basis1(double *lambda, double **gradLambda, double *nve[3], int index, double phi[2])
{
	int i,j,k;
	double c1, c2;
	if(index<0 || index>5){ 
		phi[0]=0.0; phi[1]=0.0;
	}

	if(index<3){
		i = index;
		j = (i+1)%3;
		k = (i+2)%3;
		c1 = dot_array(2, gradLambda[i], gradLambda[j])/ lpnormp_array(2, gradLambda[j], 2);
		c2 = dot_array(2, gradLambda[i], gradLambda[k])/ lpnormp_array(2, gradLambda[k], 2);
		axpbyz_array(2, 2*lambda[i], gradLambda[i], c1*(2*lambda[j]-1), gradLambda[j], phi);
		axpy_array(2, c2*(2*lambda[k]-1), gradLambda[k], phi);
	}
	else{
		i = index-3;
		c1 = -1.0 / dot_array(2, nve[i], gradLambda[i]);
		axy_array(2, c1*(2*lambda[i]-1), gradLambda[i], phi);
	}
}

/** 
 * \fn void morley_basis2(double **gradLambda, double *nve[3], int index, double phi[3])
 * \brief the second order derivative of Morley element basis function: (\partial_{xx}phi, \partial_{yy}phi, \partial_{xy}phi)
 * \param lambda the barycentric coordinate
 * \param gradLambda the gradient of  barycentric coordinate
 * \param nve[3] the unit normal vectors of three edge
 * \param index the indicator of the basis function
 * \param *phi basis function
 * \return void
 */
void morley_basis2(double **gradLambda, double *nve[3], int index, double phi[3])
{
	int i,j,k;
	double c1, c2;
	if(index<0 || index>5){ 
		phi[0]=0.0; phi[1]=0.0; phi[2]=0.0;
	}

	if(index<3){
		i = index;
		j = (i+1)%3;
		k = (i+2)%3;
		c1 = dot_array(2, gradLambda[i], gradLambda[j])/ lpnormp_array(2, gradLambda[j], 2);
		c2 = dot_array(2, gradLambda[i], gradLambda[k])/ lpnormp_array(2, gradLambda[k], 2);
		phi[0]= 2*gradLambda[i][0]*gradLambda[i][0]+2*c1*gradLambda[j][0]*gradLambda[j][0]+2*c2*gradLambda[k][0]*gradLambda[k][0];
		phi[1]= 2*gradLambda[i][1]*gradLambda[i][1]+2*c1*gradLambda[j][1]*gradLambda[j][1]+2*c2*gradLambda[k][1]*gradLambda[k][1];
		phi[2]= 2*gradLambda[i][0]*gradLambda[i][1]+2*c1*gradLambda[j][0]*gradLambda[j][1]+2*c2*gradLambda[k][0]*gradLambda[k][1];
	}
	else{
		i = index-3;
		c1 = -1.0 / dot_array(2, nve[i], gradLambda[i]);
		phi[0]= 2*c1*gradLambda[i][0]*gradLambda[i][0];
		phi[1]= 2*c1*gradLambda[i][1]*gradLambda[i][1];
		phi[2]= 2*c1*gradLambda[i][0]*gradLambda[i][1];
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
	int in, ie, ii, i1, i2, l, i;
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
			*phi = lambda[index]*(2.0*lambda[index]-1.0);
		}
		else
		{
			in = index-3;
			i1=(in+1)%3;
			i2=(in+2)%3;
			*phi = 4.0*lambda[i1]*lambda[i2];
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

	/********* Bernstein basis ***********/
	else if(dop>=6)  // it works for any dop >= 2
	{
		if(index<3)
		{
			*phi =  pow(lambda[index], dop);
		}
		else if(index < 3*dop) 
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1); // ii=0, 1 or 2 if dop=4
			i1 = (ie+1)%3;
			i2 = (ie+2)%3;
			*phi =  4.0 *pow(lambda[i1], ii+1)*pow(lambda[i2], dop-1-ii);
		}
		else
		{
			ii = index-3*dop;
			for(l=0;l<dop-2;l++)
			{
				if(ii<(l+1)*(l+2)/2) 
					break;
			}
			i=ii-(l+1)*l/2;
			*phi =  27.0*pow(lambda[0], dop-2-l)*pow(lambda[1], l-i+1)*pow(lambda[2], i+1);
		}
	}
}

/** 
 * \fn void lagrange_basis1(double *lambda, double **gradLambda, int index, int dop, double phi[2])
 * \brief the first order derivative of Lagrange element basis function: (\partial_{x}phi, \partial_{y}phi)
 * \param lambda pointer to the area coordiante
 * \param gradLambda pointer to the gradient of the barycentric coordinate
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param phi[2] the first order derivative of Morley element basis function: (\partial_{x}phi, \partial_{y}phi)
 * \return void
 */
void lagrange_basis1(double *lambda, double **gradLambda, int index, int dop, double phi[2])
{
	double s, eta[3], xi[3], c0, c1, c2;	
	int in, ie, ii, i1, i2, i3, l, i;
	s = 0.5 / (gradLambda[0][0]*gradLambda[1][1] - gradLambda[0][1]*gradLambda[1][0]);
	for(int i=0;i<3;i++){
		eta[i] = gradLambda[i][0]*2*s;
		xi[i] = -gradLambda[i][1]*2*s;
	}

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
		copy_array(2, gradLambda[index], phi);

	else if(dop==2)
	{
		if(index<3)
			axy_array(2, 4.0*lambda[index]-1.0, gradLambda[index], phi);
		else
		{
			in = index-3;
			i1=(in+1)%3;
			i2=(in+2)%3;
			axpbyz_array(2, 4.0*lambda[i2], gradLambda[i1], 4*lambda[i1], gradLambda[i2], phi);
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

	/********* Bernstein basis ***********/
	else if(dop>=6)  // it works for any dop >= 2
	{
		if(index<3)
		{
			// *phi =  pow(lambda[index], dop);
			axy_array(2, dop*pow(lambda[index], dop-1), gradLambda[index], phi);
		}
		else if(index < 3*dop) 
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1); // ii=0, 1 or 2 if dop=4
			i1 = (ie+1)%3;
			i2 = (ie+2)%3;
			// *phi =  4.0 *pow(lambda[i1], ii+1)*pow(lambda[i2], dop-1-ii);
			c1 =  4.0 *(ii+1)*pow(lambda[i1], ii)*pow(lambda[i2], dop-1-ii);
			c2 =  4.0 *(dop-1-ii)*pow(lambda[i1], ii+1)*pow(lambda[i2], dop-2-ii);
			axpbyz_array(2, c1, gradLambda[i1], c2, gradLambda[i2], phi);
		}
		else
		{
			ii = index-3*dop;
			for(l=0;l<dop-2;l++)
			{
				if(ii<(l+1)*(l+2)/2) 
					break;
			}
			i=ii-(l+1)*l/2;
			// *phi =  27.0*pow(lambda[0], dop-2-l)*pow(lambda[1], l-i+1)*pow(lambda[2], i+1);
			c0 =  27.0*(dop-2-l)*pow(lambda[0], dop-3-l)*pow(lambda[1], l-i+1)*pow(lambda[2], i+1);
			c1 =  27.0*(l-i+1)*pow(lambda[0], dop-2-l)*pow(lambda[1], l-i)*pow(lambda[2], i+1);
			axpbyz_array(2, c0, gradLambda[0], c1, gradLambda[1], phi);
	        c2 =  27.0*(i+1)*pow(lambda[0], dop-2-l)*pow(lambda[1], l-i+1)*pow(lambda[2], i);
			axpy_array(2, c2, gradLambda[2], phi);
		}
	}
}

/** 
 * \fn void lagrange_basis2(double *lambda, double **gradLambda, int index, int dop, double phi[3])
 * \brief the second order derivative of Lagrange element basis function: (\partial_{xx}phi, \partial_{yy}phi, \partial_{xy}phi)
 * \param *lambda pointer to the area coordiante
 * \param **gradLambda pointer to the gradient of the barycentric coordinate
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param phi[2] the first order derivative of Morley element basis function: (\partial_{xx}phi, \partial_{yy}phi, \partial_{xy}phi)
 * \return void
 */
void lagrange_basis2(double *lambda, double **gradLambda, int index, int dop, double phi[3])
{
	double s, eta[3], xi[3], c0, c1, c2;	
	int in, ie, ii, i1, i2, i3, l, i;
	s = 0.5 / (gradLambda[0][0]*gradLambda[1][1] - gradLambda[0][1]*gradLambda[1][0]);
	for(int i=0;i<3;i++){
		eta[i] = gradLambda[i][0]*2*s;
		xi[i] = -gradLambda[i][1]*2*s;
	}

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
			phi[0]=4*gradLambda[index][0]*gradLambda[index][0];
			phi[1]=4*gradLambda[index][1]*gradLambda[index][1];
			phi[2]=4*gradLambda[index][0]*gradLambda[index][1];
		}
		else
		{
			i1=(index+1)%3;
			i2=(index+2)%3;
			phi[0]=8*gradLambda[i1][0]*gradLambda[i2][0];
			phi[1]=8*gradLambda[i1][1]*gradLambda[i2][1];
			phi[2]=4*(gradLambda[i1][0]*gradLambda[i2][1]+gradLambda[i1][1]*gradLambda[i2][0]);
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

	/********* Bernstein basis ***********/
	else if(dop>=6)  // it works for any dop >= 2
	{
		if(index<3)
		{
			// *phi =  pow(lambda[index], dop);
			c0 = dop*(dop-1)*pow(lambda[index], dop-2);
			phi[0] = c0*gradLambda[index][0]*gradLambda[index][0];
			phi[1] = c0*gradLambda[index][1]*gradLambda[index][1];
			phi[2] = c0*gradLambda[index][0]*gradLambda[index][1];
		}
		else if(index < 3*dop) 
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1); // ii=0, 1 or 2 if dop=4
			i1 = (ie+1)%3;
			i2 = (ie+2)%3;
			// *phi =  4.0 *pow(lambda[i1], ii+1)*pow(lambda[i2], dop-1-ii);
			c0 = 4.0*(ii+1)*ii*pow(lambda[i1], ii-1)*pow(lambda[i2], dop-1-ii);
			phi[0] = c0*gradLambda[i1][0]*gradLambda[i1][0];
			phi[1] = c0*gradLambda[i1][1]*gradLambda[i1][1];
			phi[2] = c0*gradLambda[i1][0]*gradLambda[i1][1];
			c0 = 4.0*(dop-1-ii)*(dop-2-ii)*pow(lambda[i1], ii+1)*pow(lambda[i2], dop-3-ii);
			phi[0] += c0*gradLambda[i2][0]*gradLambda[i2][0];
			phi[1] += c0*gradLambda[i2][1]*gradLambda[i2][1];
			phi[2] += c0*gradLambda[i2][0]*gradLambda[i2][1];
			c0 = 4.0*(ii+1)*(dop-1-ii)*pow(lambda[i1], ii)*pow(lambda[i2], dop-2-ii);
			phi[0] += c0*2*gradLambda[i1][0]*gradLambda[i2][0];
			phi[1] += c0*2*gradLambda[i1][1]*gradLambda[i2][1];
			phi[2] += c0*(gradLambda[i1][0]*gradLambda[i2][1]+gradLambda[i1][1]*gradLambda[i2][0]);
		}
		else
		{
			ii = index-3*dop;
			for(l=0;l<dop-2;l++)
			{
				if(ii<(l+1)*(l+2)/2) 
					break;
			}
			i=ii-(l+1)*l/2;
			// *phi =  27.0*pow(lambda[0], dop-2-l)*pow(lambda[1], l-i+1)*pow(lambda[2], i+1);
			c0 = 27.0*(dop-2-l)*(dop-3-l)*pow(lambda[0], dop-4-l)*pow(lambda[1], l-i+1)*pow(lambda[2], i+1);
			phi[0] = c0*gradLambda[0][0]*gradLambda[0][0];
			phi[1] = c0*gradLambda[0][1]*gradLambda[0][1];
			phi[2] = c0*gradLambda[0][0]*gradLambda[0][1];
			c0 = 27.0*(l-i+1)*(l-i)*pow(lambda[0], dop-2-l)*pow(lambda[1], l-i-1)*pow(lambda[2], i+1);
			phi[0] += c0*gradLambda[1][0]*gradLambda[1][0];
			phi[1] += c0*gradLambda[1][1]*gradLambda[1][1];
			phi[2] += c0*gradLambda[1][0]*gradLambda[1][1];
			c0 = 27.0*(i+1)*i*pow(lambda[0], dop-2-l)*pow(lambda[1], l-i+1)*pow(lambda[2], i-1);
			phi[0] += c0*gradLambda[2][0]*gradLambda[2][0];
			phi[1] += c0*gradLambda[2][1]*gradLambda[2][1];
			phi[2] += c0*gradLambda[2][0]*gradLambda[2][1];
			c0 = 27.0*(l-i+1)*(i+1)*pow(lambda[0], dop-2-l)*pow(lambda[1], l-i)*pow(lambda[2], i);
			phi[0] += c0*2*gradLambda[1][0]*gradLambda[2][0];
			phi[1] += c0*2*gradLambda[1][1]*gradLambda[2][1];
			phi[2] += c0*(gradLambda[1][0]*gradLambda[2][1]+gradLambda[1][1]*gradLambda[2][0]);
			c0 = 27.0*(dop-2-l)*(i+1)*pow(lambda[0], dop-3-l)*pow(lambda[1], l-i+1)*pow(lambda[2], i);
			phi[0] += c0*2*gradLambda[2][0]*gradLambda[0][0];
			phi[1] += c0*2*gradLambda[2][1]*gradLambda[0][1];
			phi[2] += c0*(gradLambda[2][0]*gradLambda[0][1]+gradLambda[2][1]*gradLambda[0][0]);
			c0 = 27.0*(dop-2-l)*(l-i+1)*pow(lambda[0], dop-3-l)*pow(lambda[1], l-i)*pow(lambda[2], i+1);
			phi[0] += c0*2*gradLambda[0][0]*gradLambda[1][0];
			phi[1] += c0*2*gradLambda[0][1]*gradLambda[1][1];
			phi[2] += c0*(gradLambda[0][0]*gradLambda[1][1]+gradLambda[0][1]*gradLambda[1][0]);
		}
	}
}

/** 
 * \fn void bernstein2d_basis(double *lambda, int index, int dop, double *phi)
 * \brief Bernstein basis function of Lagrange element
 * \param *lambda pointer to the area coordiante
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param *phi basis function
 * \return void
 */
void bernstein2d_basis(double *lambda, int index, int dop, double *phi)
{
	int in, ie, ii, i1, i2, l, i;
	int dofs = (dop+1)*(dop+2)/2; // degrees of freedom
	if(index>= dofs || index<0){
		*phi=0;
		return;
	}

	if(dop==0){
		*phi = 1;
	}

	else if(dop==1){
		*phi = lambda[index];
	} // dop=1
	
	else if(dop>=2) 
	{
		if(index<3){
			*phi =  pow(lambda[index], dop);
		}
		else if(index < 3*dop){
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1); // ii=0, 1 or 2 if dop=4
			i1 = (ie+1)%3;
			i2 = (ie+2)%3;
			*phi =  4.0 *pow(lambda[i1], ii+1)*pow(lambda[i2], dop-1-ii);
		}
		else{
			ii = index-3*dop;
			for(l=0;l<dop-2;l++)
			{
				if(ii<(l+1)*(l+2)/2) 
					break;
			}
			i=ii-(l+1)*l/2;
			*phi =  27.0*pow(lambda[0], dop-2-l)*pow(lambda[1], l-i+1)*pow(lambda[2], i+1);
		}
	}
}

/** 
 * \fn void bernstein2d_basis1(double *lambda, double **gradLambda, int index, int dop, double phi[2])
 * \brief the first order derivative of Bernstein element basis function: (\partial_{x}phi, \partial_{y}phi)
 * \param *lambda pointer to the area coordiante
 * \param **gradLambda pointer to the gradient of the barycentric coordinate
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param phi[2] the first order derivative of Morley element basis function: (\partial_{x}phi, \partial_{y}phi)
 * \return void
 */
void bernstein2d_basis1(double *lambda, double **gradLambda, int index, int dop, double phi[2])
{
	double c0, c1, c2;	
	int in, ie, ii, i1, i2, i3, l, i;
	// s = 0.5 / (gradLambda[0][0]*gradLambda[1][1] - gradLambda[0][1]*gradLambda[1][0]);
	// for(int i=0;i<3;i++){
	// 	eta[i] = gradLambda[i][0]*2*s;
	// 	xi[i] = -gradLambda[i][1]*2*s;
	// }

	int dofs = (dop+1)*(dop+2)/2; // degrees of freedom
	if(index>= dofs || index<0)
	{
		phi[0]=0;
		phi[1]=0;
		return;
	}

	if(dop==0){
		phi[0]=0;
		phi[1]=0;
	} // dop=0

	else if(dop==1)
		copy_array(2, gradLambda[index], phi);

	else if(dop>=2){
		if(index<3){
			// *phi =  pow(lambda[index], dop);
			axy_array(2, dop*pow(lambda[index], dop-1), gradLambda[index], phi);
		}
		else if(index < 3*dop){
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1); // ii=0, 1 or 2 if dop=4
			i1 = (ie+1)%3;
			i2 = (ie+2)%3;
			// *phi =  4.0 *pow(lambda[i1], ii+1)*pow(lambda[i2], dop-1-ii);
			c1 =  4.0 *(ii+1)*pow(lambda[i1], ii)*pow(lambda[i2], dop-1-ii);
			c2 =  4.0 *(dop-1-ii)*pow(lambda[i1], ii+1)*pow(lambda[i2], dop-2-ii);
			axpbyz_array(2, c1, gradLambda[i1], c2, gradLambda[i2], phi);
		}
		else{
			ii = index-3*dop;
			for(l=0;l<dop-2;l++){
				if(ii<(l+1)*(l+2)/2) 
					break;
			}
			i=ii-(l+1)*l/2;
			// *phi =  27.0*pow(lambda[0], dop-2-l)*pow(lambda[1], l-i+1)*pow(lambda[2], i+1);
			c0 =  27.0*(dop-2-l)*pow(lambda[0], dop-3-l)*pow(lambda[1], l-i+1)*pow(lambda[2], i+1);
			c1 =  27.0*(l-i+1)*pow(lambda[0], dop-2-l)*pow(lambda[1], l-i)*pow(lambda[2], i+1);
			axpbyz_array(2, c0, gradLambda[0], c1, gradLambda[1], phi);
	        c2 =  27.0*(i+1)*pow(lambda[0], dop-2-l)*pow(lambda[1], l-i+1)*pow(lambda[2], i);
			axpy_array(2, c2, gradLambda[2], phi);
		}
	}
}


/** 
 * \fn void cr_basis(int n, double *lambda, int index, double *phi)
 * \brief basis function of Crouzeix–Raviart element
 * \param n the dimension of space
 * \param *lambda pointer to the area coordiante
 * \param index the indicator of the basis function
 * \param *phi basis function
 * \return void
 */
void cr_basis(int n, double *lambda, int index, double *phi)
{
	*phi=1-n*lambda[index];
}

/** 
 * \fn void cr_basis1(int n, double **gradLambda, int index, double *phi)
 * \brief the gradient of Crouzeix–Raviart element basis function
 * \param n the dimension of space
 * \param **gradLambda pointer to the gradient of the barycentric coordinate
 * \param index the indicator of the basis function
 * \param phi the gradient of Crouzeix–Raviart element basis function
 * \return void
 */
void cr_basis1(int n, double **gradLambda, int index, double *phi)
{
	axy_array(n, -1.0*n, gradLambda[index], phi);
}

/** 
 * \fn void rt_basis(double *lambda, double *height, double **tij, short *eorien, int index, int dop, double phi[2])
 * Arnold, D. N.; Falk, R. S. & Winther, R. Geometric decompositions and local bases for spaces of finite element differential forms Comput. Methods Appl. Mech. Engrg., 2009, 198, 1660-1672
 * \param lambda the barycentric coordinate
 * \param height pointer to height to three edges
 * \param tij pointer to three edge vectors
 * \param eorien pointer to orient of three edges
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param phi[2] basis function of Raviart-Thomas element: (phi1, phi2)
 * \return void
 */
void rt_basis(double *lambda, double *height, double **tij, short *eorien, int index, int dop, double phi[2])
{
	int dofs = dop *(dop + 2); // degrees of freedom
	
	init_array(2, phi, 0);
		
	if (dop<1) return;

	if (index >= dofs || index<0)	return;

	int i, i1, i2, j, k, jk[0];
	int in, ie, ii;
	double val;

	if (dop == 1){
		i = index;
		j = (i+1)%3;
		k = (i+2)%3;
		axpbyz_array(2, lambda[j], tij[k], -lambda[k], tij[j], phi);
		ax_array(2, eorien[i]/height[i], phi);
	}

	else if (dop == 2){
		if (index < 3*dop){
			in = index;
			i = in / dop;
			ii = in % dop;
			j = (i+1)%3;
			k = (i+2)%3;
			jk[0] = j;	jk[1] = k;
			axpbyz_array(2, lambda[j], tij[k], -lambda[k], tij[j], phi);
			ax_array(2, lambda[jk[ii]], phi);
			ax_array(2, eorien[i]/height[i], phi);
		}
		else{
			i = index - 3*dop;
			j = (i+1)%3;
			k = (i+2)%3;
			axpbyz_array(2, lambda[j], tij[k], -lambda[k], tij[j], phi);
			ax_array(2, lambda[i]/height[i], phi);
		}
	} // dop=2

	else if (dop == 3){
		if (index < 3*dop){
			in = index;
			i = in / dop;
			ii = in % dop;
			j = (i+1)%3;
			k = (i+2)%3;
			axpbyz_array(2, lambda[j], tij[k], -lambda[k], tij[j], phi);
			ax_array(2, pow(lambda[j], dop-1-ii)*pow(lambda[k], ii), phi);
			ax_array(2, eorien[i]/height[i], phi);
		}
		else{
			in = index - 3*dop;
			i = in / 3;
			ii = in % 3;
			j = (i+1)%3;
			k = (i+2)%3;
			axpbyz_array(2, lambda[j], tij[k], -lambda[k], tij[j], phi);
			ax_array(2, lambda[ii]*lambda[i]/height[i], phi);
		}
	} // dop=3

	else if (dop >= 4){
		if (index < 3*dop){
			in = index;
			i = in / dop;
			ii = in % dop;
			j = (i+1)%3;
			k = (i+2)%3;
			
			axpbyz_array(2, lambda[j], tij[k], -lambda[k], tij[j], phi);
			ax_array(2, pow(lambda[j], dop-1-ii)*pow(lambda[k], ii), phi);
			ax_array(2, eorien[i]/height[i], phi);
		}
		else{
			in = index - 3*dop;
			i = in / (dop*(dop-1)/2);
			ii = in % (dop*(dop-1)/2);
			j = (i+1)%3;
			k = (i+2)%3;
			axpbyz_array(2, lambda[j], tij[k], -lambda[k], tij[j], phi);
			bernstein2d_basis(lambda, ii, dop-2, &val);
			ax_array(2, val*lambda[i]/height[i], phi);
		}
	} // dop >= 4
}

/** 
 * \fn void rt_basisDIV(double *lambda, double **grd_lambda, double *height, double **tij, short *eorien, int index, int dop, double *phi)
 * \brief the divergence of Raviart-Thomas element basis function: \partial_{x}phi1 + \partial_{y}phi2
 * Arnold, D. N.; Falk, R. S. & Winther, R. Geometric decompositions and local bases for spaces of finite element differential forms Comput. Methods Appl. Mech. Engrg., 2009, 198, 1660-1672
 * \param lambda the barycentric coordinate
 * \param gradLambda pointer to the gradient of the barycentric coordinate
 * \param height pointer to height to three edges
 * \param tij pointer to three edge vectors
 * \param eorien pointer to orient of three edges
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param phi[4] the first order derivative of Raviart-Thomas element basis function: \partial_{x}phi1 + \partial_{y}phi2
 * \return void
 */
void rt_basisDIV(double *lambda, double **grd_lambda, double *height, double **tij, short *eorien, int index, int dop, double *phi)
{
	int dofs = dop *(dop + 2); // degrees of freedom
	
	*phi = 0;

	if (dop<1) return;

	if (index >= dofs || index<0)	return;

	int i, i1, i2, j, k, m, jk[0];
	int in, ie, ii;
	double c, val[2], phi1[2];

	if (dop == 1)
		*phi = eorien[index]*2.0/height[index];

	else if (dop == 2){
		if (index < 3*dop){
			in = index;
			i = in / dop;
			ii = in % dop;
			j = (i+1)%3;
			k = (i+2)%3;
			jk[0] = j;	jk[1] = k;
			*phi = eorien[i]*3.0*lambda[jk[ii]]/height[i];
		}
		else{
			i = index - 3*dop;
			// j = (i+1)%3;
			// k = (i+2)%3;
			*phi = (3.0*lambda[i]-1)/height[i];
		}
	} // dop=2

	else if (dop == 3){
		if (index < 3*dop){
			in = index;
			i = in / dop;
			ii = in % dop;
			j = (i+1)%3;
			k = (i+2)%3;
			*phi = eorien[i]*4.0*pow(lambda[j], dop-1-ii)*pow(lambda[k], ii)/height[i];
		}
		else{
			in = index - 3*dop;
			i = in / 3;
			ii = in % 3;
			j = (i+1)%3;
			k = (i+2)%3;
			*phi = lambda[ii]*(3.0*lambda[i]-1)/height[i];
			axpbyz_array(2, lambda[j], tij[k], -lambda[k], tij[j], val);
			ax_array(2, lambda[i]/height[i], val);
			*phi += dot_array(2, val, grd_lambda[ii]);
		}
	} // dop=3

	else if (dop >= 4){
		if (index < 3*dop){
			in = index;
			i = in / dop;
			ii = in % dop;
			j = (i+1)%3;
			k = (i+2)%3;
			*phi = eorien[i]*(dop+1)*pow(lambda[j], dop-1-ii)*pow(lambda[k], ii)/height[i];
		}
		else{
			in = index - 3*dop;
			i = in / (dop*(dop-1)/2);
			ii = in % (dop*(dop-1)/2);
			j = (i+1)%3;
			k = (i+2)%3;
			bernstein2d_basis(lambda, ii, dop-2, &c);
			bernstein2d_basis1(lambda, grd_lambda, ii, dop-2, val);
			*phi = c*(3.0*lambda[i]-1)/height[i];
			axpbyz_array(2, lambda[j], tij[k], -lambda[k], tij[j], phi1);
			ax_array(2, lambda[i]/height[i], phi1);
			*phi += dot_array(2, val, phi1);
		}
	} // dop >= 4
}

/** 
 * \fn void rt_basis1(double *lambda, double **grd_lambda, double *height, double **tij, short *eorien, int index, int dop, double phi[4])
 * \brief the first order derivative of Raviart-Thomas element basis function: (\partial_{x}phi1, \partial_{y}phi1, \partial_{x}phi2, \partial_{y}phi2)
 * Arnold, D. N.; Falk, R. S. & Winther, R. Geometric decompositions and local bases for spaces of finite element differential forms Comput. Methods Appl. Mech. Engrg., 2009, 198, 1660-1672
 * \param lambda the barycentric coordinate
 * \param gradLambda pointer to the gradient of the barycentric coordinate
 * \param height pointer to height to three edges
 * \param tij pointer to three edge vectors
 * \param eorien pointer to orient of three edges
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param phi[4] the first order derivative of Raviart-Thomas element basis function: (\partial_{x}phi1, \partial_{y}phi1, \partial_{x}phi2, \partial_{y}phi2)
 * \return void
 */
void rt_basis1(double *lambda, double **grd_lambda, double *height, double **tij, short *eorien, int index, int dop, double phi[4])
{
	double phi1[4], phi2[4];

	int dofs = dop *(dop + 2); // degrees of freedom
	
	init_array(4, phi, 0);

	if (dop<1) return;

	if (index >= dofs || index<0)	return;

	int i, i1, i2, j, k, m, jk[0];
	int in, ie, ii;
	double c, val[2];

	if (dop == 1){
		i = index;
		j = (i+1)%3;
		k = (i+2)%3;
		for(m=0;m<2;m++)
			axpbyz_array(2, tij[k][m], grd_lambda[j], -tij[j][m], grd_lambda[k], phi+2*m);
		ax_array(4, eorien[i]/height[i], phi);
	}

	else if (dop == 2){
		if (index < 3*dop){
			in = index;
			i = in / dop;
			ii = in % dop;
			j = (i+1)%3;
			k = (i+2)%3;
			jk[0] = j;	jk[1] = k;
			for(m=0;m<2;m++){
				axpbyz_array(2, tij[k][m]*lambda[jk[ii]], grd_lambda[j], tij[k][m]*lambda[j], grd_lambda[jk[ii]], phi1+2*m);
				axpbyz_array(2, tij[j][m]*lambda[jk[ii]], grd_lambda[k], tij[j][m]*lambda[k], grd_lambda[jk[ii]], phi2+2*m);
			}
			axpyz_array(4, -1, phi2, phi1, phi);
			ax_array(4, eorien[i]/height[i], phi);
		}
		else{
			i = index - 3*dop;
			j = (i+1)%3;
			k = (i+2)%3;
			for(m=0;m<2;m++){
				axpbyz_array(2, tij[k][m]*lambda[i], grd_lambda[j], tij[k][m]*lambda[j], grd_lambda[i], phi1+2*m);
				axpbyz_array(2, tij[j][m]*lambda[i], grd_lambda[k], tij[j][m]*lambda[k], grd_lambda[i], phi2+2*m);
			}
			axpyz_array(4, -1, phi2, phi1, phi);
			ax_array(4, 1./height[i], phi);
		}
	} // dop=2

	else if (dop == 3){
		if (index < 3*dop){
			in = index;
			i = in / dop;
			ii = in % dop;
			j = (i+1)%3;
			k = (i+2)%3;
			if(ii == 0){
				for(m=0;m<2;m++){
					axpbyz_array(2, -tij[j][m]*lambda[j]*lambda[j], grd_lambda[k], -tij[j][m]*lambda[k]*lambda[j]*2, grd_lambda[j], phi+2*m);
					axpy_array(2, tij[k][m]*lambda[j]*lambda[j]*3, grd_lambda[j], phi+2*m);
				}
			}
			else if(ii == dop -1){
				for(m=0;m<2;m++){
					axpbyz_array(2, tij[k][m]*lambda[k]*lambda[k], grd_lambda[j], tij[k][m]*lambda[j]*lambda[k]*2, grd_lambda[k], phi+2*m);
					axpy_array(2, -tij[j][m]*lambda[k]*lambda[k]*3, grd_lambda[k], phi+2*m);
				}
			}
			else{
				for(m=0;m<2;m++){
					axpbyz_array(2, tij[k][m]*lambda[k]*lambda[j]*2, grd_lambda[j], tij[k][m]*lambda[j]*lambda[j], grd_lambda[k], phi1+2*m);
					axpbyz_array(2, -tij[j][m]*lambda[j]*lambda[k]*2, grd_lambda[k], -tij[j][m]*lambda[k]*lambda[k], grd_lambda[j], phi2+2*m);
				}
				axpyz_array(4, -1, phi2, phi1, phi);
			}
			ax_array(4, eorien[i]/height[i], phi);
		}
		else{
			in = index - 3*dop;
			i = in / 3;
			ii = in % 3;
			j = (i+1)%3;
			k = (i+2)%3;
			for(m=0;m<2;m++){
				axpbyz_array(2, tij[k][m]*lambda[ii]*lambda[i], grd_lambda[j], tij[k][m]*lambda[ii]*lambda[j], grd_lambda[i], phi1+2*m);
				axpy_array(2, tij[k][m]*lambda[i]*lambda[j], grd_lambda[ii], phi1+2*m);
				axpbyz_array(2, tij[j][m]*lambda[ii]*lambda[i], grd_lambda[k], tij[j][m]*lambda[ii]*lambda[k], grd_lambda[i], phi2+2*m);
				axpy_array(2, tij[j][m]*lambda[i]*lambda[k], grd_lambda[ii], phi2+2*m);
			}
			axpyz_array(4, -1, phi2, phi1, phi);
			ax_array(4, 1.0/height[i], phi);
		}
	} // dop=3

	else if (dop >= 4){
		if (index < 3*dop){
			in = index;
			i = in / dop;
			ii = in % dop;
			j = (i+1)%3;
			k = (i+2)%3;
			if(ii == 0){
				for(m=0;m<2;m++){
					axpbyz_array(2, -tij[j][m]*pow(lambda[j], dop-1), grd_lambda[k], -tij[j][m]*lambda[k]*pow(lambda[j], dop-2)*(dop-1), grd_lambda[j], phi+2*m);
					axpy_array(2, tij[k][m]*pow(lambda[j], dop-1)*dop, grd_lambda[j], phi+2*m);
				}
			}
			else if(ii == dop -1){
				for(m=0;m<2;m++){
					axpbyz_array(2, tij[k][m]*pow(lambda[k], dop-1), grd_lambda[j], tij[k][m]*lambda[j]*pow(lambda[k], dop-2)*(dop-1), grd_lambda[k], phi+2*m);
					axpy_array(2, -tij[j][m]*pow(lambda[k], dop-1)*dop, grd_lambda[k], phi+2*m);
				}
			}
			else{
				for(m=0;m<2;m++){
					axpbyz_array(2, tij[k][m]*pow(lambda[k], ii)*pow(lambda[j], dop-1-ii)*(dop-ii), grd_lambda[j], tij[k][m]*pow(lambda[j], dop-ii)*pow(lambda[k], ii-1)*ii, grd_lambda[k], phi1+2*m);
					axpbyz_array(2, -tij[j][m]*pow(lambda[j], dop-1-ii)*pow(lambda[k], ii)*(ii+1), grd_lambda[k], -tij[j][m]*pow(lambda[k], ii+1)*pow(lambda[j], dop-2-ii)*(dop-1-ii), grd_lambda[j], phi2+2*m);
				}
				axpyz_array(4, -1, phi2, phi1, phi);
			}
			ax_array(4, eorien[i]/height[i], phi);
		}
		else{
			in = index - 3*dop;
			i = in / (dop*(dop-1)/2);
			ii = in % (dop*(dop-1)/2);
			j = (i+1)%3;
			k = (i+2)%3;
			bernstein2d_basis(lambda, ii, dop-2, &c);
			bernstein2d_basis1(lambda, grd_lambda, ii, dop-2, val);
			for(m=0;m<2;m++){
				axpbyz_array(2, tij[k][m]*c*lambda[i], grd_lambda[j], tij[k][m]*c*lambda[j], grd_lambda[i], phi1+2*m);				
				axpy_array(2, tij[k][m]*lambda[i]*lambda[j], val, phi1+2*m);
				axpbyz_array(2, tij[j][m]*c*lambda[i], grd_lambda[k], tij[j][m]*c*lambda[k], grd_lambda[i], phi2+2*m);
				axpy_array(2, tij[j][m]*lambda[i]*lambda[k], val, phi2+2*m);
			}
			axpyz_array(4, -1, phi2, phi1, phi);
			ax_array(4, 1.0/height[i], phi);
		}
	} // dop >= 4
}

/** 
 * \fn void bdm_basis(double *lambda, double *height, double **tij, double **nv, double **tv, short *eorien, int index, int dop, double phi[2])
 * \param lambda the barycentric coordinate
 * \param height pointer to height to three edges
 * \param tij pointer to three edge vectors
 * \param eorien pointer to orient of three edges
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param phi[2] basis function of Brezzi–Douglas–Marini element: (phi1, phi2)
 * \return void
 */
void bdm_basis(double *lambda, double *height, double **tij, double **nv, double **tv, short *eorien, int index, int dop, double phi[2])
{
	int dofs = (dop + 1)*(dop + 2); // degrees of freedom
	
	init_array(2, phi, 0);
		
	if (dop<1) return;

	if (index >= dofs || index<0)	return;

	int i, i1, i2, j, k;
	int in, ie, ii;
	double val;

	if (index < 3*(dop+1)){
		in = index;
		i = in / (dop+1);
		ii = in % (dop+1);
		j = (i+1)%3;
		k = (i+2)%3;
		if(ii==0){
			lagrange_basis(lambda, j, dop, &val);
			axy_array(2, eorien[i]*val/height[i], tij[k], phi);
		}
		else if(ii==dop){
			lagrange_basis(lambda, k, dop, &val);
			axy_array(2, -eorien[i]*val/height[i], tij[j], phi);
		}
		else{
			lagrange_basis(lambda, 3 + i*(dop-1) + ii -1, dop, &val);
			axy_array(2, eorien[i]*val, nv[i], phi);
		}
	}
	else if(index < 6*dop){
		in = index - 3*(dop+1);
		i = in / (dop-1);
		ii = in % (dop-1);
		lagrange_basis(lambda, 3 + i*(dop-1) + ii, dop, &val);
		axy_array(2, val, tv[i], phi);
	}
	else{
		in = index - 6*dop;
		i = in / ((dop-1)*(dop-2)/2);
		ii = in % ((dop-1)*(dop-2)/2);
		lagrange_basis(lambda, 3*dop + ii, dop, &phi[i]);
	}
}

/** 
 * \fn void bdm_basisDIV(double *lambda, double **gradLambda, double *height, double **tij, double **nv, double **tv, short *eorien, int index, int dop, double *phi)
 * \param lambda the barycentric coordinate
 * \param height pointer to height to three edges
 * \param tij pointer to three edge vectors
 * \param eorien pointer to orient of three edges
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param phi[2] basis function of Brezzi–Douglas–Marini element: (phi1, phi2)
 * \return void
 */
void bdm_basisDIV(double *lambda, double **gradLambda, double *height, double **tij, double **nv, double **tv, short *eorien, int index, int dop, double *phi)
{
	int dofs = (dop + 1)*(dop + 2); // degrees of freedom
	
	*phi = 0;
		
	if (dop<1) return;

	if (index >= dofs || index<0)	return;

	int i, i1, i2, j, k;
	int in, ie, ii;
	double val[2];

	if (index < 3*(dop+1)){
		in = index;
		i = in / (dop+1);
		ii = in % (dop+1);
		j = (i+1)%3;
		k = (i+2)%3;
		if(ii==0){
			lagrange_basis1(lambda, gradLambda, j, dop, val);
			*phi = eorien[i]*dot_array(2, val, tij[k])/height[i];
		}
		else if(ii==dop){
			lagrange_basis1(lambda, gradLambda, k, dop, val);
			*phi = -eorien[i]*dot_array(2, val, tij[j])/height[i];
		}
		else{
			lagrange_basis1(lambda, gradLambda, 3 + i*(dop-1) + ii -1, dop, val);
			*phi = eorien[i]*dot_array(2, val, nv[i]);
		}
	}
	else if(index < 6*dop){
		in = index - 3*(dop+1);
		i = in / (dop-1);
		ii = in % (dop-1);
		lagrange_basis1(lambda, gradLambda, 3 + i*(dop-1) + ii, dop, val);
		*phi = dot_array(2, val, tv[i]);
	}
	else{
		in = index - 6*dop;
		i = in / ((dop-1)*(dop-2)/2);
		ii = in % ((dop-1)*(dop-2)/2);
		lagrange_basis1(lambda, gradLambda, 3*dop + ii, dop, val);
		*phi = val[i];
	}
}

/** 
 * \fn void bdm_basis1(double *lambda, double **gradLambda, double *height, double **tij, double **nv, double **tv, short *eorien, int index, int dop, double phi[4])
 * \param lambda the barycentric coordinate
 * \param height pointer to height to three edges
 * \param tij pointer to three edge vectors
 * \param eorien pointer to orient of three edges
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param phi[2] basis function of Brezzi–Douglas–Marini element: (phi1, phi2)
 * \return void
 */
void bdm_basis1(double *lambda, double **gradLambda, double *height, double **tij, double **nv, double **tv, short *eorien, int index, int dop, double phi[4])
{
	int dofs = (dop + 1)*(dop + 2); // degrees of freedom
	
	init_array(4, phi, 0);
	
	if (dop<1) return;

	if (index >= dofs || index<0)	return;

	int i, i1, i2, j, k, m;
	int in, ie, ii;
	double val[2];

	if (index < 3*(dop+1)){
		in = index;
		i = in / (dop+1);
		ii = in % (dop+1);
		j = (i+1)%3;
		k = (i+2)%3;
		if(ii==0){
			lagrange_basis1(lambda, gradLambda, j, dop, val);
			for(m=0;m<2;m++)
				axy_array(2, tij[k][m], val, phi+2*m);
			ax_array(4, eorien[i]/height[i], phi);
		}
		else if(ii==dop){
			lagrange_basis1(lambda, gradLambda, k, dop, val);
			for(m=0;m<2;m++)
				axy_array(2, tij[j][m], val, phi+2*m);
			ax_array(4, -eorien[i]/height[i], phi);
		}
		else{
			lagrange_basis1(lambda, gradLambda, 3 + i*(dop-1) + ii -1, dop, val);
			for(m=0;m<2;m++)
				axy_array(2, nv[i][m], val, phi+2*m);
			ax_array(4, eorien[i], phi);
		}
	}
	else if(index < 6*dop){
		in = index - 3*(dop+1);
		i = in / (dop-1);
		ii = in % (dop-1);
		lagrange_basis1(lambda, gradLambda, 3 + i*(dop-1) + ii, dop, val);
		for(m=0;m<2;m++)
			axy_array(2, tv[i][m], val, phi+2*m);
	}
	else{
		in = index - 6*dop;
		i = in / ((dop-1)*(dop-2)/2);
		ii = in % ((dop-1)*(dop-2)/2);
		lagrange_basis1(lambda, gradLambda, 3*dop + ii, dop, val);
		copy_array(2, val, phi+2*i);
	}
}

/**
* \fn void bdm_basis(double *lambda, double s, double eta[3], double xi[3], double **nv, double **nve, int index, int dop, double phi[2])
* \brief basis function of Brezzi-Douglas-Marini element: (phi1, phi2)
* Arnold, D. N.; Falk, R. S. & Winther, R. Geometric decompositions and local bases for spaces of finite element differential forms Comput. Methods Appl. Mech. Engrg., 2009, 198, 1660-1672
* \param x the x-axis coordiante
* \param y the y-axis coordiante
* \param (*T)[2] point the coordiantes of all vertices of current element
* \param s the area of the triangule
* \param eta[3] some auxiliary parameter
* \param xi[3] some auxiliary parameter
* \param orient[3] some auxiliary parameter
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param phi[2] basis function of Raviart-Thomas element: (phi1, phi2)
* \return void
*/
void bdm_basisOld(double *lambda, double s, double eta[3], double xi[3], double **nv, double **nve, int index, int dop, double phi[2])
{
	if (dop<1)
	{
		phi[0] = 0;
		phi[1] = 0;
		return;
	}

	if (index >= (dop + 1)*(dop + 2) || index<0)
	{
		phi[0] = 0;
		phi[1] = 0;
		return;
	}

	int i, i1, i2, j, k;
	int in, ie, ii;
	double orient[3], PHI[3][2], PSI[2][2];
	for (i = 0; i < 3; i++)
		orient[i] = nv[i][0] * nve[i][0] + nv[i][1] * nve[i][1];

	if (dop == 1)
	{
		in = index;
		ie = in / 2;
		ii = in % 2;
		i = (ie + 1) % 3;
		j = (ie + 2) % 3;
		if (ii == 0)
		{
			phi[0] = orient[ie] * lambda[i] * xi[j] / (2 * s);
			phi[1] = orient[ie] * lambda[i] * eta[j] / (2 * s);
		}
		else
		{
			phi[0] = -orient[ie] * lambda[j] * xi[i] / (2 * s);
			phi[1] = -orient[ie] * lambda[j] * eta[i] / (2 * s);
		}
	}

	else if (dop == 2)
	{
		if (index<9)
		{
			in = index;
			ie = in / 3;
			ii = in % 3;
			i = (ie + 1) % 3;
			j = (ie + 2) % 3;
			if (ii == 0)
			{
				phi[0] = orient[ie] * lambda[i] * lambda[i] * xi[j] / (2 * s);
				phi[1] = orient[ie] * lambda[i] * lambda[i] * eta[j] / (2 * s);
			}
			else if (ii == 1)
			{
				phi[0] = orient[ie] * lambda[i] * lambda[j] * (xi[j] - xi[i]) / (2 * s);
				phi[1] = orient[ie] * lambda[i] * lambda[j] * (eta[j] - eta[i]) / (2 * s);
			}
			else
			{
				phi[0] = -orient[ie] * lambda[j] * lambda[j] * xi[i] / (2 * s);
				phi[1] = -orient[ie] * lambda[j] * lambda[j] * eta[i] / (2 * s);
			}
		}
		else
		{
			i = index - 9;
			j = (i + 1) % 3;
			k = (i + 2) % 3;
			phi[0] = lambda[j] * lambda[k] * xi[i] / (2 * s);
			phi[1] = lambda[j] * lambda[k] * eta[i] / (2 * s);
		}
	} // dop=2

	else if (dop == 3)
	{
		if (index<12)
		{
			in = index;
			ie = in / 4;
			ii = in % 4;
			i = (ie + 1) % 3;
			j = (ie + 2) % 3;
			if (ii == 0)
			{
				phi[0] = orient[ie] * lambda[i] * lambda[i] * lambda[i] * xi[j] / (2 * s);
				phi[1] = orient[ie] * lambda[i] * lambda[i] * lambda[i] * eta[j] / (2 * s);
			}
			else if (ii == 1)
			{
				phi[0] = orient[ie] * lambda[i] * lambda[i] * lambda[j] * (2 * xi[j] - xi[i]) / (2 * s);
				phi[1] = orient[ie] * lambda[i] * lambda[i] * lambda[j] * (2 * eta[j] - eta[i]) / (2 * s);
			}
			else if (ii == 2)
			{
				phi[0] = orient[ie] * lambda[i] * lambda[j] * lambda[j] * (xi[j] - 2 * xi[i]) / (2 * s);
				phi[1] = orient[ie] * lambda[i] * lambda[j] * lambda[j] * (eta[j] - 2 * eta[i]) / (2 * s);
			}
			else
			{
				phi[0] = -orient[ie] * lambda[j] * lambda[j] * lambda[j] * xi[i] / (2 * s);
				phi[1] = -orient[ie] * lambda[j] * lambda[j] * lambda[j] * eta[i] / (2 * s);
			}
		}
		else
		{
			in = index - 12;
			ie = in / 3;
			ii = in % 3;
			i = ie;
			j = (i + 1) % 3;
			k = (i + 2) % 3;
			phi[0] = lambda[ii] * lambda[j] * lambda[k] * xi[i] / (2 * s);
			phi[1] = lambda[ii] * lambda[j] * lambda[k] * eta[i] / (2 * s);
		}
	} // dop=2

	else
	{
		phi[0] = 0;
		phi[1] = 0;
	}
}

/**
* \fn void bdm_basis1(double *lambda, double s, double eta[3], double xi[3], double **nv, double **nve, int index, int dop, double phi[4])
* \brief the first order derivative of Brezzi-Douglas-Marini element basis function: (\partial_{x}phi1, \partial_{y}phi1, \partial_{x}phi2, \partial_{y}phi2)
* Arnold, D. N.; Falk, R. S. & Winther, R. Geometric decompositions and local bases for spaces of finite element differential forms Comput. Methods Appl. Mech. Engrg., 2009, 198, 1660-1672
* \param x the x-axis coordiante
* \param y the y-axis coordiante
* \param (*T)[2] point the coordiantes of all vertices of current element
* \param s the area of the triangule
* \param eta[3] some auxiliary parameter
* \param xi[3] some auxiliary parameter
* \param orient[3] some auxiliary parameter
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param phi[4] the first order derivative of Raviart-Thomas element basis function: (\partial_{x}phi1, \partial_{y}phi1, \partial_{x}phi2, \partial_{y}phi2)
* \return void
*/
void bdm_basis1Old(double *lambda, double s, double eta[3], double xi[3], double **nv, double **nve, int index, int dop, double phi[4])
{
	if (dop<1)
	{
		phi[0] = 0;
		phi[1] = 0;
		phi[2] = 0;
		phi[3] = 0;
		return;
	}

	if (index >= (dop + 1)*(dop + 2) || index<0)
	{
		phi[0] = 0;
		phi[1] = 0;
		phi[2] = 0;
		phi[3] = 0;
		return;
	}

	int i, i1, i2, j, k;
	int in, ie, ii;
	double orient[3], PHI[3][2], PSI[2][2], PHI1[3][4], PSI1[2][4], grad[2];
	for (i = 0; i < 3; i++)
		orient[i] = nv[i][0] * nve[i][0] + nv[i][1] * nve[i][1];

	if (dop == 1)
	{
		in = index;
		ie = in / 2;
		ii = in % 2;
		i = (ie + 1) % 3;
		j = (ie + 2) % 3;
		if (ii == 0)
		{
			grad[0] = eta[i] / (2 * s);
			grad[1] = -xi[i] / (2 * s);
			phi[0] = orient[ie] * grad[0] * xi[j] / (2 * s);
			phi[1] = orient[ie] * grad[1] * xi[j] / (2 * s);
			phi[2] = orient[ie] * grad[0] * eta[j] / (2 * s);
			phi[3] = orient[ie] * grad[1] * eta[j] / (2 * s);
		}
		else
		{
			grad[0] = eta[j] / (2 * s);
			grad[1] = -xi[j] / (2 * s);
			phi[0] = -orient[ie] * grad[0] * xi[i] / (2 * s);
			phi[1] = -orient[ie] * grad[1] * xi[i] / (2 * s);
			phi[2] = -orient[ie] * grad[0] * eta[i] / (2 * s);
			phi[3] = -orient[ie] * grad[1] * eta[i] / (2 * s);
		}
	}

	else if (dop == 2)
	{
		if (index<9)
		{
			in = index;
			ie = in / 3;
			ii = in % 3;
			i = (ie + 1) % 3;
			j = (ie + 2) % 3;
			if (ii == 0)
			{
				grad[0] = 2 * lambda[i] * eta[i] / (2 * s);
				grad[1] = -2 * lambda[i] * xi[i] / (2 * s);
				phi[0] = orient[ie] * grad[0] * xi[j] / (2 * s);
				phi[1] = orient[ie] * grad[1] * xi[j] / (2 * s);
				phi[2] = orient[ie] * grad[0] * eta[j] / (2 * s);
				phi[3] = orient[ie] * grad[1] * eta[j] / (2 * s);
			}
			else if (ii == 1)
			{
				grad[0] = (lambda[i] * eta[j] + lambda[j] * eta[i]) / (2 * s);
				grad[1] = -(lambda[i] * xi[j] + lambda[j] * xi[i]) / (2 * s);
				phi[0] = orient[ie] * grad[0] * (xi[j] - xi[i]) / (2 * s);
				phi[1] = orient[ie] * grad[1] * (xi[j] - xi[i]) / (2 * s);
				phi[2] = orient[ie] * grad[0] * (eta[j] - eta[i]) / (2 * s);
				phi[3] = orient[ie] * grad[1] * (eta[j] - eta[i]) / (2 * s);
			}
			else
			{
				grad[0] = 2 * lambda[j] * eta[j] / (2 * s);
				grad[1] = -2 * lambda[j] * xi[j] / (2 * s);
				phi[0] = -orient[ie] * grad[0] * xi[i] / (2 * s);
				phi[1] = -orient[ie] * grad[1] * xi[i] / (2 * s);
				phi[2] = -orient[ie] * grad[0] * eta[i] / (2 * s);
				phi[3] = -orient[ie] * grad[1] * eta[i] / (2 * s);
			}
		}
		else
		{
			i = index - 9;
			j = (i + 1) % 3;
			k = (i + 2) % 3;
			grad[0] = (lambda[k] * eta[j] + lambda[j] * eta[k]) / (2 * s);
			grad[1] = -(lambda[k] * xi[j] + lambda[j] * xi[k]) / (2 * s);
			phi[0] = grad[0] * xi[i] / (2 * s);
			phi[1] = grad[1] * xi[i] / (2 * s);
			phi[2] = grad[0] * eta[i] / (2 * s);
			phi[3] = grad[1] * eta[i] / (2 * s);
		}
	} // dop=2

	else if (dop == 3)
	{
		if (index<12)
		{
			in = index;
			ie = in / 4;
			ii = in % 4;
			i = (ie + 1) % 3;
			j = (ie + 2) % 3;
			if (ii == 0)
			{
				grad[0] = 3 * lambda[i] * lambda[i] * eta[i] / (2 * s);
				grad[1] = -3 * lambda[i] * lambda[i] * xi[i] / (2 * s);
				phi[0] = orient[ie] * grad[0] * xi[j] / (2 * s);
				phi[1] = orient[ie] * grad[1] * xi[j] / (2 * s);
				phi[2] = orient[ie] * grad[0] * eta[j] / (2 * s);
				phi[3] = orient[ie] * grad[1] * eta[j] / (2 * s);
			}
			else if (ii == 1)
			{
				grad[0] = (2 * lambda[i] * lambda[j] * eta[i] + lambda[i] * lambda[i] * eta[j]) / (2 * s);
				grad[1] = -(2 * lambda[i] * lambda[j] * xi[i] + lambda[i] * lambda[i] * xi[j]) / (2 * s);
				phi[0] = orient[ie] * grad[0] * (2 * xi[j] - xi[i]) / (2 * s);
				phi[1] = orient[ie] * grad[1] * (2 * xi[j] - xi[i]) / (2 * s);
				phi[2] = orient[ie] * grad[0] * (2 * eta[j] - eta[i]) / (2 * s);
				phi[3] = orient[ie] * grad[1] * (2 * eta[j] - eta[i]) / (2 * s);
			}
			else if (ii == 2)
			{
				grad[0] = (2 * lambda[i] * lambda[j] * eta[j] + lambda[j] * lambda[j] * eta[i]) / (2 * s);
				grad[1] = -(2 * lambda[i] * lambda[j] * xi[j] + lambda[j] * lambda[j] * xi[i]) / (2 * s);
				phi[0] = orient[ie] * grad[0] * (xi[j] - 2 * xi[i]) / (2 * s);
				phi[1] = orient[ie] * grad[1] * (xi[j] - 2 * xi[i]) / (2 * s);
				phi[2] = orient[ie] * grad[0] * (eta[j] - 2 * eta[i]) / (2 * s);
				phi[3] = orient[ie] * grad[1] * (eta[j] - 2 * eta[i]) / (2 * s);
			}
			else
			{
				grad[0] = 3 * lambda[j] * lambda[j] * eta[j] / (2 * s);
				grad[1] = -3 * lambda[j] * lambda[j] * xi[j] / (2 * s);
				phi[0] = -orient[ie] * grad[0] * xi[i] / (2 * s);
				phi[1] = -orient[ie] * grad[1] * xi[i] / (2 * s);
				phi[2] = -orient[ie] * grad[0] * eta[i] / (2 * s);
				phi[3] = -orient[ie] * grad[1] * eta[i] / (2 * s);
			}
		}
		else
		{
			in = index - 12;
			ie = in / 3;
			ii = in % 3;
			i = ie;
			j = (i + 1) % 3;
			k = (i + 2) % 3;
			grad[0] = (eta[ii] * lambda[j] * lambda[k] + lambda[ii] * eta[j] * lambda[k] + lambda[ii] * lambda[j] * eta[k]) / (2 * s);
			grad[1] = -(xi[ii] * lambda[j] * lambda[k] + lambda[ii] * xi[j] * lambda[k] + lambda[ii] * lambda[j] * xi[k]) / (2 * s);
			phi[0] = grad[0] * xi[i] / (2 * s);
			phi[1] = grad[1] * xi[i] / (2 * s);
			phi[2] = grad[0] * eta[i] / (2 * s);
			phi[3] = grad[1] * eta[i] / (2 * s);
		}
	} // dop=2

	else
	{
		phi[0] = 0;
		phi[1] = 0;
		phi[2] = 0;
		phi[3] = 0;
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
* \param gradLambda gradient of the barycentric coordinate
* \param **nv the unit normal vectors of the three edges
* \param **tv the unit tangential vectors of the three edges
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param double(*phi)[2] the first order derivative of basis function of Hu-Zhang element
* \return void
*/
void huzhang_basis1(double *lambda, double **gradLambda, double **nv, double **tv, int index, int dop, double(*phi)[2])
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
			lagrange_basis1(lambda, gradLambda, index % 3, dop, val);
			phi[index / 3][0] = val[0]; 
			phi[index / 3][1] = val[1];
		}
		else if (index<9 + (dop - 1) * 3)
		{
			i = index - 9;
			lagrange_basis1(lambda, gradLambda, 3 + i, dop, val);
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
			lagrange_basis1(lambda, gradLambda, 3 + i, dop, val);
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
			lagrange_basis1(lambda, gradLambda, 3 + i, dop, val);
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
			lagrange_basis1(lambda, gradLambda, dop * 3 + i % (dofs - dop * 3), dop, val);
			phi[i / (dofs - dop * 3)][0] = val[0];
			phi[i / (dofs - dop * 3)][1] = val[1];
		}
	}

}

/** 
 * \fn void huzhang_basisDIV(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, int index, int dop, double phi[2])
 * \brief divergence of basis function of Hu-Zhang element
 * \param *lambda pointer to the area coordiante
 * \param gradLambda gradient of the barycentric coordinate
 * \param **nv the unit normal vectors of the three edges
 * \param **tv the unit tangential vectors of the three edges
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param *phi divergence of basis function
 * \return void
 */
void huzhang_basisDIV(double *lambda, double **gradLambda, double **nv, double **tv, int index, int dop, double phi[2])
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
			lagrange_basis1(lambda, gradLambda, index, dop, val);
			phi[0] = val[0];
			phi[1] = 0;
		}
		else if(index<6)
		{
			lagrange_basis1(lambda, gradLambda, index%3, dop, val);
			phi[0] = 0;
			phi[1] = val[1];
		}
		else if(index<9)
		{
			lagrange_basis1(lambda, gradLambda, index%3, dop, val);
			phi[0] = val[1];
			phi[1] = val[0];
		}
		else if(index<9+(dop-1)*3)
		{
			i=index-9;
			lagrange_basis1(lambda, gradLambda, 3+i, dop, val);
			phi[0]=val[0]*nn[i/(dop-1)][0]+val[1]*nn[i/(dop-1)][2];
			phi[1]=val[0]*nn[i/(dop-1)][2]+val[1]*nn[i/(dop-1)][1];
		}
		else if(index<9+(dop-1)*6)
		{
			i=index-9-(dop-1)*3;
			lagrange_basis1(lambda, gradLambda, 3+i, dop, val);
			phi[0]=val[0]*nt[i/(dop-1)][0]+val[1]*nt[i/(dop-1)][2];
			phi[1]=val[0]*nt[i/(dop-1)][2]+val[1]*nt[i/(dop-1)][1];
		}
		else if(index<9+(dop-1)*9)
		{
			i=index-9-(dop-1)*6;
			lagrange_basis1(lambda, gradLambda, 3+i, dop, val);
			phi[0]=val[0]*tt[i/(dop-1)][0]+val[1]*tt[i/(dop-1)][2];
			phi[1]=val[0]*tt[i/(dop-1)][2]+val[1]*tt[i/(dop-1)][1];
		}
		else if(index<dop*9+(dofs-dop*3))
		{
			i=index-dop*9;
			lagrange_basis1(lambda, gradLambda, dop*3+i, dop, val);
			phi[0] = val[0];
			phi[1] = 0;
		}
		else if(index<dop*9+(dofs-dop*3)*2)
		{
			i=index-dop*9-(dofs-dop*3);
			lagrange_basis1(lambda, gradLambda, dop*3+i, dop, val);
			phi[0] = 0;
			phi[1] = val[1];
		}
		else
		{
			i=index-dop*9-(dofs-dop*3)*2;
			lagrange_basis1(lambda, gradLambda, dop*3+i, dop, val);
			phi[0] = val[1];
			phi[1] = val[0];
		}
	} 
}

/**
* \fn void huzhang_basisROT(double *lambda, double **gradLambda, double **nv, double **tv, int index, int dop, double phi[2])
* \brief rotation of basis function of Hu-Zhang element
* \param *lambda pointer to the area coordiante
* \param gradLambda gradient of the barycentric coordinate
* \param **nv the unit normal vectors of the three edges
* \param **tv the unit tangential vectors of the three edges
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi divergence of basis function
* \return void
*/
void huzhang_basisROT(double *lambda, double **gradLambda, double **nv, double **tv, int index, int dop, double phi[2])
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
			lagrange_basis1(lambda, gradLambda, index, dop, val);
			phi[0] = -val[1];
			phi[1] = 0;
		}
		else if (index<6)
		{
			lagrange_basis1(lambda, gradLambda, index % 3, dop, val);
			phi[0] = 0;
			phi[1] = val[0];
		}
		else if (index<9)
		{
			lagrange_basis1(lambda, gradLambda, index % 3, dop, val);
			phi[0] = val[0];
			phi[1] = -val[1];
		}
		else if (index<9 + (dop - 1) * 3)
		{
			i = index - 9;
			lagrange_basis1(lambda, gradLambda, 3 + i, dop, val);
			phi[0] = -val[1] * nn[i / (dop - 1)][0] + val[0] * nn[i / (dop - 1)][2];
			phi[1] = -val[1] * nn[i / (dop - 1)][2] + val[0] * nn[i / (dop - 1)][1];
		}
		else if (index<9 + (dop - 1) * 6)
		{
			i = index - 9 - (dop - 1) * 3;
			lagrange_basis1(lambda, gradLambda, 3 + i, dop, val);
			phi[0] = -val[1] * nt[i / (dop - 1)][0] + val[0] * nt[i / (dop - 1)][2];
			phi[1] = -val[1] * nt[i / (dop - 1)][2] + val[0] * nt[i / (dop - 1)][1];
		}
		else if (index<9 + (dop - 1) * 9)
		{
			i = index - 9 - (dop - 1) * 6;
			lagrange_basis1(lambda, gradLambda, 3 + i, dop, val);
			phi[0] = -val[1] * tt[i / (dop - 1)][0] + val[0] * tt[i / (dop - 1)][2];
			phi[1] = -val[1] * tt[i / (dop - 1)][2] + val[0] * tt[i / (dop - 1)][1];
		}
		else if (index<dop * 9 + (dofs - dop * 3))
		{
			i = index - dop * 9;
			lagrange_basis1(lambda, gradLambda, dop * 3 + i, dop, val);
			phi[0] = -val[1];
			phi[1] = 0;
		}
		else if (index<dop * 9 + (dofs - dop * 3) * 2)
		{
			i = index - dop * 9 - (dofs - dop * 3);
			lagrange_basis1(lambda, gradLambda, dop * 3 + i, dop, val);
			phi[0] = 0;
			phi[1] = val[0];
		}
		else
		{
			i = index - dop * 9 - (dofs - dop * 3) * 2;
			lagrange_basis1(lambda, gradLambda, dop * 3 + i, dop, val);
			phi[0] = val[0];
			phi[1] = -val[1];
		}
	}
}

/**
* \fn void huzhang_basisCurlTrace(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, int index, int dop, double phi[2])
* \brief curl of the trace of basis function of Hu-Zhang element
* \param *lambda pointer to the area coordiante
* \param gradLambda gradient of the barycentric coordinate
* \param **nv the unit normal vectors of the three edges
* \param **tv the unit tangential vectors of the three edges
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi divergence of basis function
* \return void
*/
void huzhang_basisCurlTrace(double *lambda, double **gradLambda, double **nv, double **tv, int index, int dop, double phi[2])
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
			lagrange_basis1(lambda, gradLambda, index, dop, val);
			phi[0] = -val[1];
			phi[1] = val[0];
		}
		else if (index<6)
		{
			lagrange_basis1(lambda, gradLambda, index % 3, dop, val);
			phi[0] = -val[1];
			phi[1] = val[0];
		}
		else if (index<9)
		{
			lagrange_basis1(lambda, gradLambda, index % 3, dop, val);
			phi[0] = 0;
			phi[1] = 0;
		}
		else if (index<9 + (dop - 1) * 3)
		{
			i = index - 9;
			lagrange_basis1(lambda, gradLambda, 3 + i, dop, val);
			phi[0] = -val[1] * (nn[i / (dop - 1)][0] + nn[i / (dop - 1)][1]);
			phi[1] = val[0] * (nn[i / (dop - 1)][0] + nn[i / (dop - 1)][1]);
		}
		else if (index<9 + (dop - 1) * 6)
		{
			i = index - 9 - (dop - 1) * 3;
			lagrange_basis1(lambda, gradLambda, 3 + i, dop, val);
			phi[0] = -val[1] * (nt[i / (dop - 1)][0] + nt[i / (dop - 1)][1]);
			phi[1] = val[0] * (nt[i / (dop - 1)][0] + nt[i / (dop - 1)][1]);
		}
		else if (index<9 + (dop - 1) * 9)
		{
			i = index - 9 - (dop - 1) * 6;
			lagrange_basis1(lambda, gradLambda, 3 + i, dop, val);
			phi[0] = -val[1] * (tt[i / (dop - 1)][0] + tt[i / (dop - 1)][1]);
			phi[1] = val[0] * (tt[i / (dop - 1)][0] + tt[i / (dop - 1)][1]);
		}
		else if (index<dop * 9 + (dofs - dop * 3))
		{
			i = index - dop * 9;
			lagrange_basis1(lambda, gradLambda, dop * 3 + i, dop, val);
			phi[0] = -val[1];
			phi[1] = val[0];
		}
		else if (index<dop * 9 + (dofs - dop * 3) * 2)
		{
			i = index - dop * 9 - (dofs - dop * 3);
			lagrange_basis1(lambda, gradLambda, dop * 3 + i, dop, val);
			phi[0] = -val[1];
			phi[1] = val[0];
		}
		else
		{
			i = index - dop * 9 - (dofs - dop * 3) * 2;
			lagrange_basis1(lambda, gradLambda, dop * 3 + i, dop, val);
			phi[0] = 0;
			phi[1] = 0;
		}
	}
}

/**
* \fn void huzhang_basisROTROT(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, int index, int dop, double *phi)
* \brief rotrot of basis function of Hu-Zhang element
* \param *lambda pointer to the area coordiante
* \param gradLambda gradient of the barycentric coordinate
* \param **nv the unit normal vectors of the three edges
* \param **tv the unit tangential vectors of the three edges
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi divergence of basis function
* \return void
*/
void huzhang_basisROTROT(double *lambda, double **gradLambda, double **nv, double **tv, int index, int dop, double *phi)
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
			lagrange_basis2(lambda, gradLambda, index, dop, val);
			*phi = val[1];
		}
		else if (index<6)
		{
			lagrange_basis2(lambda, gradLambda, index % 3, dop, val);
			*phi = val[0];
		}
		else if (index<9)
		{
			lagrange_basis2(lambda, gradLambda, index % 3, dop, val);
			*phi = -2 * val[2];
		}
		else if (index<9 + (dop - 1) * 3)
		{
			i = index - 9;
			lagrange_basis2(lambda, gradLambda, 3 + i, dop, val);
			*phi = val[1] * nn[i / (dop - 1)][0] + val[0] * nn[i / (dop - 1)][1] - 2 * val[2] * nn[i / (dop - 1)][2];
		}
		else if (index<9 + (dop - 1) * 6)
		{
			i = index - 9 - (dop - 1) * 3;
			lagrange_basis2(lambda, gradLambda, 3 + i, dop, val);
			*phi = val[1] * nt[i / (dop - 1)][0] + val[0] * nt[i / (dop - 1)][1] - 2 * val[2] * nt[i / (dop - 1)][2];
		}
		else if (index<9 + (dop - 1) * 9)
		{
			i = index - 9 - (dop - 1) * 6;
			lagrange_basis2(lambda, gradLambda, 3 + i, dop, val);
			*phi = val[1] * tt[i / (dop - 1)][0] + val[0] * tt[i / (dop - 1)][1] - 2 * val[2] * tt[i / (dop - 1)][2];
		}
		else if (index<dop * 9 + (dofs - dop * 3))
		{
			i = index - dop * 9;
			lagrange_basis2(lambda, gradLambda, dop * 3 + i, dop, val);
			*phi = val[1];
		}
		else if (index<dop * 9 + (dofs - dop * 3) * 2)
		{
			i = index - dop * 9 - (dofs - dop * 3);
			lagrange_basis2(lambda, gradLambda, dop * 3 + i, dop, val);
			*phi = val[0];
		}
		else
		{
			i = index - dop * 9 - (dofs - dop * 3) * 2;
			lagrange_basis2(lambda, gradLambda, dop * 3 + i, dop, val);
			*phi = -2 * val[2];
		}
	}
}

/**
* \fn void huzhang_basisLaplaceTrace(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, int index, int dop, double *phi)
* \brief Laplace of the trace of basis function of Hu-Zhang element
* \param *lambda pointer to the area coordiante
* \param gradLambda gradient of the barycentric coordinate
* \param **nv the unit normal vectors of the three edges
* \param **tv the unit tangential vectors of the three edges
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi divergence of basis function
* \return void
*/
void huzhang_basisLaplaceTrace(double *lambda, double **gradLambda, double **nv, double **tv, int index, int dop, double *phi)
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
			lagrange_basis2(lambda, gradLambda, index, dop, val);
			*phi = val[0] + val[1];
		}
		else if (index<6)
		{
			lagrange_basis2(lambda, gradLambda, index % 3, dop, val);
			*phi = val[0] + val[1];
		}
		else if (index<9)
		{
			lagrange_basis2(lambda, gradLambda, index % 3, dop, val);
			*phi = 0;
		}
		else if (index<9 + (dop - 1) * 3)
		{
			i = index - 9;
			lagrange_basis2(lambda, gradLambda, 3 + i, dop, val);
			*phi = (val[0] + val[1]) * (nn[i / (dop - 1)][0] + nn[i / (dop - 1)][1]);
		}
		else if (index<9 + (dop - 1) * 6)
		{
			i = index - 9 - (dop - 1) * 3;
			lagrange_basis2(lambda, gradLambda, 3 + i, dop, val);
			*phi = (val[0] + val[1]) * (nt[i / (dop - 1)][0] + nt[i / (dop - 1)][1]);
		}
		else if (index<9 + (dop - 1) * 9)
		{
			i = index - 9 - (dop - 1) * 6;
			lagrange_basis2(lambda, gradLambda, 3 + i, dop, val);
			*phi = (val[0] + val[1]) * (tt[i / (dop - 1)][0] + tt[i / (dop - 1)][1]);
		}
		else if (index<dop * 9 + (dofs - dop * 3))
		{
			i = index - dop * 9;
			lagrange_basis2(lambda, gradLambda, dop * 3 + i, dop, val);
			*phi = val[0] + val[1];
		}
		else if (index<dop * 9 + (dofs - dop * 3) * 2)
		{
			i = index - dop * 9 - (dofs - dop * 3);
			lagrange_basis2(lambda, gradLambda, dop * 3 + i, dop, val);
			*phi = val[0] + val[1];
		}
		else
		{
			i = index - dop * 9 - (dofs - dop * 3) * 2;
			lagrange_basis2(lambda, gradLambda, dop * 3 + i, dop, val);
			*phi = 0;
		}
	}
}

/** 
 * \fn void divS_huangzhou_prebasis(double *lambda, double **tij, int i, double phi[3])
 * \brief preparation for basis function of Huang-Zhou element
 * \param *lambda pointer to the area coordiante
 * \param **nv the unit normal vectors of the three edges
 * \param **tv the unit tangential vectors of the three edges
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param *phi basis function
 * \return void
 */
void divS_huangzhou_prebasis(double *lambda, double **tij, int i, double phi[3])
{
	double val, c1, c2;

	phi[0]=0;
	phi[1]=0;
	phi[2]=0;
	
	int j = (i+1)%3; 
	int k = (i+2)%3;

	c1 = dot_array(2, tij[i], tij[k]);
	c2 = dot_array(2, tij[i], tij[j]);

	double T[3][3];
	for(int i1=0;i1<3;i1++){
		T[i1][0]=tij[i1][0]*tij[i1][0]; 
		T[i1][1]=tij[i1][1]*tij[i1][1]; 
		T[i1][2]=tij[i1][0]*tij[i1][1];
	}

	val = -10*c1*lambda[j]*lambda[k]*lambda[k];
	axpy_array(3, val, T[j], phi);

	val = -10*c2*(lambda[k]+lambda[i])*lambda[j]*lambda[k];
	axpy_array(3, val, T[k], phi);

	val = ((6*c2-9*c1)*lambda[j] + (9*c2-6*c1)*lambda[k])*lambda[j]*lambda[k];
	axpy_array(3, val, T[i], phi);

	val = ((6*c2-9*c1)*lambda[j] + (2*c1-3*c2)*lambda[k])*lambda[k]*lambda[i];
	axpy_array(3, val, T[j], phi);

	val = ((3*c1-2*c2)*lambda[j] + (9*c2-6*c1)*lambda[k])*lambda[i]*lambda[j];
	axpy_array(3, val, T[k], phi);
}

/** 
 * \fn void divS_huangzhou_prebasisDIV(double *lambda, double **tij, int i, double phi[2])
 * \brief preparation for basis function of Huang-Zhou element
 * \param *lambda pointer to the area coordiante
 * \param **nv the unit normal vectors of the three edges
 * \param **tv the unit tangential vectors of the three edges
 * \param **tij the tangential vectors from vertex i to vertex j
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param *phi basis function
 * \return void
 */
void divS_huangzhou_prebasisDIV(double *lambda, double **tij, int i, double phi[2])
{
	double val, c1, c2;

	phi[0]=0;
	phi[1]=0;
	
	int j = (i+1)%3; 
	int k = (i+2)%3;

	c1 = dot_array(2, tij[i], tij[k]);
	c2 = dot_array(2, tij[i], tij[j]);

	val = 20*c1*lambda[j]*lambda[k];
	axpy_array(2, val, tij[j], phi);

	val = 10*c2*lambda[k]*(2*lambda[j]-1);
	axpy_array(2, val, tij[k], phi);

	val = ((6*c2-9*c1)*lambda[j] + (9*c2-6*c1)*lambda[k])*(lambda[j]-lambda[k]);
	val += 3*(c1+c2)*lambda[j]*lambda[k];
	axpy_array(2, val, tij[i], phi);

	val = ((6*c2-9*c1)*lambda[j] + (2*c1-3*c2)*lambda[k])*(lambda[k]-lambda[i]);
	val += (3*c2-2*c1)*lambda[k]*lambda[i];
	axpy_array(2, val, tij[j], phi);

	val = ((3*c1-2*c2)*lambda[j] + (9*c2-6*c1)*lambda[k])*(lambda[i]-lambda[j]);
	val += (3*c1-2*c2)*lambda[i]*lambda[j];
	axpy_array(2, val, tij[k], phi);
}

/** 
 * \fn void divS_huangzhou_basis(double *lambda, double s, double **nv, double **tv, double **tij, int index, double phi[3])
 * \brief basis function of Huang-Zhou element
 * \param *lambda pointer to the area coordiante
 * \param s pointer to the area of the triangle
 * \param **nv the unit normal vectors of the three edges
 * \param **tv the unit tangential vectors of the three edges
 * \param **tij the tangential vectors from vertex i to vertex j
 * \param index the indicator of the basis function
 * \param *phi basis function
 * \return void
 */
void divS_huangzhou_basis(double *lambda, double s, double **nv, double **tv, double **tij, int index, double phi[3])
{
	phi[0]=0;
	phi[1]=0;
	phi[2]=0;
	if(index>= 21 || index<0)
		return;

	double val, c1[2], c2[2], phi0[3];
	double nn[3][3], nt[3][3], tt[3][3];
	int i, j, k, ii;
	for(i=0;i<3;i++){
		nn[i][0]=nv[i][0]*nv[i][0]; nn[i][1]=nv[i][1]*nv[i][1]; nn[i][2]=nv[i][0]*nv[i][1];
		nt[i][0]=nv[i][0]*tv[i][0]; nt[i][1]=nv[i][1]*tv[i][1]; nt[i][2]=(nv[i][0]*tv[i][1]+nv[i][1]*tv[i][0])/2;
		tt[i][0]=tv[i][0]*tv[i][0]; tt[i][1]=tv[i][1]*tv[i][1]; tt[i][2]=tv[i][0]*tv[i][1];
	}
	
	c1[0] = 36.0; c2[0] = -3.0/(2.0*s*s);
	c1[1] = -24.0; c2[1] = 3.0/(2.0*s*s);

	if(index<9){
		i = index%3;
		// val = lambda[i]*lambda[i];
		lagrange_basis(lambda, i, 2, &val);
		phi[index/3] = val;
	}
	else if(index<15){
		// i=index-9;
		ii = (index-9)%2;
		i = (index-9)/2;
		j = (i+1)%3; k = (i+2)%3;

		val = lambda[j]*lambda[k];
		axy_array(3, c1[ii]*val, nn[i], phi);
		
		divS_huangzhou_prebasis(lambda, tij, i, phi0);
		axpy_array(3, c2[ii], phi0, phi);
	}
	else if(index<18){
		i=index-15;
		// j = (i+1)%3; k = (i+2)%3;
		// val = lambda[j]*lambda[k];
		lagrange_basis(lambda, 3+i, 2, &val);
		axy_array(3, val, nt[i], phi);
	}
	else{
		i=index-18;
		// j = (i+1)%3; k = (i+2)%3;
		// val = lambda[j]*lambda[k];
		lagrange_basis(lambda, 3+i, 2, &val);
		axy_array(3, val, tt[i], phi);
	}
}

/** 
 * \fn void divS_huangzhou_basisDIV(double *lambda, double **gradLambda, double s, double **nv, double **tv, double **tij, int index, double phi[2])
 * \brief divergence of basis function of Huang-Zhou element
 * \param *lambda pointer to the area coordiante
 * \param gradLambda gradient of the barycentric coordinate
 * \param **nv the unit normal vectors of the three edges
 * \param **tv the unit tangential vectors of the three edges
 * \param index the indicator of the basis function
 * \param *phi divergence of basis function
 * \return void
 */
void divS_huangzhou_basisDIV(double *lambda, double **gradLambda, double s, double **nv, double **tv, double **tij, int index, double phi[2])
{
	phi[0]=0;
	phi[1]=0;
	if(index>= 21 || index<0)
		return;

	double val[2], c1[2], c2[2], phi0[2];
	double nn[3][3], nt[3][3], tt[3][3];
	int i, j, k, ii;
	for(i=0;i<3;i++){
		nn[i][0]=nv[i][0]*nv[i][0]; nn[i][1]=nv[i][1]*nv[i][1]; nn[i][2]=nv[i][0]*nv[i][1];
		nt[i][0]=nv[i][0]*tv[i][0]; nt[i][1]=nv[i][1]*tv[i][1]; nt[i][2]=(nv[i][0]*tv[i][1]+nv[i][1]*tv[i][0])/2;
		tt[i][0]=tv[i][0]*tv[i][0]; tt[i][1]=tv[i][1]*tv[i][1]; tt[i][2]=tv[i][0]*tv[i][1];
	}

	c1[0] = 36.0; c2[0] = -3.0/(2.0*s*s);
	c1[1] = -24.0; c2[1] = 3.0/(2.0*s*s);

	if(index<3){
		i = index;
		// axy_array(2, 2.0*lambda[i], gradLambda[i], val);
		lagrange_basis1(lambda, gradLambda, i, 2, val);
		phi[0] = val[0];
		phi[1] = 0;
	}
	else if(index<6){
		i = index%3;
		// axy_array(2, 2.0*lambda[i], gradLambda[i], val);
		lagrange_basis1(lambda, gradLambda, i, 2, val);
		phi[0] = 0;
		phi[1] = val[1];
	}
	else if(index<9){
		i = index%3;
		// axy_array(2, 2.0*lambda[i], gradLambda[i], val);
		lagrange_basis1(lambda, gradLambda, i, 2, val);
		phi[0] = val[1];
		phi[1] = val[0];
	}
	else if(index<15){
		// i=index-9;
		ii = (index-9)%2;
		i = (index-9)/2;
		j = (i+1)%3; k = (i+2)%3;

		axpbyz_array(2, lambda[j], gradLambda[k], lambda[k], gradLambda[j], val);
		phi[0]= (val[0]*nn[i][0]+val[1]*nn[i][2]) * c1[ii];
		phi[1]= (val[0]*nn[i][2]+val[1]*nn[i][1]) * c1[ii];

		divS_huangzhou_prebasisDIV(lambda, tij, i, phi0);
		axpy_array(2, c2[ii], phi0, phi);
	}
	else if(index<18){
		i=index-15;
		// j = (i+1)%3; k = (i+2)%3;
		// axpbyz_array(2, lambda[j], gradLambda[k], lambda[k], gradLambda[j], val);
		lagrange_basis1(lambda, gradLambda, 3+i, 2, val);
		phi[0]=val[0]*nt[i][0]+val[1]*nt[i][2];
		phi[1]=val[0]*nt[i][2]+val[1]*nt[i][1];
	}
	else{
		i=index-18;
		// j = (i+1)%3; k = (i+2)%3;
		// axpbyz_array(2, lambda[j], gradLambda[k], lambda[k], gradLambda[j], val);
		lagrange_basis1(lambda, gradLambda, 3+i, 2, val);
		phi[0]=val[0]*tt[i][0]+val[1]*tt[i][2];
		phi[1]=val[0]*tt[i][2]+val[1]*tt[i][1];
	}
}


/**
* \fn void mini_basis(double *lambda, int index, double *phi)
* \brief basis function of MINI element for Stokes equation
* \param *lambda pointer to the area coordiante
* \param index the indicator of the basis function
* \param *phi basis function
* \return void
*/
void mini_basis(double *lambda, int index, double *phi)
{
	*phi = 0;
	if (index >= 4 || index<0)
		return;

	// the last dof is the function value at the centroid of the triangle
	if (index < 3)
		*phi = lambda[index] - 9.0*lambda[0] * lambda[1] * lambda[2];
	else
		*phi = 27.0*lambda[0] * lambda[1] * lambda[2];
	/********* the last dof is the mean value of the function on triangle
	if (index < 3)
		*phi = lambda[index] - 20.0*lambda[0] * lambda[1] * lambda[2];
	else
		*phi = 60.0*lambda[0] * lambda[1] * lambda[2];
		*/
}

/**
* \fn void mini_basis1(double *lambda, double **gradLambda, int index, double phi[2])
* \brief the first order derivative of MINI element basis function: (\partial_{x}phi, \partial_{y}phi)
* \param *lambda pointer to the area coordiante
* \param gradLambda gradient of the barycentric coordinate
* \param index the indicator of the basis function
* \param phi[2] the first order derivative of Morley element basis function: (\partial_{x}phi, \partial_{y}phi)
* \return void
*/
void mini_basis1(double *lambda, double **gradLambda, int index, double phi[2])
{
	phi[0] = 0;
	phi[1] = 0;
	if (index >= 4 || index<0)
		return;

	// the last dof is the function value at the centroid of the triangle
	axpbyz_array(2, lambda[1]*lambda[2], gradLambda[0], lambda[2]*lambda[0], gradLambda[1], phi);
	axpy_array(2, lambda[0]*lambda[1], gradLambda[2], phi);	
	if (index < 3){
		ax_array(2, -9, phi);
		axpy_array(2, 1.0, gradLambda[index], phi);
	}
	else
		ax_array(2, 27, phi);

	/********* the last dof is the mean value of the function on triangle

	axpbyz_array(2, lambda[1]*lambda[2], gradLambda[0], lambda[2]*lambda[0], gradLambda[1], phi);
	axpy_array(2, lambda[0]*lambda[1], gradLambda[2], phi);	
	if (index < 3){
		ax_array(2, -20, phi);
		axpy_array(2, 1.0, gradLambda[index], phi);
	}
	else
		ax_array(2, 60, phi);*/
}


// /** 
//  * \fn void H1S_anhuang_basis0(double *lambda, int index, double phi[3])
//  * \brief basis function of An-Huang element
//  * \param *lambda pointer to the area coordiante
//  * \param index the indicator of the basis function
//  * \param *phi basis function
//  * \return void
//  */
// void H1S_anhuang_basis0(double *lambda, int index, double phi[3])
// {
// 	phi[0]=0;
// 	phi[1]=0;
// 	phi[2]=0;
// 	if(index>= 12 || index<0)
// 		return;

// 	double val;
// 	int i;

// 	if(index<9){
// 		i = index%3;
// 		val = lambda[i];
// 		phi[index/3] = val;
// 	}
// 	else{
// 		i=index-9;
// 		val = lambda[0]*lambda[1]*lambda[2];
// 		phi[i] = val;
// 	} 
// }

// /** 
//  * \fn void H1S_anhuang_basis1(double *lambda, double **gradLambda, int index, double(*phi)[2])
//  * \brief the first order derivative of basis function of An-Huang element
//  * \param *lambda pointer to the area coordiante
//  * \param gradLambda gradient of the barycentric coordinate
//  * \param index the indicator of the basis function
//  * \param double(*phi)[2] the first order derivative of basis function of An-Huhang element
//  * \return void
//  */
// void H1S_anhuang_basis10(double *lambda, double **gradLambda, int index, double(*phi)[2])
// {
// 	phi[0][0] = 0; phi[0][1] = 0;
// 	phi[1][0] = 0; phi[1][1] = 0;
// 	phi[2][0] = 0; phi[2][1] = 0;
// 	if(index>= 12 || index<0)
// 		return;

// 	double val[2];
// 	int i, j;

// 	if(index<9){
// 		i = index%3;
// 		phi[index/3][0] = gradLambda[i][0];
// 		phi[index/3][1] = gradLambda[i][1];
// 	}
// 	else{
// 		i=index-9;
// 		val = lambda[0]*lambda[1]*lambda[2];
// 		axpbyz_array(2, lambda[1]*lambda[2], gradLambda[0], lambda[2]*lambda[0], gradLambda[1], val);
// 		axpy_array(2, lambda[0]*lambda[1], gradLambda[2], val);
// 		phi[i][0] = val[0];
// 		phi[i][1] = val[1];
// 	} 
// }

/**
 * \fn double area(double **tri)
 * \brief get area for triangle p1(x1,y1),p2(x2,y2),p3(x3,y3)
 * area = det([1 x1 y1;
               1 x2 y2;
			   1 x3 y3])
 * \param (*tru)[3] the axis value of the three vertices
 * \return area of the trianle
 */
double area(double **tri)
{
	double val;
	val = ((tri[1][0]-tri[0][0])*(tri[2][1]-tri[0][1])-(tri[1][1]-tri[0][1])*(tri[2][0]-tri[0][0]))/2;
	
	if(val<0) val *= -1;

	return val;
}

/**
 * \fn double area0(double *v0, double *v1, double *v2)
 * \brief get area for triangle p0(x0,y0),p1(x1,y1),p2(x2,y2)
 * \param v0 the coordinate of vertex p0
 * \param v1 the coordinate of vertex p1
 * \param v2 the coordinate of vertex p2
 * \return area of the trianle
 */
double area0(double *v0, double *v1, double *v2)
{
	double val;
	val = ((v1[0]-v0[0])*(v2[1]-v0[1])-(v1[1]-v0[1])*(v2[0]-v0[0]))/2;
	
	if(val<0) val *= -1;

	return val;
}
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
#include "matvec.h"

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
 * \fn void lagrange3d_basis(double *lambda, int index, int dop, double *phi)
 * \brief basis function of Lagrange element
 * \param *lambda pointer to the volume coordiante
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param *phi basis function
 * \return void
 */
void lagrange3d_basis(double *lambda, int index, int dop, double *phi)
{
	int fi, ei[2], ii, vi[3], fl, l, i;
	int dofs = (dop+1)*(dop+2)*(dop+3)/6; // degrees of freedom
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
		if(index<4)
		{
			*phi = lambda[index]*(2*lambda[index]-1);
		}
		else
		{
			// edge: 01, 12, 23, 30, 02, 13 
		//	ei[0] = (index-4)%4;
		//	ei[1] = (ei[0]+index/4)%4;
			edge2vv3d(index-4, ei);
			*phi = 4.0*lambda[ei[0]]*lambda[ei[1]];
		}
	} // dop=2
	
	else if(dop==3)
	{
		if(index<4)
		{
			*phi = lambda[index]*(3*lambda[index]-1)*(3*lambda[index]-2)/2.0;
		}
		else if(index < 4+6*(dop-1)) 
		{
			// edge: 01, 12, 23, 30, 02, 13 
//			ei[0] = (index-4)/(dop-1)%4;
//			ei[1] = (ei[0] + (index-4)/(dop-1)/4 +1)%4;
			edge2vv3d((index-4)/(dop-1), ei);
			ii = (index-4)%(dop-1); // ii=0 or 1 if dop=3
			*phi = lambda[ei[0]]*lambda[ei[1]]*(3*lambda[ei[ii]]-1)*9.0/2.0;
		}
		else
		{
			fi=index-16;
			*phi = 27.0*lambda[(fi+1)%4]*lambda[(fi+2)%4]*lambda[(fi+3)%4];
		}
	} // dop=3

	else if(dop==4)
	{
		if(index<4)
		{
			*phi = lambda[index]*(4*lambda[index]-1)*(4*lambda[index]-2)*(4*lambda[index]-3)/6.0;
		}
		else if(index < 4+6*(dop-1)) 
		{
			// edge: 01, 12, 23, 30, 02, 13 
			edge2vv3d((index-4)/(dop-1), ei);
			ii = (index-4)%(dop-1); // ii=0, 1 or 2 if dop=4
			switch(ii)
			{
			case 0:
				*phi = lambda[ei[0]]*lambda[ei[1]]*(4*lambda[ei[0]]-1)*(4*lambda[ei[0]]-2)*8.0/3.0;
				break; 
			case 1:
				*phi = lambda[ei[0]]*lambda[ei[1]]*(4*lambda[ei[0]]-1)*(4*lambda[ei[1]]-1)*4.0;
				break;
			case 2:
				*phi = lambda[ei[0]]*lambda[ei[1]]*(4*lambda[ei[1]]-1)*(4*lambda[ei[1]]-2)*8.0/3.0;
				break;
			default:
				*phi = 0;
			}
		}
		else if(index < 2*(dop*dop+1)) 
		{
			fi=(index-(6*dop-2))/((dop-1)*(dop-2)/2);
			face2vertices3d(fi, vi);
			ii=(index-(6*dop-2))%((dop-1)*(dop-2)/2);
			*phi = 32.0*lambda[vi[0]]*lambda[vi[1]]*lambda[vi[2]]*(4*lambda[vi[ii]]-1);
		}
		else
		{
			*phi = 256.0*lambda[0]*lambda[1]*lambda[2]*lambda[3];
		}
	} // dop=4

	/********* Bernstein basis ***********/
	else if(dop>=5)  // it works for any dop >= 1
	{
		if(index<4)
		{
			*phi =  pow(lambda[index], dop);
		}
		else if(index < 4+6*(dop-1)) 
		{
			// edge: 01, 12, 23, 30, 02, 13 
			edge2vv3d((index-4)/(dop-1), ei);
			ii = (index-4)%(dop-1); // ii=0, 1 or 2 if dop=4
			*phi =  4.0 *pow(lambda[ei[0]], ii+1)*pow(lambda[ei[1]], dop-1-ii);
		}
		else if(index < 2*(dop*dop+1)) 
		{
			fi=(index-(6*dop-2))/((dop-1)*(dop-2)/2);
			face2vertices3d(fi, vi);
			ii=(index-(6*dop-2))%((dop-1)*(dop-2)/2);
			for(l=0;l<dop-2;l++)
			{
				if(ii<(l+1)*(l+2)/2) 
					break;
			}
			i=ii-(l+1)*l/2;
			*phi =  27.0*pow(lambda[vi[0]], dop-2-l)*pow(lambda[vi[1]], l-i+1)*pow(lambda[vi[2]], i+1);
		}
		else
		{	
			ii=index-2*(dop*dop+1);
			for(fl=0;fl<dop-3;fl++)
			{
				if(ii<(fl+1)*(fl+2)*(fl+3)/6) 
					break;
			}
			fi=ii-(fl+1)*(fl+2)*fl/6;
			for(l=0;l<fl+1;l++)
			{
				if(fi<(l+1)*(l+2)/2) 
					break;
			}
			i=fi-(l+1)*l/2;
			*phi = 256.0*pow(lambda[0],dop-3-fl)*pow(lambda[1],fl-l+1)*pow(lambda[2],l-i+1)*pow(lambda[3],i+1);
		}
	}
	
}

/** 
 * \fn void lagrange3d_basis1(double *lambda, double **gradLambda, int index, int dop, double phi[3])
 * \brief the first order derivative of Lagrange element basis function: (\partial_{x}phi, \partial_{y}phi, \partial_{z}phi)
 * \param *lambda pointer to the volume coordiante
 * \param gradLambda gradient of the barycentric coordinate
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param phi[3] the first order derivative of lagrangian element basis function: (\partial_{x}phi, \partial_{y}phi, \partial_{z}phi)
 * \return void
 */
void lagrange3d_basis1(double *lambda, double **gradLambda, int index, int dop, double phi[3])
{
	int fl, l, i;
	int in, ei[2], fi, ii, i1, i2, i3, vi[3];
	int dofs = (dop+1)*(dop+2)*(dop+3)/6; // degrees of freedom

	if(index>= dofs || index<0)
	{
		phi[0]=0; phi[1]=0; phi[2]=0;
		return;
	}

	// gradLambda[0][0]=(eta[1]*zeta[0]-eta[0]*zeta[1])/(6*v); gradLambda[0][1]=(zeta[1]*xi[0]-zeta[0]*xi[1])/(6*v); gradLambda[0][2]=(xi[1]*eta[0]-xi[0]*eta[1])/(6*v);
	// gradLambda[1][0]=-(eta[2]*zeta[1]-eta[1]*zeta[2])/(6*v); gradLambda[1][1]=-(zeta[2]*xi[1]-zeta[1]*xi[2])/(6*v); gradLambda[1][2]=-(xi[2]*eta[1]-xi[1]*eta[2])/(6*v);
	// gradLambda[2][0]=(eta[3]*zeta[2]-eta[2]*zeta[3])/(6*v); gradLambda[2][1]=(zeta[3]*xi[2]-zeta[2]*xi[3])/(6*v); gradLambda[2][2]=(xi[3]*eta[2]-xi[2]*eta[3])/(6*v);
	// gradLambda[3][0]=-(eta[0]*zeta[3]-eta[3]*zeta[0])/(6*v); gradLambda[3][1]=-(zeta[0]*xi[3]-zeta[3]*xi[0])/(6*v); gradLambda[3][2]=-(xi[0]*eta[3]-xi[3]*eta[0])/(6*v);

	if(dop==0)
	{
		phi[0]=0; phi[1]=0; phi[2]=0;
	} // dop=0

	else if(dop==1)
	{
		copy_array(3, gradLambda[index], phi);
	} // dop=1

	else if(dop==2)
	{
		if(index<4)
		{
			axy_array(3, 4*lambda[index]-1, gradLambda[index], phi);
		}
		else
		{
			// edge: 01, 12, 23, 30, 02, 13 
		//	ei[0] = (index-4)%4;
		//	ei[1] = (ei[0]+index/4)%4;
			edge2vv3d(index-4, ei);
			axy_array(3, 4.0*lambda[ei[1]], gradLambda[ei[0]], phi);
			axpy_array(3, 4.0*lambda[ei[0]], gradLambda[ei[1]], phi);
		}
	} // dop=2
	
	else if(dop==3)
	{
		if(index<4)
		{
			axy_array(3, ((3*lambda[index]-1)*(3*lambda[index]-2) + 3*lambda[index]*(3*lambda[index]-2) + 3*lambda[index]*(3*lambda[index]-1))/2.0, gradLambda[index], phi);
		}
		else if(index < 4+6*(dop-1)) 
		{
			// edge: 01, 12, 23, 30, 02, 13 
//			ei[0] = (index-4)/(dop-1)%4;
//			ei[1] = (ei[0] + (index-4)/(dop-1)/4 +1)%4;
			edge2vv3d((index-4)/(dop-1), ei);
			ii = (index-4)%(dop-1); // ii=0 or 1 if dop=3
			axy_array(3, lambda[ei[1]]*(3*lambda[ei[ii]]-1)*9.0/2.0, gradLambda[ei[0]], phi);
			axpy_array(3, lambda[ei[0]]*(3*lambda[ei[ii]]-1)*9.0/2.0, gradLambda[ei[1]], phi);
			axpy_array(3, lambda[ei[0]]*lambda[ei[1]]*3*9.0/2.0, gradLambda[ei[ii]], phi);
		}
		else
		{
			fi=index-16;
			axy_array(3, 27.0*lambda[(fi+2)%4]*lambda[(fi+3)%4], gradLambda[(fi+1)%4], phi);
			axpy_array(3, 27.0*lambda[(fi+1)%4]*lambda[(fi+3)%4], gradLambda[(fi+2)%4], phi);
			axpy_array(3, 27.0*lambda[(fi+1)%4]*lambda[(fi+2)%4], gradLambda[(fi+3)%4], phi);

		}
	} // dop=3

	else if(dop==4)
	{
		if(index<4)
		{
			axy_array(3, (8*lambda[index]-1)*(4*lambda[index]-2)*(4*lambda[index]-3)/6.0, gradLambda[index], phi);
			axpy_array(3, lambda[index]*(4*lambda[index]-1)*(32*lambda[index]-20)/6.0, gradLambda[index], phi);
		}
		else if(index < 4+6*(dop-1)) 
		{
			// edge: 01, 12, 23, 30, 02, 13 
			edge2vv3d((index-4)/(dop-1), ei);
			ii = (index-4)%(dop-1); // ii=0, 1 or 2 if dop=4
			switch(ii)
			{
			case 0:
				axy_array(3, lambda[ei[1]]*(4*lambda[ei[0]]-1)*(4*lambda[ei[0]]-2)*8.0/3.0, gradLambda[ei[0]], phi);
				axpy_array(3, lambda[ei[0]]*(4*lambda[ei[0]]-1)*(4*lambda[ei[0]]-2)*8.0/3.0, gradLambda[ei[1]], phi);
				axpy_array(3, lambda[ei[0]]*lambda[ei[1]]*(32*lambda[ei[0]]-12)*8.0/3.0, gradLambda[ei[0]], phi);
				break; 
			case 1:
				axy_array(3, (8*lambda[ei[0]]-1)*lambda[ei[1]]*(4*lambda[ei[1]]-1)*4.0, gradLambda[ei[0]], phi);
				axpy_array(3, lambda[ei[0]]*(4*lambda[ei[0]]-1)*(8*lambda[ei[1]]-1)*4.0, gradLambda[ei[1]], phi);
				break;
			case 2:
				axy_array(3, lambda[ei[1]]*(4*lambda[ei[1]]-1)*(4*lambda[ei[1]]-2)*8.0/3.0, gradLambda[ei[0]], phi);
				axpy_array(3, lambda[ei[0]]*(4*lambda[ei[1]]-1)*(4*lambda[ei[1]]-2)*8.0/3.0, gradLambda[ei[1]], phi);
				axpy_array(3, lambda[ei[0]]*lambda[ei[1]]*(32*lambda[ei[1]]-12)*8.0/3.0, gradLambda[ei[1]], phi);
				break;
			default:
				phi[0] = 0; phi[1] = 0; phi[2] = 0;
			}
		}
		else if(index < 2*(dop*dop+1)) 
		{
			fi=(index-(6*dop-2))/((dop-1)*(dop-2)/2);
			face2vertices3d(fi, vi);
			ii=(index-(6*dop-2))%((dop-1)*(dop-2)/2);
			axy_array(3, 32.0*lambda[vi[1]]*lambda[vi[2]]*(4*lambda[vi[ii]]-1), gradLambda[vi[0]], phi);
			axpy_array(3, 32.0*lambda[vi[0]]*lambda[vi[2]]*(4*lambda[vi[ii]]-1), gradLambda[vi[1]], phi);
			axpy_array(3, 32.0*lambda[vi[0]]*lambda[vi[1]]*(4*lambda[vi[ii]]-1), gradLambda[vi[2]], phi);
			axpy_array(3, 32.0*lambda[vi[0]]*lambda[vi[1]]*lambda[vi[2]]*4, gradLambda[vi[ii]], phi);
		}
		else
		{
			axy_array(3, 256.0*lambda[1]*lambda[2]*lambda[3], gradLambda[0], phi);
			axpy_array(3, 256.0*lambda[0]*lambda[2]*lambda[3], gradLambda[1], phi);
			axpy_array(3, 256.0*lambda[0]*lambda[1]*lambda[3], gradLambda[2], phi);
			axpy_array(3, 256.0*lambda[0]*lambda[1]*lambda[2], gradLambda[3], phi);
		}
	} // dop=4

		/********* Bernstein basis ***********/
	else if(dop>=5) // it works for any dop >= 1
	{
		if(index<4)
		{
			// *phi =  pow(lambda[index], dop);
			axy_array(3, dop*pow(lambda[index], dop-1), gradLambda[index], phi);
		}
		else if(index < 4+6*(dop-1)) 
		{
			// edge: 01, 12, 23, 30, 02, 13 
			edge2vv3d((index-4)/(dop-1), ei);
			ii = (index-4)%(dop-1); // ii=0, 1 or 2 if dop=4
			// *phi =  4.0 *pow(lambda[ei[0]], ii+1)*pow(lambda[ei[1]], dop-1-ii);
			axy_array(3, 4.0 * (ii+1) * pow(lambda[ei[0]], ii)*pow(lambda[ei[1]], dop-1-ii), gradLambda[ei[0]], phi);
			axpy_array(3, 4.0 *(dop-1-ii)*pow(lambda[ei[0]], ii+1)*pow(lambda[ei[1]], dop-2-ii), gradLambda[ei[1]], phi);
		}
		else if(index < 2*(dop*dop+1)) 
		{
			fi=(index-(6*dop-2))/((dop-1)*(dop-2)/2);
			face2vertices3d(fi, vi);
			ii=(index-(6*dop-2))%((dop-1)*(dop-2)/2);
			for(l=0;l<dop-2;l++)
			{
				if(ii<(l+1)*(l+2)/2) 
					break;
			}
			i=ii-(l+1)*l/2;
			// *phi =  27.0*pow(lambda[vi[0]], dop-2-l)*pow(lambda[vi[1]], l-i+1)*pow(lambda[vi[2]], i+1);
			axy_array(3, 27.0*(dop-2-l)*pow(lambda[vi[0]], dop-3-l)*pow(lambda[vi[1]], l-i+1)*pow(lambda[vi[2]], i+1), gradLambda[vi[0]], phi);
			axpy_array(3, 27.0*(l-i+1)*pow(lambda[vi[0]], dop-2-l)*pow(lambda[vi[1]], l-i)*pow(lambda[vi[2]], i+1), gradLambda[vi[1]], phi);
			axpy_array(3, 27.0*(i+1)*pow(lambda[vi[0]], dop-2-l)*pow(lambda[vi[1]], l-i+1)*pow(lambda[vi[2]], i), gradLambda[vi[2]], phi);
		}
		else
		{	
			ii=index-2*(dop*dop+1);
			for(fl=0;fl<dop-3;fl++)
			{
				if(ii<(fl+1)*(fl+2)*(fl+3)/6) 
					break;
			}
			fi=ii-(fl+1)*(fl+2)*fl/6;
			for(l=0;l<fl+1;l++)
			{
				if(fi<(l+1)*(l+2)/2) 
					break;
			}
			i=fi-(l+1)*l/2;
			// *phi = 256.0*pow(lambda[0],dop-3-fl)*pow(lambda[1],fl-l+1)*pow(lambda[2],l-i+1)*pow(lambda[3],i+1);
			axy_array(3, 256.0*(dop-3-fl)*pow(lambda[0],dop-4-fl)*pow(lambda[1],fl-l+1)*pow(lambda[2],l-i+1)*pow(lambda[3],i+1), gradLambda[0], phi);
			axpy_array(3, 256.0*(fl-l+1)*pow(lambda[0],dop-3-fl)*pow(lambda[1],fl-l)*pow(lambda[2],l-i+1)*pow(lambda[3],i+1), gradLambda[1], phi);
			axpy_array(3, 256.0*(l-i+1)*pow(lambda[0],dop-3-fl)*pow(lambda[1],fl-l+1)*pow(lambda[2],l-i)*pow(lambda[3],i+1), gradLambda[2], phi);
			axpy_array(3, 256.0*(i+1)*pow(lambda[0],dop-3-fl)*pow(lambda[1],fl-l+1)*pow(lambda[2],l-i+1)*pow(lambda[3],i), gradLambda[3], phi);
		}
	}
}

/**
* \fn void ncp13d_basis(double *lambda, int index, double *phi)
* \brief basis function of nonconforming P1 element
* \param *lambda pointer to the barycentric coordinates
* \param index the indicator of the basis function
* \param *phi basis function
* \return void
*/
void ncp13d_basis(double *lambda, int index, double *phi)
{
	if (index >= 4 || index<0)
	{
		*phi = 0;
		return;
	}

	*phi = 1 - 3 * lambda[index];
}

/**
* \fn void ncp13d_basis1(double **grd_lambda, int index, double *phi)
* \brief basis function of nonconforming P1 element in three dimensions
* \param *lambda pointer to the barycentric coordinates
* \param index the indicator of the basis function
* \param *phi basis function
* \return void
*/
void ncp13d_basis1(double **grd_lambda, int index, double *phi)
{
	if (index >= 4 || index<0)
	{
		phi[0] = 0; phi[1] = 0; phi[2] = 0;
		return;
	}

	axy_array(3, -3, grd_lambda[index], phi);
}

/**
* \fn void morley3d_basis(double *lambda, double v, double **grd_lambda, double **nv, double **nvf, int index, double *phi)
* \brief basis function of Morley-Wang-Xu element in three dimensions
* \param *lambda pointer to the volume coordiante
* \param **grd_lambda pointer to the gradient of the volume coordiante lambda
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi basis function
* \return void
*/
void morley3d_basis(double *lambda, double v, double **grd_lambda, double **nv, double **nvf, int index, double *phi)
{
	int i, fi, ei[2], ii, vi[3];
	int dofs = 10; // degrees of freedom
	double grdlen2[4], grdip;
	
	if (index >= dofs || index<0)
	{
		*phi = 0;
		return;
	}

	double orient[4];
	for (i = 0; i < 4; i++)
		orient[i] = nv[i][0] * nvf[i][0] + nv[i][1] * nvf[i][1] + nv[i][2] * nvf[i][2];

	for (i = 0; i < 4; i++)
		grdlen2[i] = grd_lambda[i][0] * grd_lambda[i][0] + grd_lambda[i][1] * grd_lambda[i][1] + grd_lambda[i][2] * grd_lambda[i][2];

	if (index<6)
	{
		// edge: 01, 12, 23, 30, 02, 13 
		//	ei[0] = (index-4)%4;
		//	ei[1] = (ei[0]+index/4)%4;
		edge2vv3d(index, ei);
		ii = (ei[0] + 1) % 4;
		if (ii == ei[1])
			ii = (ei[0] + 2) % 4;
		ei[0] = 6 - ei[0] - ei[1] - ii;
		ei[1] = ii;
		grdip = grd_lambda[ei[0]][0] * grd_lambda[ei[1]][0] + grd_lambda[ei[0]][1] * grd_lambda[ei[1]][1] + grd_lambda[ei[0]][2] * grd_lambda[ei[1]][2];
		*phi = 1 - 2.0*(lambda[ei[0]] + lambda[ei[1]]) + 6.0*lambda[ei[0]] * lambda[ei[1]] - grdip*(lambda[ei[0]] * (3 * lambda[ei[0]] - 2) / grdlen2[ei[0]] + lambda[ei[1]] * (3 * lambda[ei[1]] - 2) / grdlen2[ei[1]]);
	}
	else
	{
		fi = index - 6;
		*phi = lambda[fi] * (3 * lambda[fi] - 2) / (2 * sqrt(grdlen2[fi]))*orient[fi];
	}
}

/**
* \fn void morley3d_basis1(double *lambda, double v, double **grd_lambda, double **nv, double **nvf, int index, double phi[3])
* \brief the first order derivative of Morley-Wang-Xu element basis function: (\partial_{x}phi, \partial_{y}phi, \partial_{z}phi)
* \param *lambda pointer to the volume coordiante
* \param v the volume of the tetrahedron
* \param **grd_lambda pointer to the gradient of the volume coordiante lambda
* \param **nv the unit normal vectors of the fourth faces
* \param **nve the unit normal vectors of the fourth faces (fixed for each face)
* \param index the indicator of the basis function
* \param phi[3] the first order derivative of Morley-Wang-Xu element basis function: (\partial_{x}phi, \partial_{y}phi, \partial_{z}phi)
* \return void
*/
void morley3d_basis1(double *lambda, double v, double **grd_lambda, double **nv, double **nvf, int index, double phi[3])
{
	int i, fi, ei[2], ii, vi[3];
	int dofs = 10; // degrees of freedom
	double grdlen2[4], grdip;

	if (index >= dofs || index<0)
	{
		phi[0] = 0; phi[1] = 0; phi[2] = 0;
		return;
	}

	double orient[4];
	for (i = 0; i < 4; i++)
		orient[i] = nv[i][0] * nvf[i][0] + nv[i][1] * nvf[i][1] + nv[i][2] * nvf[i][2];

	for (i = 0; i < 4; i++)
		grdlen2[i] = grd_lambda[i][0] * grd_lambda[i][0] + grd_lambda[i][1] * grd_lambda[i][1] + grd_lambda[i][2] * grd_lambda[i][2];

	if (index<6)
	{
		// edge: 01, 12, 23, 30, 02, 13 
		//	ei[0] = (index-4)%4;
		//	ei[1] = (ei[0]+index/4)%4;
		edge2vv3d(index, ei);
		ii = (ei[0] + 1) % 4;
		if (ii == ei[1])
			ii = (ei[0] + 2) % 4;
		ei[0] = 6 - ei[0] - ei[1] - ii;
		ei[1] = ii;
		grdip = grd_lambda[ei[0]][0] * grd_lambda[ei[1]][0] + grd_lambda[ei[0]][1] * grd_lambda[ei[1]][1] + grd_lambda[ei[0]][2] * grd_lambda[ei[1]][2];
		phi[0] = -2.0*(grd_lambda[ei[0]][0] + grd_lambda[ei[1]][0]) + 6.0*(grd_lambda[ei[0]][0] * lambda[ei[1]] + lambda[ei[0]] * grd_lambda[ei[1]][0]) - grdip*(grd_lambda[ei[0]][0] * (6 * lambda[ei[0]] - 2) / grdlen2[ei[0]] + grd_lambda[ei[1]][0] * (6 * lambda[ei[1]] - 2) / grdlen2[ei[1]]);
		phi[1] = -2.0*(grd_lambda[ei[0]][1] + grd_lambda[ei[1]][1]) + 6.0*(grd_lambda[ei[0]][1] * lambda[ei[1]] + lambda[ei[0]] * grd_lambda[ei[1]][1]) - grdip*(grd_lambda[ei[0]][1] * (6 * lambda[ei[0]] - 2) / grdlen2[ei[0]] + grd_lambda[ei[1]][1] * (6 * lambda[ei[1]] - 2) / grdlen2[ei[1]]);
		phi[2] = -2.0*(grd_lambda[ei[0]][2] + grd_lambda[ei[1]][2]) + 6.0*(grd_lambda[ei[0]][2] * lambda[ei[1]] + lambda[ei[0]] * grd_lambda[ei[1]][2]) - grdip*(grd_lambda[ei[0]][2] * (6 * lambda[ei[0]] - 2) / grdlen2[ei[0]] + grd_lambda[ei[1]][2] * (6 * lambda[ei[1]] - 2) / grdlen2[ei[1]]);
	}
	else
	{
		fi = index - 6;
		phi[0] = grd_lambda[fi][0] * (6 * lambda[fi] - 2) / (2 * sqrt(grdlen2[fi]))*orient[fi];
		phi[1] = grd_lambda[fi][1] * (6 * lambda[fi] - 2) / (2 * sqrt(grdlen2[fi]))*orient[fi];
		phi[2] = grd_lambda[fi][2] * (6 * lambda[fi] - 2) / (2 * sqrt(grdlen2[fi]))*orient[fi];
	}
}

/**
* \fn void morley3d_basis2(double *lambda, double v, double **grd_lambda, double **nv, double **nvf, int index, double phi[6])
* \brief the second order derivative of Morley-Wang-Xu element basis function: (\partial_{xx}phi, \partial_{yy}phi, \partial_{zz}phi, \partial_{yz}phi, \partial_{zx}phi, \partial_{xy}phi)
* \param *lambda pointer to the volume coordiante
* \param v the volume of the tetrahedron
* \param **grd_lambda pointer to the gradient of the volume coordiante lambda
* \param **nv the unit normal vectors of the fourth faces
* \param **nve the unit normal vectors of the fourth faces (fixed for each face)
* \param index the indicator of the basis function
* \param phi[6] the second order derivative of Morley-Wang-Xu element basis function: (\partial_{xx}phi, \partial_{yy}phi, \partial_{zz}phi, \partial_{yz}phi, \partial_{zx}phi, \partial_{xy}phi)
* \return void
*/
void morley3d_basis2(double *lambda, double v, double **grd_lambda, double **nv, double **nvf, int index, double phi[6])
{
	int i, fi, ei[2], ii, vi[3];
	int dofs = 10; // degrees of freedom
	double grdlen2[4], grdip;

	if (index >= dofs || index<0)
	{
		phi[0] = 0; phi[1] = 0; phi[2] = 0;
		phi[3] = 0; phi[4] = 0; phi[5] = 0;
		return;
	}

	double orient[4];
	for (i = 0; i < 4; i++)
		orient[i] = nv[i][0] * nvf[i][0] + nv[i][1] * nvf[i][1] + nv[i][2] * nvf[i][2];

	for (i = 0; i < 4; i++)
		grdlen2[i] = grd_lambda[i][0] * grd_lambda[i][0] + grd_lambda[i][1] * grd_lambda[i][1] + grd_lambda[i][2] * grd_lambda[i][2];

	if (index<6)
	{
		// edge: 01, 12, 23, 30, 02, 13 
		//	ei[0] = (index-4)%4;
		//	ei[1] = (ei[0]+index/4)%4;
		edge2vv3d(index, ei);
		ii = (ei[0] + 1) % 4;
		if (ii == ei[1])
			ii = (ei[0] + 2) % 4;
		ei[0] = 6 - ei[0] - ei[1] - ii;
		ei[1] = ii;
		grdip = grd_lambda[ei[0]][0] * grd_lambda[ei[1]][0] + grd_lambda[ei[0]][1] * grd_lambda[ei[1]][1] + grd_lambda[ei[0]][2] * grd_lambda[ei[1]][2];
		phi[0] = 12.0*grd_lambda[ei[0]][0] * grd_lambda[ei[1]][0] - 6 * grdip*(grd_lambda[ei[0]][0] * grd_lambda[ei[0]][0] / grdlen2[ei[0]] + grd_lambda[ei[1]][0] * grd_lambda[ei[1]][0] / grdlen2[ei[1]]);
		phi[1] = 12.0*grd_lambda[ei[0]][1] * grd_lambda[ei[1]][1] - 6 * grdip*(grd_lambda[ei[0]][1] * grd_lambda[ei[0]][1] / grdlen2[ei[0]] + grd_lambda[ei[1]][1] * grd_lambda[ei[1]][1] / grdlen2[ei[1]]);
		phi[2] = 12.0*grd_lambda[ei[0]][2] * grd_lambda[ei[1]][2] - 6 * grdip*(grd_lambda[ei[0]][2] * grd_lambda[ei[0]][2] / grdlen2[ei[0]] + grd_lambda[ei[1]][2] * grd_lambda[ei[1]][2] / grdlen2[ei[1]]);
		phi[3] = 6.0*(grd_lambda[ei[0]][1] * grd_lambda[ei[1]][2] + grd_lambda[ei[0]][2] * grd_lambda[ei[1]][1]) - 6 * grdip*(grd_lambda[ei[0]][1] * grd_lambda[ei[0]][2] / grdlen2[ei[0]] + grd_lambda[ei[1]][1] * grd_lambda[ei[1]][2] / grdlen2[ei[1]]);
		phi[4] = 6.0*(grd_lambda[ei[0]][2] * grd_lambda[ei[1]][0] + grd_lambda[ei[0]][0] * grd_lambda[ei[1]][2]) - 6 * grdip*(grd_lambda[ei[0]][2] * grd_lambda[ei[0]][0] / grdlen2[ei[0]] + grd_lambda[ei[1]][2] * grd_lambda[ei[1]][0] / grdlen2[ei[1]]);
		phi[5] = 6.0*(grd_lambda[ei[0]][0] * grd_lambda[ei[1]][1] + grd_lambda[ei[0]][1] * grd_lambda[ei[1]][0]) - 6 * grdip*(grd_lambda[ei[0]][0] * grd_lambda[ei[0]][1] / grdlen2[ei[0]] + grd_lambda[ei[1]][0] * grd_lambda[ei[1]][1] / grdlen2[ei[1]]);
	}
	else
	{
		fi = index - 6;
		phi[0] = 3 * grd_lambda[fi][0] * grd_lambda[fi][0] / sqrt(grdlen2[fi])*orient[fi];
		phi[1] = 3 * grd_lambda[fi][1] * grd_lambda[fi][1] / sqrt(grdlen2[fi])*orient[fi];
		phi[2] = 3 * grd_lambda[fi][2] * grd_lambda[fi][2] / sqrt(grdlen2[fi])*orient[fi];
		phi[3] = 3 * grd_lambda[fi][1] * grd_lambda[fi][2] / sqrt(grdlen2[fi])*orient[fi];
		phi[4] = 3 * grd_lambda[fi][2] * grd_lambda[fi][0] / sqrt(grdlen2[fi])*orient[fi];
		phi[5] = 3 * grd_lambda[fi][0] * grd_lambda[fi][1] / sqrt(grdlen2[fi])*orient[fi];
	}
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
 * \brief basis function of Hu-Zhang element in two dimesions
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
 * \fn void huzhang_basisDIV(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, int index, int dop, double phi[2])
 * \brief divergence of basis function of Hu-Zhang element in two dimesions
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
 * \fn void huzhang3d_basis(double *lambda, double (*fnn)[6], double (*fnt1)[6], double (*fnt2)[6], double (*ft1t1)[6], double (*ft2t2)[6], double (*ft1t2)[6], double (*ett)[6], double (*etn1)[6], double (*etn2)[6], double (*en1n1)[6], double (*en2n2)[6], double (*en1n2)[6], int index, int dop, double phi[6])
 * \brief basis function of Hu-Zhang element in three dimesions
 * \param *lambda pointer to the area coordiante
 * \param **fnv the unit normal vectors of the four faces
 * \param **ft1v the unit tangential vectors of the four faces
 * \param **ft2v the unit tangential vectors of the four faces
 * \param **etv the unit tangential vectors of the six edges
 * \param **en1v the unit normal vectors of the six edges
 * \param **en2v the unit normal vectors of the six edges
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param *phi basis function
 * \return void
 */
void huzhang3d_basis(double *lambda, double (*fnn)[6], double (*fnt1)[6], double (*fnt2)[6], double (*ft1t1)[6], double (*ft2t2)[6], double (*ft1t2)[6], double (*ett)[6], double (*etn1)[6], double (*etn2)[6], double (*en1n1)[6], double (*en2n2)[6], double (*en1n2)[6], int index, int dop, double phi[6])
{
	int dofs = (dop+1)*(dop+2)*(dop+3)/6; // degrees of freedom

	int i,j;
	for(i=0;i<6;i++)
		phi[i]=0;

	if(index>= dofs*6 || index<0)
		return;

	double val;
	
	if(dop==0)
	{
		phi[index] = 1;
	}
	
	else // dop>=1
	{
		if(index<24)
		{
			lagrange3d_basis(lambda, index%4, dop, &val);
			phi[index/4] = val;
		}
		else if(index<24+(dop-1)*6)
		{
			i=index-24;
			lagrange3d_basis(lambda, 4+i, dop, &val);
			for(j=0;j<6;j++)
				phi[j]=val*en1n1[i/(dop-1)][j];
		}
		else if(index<24+(dop-1)*12)
		{
			i=index-24-(dop-1)*6;
			lagrange3d_basis(lambda, 4+i, dop, &val);
			for(j=0;j<6;j++)
				phi[j]=val*en2n2[i/(dop-1)][j];
		}
		else if(index<24+(dop-1)*18)
		{
			i=index-24-(dop-1)*12;
			lagrange3d_basis(lambda, 4+i, dop, &val);
			for(j=0;j<6;j++)
				phi[j]=val*en1n2[i/(dop-1)][j];
		}
		else if(index<24+(dop-1)*24)
		{
			i=index-24-(dop-1)*18;
			lagrange3d_basis(lambda, 4+i, dop, &val);
			for(j=0;j<6;j++)
				phi[j]=val*etn1[i/(dop-1)][j];
		}
		else if(index<24+(dop-1)*30)
		{
			i=index-24-(dop-1)*24;
			lagrange3d_basis(lambda, 4+i, dop, &val);
			for(j=0;j<6;j++)
				phi[j]=val*etn2[i/(dop-1)][j];
		}
		else if(index<24+(dop-1)*30+(dop-1)*(dop-2)/2*4)
		{
			i=index-24-(dop-1)*30;
			lagrange3d_basis(lambda, 4+(dop-1)*6+i, dop, &val);
			for(j=0;j<6;j++)
				phi[j]=val*fnn[i/((dop-1)*(dop-2)/2)][j];
		}
		else if(index<24+(dop-1)*30+(dop-1)*(dop-2)/2*8)
		{
			i=index-24-(dop-1)*30-(dop-1)*(dop-2)/2*4;
			lagrange3d_basis(lambda, 4+(dop-1)*6+i, dop, &val);
			for(j=0;j<6;j++)
				phi[j]=val*fnt1[i/((dop-1)*(dop-2)/2)][j];
		}
		else if(index<24+(dop-1)*30+(dop-1)*(dop-2)/2*12)
		{
			i=index-24-(dop-1)*30-(dop-1)*(dop-2)/2*8;
			lagrange3d_basis(lambda, 4+(dop-1)*6+i, dop, &val);
			for(j=0;j<6;j++)
				phi[j]=val*fnt2[i/((dop-1)*(dop-2)/2)][j];
		}
		else if(index<24+(dop-1)*36+(dop-1)*(dop-2)/2*12)
		{
			i=index-24-(dop-1)*30-(dop-1)*(dop-2)/2*12;
			lagrange3d_basis(lambda, 4+i, dop, &val);
			for(j=0;j<6;j++)
				phi[j]=val*ett[i/(dop-1)][j];
		}
		else if(index<24+(dop-1)*36+(dop-1)*(dop-2)/2*16)
		{
			i=index-24-(dop-1)*36-(dop-1)*(dop-2)/2*12;
			lagrange3d_basis(lambda, 4+(dop-1)*6+i, dop, &val);
			for(j=0;j<6;j++)
				phi[j]=val*ft1t1[i/((dop-1)*(dop-2)/2)][j];
		}
		else if(index<24+(dop-1)*36+(dop-1)*(dop-2)/2*20)
		{
			i=index-24-(dop-1)*36-(dop-1)*(dop-2)/2*16;
			lagrange3d_basis(lambda, 4+(dop-1)*6+i, dop, &val);
			for(j=0;j<6;j++)
				phi[j]=val*ft2t2[i/((dop-1)*(dop-2)/2)][j];
		}
		else if(index<24+(dop-1)*36+(dop-1)*(dop-2)/2*24)
		{
			i=index-24-(dop-1)*36-(dop-1)*(dop-2)/2*20;
			lagrange3d_basis(lambda, 4+(dop-1)*6+i, dop, &val);
			for(j=0;j<6;j++)
				phi[j]=val*ft1t2[i/((dop-1)*(dop-2)/2)][j];
		}
		else
		{
			// 4+(dop-1)*6+(dop-1)*(dop-2)/2*4 = 2*(dop*dop+1)
			// 24+(dop-1)*36+(dop-1)*(dop-2)/2*24 = 12*(dop*dop+1)
			i=index-12*(dop*dop+1);
			lagrange3d_basis(lambda, 2*(dop*dop+1)+i%(dofs - 2*(dop*dop+1)), dop, &val);
			phi[i/(dofs - 2*(dop*dop+1))]=val;
		}
	} 
}

/** 
 * \fn void huzhang3d_basisDIV(double *lambda, double v, double **grd_lambda, double (*fnn)[6], double (*fnt1)[6], double (*fnt2)[6], double (*ft1t1)[6], double (*ft2t2)[6], double (*ft1t2)[6], double (*ett)[6], double (*etn1)[6], double (*etn2)[6], double (*en1n1)[6], double (*en2n2)[6], double (*en1n2)[6], int index, int dop, double phi[3])
 * \brief divergence of basis function of Hu-Zhang element in three dimensions
 * \param *lambda pointer to the area coordiante
 * \param v the volume of the tetrahedron
 * \param **grd_lambda pointer to the gradient of the volume coordiante lambda
 * \param **fnv the unit normal vectors of the four faces
 * \param **ft1v the unit tangential vectors of the four faces
 * \param **ft2v the unit tangential vectors of the four faces
 * \param **etv the unit tangential vectors of the six edges
 * \param **en1v the unit normal vectors of the six edges
 * \param **en2v the unit normal vectors of the six edges
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param *phi divergence of basis function
 * \return void
 */
void huzhang3d_basisDIV(double *lambda, double v, double **grd_lambda, double (*fnn)[6], double (*fnt1)[6], double (*fnt2)[6], double (*ft1t1)[6], double (*ft2t2)[6], double (*ft1t2)[6], double (*ett)[6], double (*etn1)[6], double (*etn2)[6], double (*en1n1)[6], double (*en2n2)[6], double (*en1n2)[6], int index, int dop, double phi[3])
{
	int i;
	int dofs = (dop+1)*(dop+2)*(dop+3)/6; // degrees of freedom

	phi[0]=0;phi[1]=0;phi[2]=0;
	if(index>= dofs*6 || index<0)
		return;

	double val[3];

	if(dop==0)
	{
		return;
	}
	
	else // dop>=1
	{
		if(index<4)
		{
			lagrange3d_basis1(lambda, grd_lambda, index, dop, val);
			phi[0] = val[0]; phi[1] = 0; phi[2] = 0;
		}
		else if(index<8)
		{
			lagrange3d_basis1(lambda, grd_lambda, index%4, dop, val);
			phi[0] = 0;	phi[1] = val[1]; phi[2] = 0;
		}
		else if(index<12)
		{
			lagrange3d_basis1(lambda, grd_lambda, index%4, dop, val);
			phi[0] = 0;	phi[1] = 0; phi[2] = val[2];
		}
		else if(index<16)
		{
			lagrange3d_basis1(lambda, grd_lambda, index%4, dop, val);
			phi[0] = 0;	phi[1] = val[2]; phi[2] = val[1];
		}
		else if(index<20)
		{
			lagrange3d_basis1(lambda, grd_lambda, index%4, dop, val);
			phi[0] = val[2]; phi[1] = 0; phi[2] = val[0];
		}
		else if(index<24)
		{
			lagrange3d_basis1(lambda, grd_lambda, index%4, dop, val);
			phi[0] = val[1]; phi[1] = val[0]; phi[2] = 0;
		}
		else if(index<24+(dop-1)*6)
		{
			i=index-24;
			lagrange3d_basis1(lambda, grd_lambda, 4+i, dop, val);
			phi[0]=val[0]*en1n1[i/(dop-1)][0]+val[1]*en1n1[i/(dop-1)][5]+val[2]*en1n1[i/(dop-1)][4];
			phi[1]=val[0]*en1n1[i/(dop-1)][5]+val[1]*en1n1[i/(dop-1)][1]+val[2]*en1n1[i/(dop-1)][3];
			phi[2]=val[0]*en1n1[i/(dop-1)][4]+val[1]*en1n1[i/(dop-1)][3]+val[2]*en1n1[i/(dop-1)][2];
		}
		else if(index<24+(dop-1)*12)
		{
			i=index-24-(dop-1)*6;
			lagrange3d_basis1(lambda, grd_lambda, 4+i, dop, val);
			phi[0]=val[0]*en2n2[i/(dop-1)][0]+val[1]*en2n2[i/(dop-1)][5]+val[2]*en2n2[i/(dop-1)][4];
			phi[1]=val[0]*en2n2[i/(dop-1)][5]+val[1]*en2n2[i/(dop-1)][1]+val[2]*en2n2[i/(dop-1)][3];
			phi[2]=val[0]*en2n2[i/(dop-1)][4]+val[1]*en2n2[i/(dop-1)][3]+val[2]*en2n2[i/(dop-1)][2];
		}
		else if(index<24+(dop-1)*18)
		{
			i=index-24-(dop-1)*12;
			lagrange3d_basis1(lambda, grd_lambda, 4+i, dop, val);
			phi[0]=val[0]*en1n2[i/(dop-1)][0]+val[1]*en1n2[i/(dop-1)][5]+val[2]*en1n2[i/(dop-1)][4];
			phi[1]=val[0]*en1n2[i/(dop-1)][5]+val[1]*en1n2[i/(dop-1)][1]+val[2]*en1n2[i/(dop-1)][3];
			phi[2]=val[0]*en1n2[i/(dop-1)][4]+val[1]*en1n2[i/(dop-1)][3]+val[2]*en1n2[i/(dop-1)][2];
		}
		else if(index<24+(dop-1)*24)
		{
			i=index-24-(dop-1)*18;
			lagrange3d_basis1(lambda, grd_lambda, 4+i, dop, val);
			phi[0]=val[0]*etn1[i/(dop-1)][0]+val[1]*etn1[i/(dop-1)][5]+val[2]*etn1[i/(dop-1)][4];
			phi[1]=val[0]*etn1[i/(dop-1)][5]+val[1]*etn1[i/(dop-1)][1]+val[2]*etn1[i/(dop-1)][3];
			phi[2]=val[0]*etn1[i/(dop-1)][4]+val[1]*etn1[i/(dop-1)][3]+val[2]*etn1[i/(dop-1)][2];
		}
		else if(index<24+(dop-1)*30)
		{
			i=index-24-(dop-1)*24;
			lagrange3d_basis1(lambda, grd_lambda, 4+i, dop, val);
			phi[0]=val[0]*etn2[i/(dop-1)][0]+val[1]*etn2[i/(dop-1)][5]+val[2]*etn2[i/(dop-1)][4];
			phi[1]=val[0]*etn2[i/(dop-1)][5]+val[1]*etn2[i/(dop-1)][1]+val[2]*etn2[i/(dop-1)][3];
			phi[2]=val[0]*etn2[i/(dop-1)][4]+val[1]*etn2[i/(dop-1)][3]+val[2]*etn2[i/(dop-1)][2];
		}
		else if(index<24+(dop-1)*30+(dop-1)*(dop-2)/2*4)
		{
			i=index-24-(dop-1)*30;
			lagrange3d_basis1(lambda, grd_lambda, 4+(dop-1)*6+i, dop, val);
			phi[0]=val[0]*fnn[i/((dop-1)*(dop-2)/2)][0]+val[1]*fnn[i/((dop-1)*(dop-2)/2)][5]+val[2]*fnn[i/((dop-1)*(dop-2)/2)][4];
			phi[1]=val[0]*fnn[i/((dop-1)*(dop-2)/2)][5]+val[1]*fnn[i/((dop-1)*(dop-2)/2)][1]+val[2]*fnn[i/((dop-1)*(dop-2)/2)][3];
			phi[2]=val[0]*fnn[i/((dop-1)*(dop-2)/2)][4]+val[1]*fnn[i/((dop-1)*(dop-2)/2)][3]+val[2]*fnn[i/((dop-1)*(dop-2)/2)][2];
		}
		else if(index<24+(dop-1)*30+(dop-1)*(dop-2)/2*8)
		{
			i=index-24-(dop-1)*30-(dop-1)*(dop-2)/2*4;
			lagrange3d_basis1(lambda, grd_lambda, 4+(dop-1)*6+i, dop, val);
			phi[0]=val[0]*fnt1[i/((dop-1)*(dop-2)/2)][0]+val[1]*fnt1[i/((dop-1)*(dop-2)/2)][5]+val[2]*fnt1[i/((dop-1)*(dop-2)/2)][4];
			phi[1]=val[0]*fnt1[i/((dop-1)*(dop-2)/2)][5]+val[1]*fnt1[i/((dop-1)*(dop-2)/2)][1]+val[2]*fnt1[i/((dop-1)*(dop-2)/2)][3];
			phi[2]=val[0]*fnt1[i/((dop-1)*(dop-2)/2)][4]+val[1]*fnt1[i/((dop-1)*(dop-2)/2)][3]+val[2]*fnt1[i/((dop-1)*(dop-2)/2)][2];
		}
		else if(index<24+(dop-1)*30+(dop-1)*(dop-2)/2*12)
		{
			i=index-24-(dop-1)*30-(dop-1)*(dop-2)/2*8;
			lagrange3d_basis1(lambda, grd_lambda, 4+(dop-1)*6+i, dop, val);
			phi[0]=val[0]*fnt2[i/((dop-1)*(dop-2)/2)][0]+val[1]*fnt2[i/((dop-1)*(dop-2)/2)][5]+val[2]*fnt2[i/((dop-1)*(dop-2)/2)][4];
			phi[1]=val[0]*fnt2[i/((dop-1)*(dop-2)/2)][5]+val[1]*fnt2[i/((dop-1)*(dop-2)/2)][1]+val[2]*fnt2[i/((dop-1)*(dop-2)/2)][3];
			phi[2]=val[0]*fnt2[i/((dop-1)*(dop-2)/2)][4]+val[1]*fnt2[i/((dop-1)*(dop-2)/2)][3]+val[2]*fnt2[i/((dop-1)*(dop-2)/2)][2];
		}
		else if(index<24+(dop-1)*36+(dop-1)*(dop-2)/2*12)
		{
			i=index-24-(dop-1)*30-(dop-1)*(dop-2)/2*12;
			lagrange3d_basis1(lambda, grd_lambda, 4+i, dop, val);
			phi[0]=val[0]*ett[i/(dop-1)][0]+val[1]*ett[i/(dop-1)][5]+val[2]*ett[i/(dop-1)][4];
			phi[1]=val[0]*ett[i/(dop-1)][5]+val[1]*ett[i/(dop-1)][1]+val[2]*ett[i/(dop-1)][3];
			phi[2]=val[0]*ett[i/(dop-1)][4]+val[1]*ett[i/(dop-1)][3]+val[2]*ett[i/(dop-1)][2];
		}
		else if(index<24+(dop-1)*36+(dop-1)*(dop-2)/2*16)
		{
			i=index-24-(dop-1)*36-(dop-1)*(dop-2)/2*12;
			lagrange3d_basis1(lambda, grd_lambda, 4+(dop-1)*6+i, dop, val);
			phi[0]=val[0]*ft1t1[i/((dop-1)*(dop-2)/2)][0]+val[1]*ft1t1[i/((dop-1)*(dop-2)/2)][5]+val[2]*ft1t1[i/((dop-1)*(dop-2)/2)][4];
			phi[1]=val[0]*ft1t1[i/((dop-1)*(dop-2)/2)][5]+val[1]*ft1t1[i/((dop-1)*(dop-2)/2)][1]+val[2]*ft1t1[i/((dop-1)*(dop-2)/2)][3];
			phi[2]=val[0]*ft1t1[i/((dop-1)*(dop-2)/2)][4]+val[1]*ft1t1[i/((dop-1)*(dop-2)/2)][3]+val[2]*ft1t1[i/((dop-1)*(dop-2)/2)][2];
		}
		else if(index<24+(dop-1)*36+(dop-1)*(dop-2)/2*20)
		{
			i=index-24-(dop-1)*36-(dop-1)*(dop-2)/2*16;
			lagrange3d_basis1(lambda, grd_lambda, 4+(dop-1)*6+i, dop, val);
			phi[0]=val[0]*ft2t2[i/((dop-1)*(dop-2)/2)][0]+val[1]*ft2t2[i/((dop-1)*(dop-2)/2)][5]+val[2]*ft2t2[i/((dop-1)*(dop-2)/2)][4];
			phi[1]=val[0]*ft2t2[i/((dop-1)*(dop-2)/2)][5]+val[1]*ft2t2[i/((dop-1)*(dop-2)/2)][1]+val[2]*ft2t2[i/((dop-1)*(dop-2)/2)][3];
			phi[2]=val[0]*ft2t2[i/((dop-1)*(dop-2)/2)][4]+val[1]*ft2t2[i/((dop-1)*(dop-2)/2)][3]+val[2]*ft2t2[i/((dop-1)*(dop-2)/2)][2];
		}
		else if(index<24+(dop-1)*36+(dop-1)*(dop-2)/2*24)
		{
			i=index-24-(dop-1)*36-(dop-1)*(dop-2)/2*20;
			lagrange3d_basis1(lambda, grd_lambda, 4+(dop-1)*6+i, dop, val);
			phi[0]=val[0]*ft1t2[i/((dop-1)*(dop-2)/2)][0]+val[1]*ft1t2[i/((dop-1)*(dop-2)/2)][5]+val[2]*ft1t2[i/((dop-1)*(dop-2)/2)][4];
			phi[1]=val[0]*ft1t2[i/((dop-1)*(dop-2)/2)][5]+val[1]*ft1t2[i/((dop-1)*(dop-2)/2)][1]+val[2]*ft1t2[i/((dop-1)*(dop-2)/2)][3];
			phi[2]=val[0]*ft1t2[i/((dop-1)*(dop-2)/2)][4]+val[1]*ft1t2[i/((dop-1)*(dop-2)/2)][3]+val[2]*ft1t2[i/((dop-1)*(dop-2)/2)][2];
		}
		else if(index<12*(dop*dop+1)+(dofs-2*(dop*dop+1)))
		{
			i=index-12*(dop*dop+1);
			lagrange3d_basis1(lambda, grd_lambda, 2*(dop*dop+1)+i, dop, val);
			phi[0] = val[0]; phi[1] = 0; phi[2] = 0;
		}
		else if(index<12*(dop*dop+1)+(dofs-2*(dop*dop+1))*2)
		{
			i=index-12*(dop*dop+1)-(dofs-2*(dop*dop+1));
			lagrange3d_basis1(lambda, grd_lambda, 2*(dop*dop+1)+i, dop, val);
			phi[0] = 0;	phi[1] = val[1]; phi[2] = 0;
		}
		else if(index<12*(dop*dop+1)+(dofs-2*(dop*dop+1))*3)
		{
			i=index-12*(dop*dop+1)-(dofs-2*(dop*dop+1))*2;
			lagrange3d_basis1(lambda, grd_lambda, 2*(dop*dop+1)+i, dop, val);
			phi[0] = 0;	phi[1] = 0; phi[2] = val[2];
		}
		else if(index<12*(dop*dop+1)+(dofs-2*(dop*dop+1))*4)
		{
			i=index-12*(dop*dop+1)-(dofs-2*(dop*dop+1))*3;
			lagrange3d_basis1(lambda, grd_lambda, 2*(dop*dop+1)+i, dop, val);
			phi[0] = 0;	phi[1] = val[2]; phi[2] = val[1];
		}
		else if(index<12*(dop*dop+1)+(dofs-2*(dop*dop+1))*5)
		{
			i=index-12*(dop*dop+1)-(dofs-2*(dop*dop+1))*4;
			lagrange3d_basis1(lambda, grd_lambda, 2*(dop*dop+1)+i, dop, val);
			phi[0] = val[2]; phi[1] = 0; phi[2] = val[0];
		}
		else
		{
			i=index-12*(dop*dop+1)-(dofs-2*(dop*dop+1))*5;
			lagrange3d_basis1(lambda, grd_lambda, 2*(dop*dop+1)+i, dop, val);
			phi[0] = val[1]; phi[1] = val[0]; phi[2] = 0;
		}
	} 
}

/**
* \fn void nedelec1st3d_basis(double *lambda, double **grd_lambda, short *eorien, int **fpermi, int index, int dop, double phi[3])
* \brief basis function of the first kind Nedelec element in three dimensions
* \param *lambda pointer to the area coordiante
* \param **grd_lambda pointer to the gradient of the volume coordiante lambda
* \param **etv the unit tangential vectors of the six edges
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi basis function
* \return void
*/
void nedelec1st3d_basis(double *lambda, double **grd_lambda, short *eorien, int **fpermi, int index, int dop, double phi[3])
{
	int fi[3], ei[2], vi[3], i0, i1, i2, *permi;
	int i, j, ii, idx;
	int dofs = dop*(dop + 2)*(dop + 3) / 2; // degrees of freedom
	double c1, c2;
	
	init_array(3, phi, 0);

	if (index >= dofs || index<0)
		return;

	if (dop < 1)
		return;

	double val;
	
	if (dop == 1)
	{
		// edge: 01, 12, 23, 30, 02, 13 
		//	ei[0] = (index-4)%4; for Lagrange
		//	ei[1] = (ei[0]+index/4)%4; for Lagrange
		edge2vv3d(index, ei);
		c1 = lambda[ei[0]];
		c2 = -lambda[ei[1]];
		axpbyz_array(3, c1, grd_lambda[ei[1]], c2, grd_lambda[ei[0]], phi);
		ax_array(3, eorien[index], phi);
	}

	else if (dop == 2)
	{
		if (index < 12) // dop * 6
		{
			idx = index / dop;
			edge2vv3d(idx, ei);
			ii = index%dop; // ii=0 or 1 if dop=2
			c1 = lambda[ei[ii]] * lambda[ei[0]];
			c2 = -lambda[ei[ii]] * lambda[ei[1]];
			axpbyz_array(3, c1, grd_lambda[ei[1]], c2, grd_lambda[ei[0]], phi);
			ax_array(3, eorien[idx], phi);
			//	phi[j] = h[idx] * lambda[ei[ii]] * (lambda[ei[perm[0]]] * grd_lambda[ei[perm[1]]][j] - lambda[ei[perm[1]]] * grd_lambda[ei[perm[0]]][j]);
		}
		else if (index < 20)
		{
			index -= 12;
			idx = index / 2;
			face2vertices3d(idx, vi);
			permi = fpermi[idx];
			fi[0] = vi[permi[0]]; fi[1] = vi[permi[1]]; fi[2] = vi[permi[2]];
			ii = index % 2; // ii=0 or 1 if dop=2
			if (ii == 0)
			{
				c1 = lambda[fi[2]] * lambda[fi[0]];
				c2 = -lambda[fi[2]] * lambda[fi[1]];
				axpbyz_array(3, c1, grd_lambda[fi[1]], c2, grd_lambda[fi[0]], phi);
			}
			else
			{
				c1 = lambda[fi[1]] * lambda[fi[0]];
				c2 = -lambda[fi[1]] * lambda[fi[2]];
				axpbyz_array(3, c1, grd_lambda[fi[2]], c2, grd_lambda[fi[0]], phi);
			}
		}
	}

	else if (dop == 3)
	{
		if (index < 18) // dop * 6
		{
			idx = index / dop;
			edge2vv3d(idx, ei);
			ii = index%dop; // ii=0, 1 or 2 if dop=3
			if(ii==0)
			{
				c1 = lambda[ei[0]] * lambda[ei[0]] * lambda[ei[0]];
				c2 = -lambda[ei[0]] * lambda[ei[0]] * lambda[ei[1]];
			}
			else if(ii==1)
			{
				c1 = lambda[ei[0]] * lambda[ei[1]] * lambda[ei[0]];
				c2 = -lambda[ei[0]] * lambda[ei[1]] * lambda[ei[1]];
			}
			else
			{
				c1 = lambda[ei[1]] * lambda[ei[1]] * lambda[ei[0]];
				c2 = -lambda[ei[1]] * lambda[ei[1]] * lambda[ei[1]];
			}
			axpbyz_array(3, c1, grd_lambda[ei[1]], c2, grd_lambda[ei[0]], phi);
			ax_array(3, eorien[idx], phi);
		}
		else if (index < 42) // dop * 6 + (dop - 1)*dop * 4
		{
			index -= 18; // dop * 6
			idx = index / 6; // (dop - 1)*dop
			face2vertices3d(idx, vi);
			permi = fpermi[idx];
			fi[0] = vi[permi[0]]; fi[1] = vi[permi[1]]; fi[2] = vi[permi[2]];
			ii = index % 6; // (dop - 1)*dop
			j = ii / 3; // (dop - 1)*dop/2, j = 0 or 1
			i = ii % 3;
			if (j == 0)
			{
				c1 = lambda[vi[i]] * lambda[fi[2]] * lambda[fi[0]];
				c2 = -lambda[vi[i]] * lambda[fi[2]] * lambda[fi[1]];
				axpbyz_array(3, c1, grd_lambda[fi[1]], c2, grd_lambda[fi[0]], phi);
			}
			else
			{
				c1 = lambda[vi[i]] * lambda[fi[1]] * lambda[fi[0]];
				c2 = -lambda[vi[i]] * lambda[fi[1]] * lambda[fi[2]];
				axpbyz_array(3, c1, grd_lambda[fi[2]], c2, grd_lambda[fi[0]], phi);
			}
		}
		else
		{
			index -= 42; // dop * 6 + (dop - 1)*dop * 4
			if (index == 0)
			{
				c1 = lambda[2] * lambda[3] * lambda[0];
			}
			else if (index == 1)
			{
				c1 = lambda[1] * lambda[3] * lambda[0];
			}
			else
			{
				c1 = lambda[1] * lambda[2] * lambda[0];
			}
			c2 = -lambda[1] * lambda[2] * lambda[3];
			axpbyz_array(3, c1, grd_lambda[index+1], c2, grd_lambda[0], phi);
		}
	}
}

/**
* \fn void nedelec1st3d_basisCurl(double *lambda, double **grd_lambda, short *eorien, int **fpermi, int index, int dop, double phi[3])
* \brief basis function of the first kind Nedelec element in three dimensions
* \param *lambda pointer to the area coordiante
* \param **grd_lambda pointer to the gradient of the volume coordiante lambda
* \param **etv the unit tangential vectors of the six edges
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi basis function
* \return void
*/
void nedelec1st3d_basisCurl(double *lambda, double **grd_lambda, short *eorien, int **fpermi, int index, int dop, double phi[3])
{
	int fi[3], ei[2], vi[3], i0, i1, i2, *permi;
	int i, j, ii, idx;
	int dofs = dop*(dop + 2)*(dop + 3) / 2; // degrees of freedom
	double c1, c2, val0[3], val1[3], val2[3];
	int face;
	
	init_array(3, phi, 0);

	if (index >= dofs || index<0)
		return;

	if (dop < 1)
		return;

	if (dop == 1)
	{
		// edge: 01, 12, 23, 30, 02, 13 
		//	ei[0] = (index-4)%4;
		//	ei[1] = (ei[0]+index/4)%4;
		edge2vv3d(index, ei);
		cross_array(grd_lambda[ei[0]], grd_lambda[ei[1]], phi);
		ax_array(3, 2.0*eorien[index], phi);
	}

	else if (dop == 2)
	{
		if (index < 12)
		{
			idx = index / dop;
			edge2vv3d(idx, ei);
			ii = index%dop; // ii=0 or 1 if dop=2
			cross_array(grd_lambda[ei[0]], grd_lambda[ei[1]], phi);
			ax_array(3, 3*lambda[ei[ii]]*eorien[idx], phi);
		}
		else if (index < 20)
		{
			index -= 12;
			idx = index / 2;
			face2vertices3d(idx, vi);
			permi = fpermi[idx];
			fi[0] = vi[permi[0]]; fi[1] = vi[permi[1]]; fi[2] = vi[permi[2]];
			ii = index % 2; // ii=0 or 1 if dop=2
			if (ii == 0)
			{
				cross_array(grd_lambda[fi[2]], grd_lambda[fi[1]], val0);
				cross_array(grd_lambda[fi[0]], grd_lambda[fi[2]], val1);
				cross_array(grd_lambda[fi[0]], grd_lambda[fi[1]], val2);
				axpbyz_array(3, lambda[fi[0]], val0, lambda[fi[1]], val1, phi);
				axpy_array(3, 2*lambda[fi[2]], val2, phi);
			}
			else
			{
				cross_array(grd_lambda[fi[1]], grd_lambda[fi[2]], val0);
				cross_array(grd_lambda[fi[0]], grd_lambda[fi[1]], val2);
				cross_array(grd_lambda[fi[0]], grd_lambda[fi[2]], val1);
				axpbyz_array(3, lambda[fi[0]], val0, lambda[fi[2]], val2, phi);
				axpy_array(3, 2*lambda[fi[1]], val1, phi);	
			}
		}
	}

	else if (dop == 3)
	{
		if (index < 18) // dop * 6
		{
			idx = index / dop;
			edge2vv3d(idx, ei);
			ii = index%dop; // ii=0, 1 or 2 if dop=3
			
			cross_array(grd_lambda[ei[0]], grd_lambda[ei[1]], phi);
			if(ii==0)
			{
				ax_array(3, 4*lambda[ei[0]]*lambda[ei[0]], phi);
			}
			else if(ii==1)
			{
				ax_array(3, 4*lambda[ei[0]]*lambda[ei[1]], phi);
			}
			else
			{
				ax_array(3, 4*lambda[ei[1]]*lambda[ei[1]], phi);
			}
			ax_array(3, eorien[idx], phi);
		}
		else if (index < 42) // dop * 6 + (dop - 1)*dop * 4
		{
			index -= 18; // dop * 6
			idx = index / 6; // (dop - 1)*dop
			face2vertices3d(idx, vi);
			permi = fpermi[idx];
			fi[0] = vi[permi[0]]; fi[1] = vi[permi[1]]; fi[2] = vi[permi[2]];
			ii = index % 6; // (dop - 1)*dop
			j = ii / 3; // (dop - 1)*dop/2, j = 0 or 1
			i = ii % 3;
			if (j == 0)
			{
				axpbyz_array(3, lambda[fi[2]]*lambda[fi[0]], grd_lambda[vi[i]], lambda[vi[i]]*lambda[fi[0]], grd_lambda[fi[2]], val0);
				axpy_array(3, lambda[vi[i]]*lambda[fi[2]], grd_lambda[fi[0]], val0);
				cross_array(val0, grd_lambda[fi[1]], val1);

				axpbyz_array(3, lambda[fi[2]]*lambda[fi[1]], grd_lambda[vi[i]], lambda[vi[i]]*lambda[fi[1]], grd_lambda[fi[2]], val0);
				axpy_array(3, lambda[vi[i]]*lambda[fi[2]], grd_lambda[fi[1]], val0);
				cross_array(val0, grd_lambda[fi[0]], val2);
			}
			else
			{
				axpbyz_array(3, lambda[fi[1]]*lambda[fi[0]], grd_lambda[vi[i]], lambda[vi[i]]*lambda[fi[0]], grd_lambda[fi[1]], val0);
				axpy_array(3, lambda[vi[i]]*lambda[fi[1]], grd_lambda[fi[0]], val0);
				cross_array(val0, grd_lambda[fi[2]], val1);

				axpbyz_array(3, lambda[fi[1]]*lambda[fi[2]], grd_lambda[vi[i]], lambda[vi[i]]*lambda[fi[2]], grd_lambda[fi[1]], val0);
				axpy_array(3, lambda[vi[i]]*lambda[fi[1]], grd_lambda[fi[2]], val0);
				cross_array(val0, grd_lambda[fi[0]], val2);
			}
			axpbyz_array(3, 1.0, val1, -1.0, val2, phi);
		}
		else
		{
			index -= 42; // dop * 6 + (dop - 1)*dop * 4
			if (index == 0)
			{
				axpbyz_array(3, lambda[3]*lambda[0], grd_lambda[2], lambda[2]*lambda[0], grd_lambda[3], val0);
				axpy_array(3, lambda[2]*lambda[3], grd_lambda[0], val0);
				cross_array(val0, grd_lambda[1], val1);
			}
			else if (index == 1)
			{
				axpbyz_array(3, lambda[3]*lambda[0], grd_lambda[1], lambda[1]*lambda[0], grd_lambda[3], val0);
				axpy_array(3, lambda[1]*lambda[3], grd_lambda[0], val0);
				cross_array(val0, grd_lambda[2], val1);
			}
			else
			{
				axpbyz_array(3, lambda[2]*lambda[0], grd_lambda[1], lambda[1]*lambda[0], grd_lambda[2], val0);
				axpy_array(3, lambda[1]*lambda[2], grd_lambda[0], val0);
				cross_array(val0, grd_lambda[3], val1);
			}
			axpbyz_array(3, lambda[2]*lambda[3], grd_lambda[1], lambda[1]*lambda[3], grd_lambda[2], val0);
			axpy_array(3, lambda[1]*lambda[2], grd_lambda[3], val0);
			cross_array(val0, grd_lambda[0], val2);
			axpbyz_array(3, 1.0, val1, -1.0, val2, phi);
		}
	}
}

/**
* \fn void nedelec2nd3d_basis(double *lambda, double **grd_lambda, int **epermi, int index, int dop, double phi[3])
* \brief basis function of the second kind Nedelec element in three dimensions
* \param *lambda pointer to the area coordiante
* \param **grd_lambda pointer to the gradient of the volume coordiante lambda
* \param **etv the unit tangential vectors of the six edges
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi basis function
* \return void
*/
void nedelec2nd3d_basis(double *lambda, double **grd_lambda, int **eperm, int index, int dop, double phi[3])
{
	int fi[3], ei[2], vi[3], i0, i1, i2, *perm;
	int i, j, ii, idx;
	int dofs = (dop + 1)*(dop + 2)*(dop + 3) / 2; // degrees of freedom
	
	init_array(3, phi, 0);

	if (index >= dofs || index<0)
		return;

	if (dop > 2 || dop < 1)
		return;

	double val;

	if (dop == 1)
	{
		// edge: 01, 12, 23, 30, 02, 13 
		//	ei[0] = (index-4)%4; for Lagrange
		//	ei[1] = (ei[0]+index/4)%4; for Lagrange
		idx = index / 2;
		edge2vv3d(idx, ei);
		ii = index % 2;
		if (ii == 0)
			axy_array(3, lambda[ei[0]], grd_lambda[ei[1]], phi);
		else
			axy_array(3, lambda[ei[1]], grd_lambda[ei[0]], phi);
	}

	else if (dop == 2)
	{
		if (index < 18)
		{
			idx = index / 3;
			edge2vv3d(idx, vi);
			perm = eperm[idx];
			ei[0] = vi[0]; ei[1] = vi[1];
			ii = index % 3;
			if (ii == 0)
				axy_array(3, lambda[ei[0]] * lambda[ei[0]], grd_lambda[ei[1]], phi);
			else if (ii == 1)
			{
				ei[0] = vi[perm[0]]; ei[1] = vi[perm[1]];
				copy_array(3, grd_lambda[ei[1]], phi);
				axpy_array(3, -1.0, grd_lambda[ei[0]], phi);
				ax_array(3, lambda[ei[0]] * lambda[ei[1]], phi);
			}
			else
				axy_array(3, lambda[ei[1]] * lambda[ei[1]], grd_lambda[ei[0]], phi);
		}
		else if (index < 30)
		{
			index -= 18;
			idx = index / 3;
			face2vertices3d(idx, fi);
			ii = index % 3;
			if (ii == 0)
				axy_array(3, lambda[fi[1]] * lambda[fi[2]], grd_lambda[fi[0]], phi);
			else if (ii == 1)
				axy_array(3, lambda[fi[2]] * lambda[fi[0]], grd_lambda[fi[1]], phi);
			else if (ii == 2)
				axy_array(3, lambda[fi[0]] * lambda[fi[1]], grd_lambda[fi[2]], phi);
		}
	}
}

/**
* \fn void nedelec2nd3d_basisCurl(double *lambda, double v, double s[4], double h[6], double **grd_lambda, ELEMENT *elements, idenmat *elementFace, FACE *faces, idenmat *elementEdge, EDGE *edges, dennode *nodes, int element, int index, int dop, double phi[3])
* \brief basis function of the second kind Nedelec element in three dimensions
* \param *lambda pointer to the area coordiante
* \param **grd_lambda pointer to the gradient of the volume coordiante lambda
* \param **etv the unit tangential vectors of the six edges
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi basis function
* \return void
*/
void nedelec2nd3d_basisCurl(double *lambda, double **grd_lambda, int **eperm, int index, int dop, double phi[3])
{
	int fi[3], ei[2], vi[3], i0, i1, i2, *perm;
	int i, j, ii, idx;
	int dofs = (dop + 1)*(dop + 2)*(dop + 3) / 2; // degrees of freedom
	double val1[3], val2[3];

	init_array(3, phi, 0);

	if (index >= dofs || index<0)
		return;

	if (dop > 2 || dop < 1)
		return;

	if (dop == 1)
	{
		// edge: 01, 12, 23, 30, 02, 13 
		//	ei[0] = (index-4)%4;
		//	ei[1] = (ei[0]+index/4)%4;
		idx = index / 2;
		edge2vv3d(idx, ei);
		ii = index % 2;
		if (ii == 0)
			cross_array(grd_lambda[ei[0]], grd_lambda[ei[1]], phi);
		else
			cross_array(grd_lambda[ei[1]], grd_lambda[ei[0]], phi);
	}

	else if (dop == 2)
	{
		if (index < 18)
		{
			idx = index / 3;
			edge2vv3d(idx, vi);
			perm = eperm[idx];
			ei[0] = vi[0]; ei[1] = vi[1];
			ii = index % 3;
			if (ii == 0)
			{
				cross_array(grd_lambda[ei[0]], grd_lambda[ei[1]], phi);
				ax_array(3, 2 * lambda[ei[0]], phi);
			}
			else if (ii == 1)
			{
				ei[0] = vi[perm[0]]; ei[1] = vi[perm[1]];
				axpbyz_array(3, lambda[ei[1]], grd_lambda[ei[0]], lambda[ei[0]], grd_lambda[ei[1]], val1);
				copy_array(3, grd_lambda[ei[1]], val2);
				axpy_array(3, -1.0, grd_lambda[ei[0]], val2);
				cross_array(val1, val2, phi);
			}
			else
			{
				cross_array(grd_lambda[ei[1]], grd_lambda[ei[0]], phi);
				ax_array(3, 2 * lambda[ei[1]], phi);
			}
		}
		else if (index < 30)
		{
			index -= 18;
			idx = index / 3;
			face2vertices3d(idx, fi);
			ii = index % 3;
			if (ii == 0)
			{
				axpbyz_array(3, lambda[fi[1]], grd_lambda[fi[2]], lambda[fi[2]], grd_lambda[fi[1]], val1);
				cross_array(val1, grd_lambda[fi[0]], phi);
			}
			else if (ii == 1)
			{
				axpbyz_array(3, lambda[fi[2]], grd_lambda[fi[0]], lambda[fi[0]], grd_lambda[fi[2]], val1);
				cross_array(val1, grd_lambda[fi[1]], phi);
			}
			else
			{
				axpbyz_array(3, lambda[fi[0]], grd_lambda[fi[1]], lambda[fi[1]], grd_lambda[fi[0]], val1);
				cross_array(val1, grd_lambda[fi[2]], phi);
			}
		}
	}
}

/**
* \fn void chhcurlHermite3d_basis(double *lambda, double **grd_lambda, double **etv, int index, int dop, double phi[3])
* \brief basis function of Hermite-type Christiansen-Hu-Hu H(curl) element in three dimensions
* \param *lambda pointer to the area coordiante
* \param **grd_lambda pointer to the gradient of the volume coordiante lambda
* \param **etv the unit tangential vectors of the six edges
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi basis function
* \return void
*/
void chhcurlHermite3d_basis(double *lambda, double **grd_lambda, double **etv, int index, int dop, double phi[3])
{
	int fi, ei[2], vi[3], i0, i1, i2;
	int i, j;
	int dofs = (dop + 1)*(dop + 2)*(dop + 3) / 6; // degrees of freedom
	
	init_array(3, phi, 0);

	if (index >= dofs * 3 || index<0)
		return;

	if (dop > 2)
		return;

	double val;

	if (dop == 0)
	{
		phi[index] = 1;
	}

	else // dop>=1
	{
		if (index<12)
		{
			lagrange3d_basis(lambda, index % 4, dop, &val);
			phi[index / 4] = val;
		}
		else if (index<12 + (dop - 1) * 6)
		{
			i = index - 12;
			lagrange3d_basis(lambda, 4 + i, dop, &val);
			axy_array(3, val, etv[i / (dop - 1)], phi);
		}
		else if (index<12 + (dop - 1) * 6 + (dop - 1)*(dop + 1) * 4)
		{
			if (dop == 2)
			{
				i = index - 12 - (dop - 1) * 6;
				fi = i / ((dop - 1)*(dop + 1));
				face2vertices3d(fi, vi);
				i0 = i % ((dop - 1)*(dop + 1));
				i1 = (i0 + 1) % 3;
				i2 = (i0 + 2) % 3;
				axy_array(3, lambda[vi[i1]] * lambda[vi[i2]], grd_lambda[vi[i0]], phi);
			}
		}
		else
		{
		}
	}
}

/**
* \fn void chhcurlHermite3d_basisCurl(double *lambda, double **grd_lambda, double **etv, int index, int dop, double phi[3])
* \brief basis function of Hermite-type Christiansen-Hu-Hu H(curl) element in three dimensions
* \param *lambda pointer to the area coordiante
* \param **grd_lambda pointer to the gradient of the volume coordiante lambda
* \param **etv the unit tangential vectors of the six edges
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi basis function
* \return void
*/
void chhcurlHermite3d_basisCurl(double *lambda, double **grd_lambda, double **etv, int index, int dop, double phi[3])
{
	int fi, ei[2], vi[3], i0, i1, i2;
	int i, j;
	int dofs = (dop + 1)*(dop + 2)*(dop + 3) / 6; // degrees of freedom

	init_array(3, phi, 0);

	if (index >= dofs * 3 || index<0)
		return;

	if (dop > 2)
		return;

	double val[3];

	if (dop == 0)
	{
		return;
	}

	else // dop>=1
	{
		if (index<4)
		{
			lagrange3d_basis1(lambda, grd_lambda, index, dop, val); 
			phi[0] = 0;	phi[1] = val[2]; phi[2] = -val[1];
		}
		else if (index < 8)
		{
			lagrange3d_basis1(lambda, grd_lambda, index % 4, dop, val);
			phi[0] = -val[2]; phi[1] = 0; phi[2] = val[0];
		}
		else if (index < 12)
		{
			lagrange3d_basis1(lambda, grd_lambda, index % 4, dop, val);
			phi[0] = val[1]; phi[1] = -val[0]; phi[2] = 0;
		}
		else if (index<12 + (dop - 1) * 6)
		{
			i = index - 12;
			lagrange3d_basis1(lambda, grd_lambda, 4 + i, dop, val);
			cross_array(val, etv[i / (dop - 1)], phi);
		}
		else if (index<12 + (dop - 1) * 6 + (dop - 1)*(dop + 1) * 4)
		{
			if (dop == 2)
			{
				i = index - 12 - (dop - 1) * 6;
				fi = i / ((dop - 1)*(dop + 1));
				face2vertices3d(fi, vi);
				i0 = i % ((dop - 1)*(dop + 1));
				i1 = (i0 + 1) % 3;
				i2 = (i0 + 2) % 3;
				axpbyz_array(3, lambda[vi[i2]], grd_lambda[vi[i1]], lambda[vi[i1]], grd_lambda[vi[i2]], val);
				cross_array(val, grd_lambda[vi[i0]], phi);
			}
		}
		else
		{
		}
	}
}

/**
* \fn void huangQuadcurl3d_facebasis(double *x, double *xK, double *lambda, double **grd_lambda, double **vertices, double *nvf[4], int l, int li, double phi[3])
* \brief face basis function of the Huang element for Quad-curl probelm in three dimensions
* \param *x pointer to the Cartesian coordinates
* \param *xK pointer to the barycenter of the tetrahedron
* \param *lambda pointer to the barycentric coordinates
* \param **grd_lambda pointer to the gradient of the volume coordiante lambda
* \param **vertices pointer to the Cartesian coordinates of the four vertices of the tetrahedron
* \param *nvf[4] the unit fixed normal vectors of the fourth faces
* \param l the index of the face
* \param li the index of one vertex on the face l
* \param *phi basis function
* \return void
*/
void huangQuadcurl3d_facebasis(double *x, double *xK, double *lambda, double **grd_lambda, double **vertices, double *nvf[4], int l, int li, double phi[3])
{
	double v0[3], v1[3], y[3], c;

	axpyz_array(3, -1.0, xK, x, y); // y = x - xK
	cross_array(nvf[l], grd_lambda[li], v0);
	cross_array(y, v0, v1);
	axpyz_array(3, -1.0, xK, vertices[l], y); // y = xl - xK
	c = dot_array(3, y, nvf[l]); 
	axpbyz_array(3, 2*lambda[l]-0.75, v1, c/4., grd_lambda[li], phi);
	axpy_array(3, 1./16., nvf[l], phi);
}

/**
* \fn void huangQuadcurl3d_facebasisCurl(double *lambda, double **grd_lambda, double *nvf[4], int l, int li, double phi[3])
* \brief curl of the face basis function of the Huang element for Quad-curl probelm in three dimensions
* \param *lambda pointer to the barycentric coordinates
* \param **grd_lambda pointer to the gradient of the volume coordiante lambda
* \param *nvf[4] the unit fixed normal vectors of the fourth faces
* \param l the index of the face
* \param li the index of one vertex on the face l
* \param *phi basis function
* \return void
*/
void huangQuadcurl3d_facebasisCurl(double *lambda, double **grd_lambda, double *nvf[4], int l, int li, double phi[3])
{
	cross_array(nvf[l], grd_lambda[li], phi);
	ax_array(3, 2*(1-3*lambda[l]), phi);
}

/**
* \fn void huangQuadcurl3d_facebasisGradCurl(double **grd_lambda, double *nvf[4], int l, int li, double phi[9])
* \brief gradcurl of the face basis function of the Huang element for Quad-curl probelm in three dimensions
* \param **grd_lambda pointer to the gradient of the volume coordiante lambda
* \param *nvf[4] the unit fixed normal vectors of the fourth faces
* \param l the index of the face
* \param li the index of one vertex on the face l
* \param *phi basis function
* \return void
*/
void huangQuadcurl3d_facebasisGradCurl(double **grd_lambda, double *nvf[4], int l, int li, double phi[9])
{
	double v0[3];
	cross_array(nvf[l], grd_lambda[li], v0);
	for(int i=0;i<3;i++)
		axy_array(3, -6*v0[i], grd_lambda[l], phi+3*i);
}

/**
* \fn void huangQuadcurl3d_basis(double *x, double *xK, double *lambda, double **grd_lambda, double **vertices, double *nvf[4], short *eorien, int **fpermi, int index, double phi[3])
* \brief basis function of the Huang element for Quad-curl probelm in three dimensions
* \param *x pointer to the Cartesian coordinates
* \param *xK pointer to the barycenter of the tetrahedron
* \param *lambda pointer to the area coordiante
* \param **grd_lambda pointer to the gradient of the volume coordiante lambda
* \param **etv the unit tangential vectors of the six edges
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi basis function
* \return void
*/
void huangQuadcurl3d_basis(double *x, double *xK, double *lambda, double **grd_lambda, double **vertices, double *nvf[4], short *eorien, int **fpermi, int index, double phi[3])
{
	int fi[3], ei[2], vi[3], *permi;
	int i, l, ii, idx;
	int dofs = 14; // number of degrees of freedom
	double v0[3], v1[3], v2[3], c;
	
	init_array(3, phi, 0);

	if (index >= dofs || index<0)
		return;

	if (index < 6)
	{
		// edge: 01, 12, 23, 30, 02, 13 
		//	ei[0] = (index-4)%4; for Lagrange
		//	ei[1] = (ei[0]+index/4)%4; for Lagrange
		edge2vv3d(index, ei);
		axpbyz_array(3, lambda[ei[0]], grd_lambda[ei[1]], -lambda[ei[1]], grd_lambda[ei[0]], phi);
		cross_array(grd_lambda[ei[1]], grd_lambda[ei[0]], v0);
		for(l=0;l<4;l++)
		{
			face2vertices3d(l, vi);
			permi = fpermi[l];
			fi[0] = vi[permi[0]]; fi[1] = vi[permi[1]]; fi[2] = vi[permi[2]];
			for(i=0;i<2;i++)
			{
				cross_array(grd_lambda[fi[1-i]], nvf[l], v1);
				cross_array(nvf[l], v1, v2);
				c = dot_array(3, v0, v2) / dot_array(3, grd_lambda[fi[i]], v1);
				huangQuadcurl3d_facebasis(x, xK, lambda, grd_lambda, vertices, nvf, l, fi[i], v1);
				axpy_array(3, c, v1, phi);
			}
		}
		ax_array(3, eorien[index], phi);
	}
	else if (index < 14)
	{
		index -= 6;
		idx = index / 2;
		face2vertices3d(idx, vi);
		permi = fpermi[idx];
		fi[0] = vi[permi[0]]; fi[1] = vi[permi[1]]; fi[2] = vi[permi[2]];
		ii = index % 2; // ii=0 or 1
		huangQuadcurl3d_facebasis(x, xK, lambda, grd_lambda, vertices, nvf, idx, fi[ii], phi);
	}
}

/**
* \fn void huangQuadcurl3d_basisCurl(double *lambda, double **grd_lambda, double *nvf[4], short *eorien, int **fpermi, int index, double phi[3])
* \brief curl of the basis function of the Huang element for Quad-curl probelm in three dimensions
* \param *lambda pointer to the area coordiante
* \param **grd_lambda pointer to the gradient of the volume coordiante lambda
* \param **etv the unit tangential vectors of the six edges
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi basis function
* \return void
*/
void huangQuadcurl3d_basisCurl(double *lambda, double **grd_lambda, double *nvf[4], short *eorien, int **fpermi, int index, double phi[3])
{
	int fi[3], ei[2], vi[3], *permi;
	int i, l, ii, idx;
	int dofs = 14; // number of degrees of freedom
	double v0[3], v1[3], v2[3], c;
	
	init_array(3, phi, 0);

	if (index >= dofs || index<0)
		return;

	if (index < 6)
	{
		// edge: 01, 12, 23, 30, 02, 13 
		//	ei[0] = (index-4)%4; for Lagrange
		//	ei[1] = (ei[0]+index/4)%4; for Lagrange
		edge2vv3d(index, ei);
		cross_array(grd_lambda[ei[0]], grd_lambda[ei[1]], phi);
		ax_array(3, 2.0, phi);
		cross_array(grd_lambda[ei[1]], grd_lambda[ei[0]], v0);
		for(l=0;l<4;l++)
		{
			face2vertices3d(l, vi);
			permi = fpermi[l];
			fi[0] = vi[permi[0]]; fi[1] = vi[permi[1]]; fi[2] = vi[permi[2]];
			for(i=0;i<2;i++)
			{
				cross_array(grd_lambda[fi[1-i]], nvf[l], v1);
				cross_array(nvf[l], v1, v2);
				c = dot_array(3, v0, v2) / dot_array(3, grd_lambda[fi[i]], v1);
				huangQuadcurl3d_facebasisCurl(lambda, grd_lambda, nvf, l, fi[i], v1);
				axpy_array(3, c, v1, phi);
			}
		}
		ax_array(3, eorien[index], phi);
	}
	else if (index < 14)
	{
		index -= 6;
		idx = index / 2;
		face2vertices3d(idx, vi);
		permi = fpermi[idx];
		fi[0] = vi[permi[0]]; fi[1] = vi[permi[1]]; fi[2] = vi[permi[2]];
		ii = index % 2; // ii=0 or 1
		huangQuadcurl3d_facebasisCurl(lambda, grd_lambda, nvf, idx, fi[ii], phi);
	}
}

/**
* \fn void huangQuadcurl3d_basisGradCurl(double **grd_lambda, double *nvf[4], short *eorien, int **fpermi, int index, double phi[9])
* \brief gradcurl of the basis function of the Huang element for Quad-curl probelm in three dimensions
* \param *lambda pointer to the area coordiante
* \param **grd_lambda pointer to the gradient of the volume coordiante lambda
* \param **etv the unit tangential vectors of the six edges
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi basis function
* \return void
*/
void huangQuadcurl3d_basisGradCurl(double **grd_lambda, double *nvf[4], short *eorien, int **fpermi, int index, double phi[9])
{
	int fi[3], ei[2], vi[3], *permi;
	int i, l, ii, idx;
	int dofs = 14; // number of degrees of freedom
	double v0[9], v1[9], v2[9], c;
	
	init_array(9, phi, 0);

	if (index >= dofs || index<0)
		return;

	if (index < 6)
	{
		// edge: 01, 12, 23, 30, 02, 13 
		//	ei[0] = (index-4)%4; for Lagrange
		//	ei[1] = (ei[0]+index/4)%4; for Lagrange
		edge2vv3d(index, ei);
		cross_array(grd_lambda[ei[1]], grd_lambda[ei[0]], v0);
		for(l=0;l<4;l++)
		{
			face2vertices3d(l, vi);
			permi = fpermi[l];
			fi[0] = vi[permi[0]]; fi[1] = vi[permi[1]]; fi[2] = vi[permi[2]];
			for(i=0;i<2;i++)
			{
				cross_array(grd_lambda[fi[1-i]], nvf[l], v1);
				cross_array(nvf[l], v1, v2);
				c = dot_array(3, v0, v2) / dot_array(3, grd_lambda[fi[i]], v1);
				huangQuadcurl3d_facebasisGradCurl(grd_lambda, nvf, l, fi[i], v1);
				axpy_array(9, c, v1, phi);
			}
		}
		ax_array(9, eorien[index], phi);
	}
	else if (index < 14)
	{
		index -= 6;
		idx = index / 2;
		face2vertices3d(idx, vi);
		permi = fpermi[idx];
		fi[0] = vi[permi[0]]; fi[1] = vi[permi[1]]; fi[2] = vi[permi[2]];
		ii = index % 2; // ii=0 or 1
		huangQuadcurl3d_facebasisGradCurl(grd_lambda, nvf, idx, fi[ii], phi);
	}
}

/**
* \fn void huangzhang03d_basis(double *lambda, double **grd_lambda, int index, double phi[3])
* \brief basis function of Huang-Zhang element in three dimensions, but not dual to DoFs
* \param *lambda pointer to the area coordiante
* \param **grd_lambda pointer to the gradient of the volume coordiante lambda
* \param index the indicator of the basis function
* \param *phi basis function
* \return void
*/
void huangzhang03d_basis(double *lambda, double **grd_lambda, int index, double phi[3])
{
	int fi[3], ei[2], vi[3];
	int i, j, ii, idx;
	int dofs = 32; // degrees of freedom
	int dop = 2;
	double c1, c2;
	
	init_array(3, phi, 0);

	if (index >= dofs || index<0)
		return;

	double val;
	
	if (index < 12) // dop * 6
	{
		idx = index / dop;
		edge2vv3d(idx, ei);
		ii = index % dop; // ii=0 or 1
		c1 = lambda[ei[ii]] * lambda[ei[0]];
		c2 = -lambda[ei[ii]] * lambda[ei[1]];
		axpbyz_array(3, c1, grd_lambda[ei[1]], c2, grd_lambda[ei[0]], phi);
	}
	else if (index < 20)
	{
		index -= 12;
		idx = index / 2;
		face2vertices3d(idx, fi);
		ii = index % 2; // ii=0 or 1
		if (ii == 0)
		{
			c1 = lambda[fi[2]] * lambda[fi[0]];
			c2 = -lambda[fi[2]] * lambda[fi[1]];
			axpbyz_array(3, c1, grd_lambda[fi[1]], c2, grd_lambda[fi[0]], phi);
		}
		else
		{
			c1 = lambda[fi[1]] * lambda[fi[0]];
			c2 = -lambda[fi[1]] * lambda[fi[2]];
			axpbyz_array(3, c1, grd_lambda[fi[2]], c2, grd_lambda[fi[0]], phi);
		}
	}
	else // 20 <= index < 32
	{
		index -= 20;
		i = index / 4;
		j = index % 4;
		// b_K * lambda[j] * e_i
		phi[i] = lambda[0]*lambda[1]*lambda[2]*lambda[3]*lambda[j];
	}
}

/**
* \fn void huangzhang03d_basisCurl(double *lambda, double **grd_lambda, int index, double phi[3])
* \brief basis function of Huang-Zhang element in three dimensions
* \param *lambda pointer to the area coordiante
* \param **grd_lambda pointer to the gradient of the volume coordiante lambda
* \param index the indicator of the basis function
* \param *phi basis function
* \return void
*/
void huangzhang03d_basisCurl(double *lambda, double **grd_lambda, int index, double phi[3])
{
	int fi[3], ei[2], vi[3], i0, i1, i2;
	int i, j, ii, idx;
	int dofs = 32; // degrees of freedom
	int dop = 2;
	double val0[3], val1[3], val2[3];
	
	init_array(3, phi, 0);

	if (index >= dofs || index<0)
		return;

	if (index < 12)
	{
		idx = index / dop;
		edge2vv3d(idx, ei);
		ii = index % dop; // ii=0 or 1 if dop=2
		cross_array(grd_lambda[ei[0]], grd_lambda[ei[1]], phi);
		ax_array(3, 3*lambda[ei[ii]], phi); 
	}
	else if (index < 20)
	{
		index -= 12;
		idx = index / 2;
		face2vertices3d(idx, fi);
		ii = index % 2; // ii=0 or 1
		if (ii == 0)
		{
			cross_array(grd_lambda[fi[2]], grd_lambda[fi[1]], val0);
			cross_array(grd_lambda[fi[0]], grd_lambda[fi[2]], val1);
			cross_array(grd_lambda[fi[0]], grd_lambda[fi[1]], val2);
			axpbyz_array(3, lambda[fi[0]], val0, lambda[fi[1]], val1, phi);
			axpy_array(3, 2*lambda[fi[2]], val2, phi);
		}
		else
		{
			cross_array(grd_lambda[fi[1]], grd_lambda[fi[2]], val0);
			cross_array(grd_lambda[fi[0]], grd_lambda[fi[1]], val2);
			cross_array(grd_lambda[fi[0]], grd_lambda[fi[2]], val1);
			axpbyz_array(3, lambda[fi[0]], val0, lambda[fi[2]], val2, phi);
			axpy_array(3, 2*lambda[fi[1]], val1, phi);
		}
	}
	else // 20 <= index < 32
	{
		index -= 20;
		i = index / 4;
		j = index % 4;
		// (lambda[j] * grad b_K + b_K * grad lambda[j]) cross e_i
		axpbyz_array(3, lambda[1], grd_lambda[0], lambda[0], grd_lambda[1], val0);
		axpbyz_array(3, lambda[3], grd_lambda[2], lambda[2], grd_lambda[3], val1);
		axpbyz_array(3, lambda[2]*lambda[3], val0, lambda[0]*lambda[1], val1, val2);
		axpbyz_array(3, lambda[j], val2, lambda[0]*lambda[1]*lambda[2]*lambda[3], grd_lambda[j], val0);
		i1 = (i+1) % 3; i2 = (i+2) % 3;
		phi[i1] = val0[i2];
		phi[i2] = -val0[i1];
	}
}

/**
* \fn void huangzhang03d_basisGradCurl(double *lambda, double **grd_lambda, int index, double phi[9])
* \brief basis function of Huang-Zhang element in three dimensions
* \param *lambda pointer to the area coordiante
* \param **grd_lambda pointer to the gradient of the volume coordiante lambda
* \param index the indicator of the basis function
* \param *phi basis function
* \return void
*/
void huangzhang03d_basisGradCurl(double *lambda, double **grd_lambda, int index, double phi[9])
{
	int fi[3], ei[2], vi[3], i0, i1, i2;
	int i, j, ii, idx;
	int dofs = 32; // degrees of freedom
	int dop = 2;
	double val0[3], val1[3], val2[3];
	double v0[3];
	double gradlambdaXei[4][3], gradLambdaij[4][4][3];
	
	init_array(9, phi, 0);

	if (index >= dofs || index<0)
		return;

	if (index < 12)
	{
		idx = index / dop;
		edge2vv3d(idx, ei);
		ii = index % dop; // ii=0 or 1 if dop=2
		cross_array(grd_lambda[ei[0]], grd_lambda[ei[1]], v0);
		for(i=0;i<3;i++)
			axy_array(3, 3*v0[i], grd_lambda[ei[ii]], phi+3*i);
	}
	else if (index < 20)
	{
		index -= 12;
		idx = index / 2;
		face2vertices3d(idx, fi);
		ii = index % 2; // ii=0 or 1
		if (ii == 0)
		{
			cross_array(grd_lambda[fi[2]], grd_lambda[fi[1]], v0);
			for(i=0;i<3;i++)
				axy_array(3, v0[i], grd_lambda[fi[0]], phi+3*i);
			cross_array(grd_lambda[fi[0]], grd_lambda[fi[2]], v0);
			for(i=0;i<3;i++)
				axpy_array(3, v0[i], grd_lambda[fi[1]], phi+3*i);
			cross_array(grd_lambda[fi[0]], grd_lambda[fi[1]], v0);
			for(i=0;i<3;i++)
				axpy_array(3, 2*v0[i], grd_lambda[fi[2]], phi+3*i);
		}
		else
		{
			cross_array(grd_lambda[fi[1]], grd_lambda[fi[2]], v0);
			for(i=0;i<3;i++)
				axy_array(3, v0[i], grd_lambda[fi[0]], phi+3*i);
			cross_array(grd_lambda[fi[0]], grd_lambda[fi[1]], v0);
			for(i=0;i<3;i++)
				axpy_array(3, v0[i], grd_lambda[fi[2]], phi+3*i);
			cross_array(grd_lambda[fi[0]], grd_lambda[fi[2]], v0);
			for(i=0;i<3;i++)
				axpy_array(3, 2*v0[i], grd_lambda[fi[1]], phi+3*i);
		}
	}
	else // 20 <= index < 32
	{
		for(i=0;i<4;i++)
		{
			for(j=i;j<4;j++)
				axpbyz_array(3, lambda[i], grd_lambda[j], lambda[j], grd_lambda[i], gradLambdaij[i][j]);

			for(j=0;j<i;j++)
				copy_array(3, gradLambdaij[j][i], gradLambdaij[i][j]);
		}

		index -= 20;
		i = index / 4;
		i1 = (i+1) % 3; i2 = (i+2) % 3;
		for(j=0;j<4;j++)
		{
			gradlambdaXei[j][i] = 0;
			gradlambdaXei[j][i1] = grd_lambda[j][i2];
			gradlambdaXei[j][i2] = -grd_lambda[j][i1];
		}
		j = index % 4;
		// (grad lambda[j] cross e_i) otimes grad(lambda[0]lambda[1]lambda[2]lambda[3])
		axpbyz_array(3, lambda[2]*lambda[3], gradLambdaij[0][1], lambda[0]*lambda[1], gradLambdaij[2][3], val2);
		for(i=0;i<3;i++)
			axy_array(3, gradlambdaXei[j][i], val2, phi+3*i);
		// (grad lambda[0] cross e_i) otimes grad(lambda[j]lambda[1]lambda[2]lambda[3])
		axpbyz_array(3, lambda[2]*lambda[3], gradLambdaij[j][1], lambda[j]*lambda[1], gradLambdaij[2][3], val2);
		for(i=0;i<3;i++)
			axpy_array(3, gradlambdaXei[0][i], val2, phi+3*i);
		// (grad lambda[1] cross e_i) otimes grad(lambda[0]lambda[j]lambda[2]lambda[3])
		axpbyz_array(3, lambda[2]*lambda[3], gradLambdaij[0][j], lambda[0]*lambda[j], gradLambdaij[2][3], val2);
		for(i=0;i<3;i++)
			axpy_array(3, gradlambdaXei[1][i], val2, phi+3*i);
		// (grad lambda[2] cross e_i) otimes grad(lambda[0]lambda[1]lambda[j]lambda[3])
		axpbyz_array(3, lambda[j]*lambda[3], gradLambdaij[0][1], lambda[0]*lambda[1], gradLambdaij[j][3], val2);
		for(i=0;i<3;i++)
			axpy_array(3, gradlambdaXei[2][i], val2, phi+3*i);
		// (grad lambda[3] cross e_i) otimes grad(lambda[0]lambda[1]lambda[2]lambda[j])
		axpbyz_array(3, lambda[2]*lambda[j], gradLambdaij[0][1], lambda[0]*lambda[1], gradLambdaij[2][j], val2);
		for(i=0;i<3;i++)
			axpy_array(3, gradlambdaXei[3][i], val2, phi+3*i);
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
 * \fn double volume(double **tet)
 * \brief get volume for tetrahedron p1(x1,y1,z1),p2(x2,y2,z2),p3(x3,y3,z3),p4(x4,y4,z4)
 * volume=det([1 x1 y1 z1;
               1 x2 y2 z2;
			   1 x3 y3 z3;
			   1 x4 y4 z4])
 * \param (*tet)[3] the axis value of the four vertices
 * \return volume of the tetrahedron
 */
double volume(double **tet)
{
	double a11=tet[1][0]-tet[0][0], a12=tet[1][1]-tet[0][1], a13=tet[1][2]-tet[0][2];
	double a21=tet[2][0]-tet[0][0], a22=tet[2][1]-tet[0][1], a23=tet[2][2]-tet[0][2];
	double a31=tet[3][0]-tet[0][0], a32=tet[3][1]-tet[0][1], a33=tet[3][2]-tet[0][2];

	return (a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31)/6;
//	return abs(a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31)/6;
}
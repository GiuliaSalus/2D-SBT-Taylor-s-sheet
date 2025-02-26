#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "hsize.h"
#include "nrutil.h"
#include "int_weights.h"
#include "grid.h"
#include "printfiles.h"

//calculate free space stokeslet except singular part
void stokesletNoSingFree(double g[],double x, double y,double x0,double y0,double sing,double par[])
{
	double dx=x-x0;
 	double dy=y-y0;
	double x2=pow(dx,2);
	double y2=pow(dy,2);
	double r2=pow(dx,2)+pow(dy,2);
	double r=sqrt(r2);
	double lnr;
	double gyx,gxx,gxy,gyy;
	
	
	lnr=log(fabs(r));

	if(sing==1)
	{
		gxx=x2/(r2);
		gxy=dx*dy/(r2);
		gyx=gxy;
		gyy=y2/(r2);
	}
	else
	{
		gxx=-lnr+x2/(r2);
		gxy=dx*dy/(r2);
		gyx=gxy;
		gyy=-lnr+y2/(r2);
	}
	g[0]=gxx; 
	g[1]=gxy; 
	g[2]=gyx;  
	g[3]=gyy;
}

//calculate periodic stokeslet except singular part
void stokesletNoSingPeriodic(double g[],double x, double y,double x0,double y0,double sing,double par[])
{
	
	double wn = 1.0;
	x = x*wn;//   ! field point
  y = y*wn;
  x0 = x0*wn;//   ! singular point
  y0 = y0*wn;
	double dx=x-x0;
 	double dy=y-y0;
	double x2=pow(dx,2);
	double y2=pow(dy,2);
	double xy=dx*dy;
	double r2=x2+y2;
	
	double gyx,gxx,gxy,gyy;
	double Sxx,Sxy,Syy;
	double A,Ax,Ay;
	
	
	
	double s=0.5*log(r2);
	double chy = cosh(dy);
	double shy = sinh(dy);
	double cx  = cos(dx);
	double sx  = sin(dx);
	double d =chy-cx;
   
	A=0.5*log(2.0*d);
	Ax=0.5*sx/d;
	Ay=0.5*shy/d;
	
	Sxx=-A-dy*Ay+1.0;
	Sxy=dy*Ax;
	Syy=-A+dy*Ay;
	
	if(sing==1)
	{
	  gxx=Sxx+s;
	  gxy=Sxy;
		gyx=gxy;
		gyy=Syy+s;
	}
	else
	{
	gxx=Sxx;
	gxy=Sxy;
	gyx=gxy;
	gyy=Syy;
	}
	g[0]=gxx; 
	g[1]=gxy; 
	g[2]=gyx;  
	g[3]=gyy;
}

//calculate periodic stress tensor except singular part
void stressNoSingPeriodic(double t[],double x, double y,double x0,double y0,double sing,double par[])
{
	double wn = 1.0;
	x = x*wn;//   ! field point
  y = y*wn;
  x0 = x0*wn;//   ! singular point
  y0 = y0*wn;
	double dx=x-x0;
 	double dy=y-y0;
	double x2=pow(dx,2);
	double y2=pow(dy,2);
	
	double r2=x2+y2;
	double cf;
	double A,Ax,Ay,Axx,Ayy,Axy;

	double s=0.5*log(r2);
	double chy = cosh(dy);
  double shy = sinh(dy);
  double cx  = cos(dx);
  double sx  = sin(dx);
  double d =chy-cx;
  double d2=d*d;
  double Txxx,Txxy,Tyyx,Tyyy;
	A=0.5*log(2.0*d);
	Ax=0.5*sx/d;
	Ay=0.5*shy/d;
	
	
  Ayy = 0.5*(1.0-cx*chy)/d2;
  Axy =-0.5*sx*shy /d2;
	Axx = -Ayy;
	
	
  Txxx = -2.0*(2.0*Ax+dy*Axy);
  Txxy = -2.0*(Ay+dy*Ayy);
  Tyyx = 2.0*dy*Axy;
  Tyyy = -2.0*(Ay-dy*Ayy);
  
	if(sing==1)
	{
		cf  = 4.0/(r2*r2);
 		Txxx = Txxx+ cf*dx*dx*dx;
 	 	Txxy = Txxy+ cf*dx*dx*dy;
		Tyyx = Tyyx + cf*dx*dy*dy;
		Tyyy = Tyyy + cf*dy*dy*dy; //sing part of analytical part is zero since int 1/(x-x0) dx is zero (as long as x0 in centre of int). 
 
	}


	t[0]=Txxx*wn; 
	t[1]=Txxy*wn; 
	t[2]=Tyyx*wn;  
	t[3]=Tyyy*wn;
}



//Calculate singular part analytically
double singval(double dr,double m)
{
	double c;
	double mm=m*m;

	c=-sqrt(1.0+mm)*dr*(2.0*log(dr)-2.0+log(1.0+mm));
	return c;
} 


//create matrix 
void append_fxComplete(double **A,double d[],double par[],int N,double ds[])
{
	int i;
  // fill the last column (2Ne+1) 
  for(i=0;i<N;i++) 
  {
    A[i][2*N]=d[i];
    A[i+N][2*N]=0.0;
  }
  // fill the last row (2Ne+1) with Ne 1.0 elements and Ne+1 0.0 elements
  for(i=0;i<N;i++)
  {
  	A[2*N][i]=ds[i];
  }
}

//integrate and assemble SBTmatrix and calculate velocities
void SBTmatrix(double par[],double vec[],double nodes[],double **A,double B[])
{
  
	int i,j,k;
	int N=par[1];
	int NQ=par[3];
	double b=par[2];
	
	//submatrix initialisation
	double d[2*N];
	double **axx=NULL,**axy=NULL,**ayx=NULL,**ayy=NULL;
	double **wxx=NULL,**wxy=NULL,**wyx=NULL,**wyy=NULL;
	double **one=NULL;
	axx= dmatrix( 0, N, 0,N);
	axy= dmatrix( 0, N, 0,N);
	ayx= dmatrix( 0, N, 0,N);
	ayy= dmatrix( 0, N, 0,N);
	wxx= dmatrix( 0, N, 0,N);
	wxy= dmatrix( 0, N, 0,N);
	wyx= dmatrix( 0, N, 0,N);
	wyy= dmatrix( 0, N, 0,N);
	one= dmatrix( 0, N, 0,N);
	for(i=0;i<N;i++)
	{
		d[i]=0.0;
		d[i+N]=0.0;
		for(j=0;j<N;j++)
		{
		axx[i][j]=0.0;
		ayx[i][j]=0.0;
		axy[i][j]=0.0;
		ayy[i][j]=0.0;
		wxx[i][j]=0.0;
		wyx[i][j]=0.0;
		wxy[i][j]=0.0;
		wyy[i][j]=0.0;
		if(i==j) one[i][j]=1.0;
		else one[i][j]=0.0;
		}  
	}
 	double x0[N],y0[N],dx[N],dy[N],m[N],h[N],dr[N],wx[N],uy[N],ux[N]; //elements details
	double xN[N*NQ],yN[N*NQ],vnx[N*NQ],vny[N*NQ],hN[N*NQ]; //nodes details
	double xx,yy,hh,vx,vy;
	

	double sing; //switch for singular integration
	double W[NQ],ti[NQ];
	weights(ti,W,NQ);
	double g[4],t[4],Q;
	double cc=-0.5/PI;
	double inext=par[7];
	vectorOut(vec,x0,y0,dx,dy,m,h,uy,N);
	nodesOut(nodes,xN,yN,vnx,vny,hN,NQ,N);
	if(inext==1.0)  getvel(par,m,ux,uy,N);
	else
	{
	  for(i=0;i<N;i++) 
	  {
	    ux[i]=0.0;
	  }

	}
	
	// integrate to obtain submatrices 
	for(k=0;k<N;k++)
	{
	  
		for(i=0;i<N;i++) // loop over elements
		{
		  if(k==0) wx[i]=0.0;
		  
			sing=0.0;
			if(i==k) sing=1.0;
			
			//integration block over the element
			for(j=0;j<NQ;j++) 
			{	
				xx=xN[j+i*NQ];
				yy=yN[j+i*NQ];
	 			hh=hN[j+i*NQ];
	 			vx=vnx[j+i*NQ];
	 			vy=vny[j+i*NQ];
	 			
	 			if(par[9]==1.0) //periodic GF
	 			{
				stokesletNoSingPeriodic(g,xx,yy,x0[k],y0[k],sing,par);
				stressNoSingPeriodic(t,xx,yy,x0[k],y0[k],sing,par);
				}	
				else //free space GF
				stokesletNoSingFree(g,xx,yy,x0[k],y0[k],sing,par);
				
				
				dr[i]=0.5*sqrt(dx[i]*dx[i]);
				
				axx[k][i] = axx[k][i] + dr[i]*W[j]*cc*g[0]*hh; 
				axy[k][i] = axy[k][i] + dr[i]*W[j]*cc*g[1]*hh;
				ayx[k][i] = ayx[k][i] + dr[i]*W[j]*cc*g[2]*hh;		 
				ayy[k][i] = ayy[k][i] + dr[i]*W[j]*cc*g[3]*hh;

				wxx[k][i] = wxx[k][i] + (t[0]*vx+t[1]*vy)*(-cc)*dr[i]*W[j]*hh;
				wxy[k][i] = wxy[k][i] + (t[1]*vx+t[2]*vy)*(-cc)*dr[i]*W[j]*hh;
				wyx[k][i] = wyx[k][i] + (t[1]*vx+t[2]*vy)*(-cc)*dr[i]*W[j]*hh;
				wyy[k][i] = wyy[k][i] + (t[2]*vx+t[3]*vy)*(-cc)*dr[i]*W[j]*hh;
				
			}	
			
			//add singular part 
			if(i==k)
			{
				Q=singval(dr[i],m[i]);
				axx[i][i]=axx[i][i]+cc*Q; 
				ayy[i][i]=ayy[i][i]+cc*Q;
			}
		}
	}
	
	//Assemble A matrix
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
		A[i][j]    =axx[i][j]; 
		A[i][N+j]  =axy[i][j];
		A[i+N][j]  =ayx[i][j]; 
		A[i+N][N+j]=ayy[i][j];
		}
 	}
	
	//Asseble B vector
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
		  if(par[6]==1.0)
		  {
			B[i]=B[i]-wxy[i][j]*uy[j]; 
			B[i+N]=B[i+N]+(one[i][j]-wyy[i][j])*uy[j]; 
			d[i]=d[i]+(wxx[i][j]-one[i][j]);
			d[i+N]=d[i+N]+wyx[i][j]; 
		  }
		  else
		  {
		  B[i]=B[i]+(one[i][j])*ux[j]; 
			B[i+N]=B[i+N]+(one[i][j])*uy[j]; 
			d[i]=d[i]+(-one[i][j]); 
			d[i+N]=0.0; 
		  }
		}
	}
	B[2*N]=0.0;
	
	for(i=0;i<N;i++)
	{
	  wx[i]=0.0;
	  	for(j=0;j<NQ;j++) 
			{	
			  hh=hN[j+i*NQ];
			  dr[i]=0.5*sqrt(dx[i]*dx[i]);
			  wx[i]=dr[i]*W[j]*cc*hh+wx[i]; 
			}
	}
 	append_fxComplete(A,d,par,N,wx);
 	
  
  free_dmatrix(axx, 0, N, 0,N);
  free_dmatrix(axy, 0, N, 0,N);
  free_dmatrix(ayx, 0, N, 0,N);
  free_dmatrix(ayy, 0, N, 0,N);
  free_dmatrix(wxx, 0, N, 0,N);
  free_dmatrix(wxy, 0, N, 0,N);
  free_dmatrix(wyx, 0, N, 0,N);
  free_dmatrix(wyy, 0, N, 0,N);
  free_dmatrix(one, 0, N, 0,N);
}



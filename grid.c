#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "hsize.h"
#include "systemsolver.h"
#include "printfiles.h"
#include "elliptic.h"
#include "nrutil.h"
#include "int_weights.h"

//store variables in vector
void vectorIn(double vec[],double x0[],double y0[],double dx[],double dy[],double m[],double r[],double uy[],int N)
{
	for (int i=0;i<N;i++)
	{
		vec[i]=x0[i];
		vec[i+N]=y0[i];
		vec[i+2*N]=dx[i];
		vec[i+3*N]=dy[i];
		vec[i+4*N]=m[i];
		vec[i+5*N]=r[i];
		vec[i+6*N]=uy[i];
	}
}

//extract variables from vector
void vectorOut(double vec[],double x0[],double y0[],double dx[],double dy[],double m[],double r[],double uy[],int N)
{
	for (int i=0;i<N;i++)
	{
		x0[i]=vec[i];
		y0[i]=vec[i+N];
		dx[i]=vec[i+2*N];
		dy[i]=vec[i+3*N];
		m[i]=vec[i+4*N];
		r[i]=vec[i+5*N];
		uy[i]=vec[i+6*N];
	}
}

//store nodes variables in vector
void nodesIn(double nodes[],double xN[],double yN[],double ysN[],double ynN[],double rN[],int NQ,int N)
{
	for (int i=0;i<N*NQ;i++)
	{
		nodes[i]=xN[i];
		nodes[i+N*NQ]=yN[i];
		nodes[i+2*N*NQ]=ysN[i];
		nodes[i+3*N*NQ]=ynN[i];
		nodes[i+4*N*NQ]=rN[i];
		
	}
}

//extract nodes variables from vector
void nodesOut(double nodes[],double xN[],double yN[],double ysN[],double ynN[],double rN[],int NQ,int N)
{
	for (int i=0;i<N*NQ;i++)
	{
		xN[i]=nodes[i];
		yN[i]=nodes[i+N*NQ];
		ysN[i]=nodes[i+2*N*NQ];
		ynN[i]=nodes[i+3*N*NQ];
		rN[i]=nodes[i+4*N*NQ];
		
	}
}


double Qinextensible(double b)
{
  //returns vaule Q=1/2pi(int_[0..2Pi] (1+b^2cos(z)) dz for inextebile condition
  double a1=b*b+1.0;
  double a2=sqrt(a1);
  double aa=b/a2;
  double Fk,Ek,Q;
  double c=2.0/(PI);
  Complete_Elliptic_Integrals('k',aa,&Fk,&Ek);

  Q=c*a2*Ek;
  
  return Q;
}

//get material point velocity for inextensible sheet
void getvel(double par[],double m[],double ux[],double uy[],int N)
{
  double b=par[2];
  double Q;
  double tt[N+1],ct[N+1],st[N+1];
  int k;
	Q=Qinextensible( b);
	for( k=0;k<N;k++)
	{
	  tt[k]=atan(m[k]);
	  ct[k]=cos(tt[k]);
	  st[k]=sin(tt[k]);
	  ux[k]=1.0-Q*ct[k];
	  uy[k]=-Q*st[k];
	}
}


//get grid information for straight elements
void getGridStraightEl(double par[],double x[],double y[],double vec[],double nodes[])
{
	int N=par[1];
	double b=par[2];
	int NQ=par[3];
	double t=par[8];
	double x0[N],y0[N],dx[N],dy[N],m[N],r[N],uy[N],dr[N],uye[N+1],ux[N+1],tt[N+1],ct[N+1],st[N+1];
	double xN[N*NQ],yN[N*NQ],ysN[N*NQ],ynN[N*NQ],rN[N*NQ],vnx[N*NQ],vny[N*NQ];
	double ti[NQ],W[NQ];
	double inext=par[7];
	double Q;
	Q=Qinextensible(b);

	for(int k=0;k<N+1;k++)
	{
		uye[k]=-b*cos(x[k]-t);
	}
	for(int k=0;k<N;k++)
	{
		x0[k]=(x[k+1]+x[k])*0.5;  
		dx[k]=x[k+1]-x[k];
		dy[k]=y[k+1]-y[k]; 
		m[k]=dy[k]/dx[k];
		r[k]=sqrt(1.0+m[k]*m[k]);
		y0[k]=(y[k+1]+y[k])*0.5;
		if(inext==0.0)//extensible sheet
		{
		  uy[k]=(uye[k+1]+uye[k])*0.5;
		  ux[k]=0.0;
		}
		else        	//inextensible sheet
		{
			tt[k]=atan(m[k]);
	  	ct[k]=cos(tt[k]);
	  	st[k]=sin(tt[k]);
	  	ux[k]=1.0-Q*ct[k];
	  	uy[k]=-Q*st[k];
		}
		dr[k]=sqrt(dx[k]*dx[k]+dy[k]*dy[k]);
	}
	
	
	weights(ti,W,NQ);
	for(int i=0;i<N;i++)
	{   
		for(int j=0;j<NQ;j++)
		{
    	xN[j+i*NQ]=0.5*(x[i+1]+x[i])+0.5*(x[i+1]-x[i])*ti[j];
			yN[j+i*NQ]=m[i]*(xN[j+i*NQ]-x[i])+y[i];
			rN[j+i*NQ]=sqrt(1.0+m[i]*m[i]);
			vnx[j+i*NQ]=-m[i]/rN[j+i*NQ];
			vny[j+i*NQ]=1.0/rN[j+i*NQ];
		
   	}
 	}
	vectorIn(vec,x0,y0,dx,dy,m,r,uy,N);
	nodesIn(nodes,xN,yN,vnx,vny,rN,NQ,N);
  	
}


// calculate area below a curve
void areacurve(double x1,double dx,double f1, double f2, double f3,double f4,double *F1,double *F2,double *F3)
{
  //interpolates by interpolation polynomial area under curve f1, f2, f3, f4 where x1, x2=x1+dx, x3=x2+dx, x4=x3+dx
  //F1=int(x1,x2), F2=int(x2,x3), F3=int(x3,x2);
	*F2 = -dx*(f1-13.0*f2-13.0*f3+f4)/24.0;                         
  
	*F1 = dx*(9.0*f1+19.0*f2-5.0*f3+f4)/24.0;

	*F3 = dx*(f1-5.0*f2+19.0*f3+9.0*f4)/24.0;

}

//calculate average value of quantity
double average(double ux[],double t[],int Nel,double dt,double dd)
{
	int i;
	double F1,F2,F3;
	double Ua=0.0;
	for(i=0;i<Nel-1;i++)
	{
  	if(i==0)
  	{
    	areacurve(t[0],dt,ux[0],ux[1],ux[2],ux[3],&F1,&F2,&F3);
			Ua=Ua+F1;
    }
    else if(i==Nel-2)
    {
    	areacurve(t[Nel-4],dt,ux[Nel-4],ux[Nel-3],ux[Nel-2],ux[Nel-1],&F1,&F2,&F3);
     	Ua=Ua+F3;
    }
    else
    {
    	areacurve(t[i-1],dt,ux[i-1],ux[i],ux[i+1],ux[i+2],&F1,&F2,&F3);
     	Ua=Ua+F2;
    }
  }
   
   Ua=Ua/dd;
   return Ua;
}





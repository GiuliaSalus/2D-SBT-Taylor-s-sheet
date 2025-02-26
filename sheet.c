#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hsize.h"
#include "nrutil.h"
#include "grid.h"
#include "printfiles.h"
#include "sbtmodel.h"
#include "systemsolver.h"


/*initialise sine shape*/
void InitSin(double par[], double x[], double y[])
{


int i;
int N=par[1];
double b=par[2];
double p=par[4];
double t=par[8];

double Lx;
double dl;
if(par[9]==1.0) p=1.0;

Lx=2.0*PI*p;
dl=Lx/N;

for(i=0;i<=N;i++)
	{	
		x[i]=-p*PI+i*dl;
		y[i]=b*sin(x[i]-t);
	}
	
	
}



/*calculate normalised velocity for a fixed configuration (fixed timestep)*/
void InstantVelocity(double par[],double x[],double y[])
{
	int i,j;
	int N=par[1];
	int NQ=par[3];
	int Nout=2*N+1; 		//number of unknown is 2*N (forces x,y) + 1 (Ux)
	double b=par[2];
	double vec[N*7], nodes[5*N*NQ];
	double x0[N],y0[N],dx[N],dy[N],m[N],h[N];
	double uy[N];
	
	double B[Nout],X[Nout];
	double fx[N],fy[N];
	double U=par[5];

	double **A=NULL; 
	A= dmatrix( 0, Nout, 0,Nout);
	for(i=0;i<Nout;i++)
	{
		B[i]=0.0;
		for(j=0;j<Nout;j++)
		{
		A[i][j]=0.0;
		}
	}
	
	//integrate SBT matrix, append forces matrix, get velocity
	getGridStraightEl(par,x,y,vec,nodes);
	
	SBTmatrix(par,vec,nodes,A,B);  	//with k kernel
	
	
	//solve system
	solve_system(A,B,X,Nout);
	
	
	for(i=0;i<N;i++)
	{
	 fx[i]=X[i];
	 fy[i]=X[i+N];
	}
	vectorOut(vec,x0,y0,dx,dy,m,h,uy,N);
	//writevector("traction",x0,y0,fx,fy,N);  //print traction
	
	U=X[Nout-1];
	par[5]=U/(b*b);
	free_dmatrix(A, 0, Nout, 0,Nout);
	
}

/*Um vs t: obtain normalized velocity at different times*/
void VelVsTime(double par[],double tn[], double un[], double dt, int Nt)
{
	int N=par[1]; 	
	double x[N+1],y[N+1],xf[N+1]; 
	double t;
	double U=par[5];
	
	

	for(int i=0;i<Nt;i++)
	{
		t=0.0+dt*i;
		par[8]=t;
		
		//get shape at time t
		InitSin(par,x,y);	
		for(int i=0;i<=N;i++)
		{
			xf[i]=U*t+x[i];
		}
		printf("t= %f\n",t);
		writeshape("free-b0_5",xf,y,N+1);  //print x vs y at different times
				
 		InstantVelocity(par,x,y);
		
		U=par[5];
		
		tn[i]=t;
		un[i]=U;
		
		
		
		
	}
	//writedata("U vs t",tn,un,Nt);  //print Um vs t
}



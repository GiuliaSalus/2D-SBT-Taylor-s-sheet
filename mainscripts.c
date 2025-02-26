#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hsize.h"
#include "nrutil.h"
#include "grid.h"
#include "printfiles.h"
#include "sbtmodel.h"
#include "systemsolver.h"
#include "sheet.h"

/*collection of routines*/

	
/*calculate normalised velocity U at a fixed time step (fixed configuration)*/
void SheetVelocity(double par[])
{
	int N=par[1];
	double x[N+1],y[N+1]; //points
	double U;

	
	//get shape at time t
	InitSin(par,x,y);	
 		
 	//get velocity
	InstantVelocity(par,x,y);
		
	U=par[5];
	printf("U=%e\n",-U); 

		

}

/*check convergence of normalised velocity U at a timestep t for increasing number of elements*/
void ConvergeN(double par[])
{
	for(int i=1;i<=50;i++)
	{
		int N=10*i; 	//number of elements
		
		par[1]=N;
		double x[N+1],y[N+1]; //points
		
		//get shape at time t
		InitSin(par,x,y);	
 		InstantVelocity(par,x,y);
		double U=par[5];
		
		writevalue("U vs N",U,N);
		
	}
}



/*get normalised average (in time) velocity Um*/
double AverageVel(double par[])
{
	int Nt=51;
	double un[Nt];
	double tn[Nt];
	double dt=2*PI/(Nt-1);
	double delta=2*PI;
	double Um;
	
	VelVsTime(par,tn,un,dt,Nt);
	Um=average(un,tn,Nt,dt,delta);
	printf("Um=%f\n",Um);

	return Um;
}



/*get average normalised velocity Um vs wave amplitude b*/
void VelVsB(double par[])
{
	int Nt=20;
	int Nb=100;
	double un[Nt];
	double tn[Nt];
	double delta=2*PI;
	double dt=delta/(Nt-1);
	double b;

	double Um;
	
	for(int i=0;i<Nb;i++)
	{	
	
		b=pow(10.0,-5.+5.*i/(Nb-1)); //100 log spaced points between 10^⁻5 and 1 
		par[2]=b;
		VelVsTime(par,tn,un,dt,Nt);
		Um=average(un,tn,Nt,dt,delta);
		writevalue("Um vs b",Um,b); //print Um vs b	
	}
}

//calculate normalised average velocity Um vs p
void VelVsP(double par[])
{
	int Nt=20;
	int Np=20;
	double dN=1.0;
	
	double un[Nt];
	double tn[Nt];
	double dt=2*PI/(Nt-1);

	double Um;
	double p;
	double delta=2*PI;//time period
	int N=par[1];
	
	for(int j=1;j<=Np;j++)
	{	
		p=1.0+dN*(j-1);
		N=p*20;
		
		par[4]=p;
		par[1]=N;
	
		
		VelVsTime(par,tn,un,dt,Nt);
		Um=average(un,tn,Nt,dt,delta);

		
		writevalue("Um vs p",Um,p);
	}
}











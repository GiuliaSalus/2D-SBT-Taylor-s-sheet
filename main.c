#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hsize.h"
#include "printfiles.h"
#include "mainscripts.h"


int main()
{
//set parameters
double periods=2.0; //number of spatial periods
int N=50*periods;  //number of discretization elements
int NQ=20;					//integration points per elements (1, 2, 3, 4, 5, 6, 8, 12, 20)
double b=0.5;  		//wave amplitude, b<<1
double U=0.0;				//swimming speed initialisation
double t=0.0; 			//time step initialisation
double kkernel=0.0; //0.0 = no k kernel, 1.0 = k kernel
double inexten=1.0; //1.0 = inextensible sheet; 0.0 = extensible sheet;
double greenf=0.0;  //0.0 = free space gf; 1.0 = periodic gf

double par[50];     //array to pass all parameters

par[1]=N;
par[2]=b;
par[3]=NQ;
par[4]=periods;
par[5]=U;
par[6]=kkernel;
par[7]=inexten;
par[8]=t;
par[9]=greenf;



//SheetVelocity(par); 	//calculate normalised velocity at a fixed timpestep (fixed configuration) 
//ConvergeN(par); 			//calculate normalised velocity at a fixed timpestep (fixed configuration) vs N 
U=AverageVel(par); 		//calculate normalised average (in time) velocity 
//VelVsB(par); 					//calculate normalised average velocity vs wave amplitude b 
//VelVsP(par);					//calculate normalised average velocity vs spatial period number p



return 0;
}
















  

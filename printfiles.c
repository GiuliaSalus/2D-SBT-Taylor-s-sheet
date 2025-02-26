#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hsize.h"

/*write data file with vector vx[NN],vy[NN] and application point x[NN],y[NN], in the format x y vx vy*/
void writevector(const char* filename, double x[], double y[], double ux[],double uy[], int NN)
{
  int i;  
  FILE *fp=NULL;   /* file pointer */
   
	/* open the file */
  fp = fopen(filename,"w+") ;
  
  // write to the file: 
	for(i=0;i<NN;i++)
	{
  fprintf(fp,"%.8e   %.8e   %.8e   %.8e\n",x[i],y[i],ux[i],uy[i]);
 	}
 	
  // close file
  fclose(fp);
}


/*write data file with data x[NN], y[NN]*/
void writedata(const char* filename,double x[],double y[],int NN)
{
  int i;
  
  FILE *fp=NULL;   /* file pointer */
   
	/* open the file */
  fp = fopen(filename,"w+") ;
  
  
  // write to the file: 
	for(i=0;i<NN;i++)
	{
  fprintf(fp,"%.8e  %.10e   \n",x[i],y[i]);
 	}
 	
  // close file
  fclose(fp);
}

/*write data file with matrix*/
void writematrix(const char* filename,double **A, int NN)
{
  
  int i,j;
  
  
  FILE *fp=NULL;   /* file pointer */
   
	/* open the file */
  fp = fopen(filename,"w+") ;
  
  
  // write to the file: 
  fprintf(fp,"{");
	for(i=0;i<NN;i++)
	{
		fprintf(fp,"{");
		for(j=0;j<NN-1;j++)
		{
			fprintf(fp,"%.8f, ",A[i][j]);
		}
		fprintf(fp,"%.8f},\n",A[i][NN-1]);
 	}
  fprintf(fp,"}");
  // close file
  fclose(fp);
}


/*write data file with x value and step N at every step N*/
void writevalue(const char* filename,double x, double N)
{ 
  FILE *fp=NULL;   /* file pointer */
   
	/* open the file */
  fp = fopen(filename,"a+") ;
  
  
  // write to the file: 
	fprintf(fp,"%8e %.8e \n",N,x);
		
  // close file
  fclose(fp);
}

/*write data file with array x[NN]*/
void writearray(const char* filename,double x[],int NN)
{
  int i;
  FILE *fp=NULL;   /* file pointer */
   
	/* open the file */
  fp = fopen(filename,"w+") ;
 
  // write to the file: 
	for(i=0;i<NN;i++)
	{
  fprintf(fp,"%.8e,  \n",x[i]);
 	}
 	
  // close file
  fclose(fp);
}

/*write data file with x[NN], y[NN] at different times*/
void writeshape(const char* filename,double x[], double y[],int NN)
{
  int i;
  FILE *fp=NULL;   /* file pointer */
   
	/* open the file */
  fp = fopen(filename,"a+") ;
 
  // write to the file: 
	for(i=0;i<NN;i++)
	{
  fprintf(fp,"%.8e  %.8e   \n",x[i],y[i]);
 	}
 	
 	fprintf(fp,"\n");
 	
  // close file
  fclose(fp);
}













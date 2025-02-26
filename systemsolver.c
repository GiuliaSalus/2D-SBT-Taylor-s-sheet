#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hsize.h"
#include "nrutil.h"
#include "printfiles.h"
#include "ludcmp.h"
#include "lubksb.h"

void swap_rows(double **A, double *b, int i, int j, int size) 
{
		int k;
    double temp;
    for (k = 0; k < size; k++) 
    {    
        temp = A[i][k];
        A[i][k] = A[j][k];
        A[j][k] = temp;
    }
    temp = b[i];
    b[i] = b[j];
    b[j] = temp;
}


void solve_system(double **A, double *b, double x[], int size)
{

		int i,j,k;
		int max_row;
		double factor;
		
    for (i = 0; i < size; i++) 
    {
        // Partial Pivoting
        max_row = i;
        for (k = i + 1; k < size; k++) 
        {
            if (fabs(A[k][i]) > fabs(A[max_row][i])) 
            {
                max_row = k;
            }
        }
        if (max_row != i) 
        {
            swap_rows(A, b, i, max_row, size);
        }

        // Elimination
        for (j = i + 1; j < size; j++) 
        {
            factor = A[j][i] / A[i][i];
            for (k = i; k < size; k++) 
            {
                A[j][k] -= factor * A[i][k];
            }
            b[j] -= factor * b[i];
        }
    }

    // Back substitution
    
    for (i = size-1; i >= 0; i--) 
    {
        x[i] = b[i];
        for (j = i + 1; j < size; j++) 
        {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];

    }
}


void inverse(double **F, double *d, double **y, int row)
{
  int *indx =NULL;
  double *col ;
  int i,j ;

  col=dvector(1,row);
  //solves the linear ewuation Ax=b , returns x
	indx=ivector(1,row);
  ludcmp(F,row,indx,d) ;
  
  for(j=1;j<=row;j++)
  {
	  for(i=1;i<=row;i++)
	  {
      col[i]=0.0;
		}
    col[j]=1.0;
    lubksb(F,row,indx,col);
    for(i=1;i<=row;i++)
    {
      y[i][j]=col[i];
    }
  }
free_dvector(col,1,row);
free_ivector(indx,1,row);

}

void solveqt3(double **inv,double b[],double uu[], int row)
{
  double  sum;
  int i,j ;
  //solve from inverse matrix
	sum=0.0;    
	for(i=1;i<=row;i++)
	{
  	uu[i]=0.0 ;
    sum=0.0;
  	for(j=1;j<=row;j++)
  	{
  	  sum = (inv[i][j])*(b[j]);
      sum = (sum) +(uu[i]) ;
      uu[i]=(sum) ;
    }
  }

}

void solve_systemNR(double **A, double B[], double X[],int Nout)
{
	//change indeces in vectors and matrix to use NR solver
	double d;
	double **Aaux=NULL,**Ainv=NULL;
	double Baux[Nout+1],sol[Nout+1];
	Aaux= dmatrix( 1, Nout, 1, Nout);
  Ainv= dmatrix( 1, Nout, 1, Nout);
  
  for(int i=1;i<=Nout;i++)
	{
		Baux[i]=0.0;
		for(int j=1;j<=Nout;j++)
		{
			Aaux[i][j]=0.0;
			Ainv[i][j]=0.0;
		}
	}
  
  
  
  for(int i=1;i<=Nout;i++)
  {
  	Baux[i]=B[i-1];
  	for(int j=1;j<=Nout;j++)
  	{
  		Aaux[i][j]=A[i-1][j-1];
  	}
  }
 
	inverse(Aaux, &d, Ainv, Nout);
 	solveqt3(Ainv,Baux,sol,Nout);
 	
 	
	for(int i=1;i<=Nout;i++)
	{
		X[i-1]=sol[i];
	}
}




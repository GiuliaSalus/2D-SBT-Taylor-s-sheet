void vectorIn(double vec[],double x0[],double y0[],double dx[],double ys[],double yn[],double r[],double uy[],int N);
void vectorOut(double vec[],double x0[],double y0[],double dx[],double ys[],double yn[],double r[],double uy[],int N);

void nodesIn(double nodes[],double xN[],double yN[],double ysN[],double ynN[],double rN[],int NQ,int N);
void nodesOut(double nodes[],double xN[],double yN[],double ysN[],double ynN[],double rN[],int NQ,int N);

void getGridStraightEl(double par[],double x[],double y[],double vec[],double nodes[]);

void areacurve(double x1,double dx,double f1, double f2, double f3,double f4,double *F1,double *F2,double *F3);
double average(double ux[],double t[],int Nel,double dt, double dd);

double Qinextensible(double b);
void getvel(double par[],double m[],double ux[],double uy[],int N);

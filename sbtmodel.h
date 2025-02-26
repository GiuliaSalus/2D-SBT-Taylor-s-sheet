void stokesletNoSingFree(double g[],double x, double y,double x0,double y0,double sing,double par[]);
void stokesletNoSingPeriodic(double g[],double x, double y,double x0,double y0,double sing,double par[]);
void stressNoSingPeriodic(double t[],double x, double y,double x0,double y0,double sing,double par[]);

double singval(double dx,double m);

void SBTmatrix(double par[],double vec[],double nodes[],double **A,double B[]);

void append_fxComplete(double **A,double d[],double par[],int N,double ds[]);


//codes for complete elliptic integals of the first and second kind. Want char = 'm' for modulus= m=k^2
//i.e. E(m)=int((1-m*sinTt)^2)^(-1/2),0..Pi/2

////////////////////////////////////////////////////////////////////////////////
// File: complete_elliptic_integrals.c                                        //
// Routine(s):                                                                //
//    Complete_Elliptic_Integrals                                             //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void Complete_Elliptic_Integrals(char arg, double x, double* Fk,           //
//                                                              double* Ek)   //
//                                                                            //
//  Description:                                                              //
//     This function calculates the complete elliptic integrals of the first  //
//     and second kinds.                                                      //
//     The complete elliptic integral of the first kind is the integral from  //
//     0 to pi / 2 of the integrand                                           //
//                   dtheta / sqrt( 1 - k^2 sin^2(theta) ).                   //
//     The complete elliptic integral of the second kind is the integral from //
//     0 to pi / 2 of the integrand                                           //
//                    sqrt( 1 - k^2 sin^2(theta) ) dtheta .                   //
//     The parameter k is called the modulus.  This integral is even in k.    //
//     The modulus, k, must satisfy |k| <= 1.  If k = 0 then the integrals    //
//     can be readily evaluated.  If |k| = 1, then the integral for the       //
//     complete elliptic integral of the first kind is infinite.              //
//     Otherwise the integrals must be approximated.                          //
//                                                                            //
//     In practise the arguments the complete elliptic function of the first  //
//     kind are also given as F(pi/2 \ alpha) or F(pi/2 | m) and the complete //
//     elliptic function of the second kind are also given as E(pi/2 \ alpha) //
//     or E(pi/2 | m) where the angle alpha, called the modular angle,        //
//     satisfies k = sin(alpha) and the argument m = k^2 is simply called the //
//     parameter.                                                             //
//     In terms of these arguments K = F(pi/2 \ alpha) = F(pi/2, sin(alpha))  //
//     and K = F(pi/2 | m) = F(pi/2, sqrt(m)), where                          //
//             K = Complete_Elliptic_Integral_First_Kind( k )                 //
//     and E = E(pi/2 \ alpha) = E(pi/2, sin(alpha)) and                      //
//     E = E(pi/2 | m) = E(pi/2, sqrt(m)), where                              //
//             E = Complete_Elliptic_Integral_Second_Kind( k ).               //
//                                                                            //
//     The common mean method, sometimes called the Gauss transform, is a     //
//     variant of the descending Landen transformation in which two sequences //
//     are formed: Setting a[0] = 1 and g[0] = k', the complementary modulus, //
//     a[i] is the arithmetic average and g[i] is the geometric mean of a[i-1]//
//     and g[i-1], i.e. a[i+1] = (a[i] + g[i])/2 and g[i+1] = sqrt(a[i]*g[i]).//
//     The sequences satisfy the inequalities g[0] < g[1] < ... < a[1] < a[0].//
//     Further, lim g[n] = lim a[n] as n -> inf.                              //
//     The value of the complete elliptic integral of the first kind is       //
//     (pi/2) lim (1/G[n]) as n -> inf.                                       //
//     The value of the complete elliptic integral of the second kind is      //
//           E(k) = lim (pi/8g[n]) (4 - 2k^2 - Sum(2^j(a[j]-g[j])^2)).        //
//     where the limit is as n -> inf and the sum extends from j = 0 to n.    //
//     The sum of 2^j (a[j]^2 - g[j]^2) from j = 1 to n equals                //
//           (1/2) Sum (2^i (a[i] - g[i])^2 for i = 0,...,n-1, so that        //
//           E(k) = lim (pi/4g[n]) (2 - k^2 - Sum(2^j(a[j]^2 -g[j]^2))).      //
//                                                                            //
//  Arguments:                                                                //
//     char    arg                                                            //
//                The type of argument of the second argument of F() and E(): //
//                  If arg = 'k', then x = k, the modulus of F(pi/2,k) and    //
//                                E(pi/2,k).                                  //
//                  If arg = 'a', then x = alpha, the modular angle of        //
//                                F(pi/2 \ alpha) and E(pi/2 \ alpha), alpha  //
//                                in radians.                                 //
//                  If arg = 'm', then x = m, the parameter of F(pi/2 | m)    //
//                                and E(pi/2 | m).                            //
//                  The value of arg defaults to 'k'.                         //
//     double  x                                                              //
//                The second argument of the elliptic function F(pi/2,k),     //
//                F(pi/2 \ alpha) or F(pi/2 | m)  and E(pi/2,k),              //
//                E(pi/2 \ alpha) or E(pi/2 | m) corresponding to the value   //
//                of 'arg'.  Note that if arg = 'k', then | x | <= 1 and if   //
//                arg = 'm', then 0 <= x <= 1.                                //
//     double* Fk                                                             //
//                The complete elliptic integral of the first kind.           //
//                If |modulus| = 1, then the integral is infinite and DBL_MAX //
//                returned.                                                   //
//     double* Ek                                                             //
//                The complete elliptic integral of the second kind.          //
//                                                                            //
//  Return Value:                                                             //
//     void.  The values are returned through the argument list.              //
//                                                                            //
//  Example:                                                                  //
//     double m, k, a;                                                        //
//     double Fk, Ek;                                                         //
//                                                                            //
//     ( code to initialize a )                                               //
//                                                                            //
//     k = sin(a);                                                            //
//     m = k * k;                                                             //
//     Complete_Elliptic_Integrals( 'a', a, &Fk, &Ek );                       //
//     printf("K(alpha) = %12.6f where angle(radians) = %12.6f\n",Fk, a);     //
//     printf("E(alpha) = %12.6f where angle(radians) = %12.6f\n",Ek, a);     //
//     Complete_Elliptic_Integrals( 'k', k, &Fk, &Ek );                       //
//     printf("K(k) = %12.6f where k = %12.6f\n",Fk, k);                      //
//     printf("E(k) = %12.6f where k = %12.6f\n",Ek, k);                      //
//     Complete_Elliptic_Integrals( 'm', m, &Fk, &Ek );                       //
//     printf("K(m) = %12.6f where m = %12.6f\n",Fk, m);                      //
//     printf("E(m) = %12.6f where m = %12.6f\n",Ek, m);                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>        // required for fabs(), fabsl(), sqrtl() and M_PI_2
#include <float.h>       // required for LDBL_EPSILON, DBL_MAX
#include "hsize.h"

//static const long double PI_2 = PI*0.5;////1.5707963267948966192313216916397514L; // pi/2
//  static const long double PI_4 = PI*0.25;//0.7853981633974483096156608458198757L; // pi/4

void Complete_Elliptic_IntegralsA(double x, double* Fk, double* Ek){
double xx=fabs(x);
double k=sqrt(xx);
 double PI_2 = PI*0.5;

 double k2=k*k;
 double k4=k2*k2;
 double k6=k4*k2;
 *Ek=PI_2*(1-0.25*k2-0.04687500000*k4-0.01953125000*k6);
 *Fk=PI_2*(1+0.25*k2+0.1406250000*k4+0.09765625000*k6);
}
void Complete_Elliptic_Integrals(char arg, double x, double* Fk, double* Ek)
{
   long double k;      // modulus
   long double m;      // the parameter of the elliptic function m = modulus^2
   long double a;      // arithmetic mean
   long double g;      // geometric mean
   long double a_old;  // previous arithmetic mean
   long double g_old;  // previous geometric mean
   long double two_n;  // power of 2
   long double sum;
   long double PI_2 = PI*0.5;
    long double PI_4 = PI*0.25;
   if ( x == 0.0 ) {
      *Fk = M_PI_2;
      *Ek = M_PI_2;
      return;
   }

   switch (arg) {
      case 'k': k = fabsl((long double) x);
                m = k * k;
                break;
      case 'm': m = (long double) x;
                k = sqrtl(fabsl(m));
                break;
      case 'a': k = sinl((long double)x);
                m = k * k;
                break;
      default:  k = fabsl((long double) x);
                m = k * k;
   }

   if ( m == 1.0 ) {
      *Fk = DBL_MAX;
      *Ek = 1.0;
      return;
   }

   a = 1.0L;
   g = sqrtl(1.0L - m);
   two_n = 1.0L;
   sum = 2.0L - m;
   while (1) {
      g_old = g;
      a_old = a;
      a = 0.5L * (g_old + a_old);
      g = g_old * a_old;
      two_n += two_n;
      sum -= two_n * (a * a - g);
      if ( fabsl(a_old - g_old) <= (a_old * LDBL_EPSILON) ) break;
      g = sqrtl(g);
   } 
   *Fk = (double) (PI_2 / a);
   *Ek = (double) ((PI_4 / a) * sum);
   return;
}

/*
Amol Holkundkar
Department of Physics, BITS Pilani,
Pilani, Rajasthan - India.
e-mail: amol.h@live.com

The method is explained on following web-page: Courtesy Texas A & M University

http://www.math.tamu.edu/~mpilant/math614/Matlab/Lyapunov/LorenzSpectrum.pdf

*/

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
using namespace std;

 const int ATTRACTOR = 0 ; // 0 - Lorenz : 1 - Rossler 

 double sig,rho,bet;

 double N0 = 1e6; // Time after which Lyapunov exponent is evaluated
 double T = 1e4;  // Total simulatio time
 double dT = 1e-3;// Time Step
 

double SQ(double x)
{
 return x*x;
}

double fun_x1(double x1,double x2,double x3)
{
  if(ATTRACTOR == 0)
  {
   return sig*(x2 - x1);
  }
  else return -x2 -x3 ;
 
}
double fun_x2(double x1,double x2,double x3)
{
 if(ATTRACTOR == 0)
 {
  return rho*x1 - x1*x3 - x2;
 } 
 else return x1 + sig*x2 ; 
 
}
double fun_x3(double x1,double x2,double x3)
{
 if(ATTRACTOR == 0)
 {
  return x1*x2 - bet*x3;
 }
 else return rho + x3*(x1 - bet) ;
 
}
 

double rk4(double(*f)(double,double,double), double dt, double x1,double x2,double x3,int i)
{
	double	k1,k2,k3,k4;
	
	k1 = dt * f(x1,x2,x3);
	k2 = dt * f(x1 + 0.5*k1,x2 + 0.5*k1,x3 + 0.5*k1);
	k3 = dt * f(x1 + 0.5*k2,x2 + 0.5*k2,x3 + 0.5*k2);
	k4 = dt * f(x1 + k3,x2 + k3,x3 + k3);
	
	switch(i)
	{
	 case 0: return x1 + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;break;
	 case 1: return x2 + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;break;
	 case 2: return x3 + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;break;
        }
}    


const int dimension = 3;

struct nVector
{
 double r[dimension];
};

typedef vector <nVector> VnV;

class Utility
{
   public:
   
   nVector AddVector(nVector,nVector);
   nVector SubVector(nVector,nVector);
   nVector ScalarMul(nVector,double);
   nVector NormalizeVector(nVector) ;
   double ScalarProduct(nVector,nVector);
   double ModVector(nVector);
   VnV DoOrthogonalization(VnV);
   VnV MatrixMultiply(VnV,VnV) ;
   
   VnV Transpose(VnV x)
   {
    VnV xT;
    
    xT = x;
    
    for(int i =0;i< dimension;i++)
    {
     for(int j =0;j< dimension;j++)
     {    
      xT[j].r[i] = x[i].r[j] ;
     }
    }
     
     return xT;
   }
   
   void PrintVector(VnV x)
   {
    for(int i =0;i< dimension;i++)
    {
     for(int j =0;j< dimension;j++)
     {    
      cout << x[i].r[j]<< "    ";
     }
     cout << endl;
    }
   }	

};

int main()
{
 
 if(ATTRACTOR == 0)
 {
  cout << endl << "Lorenz Attractor : " ;
  sig = 10.0;
  rho = 28.0;
  bet = 8.0/3.0;
 }
 else
 { 
  cout << endl << "Rossler Attractor : " ;
  sig = 0.2;
  rho = 0.2;
  bet = 5.7;
 }

 
 cout << setprecision(4);
 
 Utility u;
   
 nVector dy,dummy; // dynamical variable
 nVector lyapunov; // lyapunov exponent
 nVector sumly;
 VnV v,v1,unit; // for lyapunov exponent
 VnV J; // jacabobian
 
 for(int i = 0;i< dimension;i++)
 dummy.r[i] = 0;
 
 v.assign(dimension,dummy);
 v1.assign(dimension,dummy);
 J.assign(dimension,dummy);
 
 double x1,x2,x3; 

 //initial condition..
 
 dy.r[0] = 1;dy.r[1] = 1;dy.r[2] = 1 ;
 sumly.r[0] = 0;sumly.r[1] = 0;sumly.r[2] = 0;
 
 // initial orthogonal vector
 v[0].r[0] = 1;v[0].r[1] = 0;v[0].r[2] = 0;
 v[1].r[0] = 0;v[1].r[1] = 1;v[1].r[2] = 0;
 v[2].r[0] = 0;v[2].r[1] = 0;v[2].r[2] = 1;
 
 unit = v;
      
 for(int n=1;n<=int(T/dT);n++)
 {
  
  x1 = rk4(fun_x1, dT, dy.r[0],dy.r[1],dy.r[2],0);
  x2 = rk4(fun_x2, dT, dy.r[0],dy.r[1],dy.r[2],1);
  x3 = rk4(fun_x3, dT, dy.r[0],dy.r[1],dy.r[2],2);
   
  dy.r[0] = x1;
  dy.r[1] = x2;
  dy.r[2] = x3;
  
  if(n > int(N0))
  {
  
	  if(ATTRACTOR == 0)
	  {
	   J[0].r[0] = - sig;J[0].r[1] = sig;J[0].r[2] = 0;
	   J[1].r[0] = rho - x3;J[1].r[1] = -1;J[1].r[2] = -x1;
	   J[2].r[0] = x2;J[2].r[1] = x1;J[2].r[2] = -bet;
	  }
	  else
	  {
	   J[0].r[0] = 0;J[0].r[1] = -1;J[0].r[2] = -1;
	   J[1].r[0] = 1;J[1].r[1] = sig;J[1].r[2] = 0;
	   J[2].r[0] = x3;J[2].r[1] = 0;J[2].r[2] = x1 - bet;
	  }
	  
	  J[0] = u.AddVector(unit[0],u.ScalarMul(J[0],dT) ) ;
	  J[1] = u.AddVector(unit[1],u.ScalarMul(J[1],dT) ) ;
	  J[2] = u.AddVector(unit[2],u.ScalarMul(J[2],dT) ) ;
	    
	  v1 = u.DoOrthogonalization(u.MatrixMultiply(J,u.Transpose(v)))	;
	  
	  for(int k = 0;k< dimension;k++)
	  {
	   lyapunov.r[k]  = log(u.ModVector(v1[k])) / dT;
	   v[k] = u.NormalizeVector(v1[k]);
	   sumly.r[k] += lyapunov.r[k] ;
	  }
	  
  }
    
 }
 
 cout << endl << "Time Averaged Lyapunov Exponents." << endl;
 
 cout << endl << "l1 = " << sumly.r[0]/( int(T/dT) - N0); 
 cout << endl << "l2 = " << sumly.r[1]/( int(T/dT) - N0); 
 cout << endl << "l3 = " << sumly.r[2]/( int(T/dT) - N0); 
 
 
 cout << endl ;
 return 0;
}

VnV Utility :: MatrixMultiply(VnV A,VnV B)
{
 nVector dummy ;
 
 for(int i = 0;i< dimension;i++)
 dummy.r[i] = 0;
 
 VnV C ;
 C.assign(dimension,dummy);
  
 for(int i = 0;i< dimension; i++)
 {
  for(int j = 0;j< dimension;j++)
  {
   for(int k = 0;k< dimension;k++)
   {
    C[i].r[j] += A[i].r[k] * B[k].r[j];    
   }
  }
 }
 
 return Transpose(C);
}

double Utility :: ModVector(nVector v1)
{
 double sum=0;
 
 for(int i = 0;i< dimension; i++)
 {
  sum += v1.r[i] * v1.r[i] ;
 }
 
 return sqrt(sum);
}

double Utility :: ScalarProduct(nVector v1,nVector v2)
{
 double sum = 0;
  
 for(int i =0;i< dimension;i++)
 {
   sum += v1.r[i] * v2.r[i];
 }
 return sum;
}

nVector Utility :: NormalizeVector(nVector v1)
{
 double sum = 0;
 double mul;
 
 nVector result;
 
 for(int i =0;i< dimension;i++)
 {
   sum += v1.r[i] * v1.r[i];
 }
 
 for(int i =0;i< dimension;i++)
 {
  result.r[i] = v1.r[i] / sqrt(sum);
 }
 
 return result;
}

nVector Utility :: AddVector(nVector v1,nVector v2)
{
 nVector result;
 
 for(int i = 0;i< dimension; i++)
 {
  result.r[i] = v1.r[i] + v2.r[i] ;
 }
 
 return result;
}

nVector Utility :: SubVector(nVector v1,nVector v2)
{
 nVector result;
 
 for(int i = 0;i< dimension; i++)
 {
  result.r[i] = v1.r[i] - v2.r[i] ;
 }
 
 return result;
}

nVector Utility :: ScalarMul(nVector v1,double x)
{
 nVector result;
 
 for(int i = 0;i< dimension; i++)
 {
  result.r[i] = v1.r[i] * x;
 }
 
 return result;
}

VnV Utility :: DoOrthogonalization(VnV in)
{
 VnV out,out_normalize;
 nVector sum,dummy;
 
 for(int i = 0;i< dimension;i++)
 dummy.r[i] = 0;
 
 out.assign(dimension,dummy);
 
 out[0] = in[0];
 
 for(int i = 1;i< in.size();i++)
 {
   sum = dummy;
   for(int j = 0;j<= i-1;j++)
   {
    sum = AddVector(sum,ScalarMul(out[j],ScalarProduct(out[j],in[i]) / SQ(ModVector(out[j]))  ) );
   }
  
   out[i] = SubVector(in[i],sum) ; 
   
   for(int k0=0;k0< dimension;k0++)
   {
    if(abs(out[i].r[k0]) < 1e-10)
    out[i].r[k0] = 0;
   }
 }
 
 return out;
}



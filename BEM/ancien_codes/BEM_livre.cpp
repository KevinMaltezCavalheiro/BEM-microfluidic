#include <gmsh.h>
#include <iostream>
#include <sstream>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Core> 
#include <Eigen/Sparse> 
#include <set>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <cstdio>
#include <vector>
//#include <QImage>
#include "test.h"
#include "BEM_livre.h"


using namespace std;
using namespace gmsh;
using namespace Eigen;

void bem(int argc, char **argv);

void gmatr(vector<double> &xl,vector<double> &yl,vector<double> &xm,vector<double> &ym,vector<double> &G,int N);

void rlintc(double x0,double y0,double x1,double y1,double x2,double y2,double &result);

void slintc(double x1, double y1, double x2, double y2, double &result);

void dalpha(double x0,double y0,double x1,double y1, double x2,double y2, double &result);

void hmatr(vector<double> &xl, vector<double> &yl, vector<double> &xm, vector<double> &ym, vector<double> &H, int N);

void abmatr(vector<double> &G,vector<double> &H, vector<double> &A, vector<double> &unb, vector<double> &ub, vector<int> &index, int N);

void solveq(vector<double> &A,vector<double> &unb, int N, int lsing);

void leqs(vector<double> &A,vector<double> &B,int N, int lsing);

void reorder(vector<double> &ub,vector<double> &unb,vector<int> &index,int N);

void uinter(vector<double> &xl,vector<double> &yl,vector<double> &xin,vector<double> &yin,vector<double> &ub,vector<double> &unb,vector<double> &uin,int N, double IN);

void input(vector<double> &xl,vector<double> &yl,vector<double> &xin,vector<double> &yin,vector<int> &index,vector<double> &ub,int &N,int &IN);

/*
main program: defines N and |N and calls the subroutines (N: Number of boundary elements and boundary nodes and |N: Number of internal points where the solution is computed)
input: reads the data and writes them in outputfile
gmatr: forms the matrix G
hmatr: forms the matrix H
abmatr: forms the matrices A and B
solveq: solves the system Ax=B
reorder: Rearranges the boundary values and forms the matrices UB and UNB
uinter: computes the solution at internal points
output: prints the results in outputfile
end
*/

void bem(int argc, char **argv)
{
  //N: Number of boundary elements and boundary nodes
  //IN: Number of internal points where the solution is computed
  //index: One-dimensional array in which a type of boundary condition is
  //assigned to the nodes. It is index[j] = 0 when u is prescribed,
  //while index[j] = 1 when @u=@n is prescribed
  //xl: One-dimensional array containing the x coordinates of the extreme
  //points of all the elements
  //yl: One-dimensional array containing the y coordinates of the extreme points of all the elements
  //xm: One-dimensional array containing the x coordinates of all the boundary nodes
  //ym: One-dimensional array containing the y coordinates of all the boundary nodes
  //G,H,A: square matrix explicitely described by the method
  //ub: One-dimensional array. At input, it contains the boundary values of
  //u, if INDEX(J) 5 0, or @u=@n, if INDEX(J) 5 1. At output it
  //contains all the boundary nodal values of u
  //unb: One-dimensional array containing the right-hand side vector of
  //equation AX = B given by Eq. (4.10c). At output it contains
  //the boundary nodal values of @u=@n
  //xin: One-dimensional array containing the x coordinates of the internal
  //points at which the values of u are computed
  //yin: One-dimensional array containing the y coordinates of the internal
  //points at which the values of u are computed
  //uin: One-dimensional array containing the computed values of u at the
  //internal points

  //initialisation
  int N = 16; //number of boundary elements xhich is equal to the number of boundary nodes
  int IN = 9; //number of internal points where the function u is calculated
  int lsing = 0;
  vector<int> index(N,0);
  vector<double> xl(N,0);
  vector<double> yl(N,0);
  vector<double> xm(N,0);
  vector<double> ym(N,0);
  vector<double> G(N*N,0);
  vector<double> H(N*N,0);
  vector<double> ub(N,0);
  vector<double> A(N*N,0);
  vector<double> unb(N,0);
  vector<double> xin(IN,0);
  vector<double> yin(IN,0);
  vector<double> uin(IN,0);
  /////////////////////////////////////////////////////////////////////////
  //todo input file with all the parameters for xl,yl,xin,yin,index,ub,N,IN
  /////////////////////////////////////////////////////////////////////////

  //read the input file
  //input(xl,yl,xin,yin,index,ub,N,IN);

  //compute the g matrix
  gmatr(xl,yl,xm,ym,G,N);

  //compute H matrix
  hmatr(xl,yl,xm,ym,H,N);

  //form the system of equations AX=B
  abmatr(G,H,A,unb,ub,index,N);

  //solve the system of equations
  solveq(A,unb,N,lsing);

  //form the vectors u and un of all the boundary values
  reorder(ub,unb,index,N);

  //compute the values uin of u at the internal points
  uinter(xl,yl,xin,yin,ub,unb,uin,N,IN);

  //print the results in the output file
  //output(xm,ym,ub,unb,xin,yin,uin,N,IN);

  //todo: close output file
}

//à faire
void input(vector<double> &xl,vector<double> &yl,vector<double> &xin,vector<double> &yin,vector<int> &index,vector<double> &ub,int &N,int &IN)
{
    xl = {0.0,0.25,0.5,0.75,1,1,1,1,1,0.75,0.5,0.25,0,0,0,0};
    yl = {0,0,0,0,0,0.25,0.5,0.75,1,1,1,1,1,0.75,0.5,0.25};
    index = {1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0};
    ub = {0,0,0,0,200,200,200,200,0,0,0,0,100,100,100,100};
    xin = {0.25,0.5,0.75,0.25,0.5,0.75,0.25,0.5,0.75};
    yin = {0.25,0.25,0.25,0.5,0.5,0.5,0.75,0.75,0.75};
    N = 16;
    IN = 9;
}

void gmatr(vector<double> &xl,vector<double> &yl,vector<double> &xm,vector<double> &ym,vector<double> &G,int N)
{
  //this subroutine computes the elements of the G matrix
  //compute the nodal coordinates and store them in the arrays xm and ym
  xl[N] = xl[0];
  yl[N] = yl[0];
  double result = 0;

  for(int i = 0; i <= N; i++)
  {
    xm[i] = (xl[i] + xl[i+1])/2.0;
    ym[i] = (yl[i] + yl[i+1])/2.0;
  }

  for(int i = 0; i < N; i++)
  {
    for(int j = 0; j < N; j++)
    {
      if(i != j)
      {
        rlintc(xm[i],ym[i],xl[j],yl[j],xl[j+1],yl[j+1],result);
        G[i*N+j] = result;
      }
      else if(i == j)
      {
        slintc(xl[j],yl[j],xl[j+1],yl[j+1],result);
        G[i*N+j] = result;
      }
    }
  }
}

void rlintc(double x0,double y0,double x1,double y1,double x2,double y2,double &result)
{
  //This subroutine computes the off-diagonal elements of the matrix G
  // ra = The distance of point O from the gauss integration point on the boundary element
  // wg = The weights of the Gauss integration
  // xi = The coordinates of the Gauss integration points in the interval [-1,1]
  //xc,yc = The global coordinates of the Gauss integration points

  vector<double> xc(4,0);
  vector<double> yc(4,0);
  vector<double> xi = {-0.86113631, -0.33998104, 0.33998104, 0.86113631};
  vector<double> wg = {0.34785485, 0.65214515, 0.65214515, 0.34785485};
  double ax = (x2-x1)/2.0;
  double ay = (y2-y1)/2.0;
  double bx = (x2+x1)/2.0;
  double by = (y2+y1)/2.0;

  //compute the line integral

  result = 0;

  for(int i = 0; i < 4; i++)
  {
    xc[i] = ax*xi[i] + bx;
    yc[i] = ay*xi[i] + by;
    double temp = (xc[i]-x0)*(xc[i]-x0) + (yc[i]-y0)*(yc[i]-y0);
    double ra = sqrt(temp);
    result = result + log(ra)*wg[i];
  }
  double temp = ax*ax+ay*ay;
  double sl = 2.0 * sqrt(temp);
  result = result*sl/(4*M_PI);
}

void slintc(double x1, double y1, double x2, double y2, double &result)
{
  //this subroutines computes the diagonal elements of the matrix G
  double ax = (x2-x1)/2.0;
  double ay = (y2-y1)/2.0;
  double temp = ax*ax+ay*ay;
  double sl = sqrt(temp);
  result = sl*(log(sl)-1.0)/M_PI;
}

void hmatr(vector<double> &xl, vector<double> &yl, vector<double> &xm, vector<double> &ym, vector<double> &H, int N)
{
  //This subroutine computes the elements of the H matrix
  xl[N] = xl[0];
  yl[N] = yl[0];
  double result = 0;

  for(int i = 0; i < N; i++)
  {
    xm[i] = (xl[i]+xl[i+1])/2.0;
    ym[i] = (yl[i]+yl[i+1])/2.0;
  }

  //Compute the elements of the H matrix

  for(int i = 0; i < N; i++)
  {
    for(int j = 0; j < N; j++)
    {
      if(i != j)
      {
        dalpha(xm[i],ym[i],xl[j],yl[j],xl[j+1],yl[j+1],result);
        H[i*N+j] = result;
      }
      else if(i == j)
      {
        H[i*N+j] = -0.5;
      }
    }
  }
}

void dalpha(double x0,double y0,double x1,double y1, double x2,double y2, double &result)
{
  //This subroutine computes the off diagonal elements of the H matrix
  double dy1 = y1-y0;
  double dx1 = x1-x0;
  double dy2 = y2-y0;
  double dx2 = x2-x0;
  double temp = dx1*dx1+dy1*dy1;
  double dl1 = sqrt(temp);
  double cos1 = dx1/dl1;
  double sin1 = dy1/dl1;
  double dx2r=dx2*cos1+dy2*sin1;
  double dy2r=-dx2*sin1+dy2*cos1;
  double da = atan2(dy2r,dx2r);
  result = da/(2.0*M_PI);
}

void abmatr(vector<double> &G,vector<double> &H, vector<double> &A, vector<double> &unb, vector<double> &ub, vector<int> &index, int N)
{
  //This subroutine rearranges the matrices G and H and produces the matrices A and B=UNB

  //Reorder the columns of the system of equations and store them in A
  for(int j = 0; j < N; j++)
  {
    if(index[j] == 0)
    {
      for(int i = 0; i < N; i++)
      {
        A[i*N+j] = -G[i*N+j];
      }
    }
    else if(index[j] != 0)
    {
      for(int i = 0; i < N; i++)
      {
        A[i*N+j] = H[i*N+j];
      }
    }
  }

  //compute the right hand side vector and store it in unb

  for(int i = 0; i < N; i++)
  {
    unb[i] = 0;
    for(int j = 0; j < N; j++)
    {
      if(index[j] == 0)
      {
        unb[i] = unb[i] - H[i*N+j]*ub[j];
      }
      else if(index[j] != 0)
      {
        unb[i] = unb[i] + G[i*N+j]*ub[j];
      }
    }
  }
}

void solveq(vector<double> &A,vector<double> &unb, int N, int lsing)
{
  leqs(A,unb,N,lsing);

  if(lsing == 0)
  {
    cout << "The system has been solved regularly" << endl;
  }
  else if(lsing == 1)
  {
    cout << "The system is singular" << endl;
  }
}

void leqs(vector<double> &A,vector<double> &B,int N, int lsing)
{
  //This subroutine uses Gauss elimination to solve a system of linear equations AX=B where
  //A: One-dimentional array which contains the occasional row of the two-dimentional array of the coefficients of the unknowns
  //B: One-dimentional array which contains the known coefficients
  //N: Integer denoting the number of the unknowns
  //lsing: Integer taking the values:
  //lsing = 0, if the system has been solved regularly
  //lsing = 1, if the system is singular

  //todo: corriger les if
  lsing = 0;
  double jj = -N;
  double jy = 0;
  double amax = 0;
  double ihelp = 0;
  double ij = 0;
  double imax = 0;
  double I1 = 0;
  double I2 = 0;
  double atemp = 0;
  double iqs = 0;
  double ixj = 0;
  double ijref = 0;
  double jjx = 0;
  double ny = 0;
  double NN = 0;
  double I3 = 0;

  for(int j = 0; j < N; j++)
  {
    jy = j+1;
    jj = jj + N + 1;
    amax = 0;
    ihelp = jj-j;

    for(int i = j; i < N; i++)
    {
      ij = ihelp+i;
      if((abs(amax)-abs(A[i*N+j])) < 0)
      {
        amax = A[i*N+j];
        imax = i;
      }
      else
      {
        continue;
      }
    }
    if(abs(amax) == 0)
    {
      lsing = 1;
      return;
    }
    I1 = j+N*(j-2);
    ihelp = imax-j;
    for(int k = j; k < N; k++)
    {
      I1 = I1 + N;
      I2 = I1 + ihelp;
      atemp = A[I1];
      A[I1] = A[I2];
      A[I2] = atemp;
      A[I1] = A[I1]/amax;
    }
    atemp = B[imax];
    B[imax] = B[j];
    B[j] = atemp/amax;
    if(j - N < 0 || j - N > 0)
    {
      iqs = N*(j-1);
      for(int ix = jy; ix < N; ix++)
      {
        ixj = iqs + ix;
        ihelp = j - ix;
        for(int jx = jy; jx < N; jx++)
        {
          ijref = N*(jx-1) + ix;
          jjx = ijref + ihelp;
          A[ijref] = A[ijref] - (A[ixj]*A[jjx]);
        }
        B[ix] = B[ix] - B[j]*A[ixj];
      }
    }
    if(j - N == 0)
    {
      ny = N-1;
      NN = N*N;
      for(int j = 0; j < ny; j++)
      {
        I1 = NN - j;
        I2 = N - j;
        I3 = N;
        for(int k = 0; k < N; k++)
        {
          B[I2] = B[I2] - A[I1]*B[I3];
          I1 = I1 - N;
          I3 = I3 - 1;
        }
      }
    }
  }
}


void reorder(vector<double> &ub,vector<double> &unb,vector<int> &index,int N)
{
  //This subroutine rearranges the arrays UB and UNB in such a way that all values of the function u are stored in ub while all the values of the normal derivative un are stored in unb
  double ch = 0;

  for(int i = 0; i < N; i++)
  {
    if(index[i] <= 0)
    {
      continue;
    }
    if(index[i] > 0)
    {
      ch = ub[i];
      ub[i] = unb[i];
      unb[i] = ch;
    }
  }
}

void uinter(vector<double> &xl,vector<double> &yl,vector<double> &xin,vector<double> &yin,vector<double> &ub,vector<double> &unb,vector<double> &uin,int N, double IN)
{
  //this subroutine computes the values of u at the internal points
  //Compute the values of u at the internal points

  double resh = 0;
  double resg = 0;

  for(int k = 0; k < IN; k++)
  {
    uin[k] = 0;
    for(int j = 0; j < N; j++)
    {
      double jp1 = j+1;
      dalpha(xin[k],yin[k],xl[j],yl[j],xl[jp1],yl[jp1],resh);
      rlintc(xin[k],yin[k],xl[j],yl[j],xl[jp1],yl[jp1],resg);
      uin[k] = uin[k] + resh*ub[j] - resg*unb[j];
    }
  }
}

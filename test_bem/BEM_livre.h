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
//#include <QImage>
#include <cstdio>
#include <vector>

using namespace Eigen;
using namespace std;
using namespace gmsh;

void bem(int argc, char **argv);

void gmatr(vector<double> &xl,vector<double> &yl,vector<double> &xm,vector<double> &ym,vector<double> &G,int N);

void rlintc(double x0,double y0,double x1,double y1,double x2,double y2,double &result);

void slintc(double x1, double y1, double x2, double y2, double &result);

void dalpha(double x0,double y0,double x1,double y1, double x2,double y2, double &result);

void hmatr(vector<double> &xl, vector<double> &yl, vector<double> &xm, vector<double> &ym, vector<double> &H, int N);

void abmatr(vector<double> &G,vector<double> &H, vector<double> &A, vector<double> &unb, vector<double> &ub, vector<int> &index, int N);

void solveq(vector<double> &A,vector<double> &B, int N);

void leqs(vector<double> &A,vector<double> &B,int N, int lsing);

void reorder(vector<double> &ub,vector<double> &unb,vector<int> &index,int N);

void uinter(vector<double> &xl,vector<double> &yl,vector<double> &xin,vector<double> &yin,vector<double> &ub,vector<double> &unb,vector<double> &uin,int N, double IN);


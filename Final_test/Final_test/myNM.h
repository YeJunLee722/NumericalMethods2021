#pragma once
/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Lee Ye Jun
Created          : 26-03-2018
Modified         : 18-03-2021
Language/ver     : C== in MSVS2019

Description      : myNM.h
----------------------------------------------------------------*/

#ifndef		_MY_NM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NM_H

#include "myMatrix.h"

// Matrix addition
extern	Matrix	addMat(Matrix _A, Matrix _B);

// Apply back-substitution
extern	void	backSub1(Matrix _U, Matrix _b, Matrix _x);
extern  void    fwdSub1(Matrix _U, Matrix _b, Matrix _x);

// gausElim
extern void gaussElim1(Matrix _A, Matrix _b, Matrix _U, Matrix _d);

//LU
extern void LUdecomp1(Matrix A, Matrix L, Matrix U, Matrix P);
extern void solveLU1(Matrix L, Matrix U, Matrix P, Matrix b, Matrix _x);


//Pivoting
extern int PP1(Matrix A, int i);

extern void LUdecomp_PIVOT1(Matrix _A, Matrix _L, Matrix _U, Matrix _P);

//inverse
extern void inv1(Matrix _A, Matrix _Ainv);

// Newton Rapson

extern double func1(double Tem);
extern double dfunc1(double TEM);
extern double NewtonRapson1(double  _T);


//Condition
extern double cond(Matrix _A);
extern double norm(Matrix _A);
extern double norm2(Matrix _A);
extern Matrix transpose1(Matrix _A);

//Eigen value
extern void QR_1(Matrix _A, Matrix _Q, Matrix _R, Matrix _OUTQ, Matrix _OUTR);
extern void QR_2(Matrix _A, Matrix _Q, Matrix _R, Matrix _OUTQ, Matrix _OUTR, Matrix _OUTA);
extern void eig(Matrix _A, Matrix _OUTA);

extern double Eigen_Compare(Matrix _A);
extern double Minimum(Matrix _A);

//Curve Fitting
extern Matrix	linearFit(Matrix _x, Matrix _y);
// Create a matrix from 1D-array
extern  Matrix	arr2Mat(double* _1Darray, int _rows, int _cols);
extern  Matrix linearInterp(Matrix _x, Matrix _y, Matrix _xq);

//Differrentiate
extern Matrix  gradient(Matrix _x, Matrix _y);
extern void    gradient1D(double _x[], double _y[], double _dydx[], int m);
extern Matrix  gradientFunc(double func(const double x), Matrix _x);

extern double newtonRaphsonFunc(double func(const double x), double dfunc(const double x), float _x0, float _tol);


// ODE 
extern void odeEU(double func(const double x, const double t), double y[], double t0, double tf, double h);
extern void odeEM(double func(const double x, const double t), double y[], double t0, double tf, double h);
extern void ode(double func(const double x, const double t), double y[], double t0, double tf, double h, int method);
extern void Rk3(double func(const double x, const double t), double y[], double t0, double tf, double h);
extern void Rk4(double func(const double x, const double t), double y[], double t0, double tf, double h);
extern void odeMulti(double yfunc(const double z), double zfunc(const double y, const double z, const double t), double y[][2], double t0, double tf, double h);

// inclass integral test
extern double trapz(double _x[], double _y[], int _m);
extern double integral(double func(const double x), double a, double b, int n);
extern double IntegrateRect(double _x[], double _y[], int _m);
extern double Simson83(double func(const double x), double a, double b, int n);


// inclass runge kutta

extern void odeRK2(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);

extern void odeRK4(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);

// ODE RK2:  one of 2nd order ODE <--> two of 1st order ODE
extern void sys2RK2(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init);


// Classical RK4
extern void sys2RK4(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init);
//Lagrange Cubic
extern Matrix Lagrange_cubic(double x[], double y[], double test_x[]);

#endif
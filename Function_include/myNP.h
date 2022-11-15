#pragma once

/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : HanMinung
Created          : 26-03-2018
Modified         : 10-12-2022
Language/ver     : C++ in MSVS2019

Description      : myNP.cpp
----------------------------------------------------------------*/

#define		_CRT_SECURE_NO_WARNINGS
#define		USE_MATH_DEFINES

#define		PI		(	3.14159265358979323846264338327950288419716939937510582	)

#define		FOR_LOOP(i,Initial,Final,Interval)		for(int i = Initial ; i < Final ; i = i + Interval)
#define		New_line(N)		FOR_LOOP(i,0,N,1){printf("\n");}

#define		EU		0
#define		RK2		1
#define		RK3     2
#define		RK4		3

#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*--------------------------------------------------------------------------------------------------------*/

void ErrorHandling(char* message);

void print_str(const char* message);

/*--------------------------------------------------------------------------------------------------------*/

double	Square_Number(double _x);

double  Cubic_Number(double _x);

double	Four_Square_Number(double _x);

double	Factorial(double _x);

/*--------------------------------------------------------------------------------------------------------*/
double	SinTaylor(double _x);

double	SindTaylor(double _x);

double  CosTaylor(double _x);

double  CosdTaylor(double _x);

double	SinTaylor_2(double _x);

double  ExpTaylor(double _x);

/*--------------------------------------------------------------------------------------------------------*/

double Bisection(double func(double x), double a, double b, double tol);
 
double NewtonRaphson(double func(double x), double d_func(double x), double x0, double tol);
 
double Secant(double func(double x), double x0, double x1, double tol);

/*--------------------------------------------------------------------------------------------------------*/

void gradient1D_2Points(double x[], double y[], double dydx[], int m);

void gradient1D_3Points(double x[], double y[], double dydx[], int m);

void gradient1D_4Points(double x[], double y[], double Dy_Dx[], int m);

void Function_call(double func(const double x), double X_in);

void gradientFunc(double func(const double x), double x[], double dydx[], int m);

void gradient2D(double func(const double x), double x[], double dy2dx2[], int m);

void gradient2D_array(double x[], double y[],double dy2dx2[], int m);

void printVec(double* _vec, int _row);

/*--------------------------------------------------------------------------------------------------------*/

double	Integral_Rect(double x[], double y[], int m);

double	Trapz(double x[], double y[], int m);

double  simpson13(double x[], double y[], int m);

// 13 method : I = h/3 * (x[n] + 4x[n+1] + x[n+2])
double	Integral_13(double func(const double x), double x, double y, int N);

// 38 method : I = 3h/8 * (x[n] + 3x[n+1] + 3x[n+2] + x[n+3])
double	Integral_38(double func(const double x), double x, double y, int N);

/*--------------------------------------------------------------------------------------------------------*/

double lagrange_1st(double x, double x0, double x1, double y0, double y1);

double lagrange_2nd(double x, double x0, double x1, double x2, double y0, double y1, double y2);

double lagrange_nth(double x[], double y[], double xx, int length);

double Newton_1st(double x, double x0, double x1, double y0, double y1);

double Newton_2nd(double x, double x0, double x1, double x2, double y0, double y1, double y2);

/*--------------------------------------------------------------------------------------------------------*/

void ode(double func(const double t, const double v), double t0, double tf, double h, double v0, uint8_t method);

void odeEU(double func(const double t, const double v), double t0, double tf, double h, double v0);

void odeEM(double func(const double t, const double v), double t0, double tf, double h, double v0);

void odeRK2(double func(const double t, const double v), double t0, double tf, double h, double v0);

void odeRK3(double func(const double t, const double v), double t0, double tf, double h, double v0);

void odeRK4(double func(const double t, const double v), double t0, double tf, double h, double v0);
#pragma once

/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : HanMinung
Created          : 26-03-2018
Modified         : 10-02-2022
Language/ver     : C++ in MSVS2019

Description      : myNP.cpp
----------------------------------------------------------------*/

#define		CRT_SECURE_NO_WARNINGS
#define		USE_MATH_DEFINES

#define		PI		(	3.14159265358979323846264338327950288419716939937510582	)

#define		FOR_LOOP(i,Initial,Final,Interval)		for(int i = Initial ; i < Final ; i = i + Interval)
#define		New_line(N)		FOR_LOOP(i,0,N,1){printf("\n");}


#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*--------------------------------------------------------------------------------------------------------*/

void ErrorHandling(char* message);

/*--------------------------------------------------------------------------------------------------------*/

extern double	Square_Number(double _x);

extern double   Cubic_Number(double _x);

extern double	Four_Square_Number(double _x);

extern double	Factorial(double _x);

/*--------------------------------------------------------------------------------------------------------*/
extern double	SinTaylor(double _x);

extern double	SindTaylor(double _x);

extern double	SinTaylor_2(double _x);

/*--------------------------------------------------------------------------------------------------------*/

double Bisection(double func(double x), double a, double b, double tol);

double NewtonRaphson(double func(double x), double d_func(double x), double x0, double tol);

double Secant(double func(double x), double x0, double x1, double tol);

/*--------------------------------------------------------------------------------------------------------*/

extern void	gradient1D_2Points(double x[], double y[], double dydx[], int m);

extern void gradient1D(double x[], double y[], double dydx[], int m);

extern void gradient1D_4Points(double x[], double y[], double Dy_Dx[], int m);

extern void Function_call(double func(const double x), double X_in);

extern void gradientFunc(double func(const double x), double x[], double dydx[], int m);

extern void acceleration(double func(const double x), double x[], double dydx[], int m);

extern void printVec(double* _vec, int _row);

/*--------------------------------------------------------------------------------------------------------*/

double	Integral_Rect(double x[], double y[], int m);

double	Trapz(double x[], double y[], int m);

double  simpson13(double x[], double y[], int m);

double	Integral(double func(const double x), double x, double y, int N);

double	Integral_Simpson_38(double func(const double x), double x, double y, int N);








/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : Han,Minung
Created          : 26-03-2018
Modified         : 12-12-2022
Language/ver     : C++ in MSVS2019

Description      : myMatrix.cpp
----------------------------------------------------------------*/

#define _CRT_SECURE_NO_WARNINGS

#ifndef		_MY_NP_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NP_H

#include "myMatrix.h"

// string printer
extern void print_str(const char* message);

// Dimension error check
extern void Dim_error(Matrix _A, const char* _name);

// Matrix addition
extern Matrix addMat(Matrix _A, Matrix _B);

// Matrix subtraction
extern Matrix subMat(Matrix _A, Matrix _B);

// Matrix product
extern Matrix Matproduct(Matrix _A, Matrix _B);

// Apply back-substitution
extern Matrix backSub(Matrix _U, Matrix _d);

// Apply forward-substitution
extern Matrix fwdSub(Matrix _U, Matrix _d);

// -------------------------- ODE -----------------------------------------

void ode(double func(const double t, const double v), double t0, double tf, double h, double v0, uint8_t method);

void odeEU(double func(const double t, const double v), double t0, double tf, double h, double v0);

void odeEM(double func(const double t, const double v), double t0, double tf, double h, double v0);

void odeRK2(double func(const double t, const double v), double t0, double tf, double h, double v0);

void odeRK3(double func(const double t, const double v), double t0, double tf, double h, double v0);

void odeRK4(double func(const double t, const double v), double t0, double tf, double h, double v0);

void sys2RK2(void func(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init);

void sys2RK4(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init);

#endif
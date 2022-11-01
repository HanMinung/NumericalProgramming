/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : Han,Minung
Created          : 26-03-2018
Modified         : 11-02-2022
Language/ver     : C++ in MSVS2019

Description      : myMatrix.cpp
----------------------------------------------------------------*/

#ifndef		_MY_NP_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NP_H

#include "myMatrix.h"

// Matrix addition
extern Matrix addMat(Matrix _A, Matrix _B);

// Matrix product
extern Matrix Matproduct(Matrix _A, Matrix _B);

// Apply back-substitution
extern Matrix backSub(Matrix _U, Matrix _d);

// Apply forward-substitution
extern Matrix fwdSub(Matrix _U, Matrix _d);



#endif
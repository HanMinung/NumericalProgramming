/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : Han,Minung
Created          : 26-03-2018
Modified         : 12-12-2022
Language/ver     : C++ in MSVS2019

Description      : myMatrix.cpp
----------------------------------------------------------------*/

#define _CRT_SECURE_NO_WARNINGS

#include "myNP.h"
#include "myMatrix.h"

void print_str(const char* message) {

	printf("\n****************************************************\n");
	printf("	   %s          \n", message);
	printf("****************************************************\n\n");

}


void Dim_error(Matrix _A, const char* _name) {

	if (_A.rows != _A.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at '%s' function", _name);
		printf("\n*************************************************\n");

		exit(1);

	}

}


// Matrix addition
Matrix	addMat(Matrix _A, Matrix _B){

	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addMat' function");
		printf("\n*************************************************\n");

		exit(1);
	}

	Matrix Out = zeros(_A.rows, _B.cols);

	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j] = _A.at[i][j] + _B.at[i][j];

	return Out;
}


Matrix	subMat(Matrix _A, Matrix _B){

	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addMat' function");
		printf("\n*************************************************\n");

		exit(1);
	}

	Matrix Out = zeros(_A.rows, _B.cols);

	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j] = _A.at[i][j] - _B.at[i][j];

	return Out;
}


Matrix	backSub(Matrix _U, Matrix _d)
{
	if (_U.rows != _d.rows) {

		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'backSub' function");
		printf("\n*************************************************\n");

		exit(1);
	}

	Matrix Out = zeros(_d.rows, 1);
	initMat(Out, 0);

	double sum = 0;

	FOR_LOOP_INV_INC(i,_U.rows-1, 0, 1){

		sum = 0;
		
		FOR_LOOP(j,i+1,_U.cols,1){
		
			sum = sum + (_U.at[i][j] * _d.at[j][0]);

		}

		_d.at[i][0] = (_d.at[i][0] - sum) / _U.at[i][i];
	}
	
	copyMat(_d, Out);

	return Out;
}


Matrix  fwdSub(Matrix _L, Matrix _b) {

	if (_L.rows != _b.rows) {

		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'fwdSub' function");
		printf("\n*************************************************\n");

		exit(1);
	}

	Matrix Out = zeros(_b.rows, 1);
	Matrix _copyb = zeros(_b.rows, _b.cols);

	copyMat(_b, _copyb);

	double sum = 0;

	FOR_LOOP(i,0,_copyb.rows,1){
		sum = 0;  

		FOR_LOOP(j,0,i,1) {

			sum = sum + _L.at[i][j] * _copyb.at[j][0];
		}

		_copyb.at[i][0] = (_copyb.at[i][0] - sum) / _L.at[i][i];
	}

	copyMat(_copyb, Out);
	freeMat(_copyb);

	return Out;
}


Matrix Matproduct(Matrix _A, Matrix _B) {

	if (_A.cols != _B.rows) {

		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'Matproduct' function");
		printf("\n*************************************************\n");

		exit(1);
	}

	Matrix Out = zeros(_A.rows, _B.cols);
	
	FOR_LOOP(k,0,_A.rows,1){

		FOR_LOOP(i,0,_A.cols,1){

			FOR_LOOP(j, 0, _B.cols, 1) {

				Out.at[k][j] += _A.at[k][i] * _B.at[i][j];

			}
		}
	}

	return Out;
}


/*--------------------------------------------------------------------------------------------------------*/

// select which method is used to solve ODE
// --- method ---
//	  EU	0
//	  EM	1
//	  RK2   2
//	  RK3	3
//	  RK4	4
void ode(double func(const double t, const double v), double t0, double tf, double h, double v0, uint8_t method) {

	switch (method) {

	case(0):

		odeEU(func, t0, tf, h, v0);
		break;

	case(1):

		odeEM(func, t0, tf, h, v0);
		break;

	case(2):

		odeRK2(func, t0, tf, h, v0);
		break;

	case(3):

		odeRK3(func, t0, tf, h, v0);
		break;

	case(4):

		odeRK4(func, t0, tf, h, v0);
		break;

	}
}



// 오일러 method
void odeEU(double func(const double t, const double v), double t0, double tf, double h, double v0) {

	double slope = 0;
	double N = (tf - t0) / h;

	double* t = NULL;
	double* v = NULL;

	t = (double*)malloc(sizeof(double) * (N + 1)); //t  N+1 = 101 : 0~100 ==> 101개 
	v = (double*)malloc(sizeof(double) * (N + 1));

	t[0] = t0;
	v[0] = v0;

	FOR_LOOP(i, 0, N, 1) {
		// first slope 
		slope = func(t[i], v[i]);

		// update
		v[i + 1] = v[i] + slope * h;
		t[i + 1] = t[i] + h;

	}

	print_str("odeEU");

	FOR_LOOP(i, 0, N + 1, 1) {

		printf("The odeEU v[%2d] is %lf \n", i, v[i]);
	}

	New_line(2);

	FILE* _fp;
	_fp = fopen("odeEU.txt", "wt");

	if (_fp == NULL) {
		printf("Fail");
		return;
	}

	FOR_LOOP(i, 0, N + 1, 1) {

		fprintf(_fp, "%f\n", v[i]);
	}

	fclose(_fp);

	free(t);	free(v);
}

// 오일러 modified
void odeEM(double func(const double t, const double v), double t0, double tf, double h, double v0) {

	double N = (tf - t0) / h;

	double* t = NULL;
	double* v = NULL;

	t = (double*)malloc(sizeof(double) * (N + 1));		//t  N+1 = 101 : 0~100 ==> 101개 
	v = (double*)malloc(sizeof(double) * (N + 1));

	//intial vaule
	t[0] = t0;
	v[0] = v0;

	// 1st slope & next slope : average --> Euler modified
	double slope_f = 0;
	double slope_b = 0;

	FOR_LOOP(i, 0, N, 1) {

		slope_f = func(t[i], v[i]);
		v[i + 1] = v[i] + slope_f * h;

		t[i + 1] = t[i] + h;
		slope_b = func(t[i + 1], v[i + 1]);

		// update y[i+1]
		v[i + 1] = v[i] + (slope_b + slope_f) * h / 2;

	}

	print_str("odeEM");

	FOR_LOOP(i, 0, N + 1, 1) {

		printf("The odeEM t[%2d] is %lf \n", i, v[i]);
	}

	New_line(2);

	FILE* _fp;
	_fp = fopen("odeEM.txt", "wt");

	if (_fp == NULL) {
		printf("Fail");
		return;
	}

	FOR_LOOP(i, 0, N + 1, 1) {

		fprintf(_fp, "%f\n", v[i]);
	}

	fclose(_fp);

	free(t);	free(v);

}


void odeRK2(double func(const double t, const double v), double t0, double tf, double h, double v0) {

	double N = (tf - t0) / h;

	double* t = NULL;
	double* v = NULL;

	t = (double*)malloc(sizeof(double) * (N + 1));
	v = (double*)malloc(sizeof(double) * (N + 1));

	t[0] = t0;	v[0] = v0;

	double K_1 = 0;	double K_2 = 0;

	//default
	double alpha = 1;
	double beta = alpha;

	double c2 = 1 / (2 * alpha);
	double c1 = 1 - c2;

	FOR_LOOP(i, 0, N, 1) {

		K_1 = func(t[i], v[i]);
		K_2 = func(t[i] + alpha * h, v[i] + beta * h * K_1);
		 
		v[i + 1] = v[i] + (c1 * K_1 + c2 * K_2) * h;

		t[i + 1] = t[i] + h;
	}

	print_str("RK2 method");

	FOR_LOOP(i, 0, N + 1, 1) {

		printf("The RK2 v[%2d] is %lf \n", i, v[i]);
	}

	New_line(2);

	FILE* _fp;
	_fp = fopen("odeRK2.txt", "wt");

	if (_fp == NULL) {
		printf("Fail");
		return;
	}

	FOR_LOOP(i, 0, N + 1, 1) {

		fprintf(_fp, "%f\n", v[i]);
	}

	fclose(_fp);

	free(t);	free(v);

}


void odeRK3(double func(const double t, const double v), double t0, double tf, double h, double v0) {

	double N = (tf - t0) / h;

	double* t = NULL;
	double* v = NULL;

	t = (double*)malloc(sizeof(double) * (N + 1));
	v = (double*)malloc(sizeof(double) * (N + 1));

	t[0] = t0;	v[0] = v0;

	double K_1 = 0;
	double K_2 = 0;
	double K_3 = 0;

	double alpha2 = (1.0 / 2.0);
	double alpha3 = 1.0;

	double beta21 = (1.0 / 2.0);
	double beta31 = -1.0;
	double beta32 = 2.0;

	double c1 = 1.0 / 6.0;
	double c2 = 4.0 / 6.0;
	double c3 = 1.0 / 6.0;


	FOR_LOOP(i, 0, N, 1) {

		K_1 = func(t[i], v[i]);
		K_2 = func(t[i] + alpha2 * h, v[i] + beta21 * K_1 * h);
		K_3 = func(t[i] + alpha3 * h, v[i] + beta31 * K_1 * h + beta32 * K_2 * h);

		v[i + 1] = v[i] + (c1 * K_1 + c2 * K_2 + c3 * K_3) * h;

		t[i + 1] = t[i] + h;
	}

	print_str("RK3 method");

	FOR_LOOP(i, 0, N + 1, 1) {

		printf("The RK3 v[%2d] is %lf \n", i, v[i]);
	}

	New_line(2);

	FILE* _fp;
	_fp = fopen("odeRK3.txt", "wt");

	if (_fp == NULL) {
		printf("Fail");
		return;
	}

	FOR_LOOP(i, 0, N + 1, 1) {

		fprintf(_fp, "%f\n", v[i]);
	}

	fclose(_fp);

	free(t);	free(v);

}


void odeRK4(double func(double t, double v), double t0, double tf, double h, double v0) {

	double N = (tf - t0) / h;

	double* t = NULL;
	double* v = NULL;

	t = (double*)malloc(sizeof(double) * (N + 1)); //t  N+1 = 101 : 0~100 ==> 101개 
	v = (double*)malloc(sizeof(double) * (N + 1));

	//intial vaule
	if (t != NULL && v != NULL) {
		t[0] = t0;
		v[0] = v0;
	}
	else   printf("NULL POINTER");

	double K1 = 0;
	double K2 = 0;
	double K3 = 0;
	double K4 = 0;

	double alpha2 = 1.0 / 2.0;
	double alpha3 = 1.0 / 2.0;
	double alpha4 = 1.0;
	double beta21 = 1.0 / 2.0;
	double beta31 = 0;
	double beta32 = 1.0 / 2.0;
	double beta41 = 0.0;
	double beta42 = 0.0;
	double beta43 = 1.0;

	double c1 = 1.0 / 6.0;
	double c2 = 2.0 / 6.0;
	double c3 = 2.0 / 6.0;
	double c4 = 1.0 / 6.0;

	FOR_LOOP(i, 0, N, 1) {

		K1 = func(t[i], v[i]);
		K2 = func(t[i] + 0.5 * h, v[i] + 0.5 * K1 * h);
		K3 = func(t[i] + 0.5 * h, v[i] + 0.5 * K2 * h);
		K4 = func(t[i] + 0.5 * h, v[i] + 0.5 * K3 * h);

		v[i + 1] = v[i] + (c1 * K1 + c2 * K2 + c3 * K3 + c4 * K4) * h;
		t[i + 1] = t[i] + h;
	}

	print_str("RK4 method");

	FOR_LOOP(i, 0, N + 1, 1) {

		printf("The RK4 v[%2d] is %lf \n", i, v[i]);
	}

	New_line(2);

	FILE* _fp;
	_fp = fopen("odeRK4.txt", "wt");

	if (_fp == NULL) {
		printf("Fail");
		return;
	}

	FOR_LOOP(i, 0, N + 1, 1) {

		fprintf(_fp, "%f\n", v[i]);
	}

	fclose(_fp);

	free(t);		free(v);

}

/*
	func : 2nd order differential euqationn
	y1	 : f_y ( y dot = z )
	y2   : f_z ( z dot	   )
	t0 , tf , h : initial value , final value , interval
	f_y , f_z initial value
*/
void sys2RK2(void func(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init) {

	// number of points
	double N = (tf - t0) / h + 1;

	// start point
	double ti = t0;

	// Initialization
	double K1[2] = { 0, };   // K1 = [K1_y1, K1_y2]    first slope
	double K2[2] = { 0, };   // K2 = [K2_y1, K2_y2]    second slope

	double Yin[2] = { 0, };  // u(t)
	double K1_y1 = 0;		 // F_1[t_i, y_i , z_i]	 :  dy/dt = z
	double K1_y2 = 0;		 // F_1[t_i + h , y_i + (K_y1 * h) , z_i + (K_z1 * h)]
	double K2_y1 = 0;		 // F_2[t_i, y_i , z_i]	 :  dz/dt = z dot
	double K2_y2 = 0;		 // F_2[t_i + h , y_i + (K_y1 * h) , z_i + (K_z1 * h)]

	// Initial condition
	y1[0] = y1_init;		 // y(t)
	y2[0] = y2_init;		 // z(t) = dydt(t)

	FOR_LOOP(i, 0, N - 1, 1) {

		// Slope 1 [ K1 ]
		Yin[0] = y1[i];		 // z
		Yin[1] = y2[i];		 // dzdt		

		func(ti, Yin, K1);   // Yin = Input     K1 = Output
		K1_y1 = K1[0];
		K1_y2 = K1[1];

		// Slope 2 [ K2 ]
		Yin[0] = y1[i] + K1_y1 * h;		// z              // Euler Method, get Next state(i+1) value 
		Yin[1] = y2[i] + K1_y2 * h;		// dzdt		      // Euler Method, get Next state(i+1) value

		func(ti + h, Yin, K2);		    // Yin = Input     K2 = Output
		K2_y1 = K2[0];
		K2_y2 = K2[1];


		// Update
		y1[i + 1] = y1[i] + (K1_y1 + K2_y1) * h / 2;
		y2[i + 1] = y2[i] + (K1_y2 + K2_y2) * h / 2;

		ti += h;
	}

	FOR_LOOP(i, 0, N, 1) {

		printf("The sys2RK2 v[%2d] is %lf \n", i, y1[i]);
	}

	New_line(2);

	FILE* _fp;
	_fp = fopen("odesys2RK2.txt", "wt");

	if (_fp == NULL) {
		printf("Fail");
		return;
	}

	FOR_LOOP(i, 0, N, 1) {

		fprintf(_fp, "%f\t%f\n", y1[i]);
	}

	fclose(_fp);

}


void sys2RK4(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init) {

	// number of points
	int N = (tf - t0) / h + 1;

	// statr point
	double ti = t0;

	// Initialization
	double K1[2] = { 0, };   // K1 = [K1_y1, K1_y2]
	double K2[2] = { 0, };   // K2 = [K2_y1, K2_y2]
	double K3[2] = { 0, };   // K3 = [K3_y1, K3_y2]
	double K4[2] = { 0, };   // K4 = [K4_y1, K4_y2]

	double Yin[2] = { 0 };
	double K1_y1 = 0;   // y1 :: dydt = z
	double K1_y2 = 0;   // y2 :: dzdt = zdot
	double K2_y1 = 0;
	double K2_y2 = 0;
	double K3_y1 = 0;
	double K3_y2 = 0;
	double K4_y1 = 0;
	double K4_y2 = 0;


	// Initial condition
	y1[0] = y1_init;   //  y(t)
	y2[0] = y2_init;   //  z(t) = dydt(t)

	for (int i = 0; i < N - 1; i++) {

		// Slope 1 [ K1 ]
		Yin[0] = y1[i];		// z
		Yin[1] = y2[i];		// dzdt		

		odeFunc_sys2(ti, Yin, K1);   // Yin = Input     K1 = Output
		K1_y1 = K1[0];
		K1_y2 = K1[1];


		// Slope 2 [ K2 ]
		Yin[0] = y1[i] + ((0.5) * K1_y1 * h);	     	// z              // Euler Method, get Next state(i+1) value 
		Yin[1] = y2[i] + ((0.5) * K1_y2 * h);	    	// dzdt		      // Euler Method, get Next state(i+1) value

		odeFunc_sys2(ti + (0.5 * h), Yin, K2);            // Yin = Input     K2 = Output
		K2_y1 = K2[0];
		K2_y2 = K2[1];


		// Slope 3 [ K3 ]
		Yin[0] = y1[i] + ((0.5) * K2_y1 * h);		   // z              // Euler Method, get Next state(i+1) value 
		Yin[1] = y2[i] + ((0.5) * K2_y2 * h);		   // dzdt		      // Euler Method, get Next state(i+1) value

		odeFunc_sys2(ti + (0.5 * h), Yin, K3);           // Yin = Input     K3 = Output
		K3_y1 = K3[0];
		K3_y2 = K3[1];


		// Slope 4 [ K4 ]
		Yin[0] = y1[i] + (K3_y1 * h);		          // z               // Euler Method, get Next state(i+1) value 
		Yin[1] = y2[i] + (K3_y2 * h);		          // dzdt		      // Euler Method, get Next state(i+1) value

		odeFunc_sys2(ti + h, Yin, K4);                // Yin = Input     K4 = Output
		K4_y1 = K4[0];
		K4_y2 = K4[1];

		// Update
		y1[i + 1] = y1[i] + ((K1_y1 + (2 * K2_y1) + (2 * K3_y1) + K4_y1) * (h / 6));
		y2[i + 1] = y2[i] + ((K1_y2 + (2 * K2_y2) + (2 * K3_y2) + K4_y2) * (h / 6));

		ti += h;
	}
}
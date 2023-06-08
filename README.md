# Numerical Programming

* Author : HanMinwoong
* Modified : 2022.11.12
* Mechanical and Control Engineering , Handong Global University



## myNP.c (Function definition for implementation)

-----------

* Used macro
  * for loop iteration ( FOR_LOOP )
  * newline ( New_line )

### Taylor series

*  to approximate Sin function , Cos function , Exponential function

* Basic idea : taylor series approximation of functions

  <img src="https://user-images.githubusercontent.com/99113269/197355079-e1d2453d-b1e2-4396-bb15-2ad112b3dfa5.jpg" alt="Taylor-Series" style="zoom: 67%;" />

  * Input : [ rad ]

  * Input : [ deg ]

  

### Nonlinear solver

* Bisection method
* Newton - Raphson method
* Secant method



### Differentiation

* 2 point differentiation
* 3 point differentiation
* 4 point differentiation
* Input element : array  or  user-defined function



### Integration

* Rectangular method
* Midpoint method
* Trapezoidal method
* Simpson 13 ( 2 order polynomial - integrate )
  * parameter : function, two points , N value
* Simpson 38 ( 3 order polynomial - integrate )
  * parameter : function, two points , N value



### Interpolation

* parameter : array elements , value to find 

* Lagrange 1st order polynomial
* Lagrange 2nd order polynomial
* Lagrange nth order polynomial
  * order 3 : need 3 array elements to interpolate
  * order 4 : need 4 array elements to interpolate
  * ...
* Newton 1st order polynomial
* Newton 2nd order polynomial



### ODE solver

--------

#### First order differential equation

* Basic solving algorithm

  <img src="https://user-images.githubusercontent.com/99113269/201662480-f8b4ab94-a8d1-4eba-ad2c-f91e56ca3a8c.png" alt="image" style="zoom:33%;" />

  

* Method 1 : Euler Explicit method

  * Basic idea : rectangular integral method

  <img src="https://user-images.githubusercontent.com/99113269/201665465-c8db2ab5-303c-4584-8edc-159f8fdb701a.png" alt="image" style="zoom: 50%;" />

  * Pseudo code of Euler explicit method

    ```c
    1. get initial value : t[0]  y[0]
    2. slope = func(t[0],y[0]);
    // calculate slope with current point
    3. v[i+1] = v[i] + slope * h;
       t[i+1] = t[i] + h;
    ```

  

* Method 2 : Euler Implicit method

  * Similar to above algorithm
  * Use the slope at unknown (next step)

  <img src="https://user-images.githubusercontent.com/99113269/201666944-585bb72f-4ced-4ea9-8a2c-2930d90d5ca3.png" alt="image" style="zoom:50%;" />

* Method 3 : Modified Euler's method

  * Utilize the average value of the current point slope and next point slope

    <img src="https://user-images.githubusercontent.com/99113269/201668038-150fff81-baf3-4e76-af8a-56ac30c1fabf.png" alt="image" style="zoom: 50%;" />

  * Pseudo code of modified Euler method

    ```c
    1. get initial value : t[0]  y[0]
    2. FOR_LOOP(int i = 0; i < N ; i++){
        
        slope1 = func(t[i],y[i]);
        v[i+1] = v[i] + slope * h;
        
        t[i+1] = t[i] + h;
        slope2 = func(t[i+1],y[i+1]);
        
        y[i+1] = y[i] + (slope1 + slope2)*h/2;
    }
    ```

  

* Method 4 : 2nd order Runge-Kutta method

  * Mathmatical background is as follow :

    <img src="https://user-images.githubusercontent.com/99113269/206904714-3357f2de-71bc-482a-b6c4-e32ec0e93d65.png" alt="image" style="zoom: 50%;" />

  * utilize the weight to get better performance

  * less error than above methods according to the proper weights.

  * Pseudo code of 2nd order  Runge-Kutta method

    ```c
    1. get initial value : t[0] , y[0]
    2. fix the alpha, Beta value
        --> Coefficient (c1, c2) are fixed according to these values.
    3. FOR_LOOP(int i = 0; i < N ; i++){
        
        K1 = func(t[i], y[i]);
        K2 = func(t[i] + alpha * h, v[i] + Beta * h * K_1);
        
        v[i+1] = v[i] + (c1 * K_1 * c2 * K_2) * h;
        t[i+1] = t[i] + h;
    } 
    ```

* ODE method selector : function ODE

  * Code implementation

    ```c
    void ode(double func(const double t, const double v), double t0, double tf, double h, double v0, uint8_t method) {
    	
    	switch (method) {
    		case(0) :
    			odeEU(func, t0, tf, h, v0);				// define : 0
    			break;
    
    		case(1) :
    			odeRK2(func, t0, tf, h, v0);			// define : 1
    			break;
    		
    		case(2):
    			odeRK3(func, t0, tf, h, v0);			// define : 2
    			break;
                
    		case(3):
    			odeRK4(func, t0, tf, h, v0);			// define : 3
    			break;
                
    	}
    }
    ```

  

* Code implementation of every methods are finished with file printing ( .txt ) in order to compare the result with that of MATLAB.  

  * MATLAB code (reference)

    ```matlab
    clear all; close all; clc;
    
    filepath = 'blah';
    filename = 'blahblah';
    
    addpath(filepath);
    odeEU = readmatrix('filename');
    
    a=0; b=0.1; h=0.001; 
    y0 = 0;
    
    t=a:h:b;    ytrue(1) = y0;
    t = t';
    
    [tmat,ymat] = ode45(@myFunc, [a b], y0);
    
    figure();
    plot(tmat,ymat,'-r',Linewidth=1.5);
    hold on;
    plot(t,odeEU(:),'--blue',Linewidth=1.5);
    hold off;
    xlabel('x [-]',FontSize = 14); ylabel('y [-]',FontSize =14); title("ode45 , odeEU",FontSize = 14);
    legend("ode45","odeEU",FontSize = 14);
    grid minor;
    ```



#### Second order differential equation

* sys2RK2

  * function to solve 2nd order differential equation with 2nd order Runge-Kutta method.

  * parameter : function (t, y, ydot) - differential equation , y1, y2, t0, tf , interval, initial value of z and zdot

  * include file printing process to get the y value

  * After printing, compare the result with that of MATLAB.

  * Pseudo code 

    ```c
    1. initialization :  t0 , y0 
    2. N = (tf - t0)/h
    3. z = dy/dt	
       zdot = -(k*y - c*z + F_in(t))/m;
       
    4. repetition
       for(i = 0 to N-1 , 1++){
           
           // get each K value
           t(i+1) = t(i) + h;
           K_y1 = f1(t_i, y_i, z_i);
           K_z1 = f2(t_i, y_i, z_i);
           K_y2 = f1(t_i + h, y_i + K_y1 * h, z_i + K_z1 * h);
           K_z2 = f1(t_i + h, y_i + K_y1 * h, z_i + K_z1 * h);
           
           y(i+1) = y(i) + h/2 * (K_y1 + K_y2);
           z(i+1) = z(i) + h/2 * (K_z1 + K_z2);
       
       }
    ```
    
    





---------

### myMatrix.h (Function definition for implementation)

-------------

* **<u>*BAISC STRUCTURE OF MATRIX*</u>**

  ```c
  typedef struct { 
      
  	double** at;
  	int rows, cols;
      
  }Matrix;
  
  - need double pointer to create 2nd order matrix like MxN matrix
  ```

  -------------------------------

  **<u>Basic Function lists</u>**

  

  ### Create a matrix

* Matrix createMat( int _rows, int _cols )

  * Function return type : Matrix

  * Code structure

    ```c
    Matrix Out;
    // 1. Allocate row array first
    Out.at = (double**)malloc(sizeof(double*) * _rows);
    
    // 2. Then, allocate column 
    for (int i = 0; i < _rows; i++)
    	Out.at[i] = (double*)malloc(sizeof(double) * _cols);
    	
    	// 3. Initialize row & column values of a matrix
    	Out.rows = _rows;
    	Out.cols = _cols;
    
    	return Out;
    ```

* void freeMat ( Matrix _A )

  * Free allocated column memory
  * Free allocated row memory

* void printMat ( Matrix _A , const char* _name )

  * Print matrix elements
  * Print form : name - elements

* void initMat ( Matrix _A , double _val )

  * Fill elements of Matrix _A with 'val'

* Matrix zeros ( int _rows , int _cols )

  * parameter : n ( column ) , n ( row )

  * return value : Matrix 

  * Code structure

    ```c
    Matrix	zeros(int _rows, int _cols)
    {
    	Matrix Out = createMat(_rows, _cols);
    	initMat(Out,0);
        
    	return Out;
    }
    ```

* void copyMat ( Matrix _A , Matrix _B )

  * parameter : 2 matrix

  * copy elements of matrix A to those of matrix B

  * Code structure

    ```c
    void copyMat(Matrix _A, Matrix _B) {
    
    	if (_A.cols != _B.cols || _A.rows != _B.rows) {
    		printf("\n*************************************************");
    		printf("\n  ERROR!!: dimension error - 'copyMat' function");
    		printf("\n*************************************************\n");
    		/*return createMat(0, 0);*/
    
    		exit(1);
    	}
    
    	FOR_LOOP(i, 0, _A.rows, 1) {
    		FOR_LOOP(j, 0, _A.cols, 1) {
    
    			_B.at[i][j] = _A.at[i][j];
    		}
    	}	
    }
    ```



### Transpose

* void transpose ( Matrix _A )

  * parameter : matrix _A
  * ( _A.rows , _A.cols )  --->  ( _A.cols , _A.rows )

  ```c
  Matrix transpose(Matrix _A) {
  	Matrix Out = zeros(_A.cols, _A.rows);
  
  	FOR_LOOP(i, 0, _A.rows, 1) {
  		FOR_LOOP(j, 0, _A.cols, 1) {
  
  			Out.at[j][i] = _A.at[i][j];
  		}
  	}
  	return Out;
  }
  ```

* Matrix rowExt ( Matrix _A, int in )

  * function that extracts specific row from a Matrix
  * return with matrix form

* void rowIns ( Matrix _v, int in, Matrix _A )

  * parmaters 
    * _v : a specific row wants to insert to a matrix
    * in : row index
    * _A : matrix that wants to insert a row

* void rowExchange( Matrix _A, int seq1, int seq2 )

  * with predefined functions above, we can define a function that changes two rows in a matrix.
  * Code is as follows :

  ```c
  void rowExchange(Matrix _A, int seq1, int seq2) {
  
  	Matrix row_f = rowExt(_A, seq1);
  	Matrix row_s = rowExt(_A, seq2);
      
  	rowIns(row_f, seq2, _A);
  	rowIns(row_s, seq1, _A);
  }
  ```



***<u>Main function lists</u>***

------------

### Gauss elimination

* void gaussElim ( Matrix _A , Matrix _b , Matrix _U , Matrix _d )

  * Solving a system of linear equation
  * parameter : user-defined matrix _A and _b
  * Usually, it is efficient to get matrix by typing elements in txt file and use 'txt2Mat' function.

    * It is essential to check the path of txt file !
    * It is essential to double check the size of matrices !
    * If required sizes of matrices to solve any problem are not satisfied, it is set to output an error.
  * A --> row echelon form (A) --> copy elements of A to U ( copyMat )
  
  * size of Matrices
  
    * _A : n x n
    * _b : n x 1
    * _U : n x n
    * _d : n x 1
  
    

### Back substitution method

* Matrix backSub ( Matrix _U , Matrix _d )
  * parameter
    * Matrix _U  :  row echelon form of original matrix _A
    * Matrix _b : row echelon form of original matrix _b
    
  * pseudo code of back-substitution method

    ```c
    // pseudo code
    for i = n to 1
        for j = (i+1) to n
            X_i = { b_i - a(i)(j) * X_j }
    	end
    end
    ```




### LU decomposition

* LU decomp ( Matrix _A , Matrix _L , Matrix _U )

  * parameter
    * Matrix _A  :  Matrix to decompose to L*U form
    * condition : In the case where it can be expressed in a REF form through a row contract  without a row exchange
    * Matrix _L  :  Lower triangle matrix
    * Matrix _U :  REF form

  * After LU decomposition

    * copy Matrix _A  -->  Matrix _U

    * get y = ( U * X ) with forward substitution method
    * get final X value with backward substitution method 



### LU solver


* void solveLU ( Matrix _L , Matrix _U , Matrix _b , Matrix _x )

  * parameter

    * Matrix _L : lower triangular matrix
    * Matrix _U : REF of original matrix by Gauss elimination
    * Matrix _b : input
      * Size : _L.rows x 1 
    * Matrix _x : solution to solve decomposed elements : L U
      * size : _L.rows x 1

  * Does not need to use other statement , can be made with user-predefined functions like 'backSub', 'fwdSub'.

  * Code structure :

    ```c
    void solveLU(Matrix _L, Matrix _U, Matrix _b, Matrix _x) {
    
    	Matrix _y = zeros(_b.rows,1);
    	_y = fwdSub(_L,_b);
    	_y = backSub(_U,_y);
    
    	copyMat(_y,_x);
        
    	freeMat(_y);			// free mempory allocation
    }
    ```




### Inverse matrix

* Matrix invMat ( Matrix _A )

  * return inverse matrix of _A
  * main code structure

  ```c
  Matrix invMat(Matrix _A) {
  
  	* need exception handler
  
  	uint8_t N = _A.rows;
  	
  	Matrix _L = eye(N, N);		Matrix _U = zeros(N, N);	Matrix _invA = zeros(N, N);
  	Matrix _bk = zeros(N, 1);	Matrix _x = zeros(N, 1);	Matrix _I = eye(N, N);
  
  	LUdecomp(_A, _L, _U);
  
  	FOR_LOOP(k, 0, N, 1) {
  
  		* need exception handler
  
  		copyMat(vectExt(_I, k), _bk);
  		solveLU(_L, _U, _bk, _x);
  		vectins(_x, k, _invA);	
  	}
  	return _invA;
  
  	freeMat(_L);	freeMat(_U);	freeMat(_invA);
  	freeMat(_bk);	freeMat(_x);	freeMat(_I);
  
  }
  
  ```

  

### QR householder

* void QRhouseholder( Matrix _A , Matrix _Q , Matrix _R )
  * parameter : original matrix _A , orthogonal matrix _Q , upper triangular matrix _R
  * one of the method of QR factorization
  * Necessity : a way to convert original matrix to quadrature matrix with preserving eigenvalues.
  * Needs to be predefined to implement a function that finds eigenvalues of a matrix.



### Eigenvalue

* Matrix eig ( Matrix _A )

  * parameter : A matrix to find a eigenvalues
  * Algorithm
    * represent _A with the multiplication of matrix _Q and matrix _R
    * Repeat until the sum of the elements below the diagonal elements converges under the epsilon value.
  * The eigenvalues are stored in the form of matrices.

* example code with QR household

  ```c
  Matrix _A = txt2Mat(path, "prob_matA");
  Matrix _eig = zeros(_A.rows,1);
  
  // stores eigen value of matrix _A in matrix _eig
  _eig = eig(_A);
  ```

  

### Eigenvectors

* Matrix eigVec ( Matrix _A )

  * parameter : A matrix find a eigenvectors
  * returns eigen vectors of matrix _A

* example code 

  ```
  Matrix _A = txt2Mat(path, "prob_matA");
  Matrix _eigVec = zeros( _A.rows , _A.cols );
  
  _eigVen = eigVen(_A);
  ```



### Linear regression : linear polyfit

* Mathematical basis for linear regression :

  <img src="https://user-images.githubusercontent.com/99113269/206892118-09e84eb8-8700-42c0-aa08-d15d3cd10e47.png" alt="image" style="zoom: 40%;" />

* double LinRegress (Matrix _x, Matrix _y, double find_x)

* preprocessing : convert array data to a matrix form

* Get the regressed equation from a dispersed data.

* Returns the value that user wants to predict from the regressed equation.

* Example code

  ```c
  double find_y = 0;
  double sol = 0;
  // dispersed data
  double arr_x[M] = [ - - - - - - - - - - -];
  double arr_y[M] = [ - - - - - - - - - - -];
  
  Matrix mat_x = arr2mat(arr_x, M, 1);
  Matrix mat_y = arr2mat(arr_y, M, 1);
   
  sol = LinRegress(mat_x, mat_y, find_y);
  ```



### High order curve fitting

* Mathematical basis for polyfit :

  <img src="https://user-images.githubusercontent.com/99113269/206893628-e6b94b56-6591-45a1-a69f-99c23c91134c.png" alt="image" style="zoom:40%;" />

* Return a matrix form filled with coefficients

* If 2nd order polynomial : returns coefficients for higher order to lower sequencially. (same with return value that is from MATLAB)

* parameter : Matrix _x , Matrix _y, Order of polynomial

* Example code

  ```c
  Matrix sol = zeros(n,1);
  // adjusting factor (order of polynomial)
  int P_order = k; 			
  
  // dispersed data : #m data points
  double arr_x[M] = [ - - - - - - - - - - -];
  double arr_y[M] = [ - - - - - - - - - - -];
  
  Matrix mat_x = arr2mat(arr_x, M, 1);
  Matrix mat_y = arr2mat(arr_y, M, 1);
  
  sol = polyfit(mat_x, mat_y, k);
  
  printMat(sol,"Coefficients of k_th order polynomial");
  ```




### System of Nonlinear equations

* Consists of multiple non-linear equations that need to be solved simultaneouly.

* Newton-Raphson method is applied.

* parameter : initial value, defined function, according Jacobian matrix 

* Mathmatical principle is as follow :

  <img src="https://user-images.githubusercontent.com/99113269/206907485-18f4fef7-8669-45f7-a9ad-e95f35a44157.png" alt="image" style="zoom:50%;" />

* Example code :

  ```c
  1. Precondition : 
  	predefined F(X_k) and according Jacobian matrix of F in main code
  2. main statement : 
  
  Matrix myFunc(const double _x, const double _y){
  
  	Matrix func = createMat(2, 1);
  
  	func.at[0][0] = _y - (pow(2.71828, _x / 2) + pow(2.71828, (-_x) / 2)) / 2;
  	func.at[1][0] = 9 * pow(_x, 2) + 25 * pow(_y, 2) - 225;
  
  	return func;
  }
  
  Matrix myJacob(const double _x, const double _y){
  
  	Matrix dfunc = createMat(2, 2);
  
  	dfunc.at[0][0] = -((pow(2.71828, _x / 2) - pow(2.71828, (-_x) / 2)) / 4);
  	dfunc.at[0][1] = 1;
  	dfunc.at[1][0] = 18 * _x;
  	dfunc.at[1][1] = 50 * _y;
  
  	return dfunc;
  }
  
  3. Matrix sol = zeros(2,1);
     sol = solve_nonlinear(Z, myFunc, myJacob);
  
  4. printMat(sol,"solution of non-linear system");
  ```

  

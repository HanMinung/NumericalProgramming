# Numerical Programming

* Author : HanMinung
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

  ***<u>Basic Function lists</u>**

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
    * Matrix _A  :  Matrix to decompose LU
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

  




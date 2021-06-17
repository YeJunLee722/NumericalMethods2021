/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [Ye Jun Lee]
Created          : 26-03-2018
Modified         : 18-03-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.cpp
----------------------------------------------------------------*/

#include "myNM.h"


// Matrix addition
Matrix	addMat(Matrix _A, Matrix _B)
{
	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = createMat(_A.rows, _B.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j] = _A.at[i][j] + _B.at[i][j];

	return Out;
}

//Apply forward - substitution
void fwdSub1(Matrix _A, Matrix _b, Matrix _x)
{
	copyMat(_b, _x);

	for (int i = 0; i < _A.rows; i++)
	{
		double temp = 0;
		for (int j = 0; j < i; j++)
			temp += _A.at[i][j] * _x.at[j][0];
		_x.at[i][0] = (_b.at[i][0] - temp) / _A.at[i][i];
	}
	return;
}

// Apply back-substitution
void	backSub1(Matrix _A, Matrix _b, Matrix _x)
{
	copyMat(_b, _x);
	for (int i = _A.rows; i > 0; i--)
	{
		double temp = 0;
		for (int j = i + 1; j <= _A.cols; j++)
			temp += _A.at[i - 1][j - 1] * _x.at[j - 1][0];

		_x.at[i - 1][0] = (_x.at[i - 1][0] - temp) / _A.at[i - 1][i - 1];
	}
	return;
}

void gaussElim1(Matrix _A, Matrix _b, Matrix _U, Matrix _d)
{
	int row = _A.rows;
	int col = _A.cols;
	double* m = (double*)malloc(sizeof(double) * _A.rows);

	for (int a = 0; a < row; a++)                 //처음에 _U와 _d에 _A 와 _d의 값을 저장시킨다.
	{
		_d.at[a][0] = _b.at[a][0];
		for (int b = 0; b < col; b++)
			_U.at[a][b] = _A.at[a][b];
	}
	if (row != col)                                // 정사각행렬이 아닐 경우 종료 시키는 코드
	{
		printf("This matrix is not square matrix");
		exit(0);
	}
	for (int k = 0; k < row - 1; k++)
		for (int i = k + 1; i < row; i++)
		{
			if (_U.at[k][k] == 0)			        // 나누는 수가 0일 경우 종료 시키는 코드
				exit(0);
			else
				m[i] = _U.at[i][k] / _U.at[k][k];
			for (int j = k; j < col; j++)
			{
				_U.at[i][j] = _U.at[i][j] - m[i] * _U.at[k][j];
			}

			_d.at[i][0] = _d.at[i][0] - m[i] * _d.at[k][0];
		}
}


int PP1(Matrix _A, int i)    //i는 0부터 시작 
{
	//printf("i:%d \n", i);
	int a = i;
	int h = i;
	const int o = i;

	double small;
	int large_row;
	double Pivoc;
	int row;
	int _val = 0;
	Matrix partial = createMat(1, _A.cols);
	initMat(partial, -1);
	Matrix sum = createMat(1, _A.cols);
	Matrix M = createMat(1, _A.cols);
	initMat(M, _val);

	for (; i < _A.rows; i++)
	{
		Pivoc = _A.at[i][o];
		if (Pivoc == 0)
		{
			partial.at[0][i] = 0;
			continue;
		}
		for (int j = o + 1; j < _A.cols; j++)
		{
			if (_A.at[i][j] == 0)								//열 중에 이 있으면 그 열의 값은 0으로 계산한다.
				sum.at[0][j] = 99999999999999;
			else
				sum.at[0][j] = fabs(Pivoc / _A.at[i][j]);
		}
		small = sum.at[0][o + 1];

		for (int k = o + 1; k < _A.cols - 1; k++)
		{
			if (small > sum.at[0][k + 1])
				small = sum.at[0][k + 1];
		}
		partial.at[0][i] = small;		       //a11 a12 a13 a14  a11/a12, a11/a13
		//printf("%f %d \n", partial.at[0][i],i);
	}
	double large_p = partial.at[0][o];
	large_row = o;
	for (int d = o; d < _A.rows - 1; d++)
	{
		if (large_p < partial.at[0][d + 1])
		{
			large_row = d + 1;
			large_p = partial.at[0][d + 1];
		}
	}
	//printf("largest row :%d \n", large_row);
	return large_row;
}

void LUdecomp1(Matrix _A, Matrix _L, Matrix _U, Matrix _P)
{
	int row = _A.rows;
	int col = _A.cols;
	int large;
	double* m = (double*)malloc(sizeof(double) * _A.rows);
	copyMat(_A, _U);

	if (row != col)											  // 
	{														  // 정사각행렬이 아닐 경우 종료 시키는 코드	
		printf("This matrix is not square matrix");			  //
		exit(0);
	}
	for (int k = 0; k < row - 1; k++)
	{
		for (int i = k + 1; i < row; i++)
		{
			if (_U.at[k][k] == 0)							   // 나누는 수가 0일 경우 종료 시키는 코드
				exit(0);									   //
			else
			{
				m[i] = _U.at[i][k] / _U.at[k][k];
				_L.at[i][k] = (double)m[i];							   // L을 만들어준다.
				//printMat(_L, "L");
			}
			for (int j = k; j < col; j++)
				_U.at[i][j] = _U.at[i][j] - m[i] * _U.at[k][j];
		}
	}

	//_L = addMat(_L, eye(_A.rows, _A.cols));
	//printMat(_U, "U(upper)");
	//printMat(_L, "L(lower)");
}


void LUdecomp_PIVOT1(Matrix _A, Matrix _L, Matrix _U, Matrix _P)
{
	int row = _A.rows;
	int col = _A.cols;
	int large;
	double* m = (double*)malloc(sizeof(double) * _A.rows);

	copyMat(_A, _U);

	if (row != col)								     // 정사각행렬이 아닐 경우 종료 시키는 코드
	{
		printf("This matrix is not square matrix");
		exit(0);
	}
	for (int k = 0; k < row - 1; k++)
	{
		large = PP1(_U, k);							//PP는 PARTIAL PIVOTING 함수로 S.P가 가장 높은 행을 반환한다.
		Rexchange(_L, k, large);
		Rexchange(_P, k, large);					//Rexchange는 S.P가 가장 높은 행을 PIVOT 행으로 만들어주기위해 행을 바꾸는 함수이다
		Rexchange(_U, k, large);
		for (int i = k + 1; i < row; i++)
		{
			if (_U.at[k][k] == 0)			        // 나누는 수가 0일 경우 종료 시키는 코드
				exit(0);
			else
			{
				m[i] = _U.at[i][k] / _U.at[k][k];
				_L.at[i][k] = m[i];
				//printMat(_L, "L");R
			}
			for (int j = k; j < col; j++)
				_U.at[i][j] = _U.at[i][j] - m[i] * _U.at[k][j];
		}
	}
	printMat(_P, "P");
	printMat(_U, "U(Uppter)");
	printMat(_L, "L(Lower)");
}

//
void solveLU1(Matrix _L, Matrix _U, Matrix _P, Matrix _b, Matrix _x)
{
	Matrix d = createMat(_L.rows, _b.cols);
	Matrix m_x = createMat(_L.rows, _b.cols);
	initMat(d, 0);

	d = Multi(_P, _b);
	copyMat(d, m_x);
	fwdSub1(_L, d, m_x);
	backSub1(_U, m_x, _x);
	printMat(_x, "Solution");
}

void inv1(Matrix _A, Matrix _Ainv)
{
	//LUdecomp(Matrix _A, Matrix _L, Matrix _U, Matrix _P)
	Matrix L = createMat(_A.rows, _A.cols);
	Matrix P = createMat(_A.rows, _A.cols);
	Matrix U = createMat(_A.rows, _A.cols);
	Matrix I = createMat(_A.rows, _A.cols);
	Matrix out1 = createMat(_A.rows, 1);
	Matrix out2 = createMat(_A.rows, 1);
	Matrix d = createMat(_A.rows, 1);

	I = eye(I.rows, I.cols);

	initMat(L, 0);
	initMat(P, 0);
	initMat(U, 0);
	initMat(d, 0);
	initMat(out1, 0);
	initMat(out2, 0);

	if (_A.rows != _A.cols)				//square matrix
	{
		printf("n pos n이 아님 \n");
		exit(0);
	}
	LUdecomp1(_A, L, U, P);
	for (int i = 0; i < _A.rows; i++)		//rank(A) n 인지 check		
		if (U.at[i][i] == 0)
		{
			printf("rank(A) != n");
			exit(0);
		}

	L = addMat(L, I);
	printMat(L, "L+I");

	for (int i = 0; i < _A.rows; i++)
	{
		for (int k = 0; k < _A.rows; k++)
			d.at[k][0] = I.at[k][i];
		for (int j = 0; j < _A.cols; j++)
		{
			fwdSub1(L, d, out1);
			backSub1(U, out1, out2);
		}
		for (int k = 0; k < _A.rows; k++)
			_Ainv.at[i][k] = out2.at[k][0];
	}
	printMat(_Ainv, "Ainv");
}



void QR_1(Matrix _A, Matrix _Q, Matrix _R, Matrix _OUTQ, Matrix _OUTR)
{
	Matrix H = createMat(_A.rows, _A.cols);
	Matrix v = createMat(_A.rows, 1);
	Matrix c = createMat(_A.rows, 1);
	Matrix e = createMat(_A.rows, 1);
	Matrix I = eye(_A.rows, _A.rows);
	Matrix vT = createMat(1, _A.rows);
	Matrix vvT = createMat(_A.rows, _A.cols);
	double vTv;
	int sign;

	if (_A.rows != _A.cols)
	{
		printf("This matrix is not a square matrix, No eigen matrix \n\n\n");
		exit(0);
	}

	copyMat(_A, _R);
	//printMat(_R, "R");
	copyMat(I, _Q);

	for (int k = 0; k < _A.rows - 1; k++)
	{
		initMat(c, 0);

		for (int i = k; i < _A.rows; i++)
			c.at[i][0] = _R.at[i][k];

		//printMat(c, "c");
		//printMat(_R, "R");

		if (c.at[k][0] > 0)
			sign = 1;
		else if (c.at[k][0] < 0)
			sign = -1;

		for (int i = 0; i < _A.rows; i++)
			e.at[i][0] = sign * norm(c) * I.at[i][k];

		//printMat(e, "e");
		v = addMat(c, e);
		//printMat(v, "v");
		//printf("norm:%f \n", norm(c));
		vT = transpose1(v);
		vvT = Multi(v, vT);
		vTv = norm(v) * norm(v);

		for (int i = 0; i < _A.rows; i++)
			for (int j = 0; j < _A.cols; j++)
				H.at[i][j] = I.at[i][j] - (double)2 / vTv * vvT.at[i][j];

		//printMat(H, "H");
		_Q = Multi(_Q, H);
		_R = Multi(H, _R);

		//printMat(_Q, "Q");
		//printMat(_R, "R");
	}
	copyMat(_Q, _OUTQ);
	copyMat(_R, _OUTR);
}

void QR_2(Matrix _A, Matrix _Q, Matrix _R, Matrix _OUTQ, Matrix _OUTR, Matrix _OUTA)
{
	Matrix Q = createMat(_Q.rows, _Q.cols);
	Matrix OUTQ = createMat(_Q.rows, _Q.cols);
	Matrix OUTR = createMat(_Q.rows, _Q.cols);
	Matrix R = createMat(_Q.rows, _Q.cols);
	Matrix A = createMat(_A.rows, _A.cols);

	copyMat(_Q, Q);
	copyMat(_A, A);
	copyMat(_R, R);

	for (int k = 0; k < _A.rows * 50; k++)
	{
		A = Multi(R, Q);
		QR_1(A, Q, R, OUTQ, OUTR);
		copyMat(OUTQ, Q);
		copyMat(OUTR, R);
		//printMat(A, "eigenvalue");
	}
	copyMat(A, _OUTA);
	copyMat(Q, _OUTQ);
	copyMat(R, _OUTR);
}


void eig(Matrix _A, Matrix _OUTA)
{

	Matrix Q = createMat(_A.rows, _A.cols);
	Matrix R = createMat(_A.rows, _A.cols);
	Matrix OUTQ1 = createMat(_A.rows, _A.cols);
	Matrix OUTR1 = createMat(_A.rows, _A.cols);
	Matrix OUTQ2 = createMat(_A.rows, _A.cols);
	Matrix OUTR2 = createMat(_A.rows, _A.cols);
	Matrix EigenA = createMat(_A.rows, _A.cols);

	QR_1(_A, Q, R, OUTQ1, OUTR1);				    //Part1. QRdecomp 
	QR_2(_A, OUTQ1, OUTR1, OUTQ2, OUTR2, EigenA);   //Part2. Making similar matrix

	copyMat(EigenA, _OUTA);
}

double Eigen_Compare(Matrix _A)
{
	Matrix eigen = createMat(_A.rows, 1);
	initMat(eigen, 0);

	for (int i = 0; i < _A.rows; i++)
		eigen.at[i][0] = fabs(_A.at[i][i]);

	printMat(_A, "a");
	double Max = eigen.at[0][0];
	double Min = eigen.at[0][0];
	double result;

	for (int k = 0; k < _A.rows; k++)
	{
		for (int i = 0; i < _A.rows; i++)
		{
			if (Max < eigen.at[i][0])
				Max = eigen.at[i][0];
			if (Min > eigen.at[i][0] && eigen.at[i][0] != 0)
				Min = eigen.at[i][0];
		}
	}
	result = Max / Min;

	if (Min == 0 || Max == 0)
	{
		printf("Eigen value 없음 \n");
		result = 0;
	}

	printf("Maximum Eigenvalue : %f \n\n", Max);
	printf("Minimum Eigenvalue : %f \n\n", Min);

	return result;
}

double cond(Matrix A)
{
	Matrix ATA = createMat(A.cols, A.cols);

	initMat(ATA, 0);
	ATA = Multi(transpose1(A), A);

	Matrix Q = createMat(ATA.rows, ATA.cols);
	Matrix outQ = createMat(ATA.rows, ATA.cols);
	Matrix outR = createMat(ATA.rows, ATA.cols);
	Matrix finalQ = createMat(ATA.rows, ATA.cols);
	Matrix finalR = createMat(ATA.rows, ATA.cols);
	Matrix R = createMat(ATA.rows, ATA.cols);
	Matrix finalA = createMat(ATA.rows, ATA.cols);

	double Condition;

	initMat(R, 0);
	initMat(Q, 0);

	initMat(outR, 0);
	initMat(outQ, 0);

	initMat(finalR, 0);
	initMat(finalQ, 0);

	printMat(ATA, "ATA");

	QR_1(ATA, Q, R, outQ, outR);
	QR_2(ATA, outQ, outR, finalQ, finalR, finalA);

	printMat(finalQ, "Q of ATA");
	printMat(finalR, "R of ATA");
	printMat(finalA, "Eigenvalue of ATA");
	Condition = Eigen_Compare(finalA);

	return sqrt(Condition);

}


Matrix linearInterp(Matrix _x, Matrix _y, Matrix _xq)
{
	int len = _xq.rows;
	int mx = _x.rows;
	int my = _y.rows;

	int flag;

	Matrix OUT = createMat(_xq.rows, 1);
	initMat(OUT, 0);


	if (mx != my)
	{
		printf("어림도 없지 돌아가!\n");
	}
	else
	{
		for (int i = 0; i < len; i++)
		{
			int k = 0;
			double x = 0;
			for (; k < mx - 1; k++)                                            // x data
			{
				if (_xq.at[i][0] > _x.at[k][0] && _xq.at[i][0] < _x.at[k + 1][0])
				{
					printf("x0:%f     x1:%f \n", _x.at[k][0], _x.at[k + 1][0]);
					x = _xq.at[i][0];
					OUT.at[i][0] = _y.at[k][0] * (x - _x.at[k + 1][0]) / (_x.at[k][0] - _x.at[k + 1][0]) + _y.at[k + 1][0] * (x - _x.at[k][0]) / (_x.at[k + 1][0] - _x.at[k][0]);
					break;
				}
				else if (_xq.at[i][0] == _x.at[k][0])
				{
					OUT.at[i][0] = _y.at[k][0];
					break;
				}
				else if (_xq.at[i][0] == _x.at[k + 1][0])
				{
					OUT.at[i][0] = _y.at[k + 1][0];
					break;
				}
			}

		}
	}
	printMat(OUT, "OUT");

	return OUT;
}



// Returns the parameters of the linear least square function.
Matrix	linearFit(Matrix _x, Matrix _y) {
	int mx = _x.rows;
	int my = _y.rows;

	double a1 = 0;
	double a0 = 0;
	Matrix ATA = createMat(2, 2);
	Matrix OUT = createMat(2, 1);
	Matrix ATY = createMat(2, 1);
	Matrix invATA = createMat(2, 2);

	initMat(OUT, 0);
	initMat(invATA, 0);
	initMat(ATY, 0);
	initMat(ATA, 0);

	double Sxx = 0;
	double Sxy = 0;
	double Sx = 0;
	double Sy = 0;
	int    n = 0;

	if (mx != my)
	{
		printf("안돼 돌아가! \n");
		ATA.at[0][0] = 0;
		ATA.at[1][0] = 0;
	}
	else
	{

		for (int i = 0; i < mx; i++)
		{
			Sxx += _x.at[i][0] * _x.at[i][0];
			Sxy += _x.at[i][0] * _y.at[i][0];
			Sx += _x.at[i][0];
			Sy += _y.at[i][0];
			n++;
		}

		ATA.at[0][0] = Sxx;
		ATA.at[0][1] = Sx;
		ATA.at[1][0] = Sx;
		ATA.at[1][1] = n;

		ATY.at[0][0] = Sxy;
		ATY.at[1][0] = Sy;
		printMat(ATY, "ATY");
		printMat(ATA, "ATA");

		inv1(ATA, invATA);
		printMat(invATA, "invATA");
		OUT = Multi(invATA, ATY);
		a1 = OUT.at[0][0];
		a0 = OUT.at[1][0];

	}

	double z_array[] = { a1, a0 };
	return arr2Mat(z_array, 2, 1);
}

// Create a matrix from 1D-array
Matrix	arr2Mat(double* _1Darray, int _rows, int _cols)
{
	Matrix Output = createMat(_rows, _cols);

	for (int i = 0; i < _rows; i++)
		for (int j = 0; j < _cols; j++)
			Output.at[i][j] = _1Darray[i * _cols + j];

	return Output;
}


//
// Differentiate
//

void    gradient1D(double _x[], double _y[], double _dydx[], int m)
{
	if (m >= 2)
	{
		for (int i = 0; i < m; i++)
		{
			if (i == 0)
				_dydx[i] = (-3 * _y[i] + 4 * _y[i + 1] - _y[i + 2]) / (2 * (_x[i + 1] - _x[i]));  // 3 point forward
			else if (0 < i && i < m - 1)
				_dydx[i] = (_y[i + 1] - _y[i - 1]) / ((_x[i + 1] - _x[i - 1]));                         // 2 point central
			else if (i == m - 1)
				_dydx[i] = (3 * _y[i] - 4 * _y[i - 1] + _y[i - 2]) / (2 * (_x[i] - _x[i - 1]));    // 3 point backward
		}
	}
	else if (m == 2)
	{
		int i = 0;
		_dydx[i] = (_y[i + 1] - _y[i]) / (2 * (_x[i + 1] - _x[i]));
		_dydx[i] = _dydx[i];
	}

}
Matrix	gradient(Matrix _x, Matrix _y)
{
	int mx = _x.rows;
	int my = _y.rows;
	Matrix grad = createMat(mx, 1);
	initMat(grad, 0);

	double h;

	if (mx != my)
		printf("어림도 없지 돌아가! \n");
	else
	{
		if (mx >= 2)
		{
			for (int i = 0; i < mx; i++)
			{
				if (i == 0)
					grad.at[i][0] = (-3 * _y.at[i][0] + 4 * _y.at[i + 1][0] - _y.at[i + 2][0]) / (2 * (_x.at[i + 1][0] - _x.at[i][0]));  // 3 point forward
				else if (0 < i && i < mx - 1)
					grad.at[i][0] = (_y.at[i + 1][0] - _y.at[i - 1][0]) / ((_x.at[i + 1][0] - _x.at[i - 1][0]));                         // 2 point central
				else if (i == mx - 1)
					grad.at[i][0] = (3 * _y.at[i][0] - 4 * _y.at[i - 1][0] + _y.at[i - 2][0]) / (2 * (_x.at[i][0] - _x.at[i - 1][0]));    // 3 point backward
			}
		}
		else if (mx == 2)
		{
			int i = 0;
			grad.at[i][0] = (_y.at[i + 1][0] - _y.at[i][0]) / (2 * (_x.at[i + 1][0] - _x.at[i][0]));
			grad.at[i + 1][0] = grad.at[i][0];
		}
	}
	return grad;
}


Matrix	gradientFunc(double func(const double x), Matrix _x)
{
	int mx = _x.rows;
	Matrix _y = createMat(mx, 1);
	Matrix grad = createMat(mx, 1);
	initMat(_y, 0);
	initMat(grad, 0);

	for (int i = 0; i < mx; i++)
		_y.at[i][0] = func(_x.at[i][0]);

	grad = gradient(_x, _y);


	return grad;
}

double newtonRaphsonFunc(double func(const double x), double dfunc(const double x), float _x0, float _tol)
{
	int N = 0;
	int Nmax = 20;

	double ep = _tol;
	double T = _x0;
	double error = 0;
	double H = 0;

	do
	{
		error = fabs(func(T));
		H = (double)func(T) / dfunc(T);
		//printf("%f\n", H);
		T = T - H;
		printf("iteration : %d \t\t\t Tem : %f \t\t\t Tol:%f  \n", N, T, error);
		N++;

	} while (N < Nmax && fabs(error) >= ep);

	return T;
}

//ODE


void odeEU(double func(const double x, const double t), double y[], double t0, double tf, double h)
{
	int iter = tf / h;
	double t = 0;
	printf("odeEU fucntion \n");
	for (int i = 0; i < iter; i++)
	{
		y[i + 1] = y[i] + func(y[i], t) * h;
		t = t + h;

	}

}

void odeEM(double func(const double x, const double t), double y[], double t0, double tf, double h)
{
	int iter = tf / h;
	double x = 0;
	double t = 0;
	double slope1 = 0;
	double slope2 = 0;
	double y_est = 0;

	printf("odeEM fucntion \n");

	for (int i = 0; i < iter; i++)
	{
		slope1 = func(y[i], t);
		y_est = y[i] + slope1 * h;
		t = t + h;
		slope2 = func(y_est, t);
		y[i + 1] = y[i] + (slope1 + slope2) / 2 * h;
	}
}

void ode(double func(const double x, const double t), double y[], double t0, double tf, double h, int method)
{
	if (method == 0)
		odeEU(func, y, t0, tf, h);
	else if (method == 1)
		odeEM(func, y, t0, tf, h);
	else if (method == 2)
		Rk3(func, y, t0, tf, h);
	else if (method == 3)
		Rk4(func, y, t0, tf, h);
}

void Rk3(double func(const double x, const double t), double y[], double t0, double tf, double h)
{
	int iter = tf / h;
	double x = 0;
	double t = 0;
	double slope1 = 0;
	double slope2 = 0;
	double slope3 = 0;

	printf("Rk3 function \n");

	for (int i = 0; i < iter; i++)
	{
		slope1 = func(y[i], t);
		slope2 = func(y[i] + (double)1 / 2 * slope1 * h, t + (double)1 / 2 * h);
		slope3 = func(y[i] - slope1 * h + 2 * slope2 * h, t + h);
		y[i + 1] = y[i] + ((double)1 / 6 * slope1 + (double)2 / 3 * slope2 + (double)1 / 6 * slope3) * h;
		t += h;
	}
}

void Rk4(double func(const double x, const double t), double y[], double t0, double tf, double h)
{
	int iter = tf / h;
	double x = 0;
	double t = 0;
	double slope1 = 0;
	double slope2 = 0;
	double slope3 = 0;
	double slope4 = 0;

	printf("Rk4 function \n");

	for (int i = 0; i < iter; i++)
	{
		slope1 = func(y[i], t);
		slope2 = func(y[i] + (double)1 / 2 * slope1 * h, t + (double)1 / 2 * h);
		slope3 = func(y[i] + (double)1 / 2 * slope2 * h, t + (double)1 / 2 * h);
		slope4 = func(y[i] + slope3 * h, t + h);
		y[i + 1] = y[i] + (slope1 + 2 * slope2 + 2 * slope3 + slope4) * (double)1 / 6 * h;
		t += h;
	}
}


void odeMulti(double yfunc(const double z), double zfunc(const double y, const double z, const double t), double y[][2], double t0, double tf, double h)
{
	double slope_y1 = 0;
	double slope_y2 = 0;
	double slope_z1 = 0;
	double slope_z2 = 0;

	double yest = 0;
	double zest = 0;
	double time = 0;
	int    iter = (tf - t0) / h;
	
	double theta = y[0][0];
	double omega = y[0][1];
	
	for (int i = 0; i < iter; i++)
	{
		slope_y1 = yfunc(omega);             
		slope_z1 = zfunc(theta, omega, time);
		yest = theta + slope_y1 * h;
		zest = omega + slope_z1 * h;
		time += h;
		slope_y2 = yfunc(zest);
		slope_z2 = zfunc(yest, zest, time);

		y[i + 1][0] = theta + (slope_y1 + slope_y2) / 2 * h;
		y[i + 1][1] = omega + (slope_z1 + slope_z2) / 2 * h;
		omega = y[i + 1][1];
		theta = y[i + 1][0];
	}
}

//
// In class
//
double trapz(double _x[], double _y[], int _m)
{
	int N = _m - 1;
	double I = 0;
	double h = 0;

	for (int i = 0; i < N; i++)
	{
		h = (_x[i + 1] - _x[i]);
		I += (_y[i] + _y[i + 1]) *h ;
	}
	I = (double)I / 2;

	return I;
}

double integral(double func(const double x), double _a, double _b, int _n)
{    // simpson 1/3
	int N = _n;                                          //interval 개수를 의미한다
	double fx_odd = 0;
	double fx_even = 0;
	double fx = 0;
	double h = (_b-_a)/_n;                               //간격이야!

	for (int i = 1; i < N; i++)
	{
		double xi = double(_a + h * i);
		printf("xi  %f i   %d \n", xi, i);
		if (i % 2 == 0 && i != 0)
			fx_even += func(xi);
		else if (i % 2 != 0)
			fx_odd += func(xi);
		//printf("%f %f \n", fx_even, fx_odd);
	}

	fx_odd = 4 * fx_odd;
	fx_even = 2 * fx_even;
	fx = h / 3 * (func(_a) + func(_b) + fx_odd + fx_even);

	return fx;
}

double Simson83(double func(const double x), double _a, double _b, int _n) 
{
	// simpson 3/8
	int N = _n-1;                                          //총 개수를 의미한다
	double fx_odd = 0;
	double fx_even = 0;
	double fx = 0;
	double h = (_b - _a) / N;                      //간격이야!

	for (int i = 1; i < N; i++)
	{
		double xi = double(_a + h * i);
		printf("xi  %f i   %d \n", xi, i);
		if (i % 3 == 1 )
			fx_even += func(xi)+ func(xi+h);
		else if (i % 3 == 0 && i != 0)
			fx_odd += func(xi);
		//printf("%f %f \n", fx_even, fx_odd);
	}

	fx_odd = 2 * fx_odd;
	fx_even = 3 * fx_even;
	fx = h *3/ 8 * (func(_a) + func(_b) + fx_odd + fx_even);

	return fx;
}
// Integration using rectangular method for discrete data inputs
double IntegrateRect(double _x[], double _y[], int _m)
{
	int N = _m - 1;
	double I = 0;

	for (int i = 0; i < N; i++)
		I += _y[i] * (_x[i + 1] - _x[i]);

	return I;
}


void odeRK2(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
{
	double C1 = 0.5;
	double C2 = 0.5;
	double alpha = 1;
	double beta = alpha;  // alpha=beta

	int N = (tf - t0) / h + 1;
	double ti = t0;
	double y_EU;
	double K1 = 0, K2 = 0;

	// Initialization 
	y[0] = y0;

	for (int i = 0; i < N - 1; i++)
	{
		// first slope  
		K1 = odeFunc(ti, y[i]);

		// Second slope  
		y_EU = y[i] + beta * K1 * h;
		K2 = odeFunc(ti + h * alpha, y_EU);

		// Update 
		y[i + 1] = y[i] + (C1 * K1 + C2 * K2) * h;
		ti += h;
	}
}

void odeRK4(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
{
	double a = 0.5;
	double b = 1;

	int N = (tf - t0) / h + 1;
	double ti = t0;
	double y_EU;
	double K1 = 0, K2 = 0, K3 = 0, K4 = 0;

	// Initialization 
	y[0] = y0;

	for (int i = 0; i < N - 1; i++)
	{
		// first slope K1 
		K1 = odeFunc(ti, y[i]);

		// Second slope K2
		K2 = odeFunc(ti + h * a, y[i] + a * K1 * h);

		// Third slope K3
		K3 = odeFunc(ti + h * a, y[i] + a * K2 * h);

		// Fourth slope K4
		K4 = odeFunc(ti + h, y[i] + K3 * h);

		// Update 
		y[i + 1] = y[i] + (1.0 / 6.0) * (K1 + 2 * K2 + 2 * K3 + K4) * h;

		ti += h;
	}

}

// ODE RK2:  one of 2nd order ODE <--> two of 1st order ODE
void sys2RK2(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init)
{
	int N = (tf - t0) / h + 1;
	double ti = t0;

	double K1[2] = { 0 };  // K1 = {K1_y1, K1_y2}
	double K2[2] = { 0 };
	double Yin[2] = { 0 };
	double K1_y1 = 0, K1_y2 = 0, K2_y1 = 0, K2_y2 = 0;

	//Y[0]==y(t)
	//Y[1]==z(y)
	//dydt[0] == z(t) == Y[1]
	//dydt[1] =?

	// Initial condition
	y1[0] = y1_init;    //y(t)
	y2[0] = y2_init;    //z(t) = dydt(t)
	//Yin[0] = y1[0];
	//Yin[1] = y2[0];

	for (int i = 0; i < N - 1; i++) 
	{
		// Slope 1 : K1
 
		odeFunc_sys2(ti,Yin , K1) ;

		K1_y1 = K1[0] ;
		K1_y2 = K1[1] ;

		// Slope 2 : K2
		Yin[0] = y1[i] + K1_y1 * h;      // yest
		Yin[1] = y2[i] + K1_y2 * h;      // zest      
		ti += h;
		odeFunc_sys2(ti, Yin, K2);
		K2_y1 = K2[0];
		K2_y2 = K2[1];

		// Update
		y1[i + 1] = y1[i] + (double)(K1_y1 + K2_y1)/ 2 * h;
		y2[i + 1] = y2[i] + (double)(K1_y2 + K2_y2)/ 2 * h;
		printf("\n\n%f %f \n\n", y1[i], y2[i]);
		
	}
}


// Classical RK4
void sys2RK4(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init)
{
	int N = (tf - t0) / h + 1;
	double ti = t0;

	double K1[2] = { 0 };  // K1 = {K1_y1, K1_y2}
	double K2[2] = { 0 };
	double K3[2] = { 0 };
	double K4[2] = { 0 };
	double Yin[2] = { 0 };
	double K1_y1 = 0, K1_y2 = 0, K2_y1 = 0, K2_y2 = 0, K3_y1 = 0, K3_y2 = 0, K4_y1 = 0, K4_y2 = 0;

	// Initial condition
	y1[0] = y1_init;    //y(t)
	y2[0] = y2_init;    //z(t) = dydt(t)

	for (int i = 0; i < N - 1; i++) {

		// Slope 1 : K1
		Yin[0] = y1[i];      // z
		Yin[1] = y2[i];      // dzdt      
		odeFunc_sys2(ti, Yin, K1);
		K1_y1 = K1[0];
		K1_y2 = K1[1];

		// Slope 2 : K2
		Yin[0] = y1[i] + 0.5 * K1_y1 * h;      // z
		Yin[1] = y2[i] + 0.5 * K1_y2 * h;      // dzdt      
		odeFunc_sys2(ti + 0.5 * h, Yin, K2);
		K2_y1 = K2[0];
		K2_y2 = K2[1];

		// Slope 3 : K3
		Yin[0] = y1[i] + 0.5 * K2_y1 * h;      // z
		Yin[1] = y2[i] + 0.5 * K2_y2 * h;      // dzdt      
		odeFunc_sys2(ti + 0.5 * h, Yin, K3);
		K3_y1 = K3[0];
		K3_y2 = K3[1];

		// Slope 4 : K4
		Yin[0] = y1[i] + K3_y1 * h;      // z
		Yin[1] = y2[i] + K3_y2 * h;      // dzdt      
		odeFunc_sys2(ti + h, Yin, K4);
		K4_y1 = K4[0];
		K4_y2 = K4[1];

		// Update
		y1[i + 1] = y1[i] + (K1_y1 + 2 * K2_y1 + 2 * K3_y1 + K4_y1) * h / 6;
		y2[i + 1] = y2[i] + (K1_y2 + 2 * K2_y2 + 2 * K3_y2 + K4_y2) * h / 6;
		ti += h;
	}
}
//Lagrange Cubic
Matrix Lagrange_cubic(double x[], double y[], double test_x[])
{

	int mx = sizeof(x) / sizeof(double);
	int my = sizeof(y) / sizeof(double);
	int n = my;
	int tx = sizeof(test_x) / sizeof(double);

	Matrix A    = createMat(n, n);
	Matrix Ainv = createMat(n, n);

	Matrix h = createMat(n - 2, 1);
	Matrix u = createMat(n - 2, 1);
	Matrix h = createMat(n - 2, 1);
	Matrix b = createMat(n - 2, 1);
	Matrix v = createMat(n - 2, 1);
	Matrix v_mat = createMat(n, n);
	Matrix a = createMat(n, 1);

	initMat(h, 0);
	initMat(A, 0);
	initMat(u, 0);
	initMat(b, 0);
	initMat(v, 0);
	initMat(v_mat, 0);
	initMat(a, 0);
	initMat(Ainv, 0);

	A.at[0][0] = 1;

	for (int i = 0; i < n - 2; i++)
		h.at[i][0] = x[i + 1] - x[i];
	
	for (int i = 1; i < n - 1; i++)
		u.at[i][0] = 2 * (h.at[i][0] + h.at[i - 1][0]);

	for (int i = 0; i < n - 3; i++)
		b.at[i][0] = (y[i + 1] - y[i]) / h.at[i][0];

	for (int i = 1; i < n - 1; i++)
		v.at[i][0] = 6 * (b.at[i][0] - b.at[i - 1][0]);
	
	for (int i = 1; i < n-1; i++)
	{
		A.at[i][i - 1] = h.at[i - 1][0];
		A.at[i][i]     = u.at[i][0];
		A.at[i][i + 1] = h.at[i][0];
	}
	A.at[n - 1][n - 1] = 1;

	for (int i = 1; i < n - 1; i++)
		v_mat.at[i][0] = v.at[i][0];

	inv1(A, Ainv);
	a = Multi(Ainv, v_mat);
	
	for (int i = 0; i < tx; i++)
	{
		int k = 0;

		for (; k < n-1; k++)
		{
			if (x[k] <= test_x[i] && x[k + 1] > test_x[i])
				break;
		}
		a.at[k][0] / (6 * h.at[k][0]) * pow((x[k + 1] - test_x[i]), 3) + a.at[k + 1][0] / (6 * h.at[k][0]) * pow((x[k + 1] - test_x[i]), 3)+ (y[k + 1] / h.at[k][0] - a.at[k + 1][0] * h.at[k][0] / 6)*(test_x - x[k]) + (y[k] / h.at[k][0] - a.at[k][0] * h.at[k][0] / 6) * (x[k + 1] - test_x[k]);
	}
	
}
/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Jan Park
Created          : 10-05-2021
Modified         : 10-05-2021
Language/ver     : C++ in MSVS2019

Description      : [Tutorial]Differentiation_student.cpp
-------------------------------------------------------------------------------*/

#include "myNM.h"

#define Assignment	8		// enter your assignment number
#define eval		0		// set 0
#define PI          3.141592
#define Eu          0
#define Em          1
#define Rk3         2
#define Rk4         3

// Define a function that defines the target equation.
double myFunc(const double x, const double t);


//ode Function
double yFunc(const double z);
double zFunc(const double y, const double z, double const t);
void odeFunc_mck(const double t, const double Y[], double dYdt[]);

//integral function
double int_Func(const double x);


int main(int argc, char* argv[])
{
	std::string path = "C:/NM_data_2021/Assignment" + std::to_string(Assignment) + "/";

#if eval
	path += "eval/";
#endif

	//// PART 1
	//printf("\n**************************************************");
	//printf("\n|                     ode                         |");
	//printf("\n**************************************************\n");


	//double h = 0.01;
	double t0 = 0;
	double tf = 1;

	//double y[101] = { 0,0 };

	//ode(myFunc, y, t0, tf, h, 3);


	//for (int i = 0; i < 101; i++)
	//	printf("%f \n", y[i]);

	// PART 2
	printf("\n**************************************************");
	printf("\n|               Multi   ode                       |");
	printf("\n**************************************************\n");

	//double Y[101][2] = { { 0.0 },{0.0} };
	//double Yd[101] = { 0.0 };
	//double Vd[101] = { 0.0 };
	//Yd[0] = PI / 2;
	//Vd[0] = 0.2;

	//Y[0][0] = 0;
	//Y[0][1] = 0.2;
	//odeMulti(yFunc, zFunc, Y, t0, tf, h);
	//sys2RK2(odeFunc_mck, Yd, Vd, t0, tf, h, Yd[0], Vd[0]);

	//for (int i = 0; i < 101; i++)
	//{
	//	for (int j = 0; j < 2; j++)
	//		printf("%f  ", Y[i][j]);
	//	printf("\n");
	//}

	//for (int i = 0; i < 101; i++)
	//{
	//		printf("%f  %f \n", Yd[i], Vd[i]);
	//}
	//// PART 3
	//printf("\n**************************************************");
	//printf("\n|               Curve Fitting                    |");
	//printf("\n**************************************************\n");


	//double X[101][2] = { {0.0},{0.0} };
	//double dXdt[101][2] = { {0.0},{0.0} };

		 //PART 2
	//printf("\n**************************************************");
	//printf("\n|              Integral                      |");
	//printf("\n**************************************************\n");
	//
	//double a = 0;
	//double b = 1;
	//double N = 100;             //iteration
	//double h = (b - a) / 100;
	//double x[100] = { 0.0 };
	//double y[100] = { 0.0 };
	//double I = 0;
	//double S83 = 0;
	//for (int i = 0; i < N + 1; i++)
	//{
	//	x[i] = a + h * i;
	//	y[i] = int_Func(x[i]);
	//	printf("%f %f \n", x[i], y[i]);
	//}
	//I = trapz(x, y, 101);
	//S83 = Simson83(int_Func, a, b, 101);

	//printf("%f \n", I);
	//printf("%f \n", S83);
	//system("pause");
	
	printf("\n**************************************************");
	printf("\n|              Cubic spline                      |");
	printf("\n**************************************************\n");

	double x[] = { 0, 1, 2.5, 3.6, 5, 7 ,8.1, 10 };
	double y[] = {0 ,   0.8415,    0.5985, -0.4425 ,-0.9589,0.6570,    0.9699, -0.5440 };
	Matrix OUT = 0;

	return 0;

}




double yFunc(const double z)
{
	return z;
}

double zFunc(const double y, const double z, const double t)
{/*
	double c = 0.16;
	double m = 0.5;
	double L = 1.2;
	double g = 9.8;
	return -c / m * z - g / L * sin(y);*/

	double m = 1;
	double c = 7;
	double k = 6.9;
	double f = 5;
	double Fin = 2 * cos(2 * PI * f * t);
	
	return 1 / m * (-k * y - c * z + Fin);
	
}

// Define a function that defines the target equation.
double myFunc(const double _x, const double t)
{
	double tau = 1;
	double  T = 1 / tau;
	double  f = 10;
	double  w = 2 * PI * f;


	double Vm = 1;

	return  -


		T * _x + T * Vm * cos(w * t);
}

void odeFunc_mck(const double t, const double Y[], double dYdt[])
{
	double m = 1;
	double c = 7;
	double k = 6.9;
	double f = 5;
	double Fin = 2 * cos(2 * PI * f * t);
	/*double c = 0.16;
	double m = 0.5;
	double L = 1.2;
	double g = 9.8;*/

	dYdt[0] = Y[1];
	// EXERCISE: MODIFY HERE
	dYdt[1] = 1 / m * (-k * Y[0] - c * Y[1] + Fin);

}

// integral fucntion

double int_Func(const double x)
{
	double w  = pow(10,4);
	double L  = 1;
	double E  = 200*pow(10,9);
	double	I = 2.1*pow(10,-4);

	return pow((-w * (pow(x, 2) / 2)),2) / (2 * E * I);

}
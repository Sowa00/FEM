#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <crtdbg.h>
#include <iomanip>

using namespace std;

double f(double x) //funkcja 1d
{
	//return (2 * x * x + 3 * x - 8);
	return (x * x - 3 * x + 6);
}

double f(double x, double y) //funkcja 2d
{
	return (-5 * x * x * y + 2 * x * y * y + 10);
}

double kwadratury1do1(int n, double* pc, double* w) // liczenie calki od -1 do 1 1d
{
	double suma = 0;
	for (int i = 0; i < n + 1; i++)
	{
		suma += f(pc[i]) * w[i];
	}
	return suma;
}

double kwadratury1do12d(double n, double* pc, double* w) //liczenie calki od -1 do 1 2d
{
	double suma = 0;
	for (int i = 0; i < n + 1; i++)
	{
		for (int j = 0; j < n + 1; j++)
		{
			suma += f(pc[i], pc[j]) * w[i] * w[j];
		}
	}
	return suma;
}

double kwadratury(double x1, double x2, int pktc, double* pc, double* w) // liczenie calki od x1 do x2 1d
{
	double suma = 0;
	double detJ = (x2 - x1) / 2;
	for (int i = 0; i < pktc; i++)
	{
		pc[i] = ((1 - pc[i]) / 2) * x1 + ((pc[i] + 1) / 2) * x2;
		//cout << pc[i] << endl;
	}
	for (int i = 0; i < pktc + 1; i++)
	{
		suma += (f(pc[i]) * w[i]) * detJ;
		//cout << suma << endl;
	}
	return suma;
}
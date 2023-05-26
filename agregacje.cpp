#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <crtdbg.h>
#include <iomanip>

#include "czytanieZPliku.h"

using namespace std;

double* agregacjaWektorP(double** wektorP, grid grid1)
{
	double* agregacjaTab = new double[grid1.nn];
	for (int i = 0; i < grid1.nn; i++)
	{
		agregacjaTab[i] = 0;
	}
	for (int i = 0; i < grid1.ne; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			agregacjaTab[grid1.el[i].id[j] - 1] += wektorP[i][j];
		}
	}

	/*cout << "Agregacja wektora P" << endl;
	for (int i = 0; i < grid1.nn; i++)
	{
		cout << agregacjaTab[i] << " ";
	}*/

	return agregacjaTab;
}

double** agregacja(double*** macierz, grid grid1)
{
	double** agregacjaTab = new double* [grid1.nn];
	for (int i = 0; i < grid1.nn; i++)
	{
		agregacjaTab[i] = new double[grid1.nn];
	}
	for (int i = 0; i < grid1.nn; i++)
	{
		for (int j = 0; j < grid1.nn; j++)
		{
			agregacjaTab[i][j] = 0;
		}
	}
	for (int i = 0; i < grid1.ne; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int q = 0; q < 4; q++)
			{
				agregacjaTab[grid1.el[i].id[j] - 1][grid1.el[i].id[q] - 1] += macierz[i][j][q];
			}
		}
	}

	////cout << "Agregacja macierzy" << endl;
	//for (int i = 0; i < grid1.nn; i++)
	//{
	//	for (int j = 0; j < grid1.nn; j++)
	//	{
	//		cout << agregacjaTab[i][j] << " ";
	//	}
	//	cout << endl;
	//}

	return agregacjaTab;
}
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <crtdbg.h>
#include <iomanip>
#include <vector>

#include "czytanieZPliku.h"
#include "kwadratury.h"
#include "macierzH.h"
#include "agregacje.h"
#include "BokHBCwektorP.h"
#include "macierzCRownania.h"

using namespace std;

int main()
{
	cout.precision(9);
	int N = 2; //punkty calkowania
	//wczytywanie z pliku
	string file = "Test1_4_4.txt";
	//string file = "Test2_4_4_MixGrid.txt";
	//string file = "Test3_31_31_kwadrat.txt";
	//showfile(file);
	int nn = node_number(file);
	int ne = el_number(file);
	grid grid1(nn, ne);	
	GlobalData data;
	data = load_from_file_data(file);
	grid1 = load_from_file_grid(file);
	//show_grid(grid1);
	//show_globaldata(data);

	double pc[4][5] = { { -(double)(1 / sqrt(3)), (double)(1 / sqrt(3)) },
						{ -(double)(sqrt(15) * 0.2), 0, (double)(sqrt(15) * 0.2) } ,
						{ -(double)(sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0))), -(double)(sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0))), (double)(sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0))), (double)(sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)))} ,
						{ -(double)((1.0 / 3.0) * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0))), -(double)((1.0 / 3.0) * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0))), 0, (double)((1.0 / 3.0) * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0))), (double)((1.0 / 3.0) * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)))} };
	double w[4][5] = {  { 1, 1 } ,
						{ (5.0 / 9.0), (8.0 / 9.0), (5.0 / 9.0) } ,
						{ (double)((18.0 - sqrt(30.0)) / 36.0), (double)((18.0 + sqrt(30.0)) / 36.0), (double)((18.0 + sqrt(30.0)) / 36.0), (double)((18.0 - sqrt(30.0)) / 36.0)} ,
						{ (double)((322.0 - 13.0 * sqrt(70.0)) / 900.0), (double)((322.0 + 13.0 * sqrt(70.0)) / 900.0), (double)(128.0 / 225.0), (double)((322.0 + 13.0 * sqrt(70.0)) / 900.0), (double)((322.0 - 13.0 * sqrt(70.0)) / 900.0)} };


	/*int n = 4;
	double funkcja1d;
	funkcja1d = kwadratury1do1(n, pc[n - 1], w[n - 1]);
	double funkcja2d;
	funkcja2d = kwadratury1do12d(n, pc[n - 1], w[n - 1]);
	cout << setprecision(9) << funkcja1d << endl << setprecision(9) << funkcja2d << endl;

	int pktc = 3;
	double fod2do10;
	fod2do10 = kwadratury(3, 8, 3, pc[pktc-2], w[pktc-2]);
	cout << setprecision(15) << fod2do10;*/

	double eta[3][16] = {{pc[0][0], pc[0][0], pc[0][1], pc[0][1]},
						  {pc[1][0], pc[1][0], pc[1][0], pc[1][1], pc[1][1], pc[1][1], pc[1][2], pc[1][2], pc[1][2]},
						  {pc[2][0], pc[2][0], pc[2][0], pc[2][0], pc[2][1], pc[2][1], pc[2][1], pc[2][1], pc[2][2], pc[2][2], pc[2][2], pc[2][2], pc[2][3], pc[2][3], pc[2][3], pc[2][3]} };
	double ksi[3][16] = {{pc[0][0], pc[0][1], pc[0][0], pc[0][1]},
						  {pc[1][0], pc[1][1], pc[1][2], pc[1][0], pc[1][1], pc[1][2], pc[1][0], pc[1][1], pc[1][2]},
						  {pc[2][0], pc[2][1], pc[2][2], pc[2][3], pc[2][0], pc[2][1], pc[2][2], pc[2][3], pc[2][0], pc[2][1], pc[2][2], pc[2][3], pc[2][0], pc[2][1], pc[2][2], pc[2][3]} };

	//double eta9[9] = {pc[1][0], pc[1][0], pc[1][0], pc[1][1], pc[1][1], pc[1][1], pc[1][2], pc[1][2], pc[1][2]};
	//double ksi9[9] = {pc[1][0], pc[1][1], pc[1][2], pc[1][0], pc[1][1], pc[1][2], pc[1][0], pc[1][1], pc[1][2]};

	//double eta16[16] = {pc[2][0], pc[2][0], pc[2][0], pc[2][0], pc[2][1], pc[2][1], pc[2][1], pc[2][1], pc[2][2], pc[2][2], pc[2][2], pc[2][2], pc[2][3], pc[2][3], pc[2][3], pc[2][3]};
	//double ksi16[16] = {pc[2][0], pc[2][1], pc[2][2], pc[2][3], pc[2][0], pc[2][1], pc[2][2], pc[2][3], pc[2][0], pc[2][1], pc[2][2], pc[2][3], pc[2][0], pc[2][1], pc[2][2], pc[2][3]};

	/*double eta4[4] = { -(double)(1 / sqrt(3)), -(double)(1 / sqrt(3)), (double)(1 / sqrt(3)), (double)(1 / sqrt(3)) };
	double ksi4[4] = { -(double)(1 / sqrt(3)), (double)(1 / sqrt(3)), -(double)(1 / sqrt(3)), (double)(1 / sqrt(3)) };*/

	/*double eta9[9] = { -(double)(sqrt(15) * 0.2), -(double)(sqrt(15) * 0.2), -(double)(sqrt(15) * 0.2), 0, 0, 0, (double)(sqrt(15) * 0.2), (double)(sqrt(15) * 0.2), (double)(sqrt(15) * 0.2) };
	double ksi9[9] = { -(double)(sqrt(15) * 0.2), 0, (double)(sqrt(15) * 0.2), -(double)(sqrt(15) * 0.2), 0, (double)(sqrt(15) * 0.2), -(double)(sqrt(15) * 0.2), 0, (double)(sqrt(15) * 0.2) };*/

	/*ElUniversal el4(2);
	el4.pochodne(eta4, ksi4);
	el4.showtab();
	cout << endl;
	ElUniversal el9(3);
	el9.pochodne(eta9, ksi9);
	el9.showtab();
	ElUniversal el16(4);
	el16.pochodne(eta16, ksi16);
	el16.showtab();*/
	//el.zad(eta4);
	//macierzHObliczanieplik(eta16, ksi16);
	/*ElUniversal el4(N);
	el4.funkcjeKsztaltu(eta[N-2], ksi[N-2]);
	el4.showtabn();*/
	

	//macierzHObliczanieplik(pc[N-2], w[N-2], grid1, data1, N);
	//macierzHObliczanieplik(pc[N - 2], w[N - 2], grid1, data1, N);

	double*** macierzH = macierzHObliczanieplik(eta[N-2], ksi[N-2], pc[N - 2], w[N - 2], grid1, data, N);
	double*** macierzC = macierzCObliczanie(eta[N - 2], ksi[N - 2], pc[N - 2], w[N - 2], grid1, data, N);
	/*for (int i = 0; i < 9; i++)
	{
		cout << "element: " << i<<endl;
		for (int j = 0; j < 4; j++)
		{
			for (int q = 0; q < 4; q++)
			{
				cout << macierzH[i][j][q] << " ";
			}
			cout << endl;
		}
	}
	cout << endl;*/
	//agregacja(macierzH, grid1);
	//cout << endl;
	double testElX[4] = {0, 0.025, 0.025, 0};
	double testElY[4] = {0, 0, 0.025, 0.025};

	Bok **boki = new Bok*[grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		boki[i] = new Bok[4];
	}
	for (int i = 0; i < grid1.ne; i++)
	{
		//cout << "wartosci pc dla elementu: " << i+1 << endl;
		for (int j = 0; j < 4; j++)
		{
			Bok b(N);
			boki[i][j] = b;
			boki[i][j].liczWartosci(grid1, pc[N - 2], w[N - 2], i, j);
			/*cout << "wartosci pc dla boku: " << j+1 << endl;
			if (j == 0)
			{
				cout << " DOLNY" << endl;
			}
			else if (j == 1)
			{
				cout << " PRAWY" << endl;
			}
			else if (j == 2)
			{
				cout << " GORNY" << endl;
			}
			else
			{
				cout << " LEWY" << endl;
			}*/
			/*for (int r = 0; r < N; r++)
			{
				for (int q = 0; q < 4; q++)
				{
					cout << boki[i][j].wartosciPc[r][q] << " ";
				}
				cout << endl;
			}
			cout << endl;*/
		}
		//cout << endl;
	}

	double ***macierzHBC = macierzHBCObliczanieplik(grid1, data, boki, N, w[N - 2]);
	double **wektorP = obliczanieWektoraP(grid1, data, boki, N, w[N - 2]);
	
	double* agregowanyWektorP = agregacjaWektorP(wektorP, grid1);
	
	double*** macierzHHBC = new double** [grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		macierzHHBC[i] = new double* [4];
		for (int j = 0; j < 4; j++)
		{
			macierzHHBC[i][j] = new double[4];
		}
	}

	for (int i = 0; i < grid1.ne; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int q = 0; q < 4; q++)
			{
				macierzHHBC[i][j][q] = macierzH[i][j][q] + macierzHBC[i][j][q];
			}
		}
	}
	/*cout << "Macierz H + HBC:\n";
	for (int i = 0; i < grid1.ne; i++)
	{
		cout << "Element: " << i + 1 << endl;
		for (int j = 0; j < 4; j++)
		{
			for (int q = 0; q < 4; q++)
			{
				cout << macierzHHBC[i][j][q] << " ";
			}
			cout << endl;
		}
	}*/

	double* t0 = new double[grid1.nn];
	for (int i = 0; i < grid1.nn; i++)
	{
		t0[i] = grid1.nd[i].t;
	}
	//cout << "\nAgregacja macierzy H:\n";
	double** agregowanaMacierzH = agregacja(macierzH, grid1);
	//cout << "\nAgregacja macierzy H+HBC:\n";
	double** agregowanaMacierzHHBC = agregacja(macierzHHBC, grid1);
	//cout << "\nAgregacja macierzy C:\n";
	double** agregowanaMacierzC = agregacja(macierzC, grid1);
	
	//rownania(agregowanaMacierzC, agregowanaMacierzHHBC, agregowanyWektorP, grid1, data, t0, data.steptime);
	minMaxTempInEachStep(agregowanaMacierzC, agregowanaMacierzHHBC, agregowanyWektorP, grid1, data, t0, data.steptime);

	delete macierzHHBC;
	delete macierzHBC;
	delete macierzH;
}
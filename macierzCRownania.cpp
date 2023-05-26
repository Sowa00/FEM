#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <crtdbg.h>
#include <iomanip>

#include "czytanieZPliku.h"

double*** macierzCObliczanie(double* Eta, double* Ksi, double* pctab, double* wtab, grid grid1, GlobalData data, int pc) //globalna
{
	ElUniversal el4(pc);
	el4.pochodne(Eta, Ksi);
	el4.funkcjeKsztaltu(Eta, Ksi);
	//el4.showtabn();
	//el4.showtabdn();


	double**** macierzJakobianpc = new double*** [grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		macierzJakobianpc[i] = new double** [pow(pc, 2)];
		for (int j = 0; j < pow(pc, 2); j++)
		{
			macierzJakobianpc[i][j] = new double* [2];
			for (int q = 0; q < 2; q++)
			{
				macierzJakobianpc[i][j][q] = new double[2];
			}
		}
	}
	for (int i = 0; i < grid1.ne; i++)
	{
		for (int j = 0; j < pow(pc, 2); j++)
		{
			for (int q = 0; q < 2; q++)
			{
				for (int r = 0; r < 2; r++)
				{
					macierzJakobianpc[i][j][q][r] = 0;
				}
			}
		}
	}

	cout << endl;

	//Liczenie macierzy jakobiego
	for (int q = 0; q < grid1.ne; q++)
	{
		//cout << "nowy element " << q+1 << endl;
		for (int i = 0; i < pow(pc, 2); i++)
		{
			//cout << "nowy pc " << i+1 << endl;
			for (int j = 0; j < 4; j++)
			{
				macierzJakobianpc[q][i][0][0] += el4.dnksi[i][j] * grid1.nd[grid1.el[q].id[j] - 1].x;
				macierzJakobianpc[q][i][0][1] += el4.dnksi[i][j] * grid1.nd[grid1.el[q].id[j] - 1].y;
				macierzJakobianpc[q][i][1][0] += el4.dneta[i][j] * grid1.nd[grid1.el[q].id[j] - 1].x;
				macierzJakobianpc[q][i][1][1] += el4.dneta[i][j] * grid1.nd[grid1.el[q].id[j] - 1].y;
			}
		}
	}

	/*for (int r = 0; r < grid1.ne; r++)
	{
		cout << "\nelement: " << r+1 << endl;
		for (int q = 0; q < pow(pc, 2); q++)
		{

			cout << "\nmacierz jakobiego pc" << q + 1 << ":" << endl;
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					cout << macierzJakobianpc[r][q][i][j] << " ";
				}
				cout << endl;
			}
		}
	}*/

	double** detjpc = new double* [grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		detjpc[i] = new double[pow(pc, 2)];
	}

	//liczenie wyznacznika jakobiego
	for (int i = 0; i < grid1.ne; i++)
	{
		for (int j = 0; j < pow(pc, 2); j++)
		{
			detjpc[i][j] = (macierzJakobianpc[i][j][0][0] * macierzJakobianpc[i][j][1][1] - macierzJakobianpc[i][j][0][1] * macierzJakobianpc[i][j][1][0]);
		}
	}

	double**** macierzCpc = new double*** [grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		macierzCpc[i] = new double** [pow(pc, 2)];
		for (int j = 0; j < pow(pc, 2); j++)
		{
			macierzCpc[i][j] = new double* [4];
			for (int q = 0; q < 4; q++)
			{
				macierzCpc[i][j][q] = new double[4];
			}
		}
	}

	for (int k = 0; k < grid1.ne; k++)
	{

		for (int j = 0; j < 4; j++)
		{
			for (int q = 0; q < 4; q++)
			{
				for (int r = 0; r < pow(pc, 2); r++)
				{
					macierzCpc[k][r][j][q] = el4.Ntab[r][j] * el4.Ntab[r][q];
				}
			}
		}
		//cout << endl;
	}


	

	for (int i = 0; i < grid1.ne; i++)
	{
		//cout << "element: " << i << endl;
		for (int j = 0; j < pow(pc, 2); j++)
		{
			//cout << "pc: " << j << endl;
			for (int q = 0; q < 4; q++)
			{
				for (int r = 0; r < 4; r++)
				{
					macierzCpc[i][j][q][r] = macierzCpc[i][j][q][r] * data.specific_heat * data.density * detjpc[i][j];

				}
			}
		}
	}

	/*cout << "macierzc1pc" << endl;
	for (int i = 0; i < grid1.ne; i++)
	{
		cout << "element: " << i << endl;
		for (int j = 0; j < pow(pc, 2); j++)
		{
			cout << "pc: " << j << endl;
			for (int q = 0; q < 4; q++)
			{
				for (int r = 0; r < 4; r++)
				{
					cout << macierzCpc[i][j][q][r] << " ";
				}
				cout << endl;
			}
		}
	}*/

	double*** macierzC = new double** [grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		macierzC[i] = new double* [pow(pc, 2)];
		for (int j = 0; j < pow(pc, 2); j++)
		{
			macierzC[i][j] = new double[4];
		}
	}
	for (int i = 0; i < grid1.ne; i++)
	{
		for (int j = 0; j < pow(pc, 2); j++)
		{
			for (int q = 0; q < 4; q++)
			{
				macierzC[i][j][q] = 0;
			}
		}
	}
	//liczenie macierzy C
	if (pc == 2)
	{
		for (int q = 0; q < grid1.ne; q++)
		{
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					macierzC[q][i][j] += macierzCpc[q][0][i][j] * wtab[0] * wtab[0] + macierzCpc[q][1][i][j] * wtab[1] * wtab[0] +
						macierzCpc[q][2][i][j] * wtab[0] * wtab[1] + macierzCpc[q][3][i][j] * wtab[1] * wtab[1];
				}
			}
		}

	}
	else if (pc == 3)
	{
		for (int r = 0; r < grid1.ne; r++)
		{
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					macierzC[r][i][j] = macierzCpc[r][0][i][j] * wtab[0] * wtab[0] + macierzCpc[r][1][i][j] * wtab[1] * wtab[0] + macierzCpc[r][2][i][j] * wtab[2] * wtab[0] +
										macierzCpc[r][3][i][j] * wtab[0] * wtab[1] + macierzCpc[r][4][i][j] * wtab[1] * wtab[1] + macierzCpc[r][5][i][j] * wtab[2] * wtab[1] +
										macierzCpc[r][6][i][j] * wtab[0] * wtab[2] + macierzCpc[r][7][i][j] * wtab[2] * wtab[1] + macierzCpc[r][8][i][j] * wtab[2] * wtab[2];
				}
			}
		}
	}
	else if (pc == 4)
	{
		for (int r = 0; r < grid1.ne; r++)
		{
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					macierzC[r][i][j] = macierzCpc[r][0][i][j] * wtab[0] * wtab[0] + macierzCpc[r][1][i][j] * wtab[1] * wtab[0] + macierzCpc[r][2][i][j] * wtab[2] * wtab[0] + macierzCpc[r][3][i][j] * wtab[3] * wtab[0] +
						macierzCpc[r][4][i][j] * wtab[0] * wtab[1] + macierzCpc[r][5][i][j] * wtab[1] * wtab[1] + macierzCpc[r][6][i][j] * wtab[2] * wtab[1] + macierzCpc[r][7][i][j] * wtab[3] * wtab[1] +
						macierzCpc[r][8][i][j] * wtab[0] * wtab[2] + macierzCpc[r][9][i][j] * wtab[1] * wtab[2] + macierzCpc[r][10][i][j] * wtab[2] * wtab[2] + macierzCpc[r][11][i][j] * wtab[3] * wtab[2] +
						macierzCpc[r][12][i][j] * wtab[0] * wtab[3] + macierzCpc[r][13][i][j] * wtab[1] * wtab[3] + macierzCpc[r][14][i][j] * wtab[2] * wtab[3] + macierzCpc[r][15][i][j] * wtab[3] * wtab[3];
				}
			}
		}
	}
	/*for (int q = 0; q < grid1.ne; q++)
	{
		cout << "\nmacierz C dla el: " << q << endl;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << macierzC[q][i][j] << " ";
			}
			cout << endl;
		}
	}*/

	return macierzC;
}

double* rownania(double** agregowanaMacierzC, double** agregowanaMacierzHHBC, double* agregowanyWektorP, grid grid1, GlobalData data, double* t0)
{
	double** macierz = new double* [grid1.nn];
	for (int i = 0; i < grid1.nn; i++)
	{
		macierz[i] = new double[grid1.nn + 1];
	}
	double** mx = new double* [grid1.nn];
	for (int i = 0; i < grid1.nn; i++)
	{
		mx[i] = new double[grid1.nn];
	}
	double* t1 = new double[grid1.nn];

	double* wektorc = new double[grid1.nn];
	for (int i = 0; i < grid1.nn; i++)
	{
		wektorc[i] = 0;
	}

	for (int i = 0; i < grid1.nn; i++)
	{
		for (int j = 0; j < grid1.nn; j++)
		{
			mx[i][j] = agregowanaMacierzC[i][j] / data.steptime;
		}
	}
	for (int i = 0; i < grid1.nn; i++)
	{
		for (int j = 0; j < grid1.nn; j++)
		{
			wektorc[i] += mx[i][j] * t0[j];
			//cout << mx[i][j] << " * " << t0[i] << " = " << mx[i][j] * t0[i] << endl;
		}
		//cout << "== " << wektorc[i] << endl;
	}

	for (int n = 0; n < grid1.nn; n++)
	{
		for (int j = 0; j < grid1.nn; j++)
		{
			macierz[n][j] = agregowanaMacierzHHBC[n][j] + mx[n][j];
		}
		macierz[n][grid1.nn] = agregowanyWektorP[n] + wektorc[n];
		//cout << macierz[n][grid1.nn] << " = " << agregowanyWektorP[n] << " + " << wektorc[n] << endl;
	}

	//cout << endl << "Macierz rozszerzona: \n";
	//for (int i = 0; i < grid1.nn; i++)
	//{
	//	for (int j = 0; j < grid1.nn + 1; j++)
	//	{
	//		cout << macierz[i][j] << " ";
	//	}
	//	cout << endl;
	//}
	//cout << endl;

	/*cout << endl << "macierz[[h]+[c]/dt]: \n";
	for (int i = 0; i < grid1.nn; i++)
	{
		for (int j = 0; j < grid1.nn + 1; j++)
		{
			cout << macierz[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	

	cout << "Wektor{{P}+{[C]/dT}*{T0}}\n";
	for (int i = 0; i < grid1.nn; i++)
	{
		cout << macierz[i][grid1.nn] << " ";
	}
	cout << endl;*/
	
	

	int N = grid1.nn;

	double temp, s;
	//Sprowadzanie do macierzy trojkatnej
	for (int j = 0; j < N - 1; j++)
	{
		for (int i = j + 1; i < N; i++)
		{
			temp = macierz[i][j] / macierz[j][j];

			for (int k = 0; k < N + 1; k++)
				macierz[i][k] -= macierz[j][k] * temp;
		}
	}
	//Wyliczanie niewiadomych t1 od do³u macierzy, wyliczajac t1[15] mozna wyliczyc t1[14] i tak az do wyliczenia ostatniej niewiadomej t1[0]
	for (int i = N - 1; i >= 0; i--)
	{
		s = 0;
		for (int j = i + 1; j < N; j++)
			s += macierz[i][j] * t1[j];
		t1[i] = (macierz[i][N] - s) / macierz[i][i];
	}

	/*for (int i = 0; i < N; i++)
		cout << "t1[" << i << "]=" << t1[i] << endl;*/

	return t1;
}

void minMaxTempInEachStep(double** agregowanaMacierzC, double** agregowanaMacierzHHBC, double* agregowanyWektorP, grid grid1, GlobalData data, double* t0, int stept)
{
	for (int i = data.steptime; i <= data.time; i += data.steptime)
	{
		//cout << "\nIteracja nr: " << (i / 50)-1 << endl;
		t0 = rownania(agregowanaMacierzC, agregowanaMacierzHHBC, agregowanyWektorP, grid1, data, t0);
		double min = t0[0];
		double max = t0[0];
		for (int j = 0; j < grid1.nn; j++)
		{
			if (t0[j] < min)
			{
				min = t0[j];
			}
			if (t0[j] > max)
			{
				max = t0[j];
			}
		}
		//cout << endl;
		cout << "Time: " << i << " Min: " << min << " Max: " << max << endl;
		
	}
}
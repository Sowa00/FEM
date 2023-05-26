#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <crtdbg.h>
#include <iomanip>

#include "czytanieZPliku.h"

using namespace std;

struct Bok //sciana elementu skonczonego
{
	int pc;//ilosc punktow calkowania
	double pcx = 0;//wspolrzedna punktu calkowanie x
	double pcy = 0;//wspolrzedna punktu calkowania y
	double** wartosciPc;//wartosci funkcji ksztaltu w tabeli, ilosc wierszy = pc
	Bok()
	{
	}
	Bok(int pc)
	{
		this->pc = pc;
		wartosciPc = new double* [pc];
		for (int i = 0; i < pc; i++)
		{
			wartosciPc[i] = new double[4];
		}
		for (int i = 0; i < pc; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				wartosciPc[i][j] = 0;
			}
		}
	}
	~Bok()
	{
		//delete[] wartosciPc;
	}
	//Bok** liczWartosci(grid grid1, GlobalData data, int pc, double* eta, double* ksi)
	void liczWartosci(grid grid1, double* pctab, double* wtab, int tempi, int tempj) //liczenie funkcji ksztaltu dla bokow
	{
		if (tempj == 0)
		{
			if (grid1.nd[grid1.el[tempi].id[0] - 1].bc * grid1.nd[grid1.el[tempi].id[1] - 1].bc == 1)//dolny bok
			{
				//cout << "DOLNY SPRAWDZENIA BC: " << grid1.nd[grid1.el[tempi].id[0] - 1].bc << " * " << grid1.nd[grid1.el[tempi].id[1] - 1].bc << endl;
				for (int j = 0; j < pc; j++)
				{
					wartosciPc[j][0] = 0.25 * (1 - pctab[j]) * (1 - (-1));
					wartosciPc[j][1] = 0.25 * (1 + pctab[j]) * (1 - (-1));
					wartosciPc[j][2] = 0.25 * (1 + pctab[j]) * (1 + (-1));
					wartosciPc[j][3] = 0.25 * (1 - pctab[j]) * (1 + (-1));
				}
				//cout << "\nDOLNY BOK LICZENIE\n";
			}
			else
			{
				for (int j = 0; j < pc; j++)
				{
					wartosciPc[j][0] = 0;
					wartosciPc[j][1] = 0;
					wartosciPc[j][2] = 0;
					wartosciPc[j][3] = 0;
				}
				//cout << "brak bc dla elementu: " << tempi + 1 << " dla boku: " << tempj + 1 << endl;
			}
		}
		if (tempj == 1)
		{
			if (grid1.nd[grid1.el[tempi].id[1] - 1].bc * grid1.nd[grid1.el[tempi].id[2] - 1].bc == 1)//prawy bok
			{
				//cout << "PRAWY SPRAWDZENIA BC: " << grid1.nd[grid1.el[tempi].id[1] - 1].bc << " * " << grid1.nd[grid1.el[tempi].id[2] - 1].bc << endl;
				for (int j = 0; j < pc; j++)
				{
					wartosciPc[j][0] = 0.25 * (1 - 1) * (1 - pctab[j]);
					wartosciPc[j][1] = 0.25 * (1 + 1) * (1 - pctab[j]);
					wartosciPc[j][2] = 0.25 * (1 + 1) * (1 + pctab[j]);
					wartosciPc[j][3] = 0.25 * (1 - 1) * (1 + pctab[j]);
				}
				//cout << "\nPRAWY BOK LICZENIE\n";
			}
			else
			{
				for (int j = 0; j < pc; j++)
				{
					wartosciPc[j][0] = 0;
					wartosciPc[j][1] = 0;
					wartosciPc[j][2] = 0;
					wartosciPc[j][3] = 0;
				}
				//cout << "brak bc dla elementu: " << tempi + 1 << " dla boku: " << tempj + 1 << endl;
			}
		}
		if (tempj == 2)
		{
			if (grid1.nd[grid1.el[tempi].id[2] - 1].bc * grid1.nd[grid1.el[tempi].id[3] - 1].bc == 1)//gorny bok
			{
				//cout << "GORNY SPRAWDZENIA BC: " << grid1.nd[grid1.el[tempi].id[2] - 1].bc << " * " << grid1.nd[grid1.el[tempi].id[3] - 1].bc << endl;
				int temp = pc - 1;
				for (int j = 0; j < pc; j++)
				{
					wartosciPc[j][0] = 0.25 * (1 - pctab[temp]) * (1 - 1);
					wartosciPc[j][1] = 0.25 * (1 + pctab[temp]) * (1 - 1);
					wartosciPc[j][2] = 0.25 * (1 + pctab[temp]) * (1 + 1);
					wartosciPc[j][3] = 0.25 * (1 - pctab[temp]) * (1 + 1);
					temp--;
				}
				//cout << "\nGORNY BOK LICZENIE\n";
			}
			else
			{
				for (int j = 0; j < pc; j++)
				{
					wartosciPc[j][0] = 0;
					wartosciPc[j][1] = 0;
					wartosciPc[j][2] = 0;
					wartosciPc[j][3] = 0;
				}
				//cout << "brak bc dla elementu: " << tempi + 1 << " dla boku: " << tempj + 1 << endl;
			}
		}
		if (tempj == 3)
		{
			if (grid1.nd[grid1.el[tempi].id[3] - 1].bc * grid1.nd[grid1.el[tempi].id[0] - 1].bc == 1)//lewy bok
			{
				//cout << "LEWY SPRAWDZENIA BC: " << grid1.nd[grid1.el[tempi].id[3] - 1].bc << " * " << grid1.nd[grid1.el[tempi].id[0] - 1].bc << endl;
				int temp = pc - 1;
				for (int j = 0; j < pc; j++)
				{
					wartosciPc[j][0] = 0.25 * (1 - (-1)) * (1 - pctab[temp]);
					wartosciPc[j][1] = 0.25 * (1 + (-1)) * (1 - pctab[temp]);
					wartosciPc[j][2] = 0.25 * (1 + (-1)) * (1 + pctab[temp]);
					wartosciPc[j][3] = 0.25 * (1 - (-1)) * (1 + pctab[temp]);
					temp--;
				}
				//cout << "\nLEWY BOK LICZENIE\n";
			}
			else
			{
				for (int j = 0; j < pc; j++)
				{
					wartosciPc[j][0] = 0;
					wartosciPc[j][1] = 0;
					wartosciPc[j][2] = 0;
					wartosciPc[j][3] = 0;
				}
				//cout << "brak bc dla elementu: " << tempi + 1 << " dla boku: " << tempj + 1 << endl;
			}
		}
		/*cout << "element: " << tempi << endl;
		cout << "bok: " << tempj << endl;
		for (int j = 0; j < pc; j++)
		{
			cout << "pc: " << j << endl;
			for (int q = 0; q < 4; q++)
			{
				cout << wartosciPc[j][q] << " ";
			}
			cout << endl;
		}*/

	}
};

double*** macierzHBCObliczanieplik(grid grid1, GlobalData data, Bok** boki, int pc, double* wtab)
{
	//cout << "macierzHBCObliczanie" << endl;
	/*double** macierzHBC = new double* [grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		macierzHBC[i] = new double[pc];
	}*/

	double***** macierzpc = new double**** [grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		macierzpc[i] = new double*** [4];
		for (int j = 0; j < 4; j++)
		{
			macierzpc[i][j] = new double** [pc];
			for (int q = 0; q < pc; q++)
			{
				macierzpc[i][j][q] = new double* [4];
				for (int r = 0; r < 4; r++)
				{
					macierzpc[i][j][q][r] = new double[4];
				}
			}
		}
	}
	for (int i = 0; i < grid1.ne; i++)
	{
		for (int p = 0; p < 4; p++)
		{
			for (int j = 0; j < 4; j++)
			{
				for (int q = 0; q < 4; q++)
				{
					for (int r = 0; r < pc; r++)
					{
						macierzpc[i][p][r][j][q] = data.alfa * (boki[i][p].wartosciPc[r][j] * boki[i][p].wartosciPc[r][q]);
					}
				}
			}
		}
	}
	/*cout << "macierzepc dla hbc" << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << "bok: " << i << endl;
		for (int j = 0; j < pc; j++)
		{
			cout << "pc: " << j << endl;
			for (int q = 0; q < 4; q++)
			{
				for (int r = 0; r < 4; r++)
				{
					cout << macierzpc[i][j][q][r] << " ";
				}
				cout << endl;
			}
		}
	}*/

	double**** macierzHBC = new double*** [grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		macierzHBC[i] = new double** [4];
		for (int j = 0; j < 4; j++)
		{
			macierzHBC[i][j] = new double* [4];
			for (int q = 0; q < 4; q++)
			{
				macierzHBC[i][j][q] = new double[4];
			}
		}
	}

	for (int i = 0; i < grid1.ne; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int q = 0; q < 4; q++)
			{
				for (int r = 0; r < 4; r++)
				{
					macierzHBC[i][j][q][r] = 0;
				}
			}
		}

	}

	double** det = new double* [grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		det[i] = new double[4];
	}

	for (int i = 0; i < grid1.ne; i++)
	{
		det[i][0] = (sqrt(pow(grid1.nd[grid1.el[i].id[1] - 1].x - grid1.nd[grid1.el[i].id[0] - 1].x, 2) + pow(grid1.nd[grid1.el[i].id[1] - 1].y - grid1.nd[grid1.el[i].id[0] - 1].y, 2))) * 0.5;
		det[i][1] = (sqrt(pow(grid1.nd[grid1.el[i].id[2] - 1].x - grid1.nd[grid1.el[i].id[1] - 1].x, 2) + pow(grid1.nd[grid1.el[i].id[2] - 1].y - grid1.nd[grid1.el[i].id[1] - 1].y, 2))) * 0.5;
		det[i][2] = (sqrt(pow(grid1.nd[grid1.el[i].id[3] - 1].x - grid1.nd[grid1.el[i].id[2] - 1].x, 2) + pow(grid1.nd[grid1.el[i].id[3] - 1].y - grid1.nd[grid1.el[i].id[2] - 1].y, 2))) * 0.5;
		det[i][3] = (sqrt(pow(grid1.nd[grid1.el[i].id[0] - 1].x - grid1.nd[grid1.el[i].id[3] - 1].x, 2) + pow(grid1.nd[grid1.el[i].id[0] - 1].y - grid1.nd[grid1.el[i].id[3] - 1].y, 2))) * 0.5;
	}

	/*for (int i = 0; i < grid1.ne; i++)
	{
		cout << "dla elementu " << i+1 << endl;
		for (int j = 0; j < 4; j++)
		{
			cout << "dla boku " << j+1;
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
			}
			cout << "wyznacznik det: " << det[i][j] << endl;
		}
	}*/

	for (int e = 0; e < grid1.ne; e++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int q = 0; q < pc; q++)
			{
				for (int r = 0; r < 4; r++)
				{
					for (int j = 0; j < 4; j++)
					{
						macierzHBC[e][i][r][j] += macierzpc[e][i][q][r][j] * wtab[q];
					}
				}
			}
		}
	}
	for (int e = 0; e < grid1.ne; e++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int q = 0; q < 4; q++)
			{
				for (int r = 0; r < 4; r++)
				{
					macierzHBC[e][i][q][r] *= det[e][i];
					
				}
			}
		}
	}
	

	/*cout << "macierzepc hbc" << endl;
	for (int e = 0; e < grid1.ne; e++)
	{
		cout << "element: " << e+1 << endl;
		for (int i = 0; i < 4; i++)
		{
			cout << "bok: " << i+1;
			if (i == 0)
			{
				cout << " DOLNY" << endl;
			}
			else if (i == 1)
			{
				cout << " PRAWY" << endl;
			}
			else if (i == 2)
			{
				cout << " GORNY" << endl;
			}
			else
			{
				cout << " LEWY" << endl;
			}
			for (int q = 0; q < 4; q++)
			{
				for (int r = 0; r < 4; r++)
				{
					cout << macierzHBC[e][i][q][r] << " ";
				}
				cout << endl;
			}
		}
	}*/
	double*** macierzFinalHBC = new double** [grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		macierzFinalHBC[i] = new double* [4];
		for (int j = 0; j < 4; j++)
		{
			macierzFinalHBC[i][j] = new double[4];
		}
	}
	for (int i = 0; i < grid1.ne; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int q = 0; q < 4; q++)
			{
				macierzFinalHBC[i][j][q] = 0;
			}
		}
	}
	for (int i = 0; i < grid1.ne; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int q = 0; q < 4; q++)
			{
				for (int r = 0; r < 4; r++)
				{
					macierzFinalHBC[i][j][q] += macierzHBC[i][r][j][q];
				}
			}
		}
	}
	/*cout << "Macierz Final HBC: \n";
	for (int i = 0; i < grid1.ne; i++)
	{
		cout << "Element: " << i + 1 << endl;
		for (int j = 0; j < 4; j++)
		{
			for (int q = 0; q < 4; q++)
			{
				cout << macierzFinalHBC[i][j][q]<<" ";
			}
			cout << endl;
		}
	}*/

	return macierzFinalHBC;
}




double** obliczanieWektoraP(grid grid1, GlobalData data, Bok** boki, int pc, double* wtab)
{
	double*** wektorP = new double** [grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		wektorP[i] = new double* [4];
		for (int j = 0; j < 4; j++)
		{
			wektorP[i][j] = new double[4];
		}
	}
	for (int i = 0; i < grid1.ne; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int q = 0; q < 4; q++)
			{
				wektorP[i][j][q] = 0;
			}
		}
	}

	double** detj = new double* [grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		detj[i] = new double[4];
	}

	for (int i = 0; i < grid1.ne; i++)
	{
		detj[i][0] = (sqrt(pow(grid1.nd[grid1.el[i].id[1] - 1].x - grid1.nd[grid1.el[i].id[0] - 1].x, 2) + pow(grid1.nd[grid1.el[i].id[1] - 1].y - grid1.nd[grid1.el[i].id[0] - 1].y, 2))) * 0.5;
		detj[i][1] = (sqrt(pow(grid1.nd[grid1.el[i].id[2] - 1].x - grid1.nd[grid1.el[i].id[1] - 1].x, 2) + pow(grid1.nd[grid1.el[i].id[2] - 1].y - grid1.nd[grid1.el[i].id[1] - 1].y, 2))) * 0.5;
		detj[i][2] = (sqrt(pow(grid1.nd[grid1.el[i].id[3] - 1].x - grid1.nd[grid1.el[i].id[2] - 1].x, 2) + pow(grid1.nd[grid1.el[i].id[3] - 1].y - grid1.nd[grid1.el[i].id[2] - 1].y, 2))) * 0.5;
		detj[i][3] = (sqrt(pow(grid1.nd[grid1.el[i].id[0] - 1].x - grid1.nd[grid1.el[i].id[3] - 1].x, 2) + pow(grid1.nd[grid1.el[i].id[0] - 1].y - grid1.nd[grid1.el[i].id[3] - 1].y, 2))) * 0.5;
	}

	for (int i = 0; i < grid1.ne; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int q = 0; q < 4; q++)
			{
				for (int p = 0; p < pc; p++)
				{
					wektorP[i][j][q] += data.alfa * (wtab[p] * (boki[i][j].wartosciPc[p][q] * data.tot)) * detj[i][j];
				}
			}
		}
	}
	/*cout << "\nWEKTOR P:\n\n";
	for (int i = 0; i < grid1.ne; i++)
	{
		cout << "Element: " << i + 1 << endl;
		for (int j = 0; j < 4; j++)
		{
			cout << "Bok: " << j + 1 << endl;
			for (int q = 0; q < 4; q++)
			{
				cout << wektorP[i][j][q] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}*/

	double** wektorFinalP = new double* [grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		wektorFinalP[i] = new double[4];
	}
	for (int i = 0; i < grid1.ne; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			wektorFinalP[i][j] = 0;
		}
	}

	for (int i = 0; i < grid1.ne; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int q = 0; q < 4; q++)
			{
				wektorFinalP[i][j] += wektorP[i][q][j];
			}
		}
	}

	/*cout << "WEKTOR FINAL P: " << endl;
	for (int i = 0; i < grid1.ne; i++)
	{
		cout << "Element: " << i+1 << endl;
		for (int j = 0; j < 4; j++)
		{
			cout << wektorFinalP[i][j] << " ";
		}
		cout << endl;
	}*/

	return wektorFinalP;

}
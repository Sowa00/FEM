#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <crtdbg.h>
#include <iomanip>

#include "czytanieZPliku.h"

using namespace std;

double*** macierzHObliczanieplik(double* Eta, double* Ksi, double* pctab, double* wtab, grid grid1, GlobalData data, int pc) //obliczanie macierzy H lokalna
{
		
	double* tabEta = new double[pow(pc, 2)];
	double* tabKsi = new double[pow(pc, 2)];
	for (int i = 0; i < pow(pc, 2); i++)
	{
		tabEta[i] = Eta[i];
		tabKsi[i] = Ksi[i];
	}
	ElUniversal el4(pc);
	el4.pochodne(tabEta, tabKsi);
	//el4.showtab();

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
	double** detjodwpc = new double* [grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		detjpc[i] = new double[pow(pc, 2)];
		detjodwpc[i] = new double[pow(pc, 2)];
	}

	double**** macierzJakobianOdwrotnapc = new double*** [grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		macierzJakobianOdwrotnapc[i] = new double** [pow(pc, 2)];
		for (int j = 0; j < pow(pc, 2); j++)
		{
			macierzJakobianOdwrotnapc[i][j] = new double* [2];
			for (int q = 0; q < 2; q++)
			{
				macierzJakobianOdwrotnapc[i][j][q] = new double[2];
			}
		}
	}

	//liczenie wyznacznika, odwrotnego wyznacznika i odwrotnej macierzy jakobiego
	for (int i = 0; i < grid1.ne; i++)
	{
		for (int j = 0; j < pow(pc, 2); j++)
		{
			detjpc[i][j] = (macierzJakobianpc[i][j][0][0] * macierzJakobianpc[i][j][1][1] - macierzJakobianpc[i][j][0][1] * macierzJakobianpc[i][j][1][0]);
			detjodwpc[i][j] = 1 / detjpc[i][j];
			//cout << "\ndetj: " << detjpc[i][j] << " 1/detj: " << detjodwpc[i][j] << endl;
			macierzJakobianOdwrotnapc[i][j][0][0] = macierzJakobianpc[i][j][1][1];
			macierzJakobianOdwrotnapc[i][j][0][1] = -macierzJakobianpc[i][j][0][1];
			macierzJakobianOdwrotnapc[i][j][1][0] = -macierzJakobianpc[i][j][1][0];
			macierzJakobianOdwrotnapc[i][j][1][1] = macierzJakobianpc[i][j][0][0];
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
					macierzJakobianOdwrotnapc[i][j][q][r] *= detjodwpc[i][j];
				}
			}
		}
	}
	/*for (int r = 0; r < grid1.ne; r++)
	{
		cout << "\nelement: " << r + 1 << endl;
		for (int q = 0; q < pow(pc, 2); q++)
		{
			cout << "\nmacierz odwrotna jakobiego pc" << q + 1 << ":" << endl;
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					cout << macierzJakobianOdwrotnapc[r][q][i][j] << " ";
				}
				cout << endl;
			}
		}
	}*/

	double**** ndxpc = new double*** [grid1.ne];
	double**** ndypc = new double*** [grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		ndxpc[i] = new double** [pow(pc, 2)];
		ndypc[i] = new double** [pow(pc, 2)];
		for (int j = 0; j < pow(pc, 2); j++)
		{
			ndxpc[i][j] = new double* [pow(pc, 2)];
			ndypc[i][j] = new double* [pow(pc, 2)];
			for (int q = 0; q < pow(pc, 2); q++)
			{
				ndxpc[i][j][q] = new double[4];
				ndypc[i][j][q] = new double[4];
			}
		}
	}

	//liczenie dn/dx i dn/dy
	for (int i = 0; i < grid1.ne; i++)
	{
		for (int j = 0; j < pow(pc, 2); j++)
		{
			for (int q = 0; q < pow(pc, 2); q++)
			{
				for (int r = 0; r < 4; r++)
				{
					ndxpc[i][j][q][r] = macierzJakobianOdwrotnapc[i][j][0][0] * el4.dnksi[q][r] + macierzJakobianOdwrotnapc[i][j][0][1] * el4.dneta[q][r];
					ndypc[i][j][q][r] = macierzJakobianOdwrotnapc[i][j][1][0] * el4.dnksi[q][r] + macierzJakobianOdwrotnapc[i][j][1][1] * el4.dneta[q][r];
				}
			}
		}
	}
	/*cout << "ndx" << endl;
	for (int i = 0; i < grid1.ne; i++)
	{
		cout << "element: " << i << endl;
		for (int j = 0; j < pow(pc, 2); j++)
		{
			cout << "pc: " << j << endl;
			for (int q = 0; q < pow(pc, 2); q++)
			{
				for (int r = 0; r < 4; r++)
				{
					cout << ndxpc[i][j][q][r] << " ";
				}
				cout << endl;
			}
		}
	}*/
	/*cout << "ndy" << endl;
	for (int i = 0; i < grid1.ne; i++)
	{
		cout << "element: " << i << endl;
		for (int j = 0; j < pow(pc, 2); j++)
		{
			cout << "pc: " << j << endl;
			for (int q = 0; q < pow(pc, 2); q++)
			{
				for (int r = 0; r < 4; r++)
				{
					cout << ndypc[i][j][q][r] << " ";
				}
				cout << endl;
			}
		}
	}*/

	double**** macierzH1xpc = new double*** [grid1.ne];
	double**** macierzH1ypc = new double*** [grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		macierzH1xpc[i] = new double** [pow(pc, 2)];
		macierzH1ypc[i] = new double** [pow(pc, 2)];
		for (int j = 0; j < pow(pc, 2); j++)
		{
			macierzH1xpc[i][j] = new double* [4];
			macierzH1ypc[i][j] = new double* [4];
			for (int q = 0; q < 4; q++)
			{
				macierzH1xpc[i][j][q] = new double[4];
				macierzH1ypc[i][j][q] = new double[4];
			}
		}
	}

	//liczenie macierzy Hx i Hy
	for (int k = 0; k < grid1.ne; k++)
	{

		for (int j = 0; j < 4; j++)
		{
			for (int q = 0; q < 4; q++)
			{
				for (int r = 0; r < pow(pc, 2); r++)
				{
					macierzH1xpc[k][r][j][q] = ndxpc[k][r][r][j] * ndxpc[k][r][r][q];
					macierzH1ypc[k][r][j][q] = ndypc[k][r][r][j] * ndypc[k][r][r][q];
				}
			}
		}
		//cout << endl;
	}

	/*cout << "macierzh1xpc" << endl;
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
					cout << macierzH1xpc[i][j][q][r] << " ";
				}
				cout << endl;
			}
		}
	}*/

	double**** macierzHpc = new double*** [grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		macierzHpc[i] = new double** [pow(pc, 2)];
		for (int j = 0; j < pow(pc, 2); j++)
		{
			macierzHpc[i][j] = new double* [4];
			for (int q = 0; q < 4; q++)
			{
				macierzHpc[i][j][q] = new double[4];
			}
		}
	}

	//Liczenie macierzy H, ktore beda sumowane aby otrzymac macierz final H
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
					macierzHpc[i][j][q][r] = (macierzH1xpc[i][j][q][r] + macierzH1ypc[i][j][q][r]) * data.conductivity * detjpc[i][j];
				}
			}
		}
	}
	//for (int i = 0; i < grid1.ne; i++)
	//{
	//	cout << "\nmacierz do sumowania element("<< i+1 <<"): \n";
	//	for (int j = 0; j < pow(pc, 2); j++)
	//	{
	//		cout << "pc = " << j + 1 << endl;
	//		for (int q = 0; q < 4; q++)
	//		{
	//			for (int r = 0; r < 4; r++)
	//			{
	//				cout << macierzHpc[i][j][q][r] << " ";
	//			}
	//			cout << endl;
	//		}
	//		cout << endl;
	//	}
	//}


	double*** macierzH = new double** [grid1.ne];
	for (int i = 0; i < grid1.ne; i++)
	{
		macierzH[i] = new double* [pow(pc, 2)];
		for (int j = 0; j < pow(pc, 2); j++)
		{
			macierzH[i][j] = new double[4];
		}
	}
	for (int i = 0; i < grid1.ne; i++)
	{
		for (int j = 0; j < pow(pc, 2); j++)
		{
			for (int q = 0; q < 4; q++)
			{
				macierzH[i][j][q] = 0;
			}
		}
	}

	/*for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			macierzH[i][j] = macierzHpc[0][i][j] * w[0] * w[0] + macierzHpc[1][i][j] * w[1] * w[0] + macierzHpc[2][i][j] * w[0] * w[1] + macierzHpc[3][i][j] * w[1] * w[1];

		}
	}*/

	//liczenie macierzy H - macierz transportu ciepla
	if (pc == 2)
	{
		for (int q = 0; q < grid1.ne; q++)
		{
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					macierzH[q][i][j] += macierzHpc[q][0][i][j] * wtab[0] * wtab[0] + macierzHpc[q][1][i][j] * wtab[1] * wtab[0] +
										 macierzHpc[q][2][i][j] * wtab[0] * wtab[1] + macierzHpc[q][3][i][j] * wtab[1] * wtab[1];
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
					macierzH[r][i][j] = macierzHpc[r][0][i][j] * wtab[0] * wtab[0] + macierzHpc[r][1][i][j] * wtab[1] * wtab[0] + macierzHpc[r][2][i][j] * wtab[2] * wtab[0] +
										macierzHpc[r][3][i][j] * wtab[0] * wtab[1] + macierzHpc[r][4][i][j] * wtab[1] * wtab[1] + macierzHpc[r][5][i][j] * wtab[2] * wtab[1] +
										macierzHpc[r][6][i][j] * wtab[0] * wtab[2] + macierzHpc[r][7][i][j] * wtab[2] * wtab[1] + macierzHpc[r][8][i][j] * wtab[2] * wtab[2];
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
					macierzH[r][i][j] = macierzHpc[r][0][i][j]  * wtab[0] * wtab[0] + macierzHpc[r][1][i][j]  * wtab[1] * wtab[0] + macierzHpc[r][2][i][j]  * wtab[2] * wtab[0] + macierzHpc[r][3][i][j]  * wtab[3] * wtab[0] + 
										macierzHpc[r][4][i][j]  * wtab[0] * wtab[1] + macierzHpc[r][5][i][j]  * wtab[1] * wtab[1] + macierzHpc[r][6][i][j]  * wtab[2] * wtab[1] + macierzHpc[r][7][i][j]  * wtab[3] * wtab[1] + 
										macierzHpc[r][8][i][j]  * wtab[0] * wtab[2] + macierzHpc[r][9][i][j]  * wtab[1] * wtab[2] + macierzHpc[r][10][i][j] * wtab[2] * wtab[2] + macierzHpc[r][11][i][j] * wtab[3] * wtab[2] + 
										macierzHpc[r][12][i][j] * wtab[0] * wtab[3] + macierzHpc[r][13][i][j] * wtab[1] * wtab[3] + macierzHpc[r][14][i][j] * wtab[2] * wtab[3] + macierzHpc[r][15][i][j] * wtab[3] * wtab[3];
				}
			}
		}
	}
	/*for (int q = 0; q < grid1.ne; q++)
	{
		cout << "\nmacierz H dla el: " << q << endl;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << macierzH[q][i][j] << " ";
			}
			cout << endl;
		}
	}*/

	return macierzH;
}
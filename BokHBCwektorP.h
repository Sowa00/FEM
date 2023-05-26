#pragma once
struct Bok //sciana elementu skonczonego
{
	int pc;//ilosc punktow calkowania
	double pcx = 0;//wspolrzedna punktu calkowanie x
	double pcy = 0;//wspolrzedna punktu calkowania y
	double** wartosciPc;//wartosci funkcji w tabeli, ilosc wierszy = pc
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
	void liczWartosci(grid grid1, double* pctab, double* wtab, int tempi, int tempj)
	{
		/*cout << "WYSWIETLENIE BC DLA CALEGO ELEMENTU" << endl;
		for (int i = 0; i < 4; i++)
		{
			cout << "ID: " << grid1.el[tempi].id[i] << " bc: ";
			cout << grid1.nd[grid1.el[tempi].id[i]-1].bc << " ";
		}
		cout << endl<<endl;*/
		/*for (int i = 0; i < pc; i++)
		{
			cout << "pctab: " << pctab[i] << " wtab: " << wtab[i] << endl;
		}
		cout << endl;*/
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

double*** macierzHBCObliczanieplik(grid grid1, GlobalData data, Bok** boki, int pc, double* wtab);
double** obliczanieWektoraP(grid grid1, GlobalData data, Bok** boki, int pc, double* wtab);

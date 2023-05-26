#pragma once

#include <string>
using namespace std;

struct node
{
	double x;
	double y;
	double t = 100;
	int bc = 0;
};

struct element
{
	int id[4];
};

struct grid
{
	int nn;
	int ne;
	grid(int x, int y) : nn(x), ne(y) {};
	node* nd = new node[nn];
	element* el = new element[ne];
};

struct GlobalData
{
	int time;
	int steptime;
	int conductivity;
	int alfa;
	int tot;
	int density;
	int specific_heat;
};

struct ElUniversal
{

	int x;
	int y;
	double** dnksi;
	double** dneta;
	double** Ntab;
	ElUniversal(int b)
	{
		x = 4;
		y = pow(b, 2);
		dnksi = new double* [y];
		dneta = new double* [y];
		Ntab = new double* [y];
		for (int i = 0; i < y; i++)
		{
			dnksi[i] = new double[x];
			dneta[i] = new double[x];
			Ntab[i] = new double[x];
		}
	}
	~ElUniversal()
	{
		delete[] dnksi;
		delete[] dneta;
	}

	void showtabdn() //wyswietlanie dn/dksi i dn/deta
	{
		cout << "dN/dKsi:" << endl;
		for (int i = 0; i < y; i++)
		{
			for (int j = 0; j < x; j++)
			{
				cout << dnksi[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl << "dN/dEta:" << endl;
		for (int i = 0; i < y; i++)
		{
			for (int j = 0; j < x; j++)
			{
				cout << dneta[i][j] << " ";
			}
			cout << endl;
		}
	}
	void showtabn() //wyswietlanie dn/dksi i dn/deta
	{
		cout << "Funkcje ksztaltu:" << endl;
		for (int i = 0; i < y; i++)
		{
			for (int j = 0; j < x; j++)
			{
				cout << Ntab[i][j] << " ";
			}
			cout << endl;
		}
	}
	void pochodne(double* eta, double* ksi)
	{
		for (int i = 0; i < y; i++)
		{
			dnksi[i][0] = -0.25 * (1 - eta[i]);
			dnksi[i][1] = 0.25 * (1 - eta[i]);
			dnksi[i][2] = 0.25 * (1 + eta[i]);
			dnksi[i][3] = -0.25 * (1 + eta[i]);

			dneta[i][0] = -0.25 * (1 - ksi[i]);
			dneta[i][1] = -0.25 * (1 + ksi[i]);
			dneta[i][2] = 0.25 * (1 + ksi[i]);
			dneta[i][3] = 0.25 * (1 - ksi[i]);
		}
	}
	void funkcjeKsztaltu(double* pctab, double* wtab)
	{
		for (int i = 0; i < y; i++)
		{
			Ntab[i][0] = 0.25 * (1 - pctab[i]) * (1 - wtab[i]);
			Ntab[i][1] = 0.25 * (1 - pctab[i]) * (1 + wtab[i]);
			Ntab[i][2] = 0.25 * (1 + pctab[i]) * (1 + wtab[i]);
			Ntab[i][3] = 0.25 * (1 + pctab[i]) * (1 - wtab[i]);
		}
	}
};

void show_globaldata(GlobalData data);

void show_grid(grid gr);

GlobalData load_from_file_data(string file);

int node_number(string file);

int el_number(string file);

grid load_from_file_grid(string file);

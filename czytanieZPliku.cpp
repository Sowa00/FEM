#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <crtdbg.h>
#include <iomanip>

using namespace std;

struct node //struktura wezla
{
	double x;
	double y;
	double t = 100;
	int bc = 0;
};

struct element //struktura elementu zawierajacego id do wezlow
{
	int id[4];
};

struct grid 
{
	int nn;//ilosc wezlow
	int ne;//ilosc elementow
	grid(int x, int y) : nn(x), ne(y) {};
	node* nd = new node[nn]; //tablica przechowujaca wezly
	element* el = new element[ne]; //tablica przechowujaca elementy
};

struct GlobalData
{
	int time;		   //czas symulacji
	int steptime;	   //krok czasowy symulacji
	int conductivity;  //przewodnosc
	int alfa;		   //alfa z warunku brzegowego
	int tot;		   //temperatura otoczenia z warunku brzegowego
	int density;	   //gestosc
	int specific_heat; //cieplo wlasciwe
};

struct ElUniversal //struktura elementu uniwersalnego przechowujaca dn/dksi i dn/deta oraz funkcje ksztaltu
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
	void pochodne(double* eta, double* ksi)			//obliczanie dn/dksi i dn/deta po danych punktach calkowania
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
	void funkcjeKsztaltu(double* eta, double* ksi)
	{
		for (int i = 0; i < y; i++)
		{
			Ntab[i][0] = 0.25 * (1 - eta[i]) * (1 - ksi[i]);
			Ntab[i][1] = 0.25 * (1 - eta[i]) * (1 + ksi[i]);
			Ntab[i][2] = 0.25 * (1 + eta[i]) * (1 + ksi[i]);
			Ntab[i][3] = 0.25 * (1 + eta[i]) * (1 - ksi[i]);
		}
	}
};

void showfile(string file) //wyswietlenie calego pliku zczytanego
{
	ifstream load(file);
	string z;
	while (!load.eof())
	{
		getline(load, z);
		cout << z << endl;
	}
	load.close();
}

void show_globaldata(GlobalData data) //wyswietlenie danych zczytanych z pliku dla global data
{
	cout << "Global Data:\nSimulation Time: " << data.time << "\nSimulation Step Time: " << data.steptime
		<< "\nConductivity: " << data.conductivity << "\nAlfa: " << data.alfa << "\nTot: " << data.tot
		<< "\nDensity: " << data.density << "\nSpecific Heat: " << data.specific_heat << endl;
}

void show_grid(grid gr) //wyswietlenie danych zczytanych z pliku dla wezlow i elementow
{
	cout << "Grid:\n" << "Nodes: " << gr.nn << "\nElements: " << gr.ne << endl;
	cout << "Nodes data:\n";
	for (int i = 0; i < gr.nn; i++)
	{
		cout << i + 1 << ":\tx: " << gr.nd[i].x << "\ty: " << gr.nd[i].y << "\tInitial temp: " << gr.nd[i].t << "\tBC: " << gr.nd[i].bc << endl << endl;
	}
	for (int i = 0; i < gr.ne; i++)
	{
		cout << "Element [ID: " << i + 1 << "]:  " << "{ ";
		for (int j = 0; j < 4; j++)
		{
			cout << gr.el[i].id[j] << ", ";
		}
		cout << "}" << endl;
	}
}

GlobalData load_from_file_data(string file) //wczytanie danych dla struktury global data
{
	GlobalData gd;
	ifstream load(file);
	string z;
	double x;
	//x to node temp
	load >> z >> gd.time >> z >> gd.steptime >> z >> gd.conductivity >> z >> gd.alfa >> z >> gd.tot >> z >> x >> z >> gd.density >> z >> gd.specific_heat;
	load.close();
	return gd;
}
int node_number(string file) //wczytanie ilosci wezlow
{
	ifstream load(file);
	string z;
	int x;
	int nn;
	load >> z >> x >> z >> x >> z >> x >> z >> x >> z >> x >> z >> x >> z >> x >> z >> x; //globaldata
	load >> z >> z >> nn;
	load.close();
	return nn;
}
int el_number(string file) //wczytanie ilosci elementow
{
	ifstream load(file);
	string z;
	int x;
	int ne;
	load >> z >> x >> z >> x >> z >> x >> z >> x >> z >> x >> z >> x >> z >> x >> z >> x; //globaldata
	load >> z >> z >> z >> z >> z >> ne >> z;
	load.close();
	return ne;
}
grid load_from_file_grid(string file) //wczytanie danych dla wezlow i elementow
{
	//grid gr;
	ifstream load(file);
	string z;
	int x;
	int nn, ne;
	load >> z >> x >> z >> x >> z >> x >> z >> x >> z >> x >> z >> x >> z >> x >> z >> x; //globaldata
	load >> z >> z >> nn >> z >> z >> ne >> z; //nodes numer element number
	//cout << gr.nn << endl << gr.ne;
	grid gr(nn, ne);
	for (int i = 0; i < gr.nn; i++)
	{
		load >> z >> gr.nd[i].x >> z >> gr.nd[i].y;
	}
	load >> z >> z;
	for (int i = 0; i < gr.ne; i++)
	{
		load >> z >> gr.el[i].id[0] >> z >> gr.el[i].id[1] >> z >> gr.el[i].id[2] >> z >> gr.el[i].id[3];
	}
	load >> z;
	while (!load.eof())
	{
		load >> x;
		load >> z;
		gr.nd[x - 1].bc = 1;
	}

	return gr;
}
#pragma once
double*** macierzCObliczanie(double* Eta, double* Ksi, double* pctab, double* wtab, grid grid1, GlobalData data, int pc);
double* rownania(double** agregowanaMacierzC, double** agregowanaMacierzHHBC, double* agregowanyWektorP, grid grid1, GlobalData data, double* t0, int stept);
void minMaxTempInEachStep(double** agregowanaMacierzC, double** agregowanaMacierzHHBC, double* agregowanyWektorP, grid grid1, GlobalData data, double* t0, int stept);
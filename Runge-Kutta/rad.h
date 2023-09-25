#ifndef _RAD_H_
#define _RAD_H_

int ucitaj_matricu_iz_datoteke(float A[][100], int n, FILE *f);
int ucitaj_vektor_iz_datoteke(float x[], int n, FILE *f);
float racunaj(int p, int n, int k, float A[][100], float y[], float b[], float h);
int uslov(int n, float u[], float v[], float eps, int stepen2);

#endif
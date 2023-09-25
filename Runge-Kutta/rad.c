#include<stdio.h>

#include<stdlib.h>

#include<string.h>

#include<math.h>

#include"rad.h"


// funkcija ucitava matricu dimenzija nxn iz datoteke f
int ucitaj_matricu_iz_datoteke(float A[][100], int n, FILE * f) {

  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      if (fscanf(f, "%f ", & A[i][j]) != 1) // ukoliko fscanf nije uspeo da ucita broj, vraca kao signal 1
        return 1;

  return 0; // inace ako je prosao kroz celu petlju i nije vratio 1, vraca 0 kao signal da je ucitavanje uspesno
}

// funkcija ucitava vektor dimenzije nx1 iz datoteke f
int ucitaj_vektor_iz_datoteke(float x[], int n, FILE * f) {

  int i = 0;

  for (i; i < n; i++)
    if (fscanf(f, "%f ", & x[i]) != 1) // ukoliko fscanf nije uspeo da ucita broj, vraca kao signal 1
      return 1;

  return 0; // inace ako je prosao kroz celu petlju i nije vratio 1, vraca 0 kao signal da je ucitavanje uspesno

}

// funkcija racuna koeficijente k_1,2,3,4 
//za dimenziju sistema n, koordinatu k, matricu A, 
//pocetni vektr y, slobodan vektor b sa korakom h
//i vraca vrednost funkcije u u k-toj koordinati
float racunaj(int p, int n, int k, float A[][100], float y[], float b[], float h) {

  int i, j;
  float k1[100];
  float k2[100];
  float k3[100];
  float k4[100];
  float u[100];

  // racunam k1:

  for (i = 0; i < n; i++) {
    float s = 0;
    for (j = 0; j < n; j++) {
      s += A[i][j] * y[j];
    }
    s = s + b[i];
    k1[i] = h * s;
  }

  // racunam k2:

  for (i = 0; i < n; i++) {
    float s = 0;
    for (j = 0; j < n; j++) {
      s += A[i][j] * (y[j] + k1[j] / 2);
    }
    s = s + b[i];
    k2[i] = h * s;
  }

  // racunam k3:

  for (i = 0; i < n; i++) {
    float s = 0;
    for (j = 0; j < n; j++) {
      if (p == 4) { // koristi se formula za racunanje k3, za p = 4 
        s += A[i][j] * (y[j] + k2[j] / 2);
      } else if (p == 3) { // koristi se formula za racunanje k3, za p = 3 
        s += A[i][j] * (y[j] - k1[j] + 2 * k2[j]);
      }
    }
    s = s + b[i];
    k3[i] = h * s;
  }

  // racunam k4:

  for (i = 0; i < n; i++) {
    float s = 0;
    for (j = 0; j < n; j++) {
      s += A[i][j] * (y[j] + k3[j]);
    }
    s = s + b[i];
    k4[i] = h * s;
  }

  if (p == 4)
    return u[k] = y[k] + h / 6 * (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k]);
  else if (p == 3)
    return u[k] = y[k] + h / 6 * (k1[k] + 4 * k2[k] + k3[k]);
}

// funkcija uslov proverava da li funkcije u i v za zadato epsilon
// postizu odgovarajucu tacnost po svih n-koordinata
int uslov(int n, float u[], float v[], float eps, int stepen2) {

  int i;

  for (i = 0; i < n; i++)
    if (fabsf(u[i] - v[i]) / (stepen2 - 1) >= eps) // ukoliko nije ispunjen kriterijum za neko i=0,..,n-1, vraca kao signal -1
      return -1;

  return 1; // inace, ako je prosao kroz celu petlju i nije vratio -1, vraca 1 kao signal da je uslov zadovoljen

}
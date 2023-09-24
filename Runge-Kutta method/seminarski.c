#include<stdio.h>

#include<stdlib.h>

#include<string.h>

#include<math.h>

#include"rad.h"

int main(int argc, char ** argv) {

  if (argc != 3) {
    fprintf(stderr, "Neispravan unos! Program pokrenuti sa:\n./a.out ime_datoteke.txt p\n");
    exit(1);

  }

  FILE * f = fopen(argv[1], "r");

  if (f == NULL) {
    fprintf(stderr, "Neuspelo otvaranje datoteke za citanje!\n");
    exit(1);
  }

  int p = atoi(argv[2]); // red metode

  if (p != 3 && p != 4) {
    fprintf(stderr, "Neispravan red metode! Metoda mora biti reda 3 ili 4.\n");
    exit(1);
  }

  int stepen2 = pow(2, p); // racunam 2^p kako bih kasnije koristila u formuli

  int n; // dimenzija sistema
  float A[100][100]; // matrica sistema
  float b[100]; // vektor slobodnih clanova
  float x[100], u[100]; // pocetna tacka i vrednost funkcije u pocetnoj tacki
  int i, j;
  float h; // korak 0.25
  float eps; // tacnost metode

  // ucitavanje podataka iz datoteke

  fscanf(f, "%d", & n);
  if (ucitaj_matricu_iz_datoteke(A, n, f) == 1) {
    fprintf(stderr, "Neuspelo ucitavanje matrice!\n");
    exit(1);
  }
  if (ucitaj_vektor_iz_datoteke(b, n, f) == 1) {
    fprintf(stderr, "Neuspelo ucitavanje vektora!\n");
    exit(1);
  }
  if (ucitaj_vektor_iz_datoteke(x, n, f) == 1) {
    fprintf(stderr, "Neuspelo ucitavanje vektora!\n");
    exit(1);
  }
  if (ucitaj_vektor_iz_datoteke(u, n, f) == 1) {
    fprintf(stderr, "Neuspelo ucitavanje vektora!\n");
    exit(1);
  }

  fscanf(f, " %f %f ", & h, & eps);

  fclose(f);

  float v[100];

  /*

  	u - priblizno resenje u tacki x+h
  	v - priblizno resenje u istoj toj tacki samo za korak 2*h
  	
  	Ukoliko je |u-v|/(2^p-1)>=eps , potrebno je izvrsiti popravku:
  			up = u + (u-v)/(2^p-1)

  */

  // racunam vrednosti za u i v, za pocetni vektor u 
  for (i = 0; i < n; i++) {
    v[i] = racunaj(p, n, i, A, u, b, 2 * h); // prvo se racunaju vrednosti za v, pa onda za u
    u[i] = racunaj(p, n, i, A, u, b, h); // jer ukoliko se prvo izracuna u, onda se ta nova vrednost koristi u funkciji za v sto bi bila greska
  }

  // proveravam da li je ispunjen kriterijum zaustavljanja
  int usl = uslov(n, u, v, eps, stepen2);

  while (usl == -1) { // dok nije ispunjen kriterijum zaustavljanja
    h = h / 2;
    for (i = 0; i < n; i++) {
      if ((fabsf(u[i] - v[i]) / (stepen2 - 1)) >= eps) { // za one u[i] i v[i] za koje nije ispunjen uslov
        u[i] = u[i] + (u[i] - v[i]) / (stepen2 - 1); // racunam popravku
      }
    }
    for (i = 0; i < n; i++) { // nakon popravke, racunam nove u i v, za pocetni vektor u
      v[i] = racunaj(p, n, i, A, u, b, 2 * h);
      u[i] = racunaj(p, n, i, A, u, b, h);
    }
    usl = uslov(n, u, v, eps, stepen2); // proverim uslov za nove u i v 
  }

  // ispisujem resenje
  for (i = 0; i < n; i++) {
    printf("%f \n", u[i]);
  }

  return 0;

}
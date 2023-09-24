#include <stdio.h>

#include <stdlib.h>

#include <math.h>

#include <armadillo>

#include <iostream>     // std:://cout

#include <algorithm>// std::random_shuffle

#include <vector> // std::vector

#include <ctime>        // std::time

#include <cstdlib>

#include <initializer_list>

#include <stack>

#include <chrono>

#include <sstream>

#include <string>

#undef max
#undef min

#pragma warning(disable: 4996)
using namespace std;
using namespace arma;

int randint(int i) {
  return randi < int > (distr_param(0, i - 1));
}

//nasumicno bira k elemenata iz skupa v
// za k==n, ovo bude bas permutacija vektora od n elemenata 
// 
uvec randsample(uvec v, int k) {
  uvec d(k);
  if (k == 1) {
    int ind = randint(v.n_elem); // kada je k == 1, samo uzme bilo koji od elemnata iz vektora v
    d(0) = v(ind);
  } else {
    uvec v_shuffle = arma::shuffle(v); // kada je k > 1, uzmem nekako proizvoljne elemente
    d = v_shuffle(span(0, k - 1)); // i onda uzmem prvih k 
  }
  return d;
}

//vraca vektor koji predstavlja razliku vektora v i u, tj. v\u (kao skupovna razlika)
// vektor koji nastaje ovom razlikom je sortiran :)
//
uvec setdiff(uvec v, uvec u) {
  int l = v.n_elem;
  uvec s(l);
  int k = 0;
  for (int i = 0; i < l; i++) {
    if (all(u != v(i))) {
      s(k) = v(i);
      k += 1;
    }
  }
  if (k == 0) {
    return uvec();
  }
  return s(span(0, k - 1));
}

// funkcija koja vraca indeks minimalnog rastojanja od cvora do haba
// potrebna za pronalazenje pocetnog resenja
//
int min_rastojanje(int i, uvec x, mat D) {

  int n = D.n_rows;
  int p = x.n_elem;
  //cout << xx << endl;
  uvec xpom = zeros < uvec > (n);
  double dmin = std::numeric_limits < double > ::infinity();
  int index_min = 0;
  // pravim vektor koji sadrzi jedinice na pozicijama gde su habovi 
  // npr. x = [7, 2, 5]
  // xx = [2,5,7]
  // xpom = [0 0 1 0 0 1 0 1 0 0 .. 0]
  // 
  int br = 0;
  for (int k = 0; k < n; k++) {
    for (int br = 0; br < p; br++) {
      if (k == x(br)) {
        xpom(k) = 1;
        br += 1;
      }
    }
  }

  //cout << xpom << endl;
  // 
  // sada trazim minimalno rastojanje

  for (int s = 0; s < n; s++) {
    //cout << D(i, s) << endl;
    if (xpom(s) == 1 & D(i, s) != 0 & dmin > D(i, s)) {
      dmin = D(i, s);
      index_min = s;
    }
  }

  //cout << index_min << endl;
  //cout << dmin << endl; 

  return index_min;

}

// funkcija koja vraca matricu koja sadrzi 1 na pozicijama na kojima je i-ti cvor pridruzen k-tom habu 
// za inicajno resenje i shaking resenje
//
umat napravi_matricu(uvec x, mat D) {

  int n = D.n_rows;
  int p = x.n_elem;
  int pozicija;
  umat X = zeros < umat > (n, n);
  uvec xx = sort(x);
  //cout << xx << endl;
  // matrica nxn sa 1 na pozicijaam gde je cvor i pridruzen habu k 

  for (int a = 0; a < p; a++) {
    X(xx(a), xx(a)) = 1; // da bi i Xkk bili =1 tj. uslov da je hab dodeljen sam sebi
  }

  int ind = 0; // indikator da li je cvor i hab cvor
  for (int i = 0; i < n; i++) {
    for (int s = 0; s < p; s++) {
      //cout << i << xx(s) << endl;
      if (i == xx(s)) { // ukoliko je i jednak bilo kom hab cvoru, treba da preskocim tu vrstu
        ind = 1; // postavim indikator na 1
        //cout << ind << endl;
      }
    }
    if (ind == 0) { // ukoliko i nije hab cvor, treba da mu odredim najblizi hab
      pozicija = min_rastojanje(i, xx, D);
      //cout << pozicija << endl;
      X(i, pozicija) = 1; // i postavim poziciju u matrici X na 1 
    }
    ind = 0;
  }

  return X;

}

// funkcija koja nasumicno bira p od n cvorova koji ce biti habovi
// vraca vektor x koji sadrzi pozicije hab-cvorova
//
uvec nasumicno_pocetno_hab_resenje(mat D, int p) {
  int n = D.n_rows;
  uvec x = randsample(linspace < uvec > (0, n - 1, n), p); // izaberem nekih p pozicija na kojima su habovi iz niza [0,1,2,...,n-1]
  return x;
}

double funkcija_cilja(uvec x, umat X, mat D, mat C, double alpha, double beta) {

  double f = 0.0;
  int pozicija;
  int n = D.n_rows;
  int p = x.n_elem;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < p; k++) {
        for (int m = 0; m < p; m++) {
          if (X(i, x(k)) == 1 & X(j, x(m)) == 1 & (C(i, x(k)) + alpha * C(x(k), x(m)) + C(x(m), j)) <= beta) {
            //cout << X(i, x(k))  << endl;
            //cout << X(j, x(m)) << endl;
            //cout << (C(i, x(k)) + alpha * C(x(k), x(m)) + C(j, x(m))) << endl;
            //cout << beta << endl;
            f += D(i, j) / 1000;
          }
        }
      }
    }
  }

  return f;
}

// okolina za shaking:
// trenutni hab se gasi, a otvara se novi hab 
// koji je jedna od alokacija koja je do tada bila 
// pridruzena tom habu koji se zatvara
//
umat nasumicno_iz_okoline_Nk(int k, uvec & x, umat X, int n) {

  uvec k_hab_index = randsample(x, k); // izaberem k od p habova iz vektora xx koji sadrzi habove

  umat X_n = X;

  for (int i = 0; i < k; i++) {

    uvec kolona = X.col(k_hab_index(i)); // izvucem kolonu gde je hab
    kolona(k_hab_index(i)) = 0; // izbacim jedinicu sa pozicije gde je hab, pa ostaju samo 1 tamo gde su alokacije

    int indikator = kolona.is_zero(); // ako su u koloni sve nule, nemam poziciju za novi hab

    if (indikator == 0) { // ako kolona sadrzi bar jednu jedinicu
      uvec alokacije = find(kolona); // nadjem sve pozicije alokacija iz kolone
      uvec pozicija = randsample(alokacije, 1); // izaberem neku od njih nasumicno
      X_n.col(k_hab_index(i)) = X.col(pozicija(0)); // zamenim kolone
      X_n.col(pozicija(0)) = X.col(k_hab_index(i));
    }
    X = X_n; // predjem u to novo resenje
  }

  return X_n;
}

// okolina za VND:
// funkcija koja menja resenje tako sto menja tacno jednu alokaciju, proverava da li je doslo do poboljsanja 
// i vraca najbolje resenje
//
umat najbolje_resenje_iz_okoline_Nl(uvec x, umat X, mat D, mat C, double alpha, double beta, double & fb) {

  int p = x.n_elem;
  int n = D.n_rows;

  int ind = 1;
  umat X_poboljsanje = X;
  umat X_nova = X;

  while (ind == 1) {
    ind = 0; // dokle god je fb bolje od prethodno dobijene vrednosti funkcije cilja
    for (int i = 0; i < n; i++) { // prodjem kroz sve cvorove
      for (int s = 0; s < p; s++) { // i sve habove
        //uvec xpom = x(span(0, s));
        if (X(i, x(s)) == 1) { // za cvor i koji je dodeljen posmatranom habu
          for (int t = 0; t < p; t++) {
            X_nova = X; // vratim matricu na pocetnu kako bih menjala ponovo tacno jednu alokaciju
            if (x(s) != x(t)) { // ukoliko hab nije onaj koji vec posmatram
              X_nova(i, x(s)) = 0; // i vise nije dodljen tom habu
              X_nova(i, x(t)) = 1; // nego nekom drugom
              double f_novo = funkcija_cilja(x, X_nova, D, C, alpha, beta); // za takvo resenje izracunam fju cilja
              if (f_novo > fb) {
                ind = 1; // uporedim sa najboljom 
                fb = f_novo; // o ako je bolje postavim fb na tu novu vrednost
                X_poboljsanje = X_nova; // zapamtim poboljsano resenje
                X = X_poboljsanje;
              }
            }
            if (ind == 1) {
              break;
            }

          }

        }
        if (ind == 1) {
          break;
        }
      }
      if (ind == 1) {
        break;
      }

    }

    X = X_poboljsanje;

  }

  return X_poboljsanje;
}

//GVNS:
//
double GVNS(int n, int p, mat D, mat C, double alpha, double beta, int maxiter, int kmax, int lmax, double & vreme) {

  wall_clock timer;

  int k;
  int l;
  int iter = 0;

  // xb - x best
  // xsh - x shaking
  // xn - x novo
  // xv - x VND
  uvec xb, xsh, xn, xv, xvnd;
  umat X_b, X_sh, X_n, X_v, X_vnd;
  double fb, fsh, fv, fvnd, f_stari;

  chrono::steady_clock::time_point begin = chrono::steady_clock::now();
  chrono::steady_clock::time_point t_current;

  // pocetno resenje (nasumicno), a ujedno i trenutno najbolje resenje 
  xb = nasumicno_pocetno_hab_resenje(D, p);
  X_b = napravi_matricu(xb, D);
  fb = funkcija_cilja(xb, X_b, D, C, alpha, beta);

  //cout << "Pocetno fb\n" << fb << endl;

  while (iter < maxiter) {
    k = 1;
    while (k <= kmax) {
      // shaking faza
      X_sh = nasumicno_iz_okoline_Nk(k, xb, X_b, n); // vektor koji izabere nasumicno habove
      xsh = find(diagvec(X_sh));
      fsh = funkcija_cilja(xsh, X_sh, D, C, alpha, beta);
      l = 1;
      // VND faza
      while (l <= lmax) {
        f_stari = fsh; // pamtim vr fje cilja iz shakinga
        X_v = najbolje_resenje_iz_okoline_Nl(xsh, X_sh, D, C, alpha, beta, fsh); // fsh se prosledjuje preko pokazivaca, pa ce odavde biti povucena vrednost fje cilja u okviru VND faze
        xv = find(diagvec(X_v));
        if (fsh > f_stari) { // ako je tacno, uzimam ove vrednosti i resetujem brojac
          X_sh = X_v;
          // fsh ostaje promenjena jer je doslo do poboljsanja
          xsh = xv;
          l = 1;
        } else { // ako nije tacno, povecavam brojac i pretrazujem u sledecoj okolini za nepromenjene vrednosti
          fsh = f_stari; // vracam vrednost funkcije cilja na onu koja je bila bolja
          l = l + 1;
        }
      }
      // gotov je VND, pa proveravam da li je dosao do boljih rezultata
      X_vnd = X_sh;
      xvnd = xsh;
      fvnd = fsh;
      if (fvnd > fb) { // ako je tacno, onda uzimam ovo za najbolje vredosti i resetujem brojac
        xb = xvnd;
        X_b = X_vnd;
        fb = fvnd;
        k = 1;
        t_current = chrono::steady_clock::now();
        vreme = chrono::duration_cast < chrono::microseconds > (t_current - begin).count() / 1000000.0; //belezimo prvo vreme dostizanja najboljeg resenja
        //cout << "Novo fb\n" << fb << endl; 
        //cout << "Iteracija:\n" << iter << endl;

      } else { // ako ne, uvecavam brojac i idem na narednu okolinu za nepromenjene vrednosti 
        k = k + 1;
      }
    }
    iter = iter + 1;

  }

  //cout << "funkcija cilja:" << endl << fb << endl << "habovi" << endl << sort(xb) << endl;
  //cout << "resenje:" << endl << X_b << endl;
  //cout << vreme << endl; 

  return fb;
}

// glavni program

int main(int argc, char * argv[]) {
  //citanje podataka
  FILE * file;

  //format instance:
  /*
  n (br. korisnika/potencijalnih lokacija centara)
  p (br. centara koji treba otvoriti)
  d11,d12,...,d1n   (matrica rastojanja)
  d21,d22,...,d2n
  ....
  dn1,dn2,...,dnn
  c11,c12,...,c1n   (matrica troskova)
  c21,c22,...,c2n
  ....
  cn1,cn2,...,cnn
  alpha
  beta
  */
  arma_rng::set_seed_random();

  int n;
  int p;

  file = fopen("C:\\armadillo\\examples\\ulaz_USApHMCP.txt", "r");
  if (file == NULL) {
    printf("Error!");
    exit(1);
  }

  fscanf(file, "%d", & n);
  fscanf(file, "%d", & p);

  mat D = zeros < mat > (n, n);
  mat C = zeros < mat > (n, n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      fscanf(file, "%lf", & D(i, j));
    }
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      fscanf(file, "%lf", & C(i, j));
    }
  }

  double alpha;
  double beta;

  fscanf(file, "%lf", & alpha);
  fscanf(file, "%lf", & beta);
  double opt_sol;
  fscanf(file, "%lf", & opt_sol);

  //kraj ucitvanja podataka iz fajla 

  //cout << n << endl << p << endl;
  //cout << D << endl;
  //cout << C << endl;
  //cout << alpha << endl << beta << endl;
  //cout << opt_sol << endl;

  // parametri:

  int maxiter = 100;
  int kmax = p; // odnosi se na okoline Nk (k=1,..,k_max)

  if (n == 15) {
    maxiter = 300;
  }

  if (n == 20 || n == 25) {
    maxiter = 500;
  }

  int lmax = 1; //odnosi se na okoline Nl (l=1,...,l_max), u mom slucaju je 1

  // podaci za metaheuristiku

  vec sol_i(15); //najbolje resenje MA u i-tom izvrsavanju
  vec t_i(15); // vreme za koje je MA prvi put dobio najbolje resenje soli, tzv. pocetno vreme(obavezno);
  vec ttot_i(15); // ukupno vreme i - tog izvrsavanja MA, odnosno vreme rada do zadovoljenja nekog kriterijuma zaustavljanja(obavezno);
  vec gap_i(15);

  double vreme;
  for (int br = 0; br < 15; br++) {
    chrono::steady_clock::time_point begin = chrono::steady_clock::now(); //trenutno vreme (pocetak)
    sol_i(br) = GVNS(n, p, D, C, alpha, beta, maxiter, kmax, lmax, vreme);
    chrono::steady_clock::time_point end = chrono::steady_clock::now(); //trenutno vreme (kraj)
    double elapsed = chrono::duration_cast < chrono::microseconds > (end - begin).count() / 1000000.0; // ukupno vreme izvrsavanja
    ttot_i(br) = elapsed;
    t_i(br) = vreme;
  }

  //cout << ttot_i << endl;
  //cout << t_i << endl;

  double t = sum(t_i) / 15; // srednje pocetno vreme 
  double t_tot = sum(ttot_i) / 15; // srednje ukupno vreme izvrsavanja
  //double t_best = min(t_i); // najbolje vreme 

  //cout << "srednje pocetno vreme t:\n" << t << endl;
  //cout << "srednje ukupno vreme t_tot:\n" << t_tot << endl;
  //cout << "najbolje vreme t_best:\n" << t_best << endl;
  printf("%.4f\n", t);
  printf("%.4f\n", t_tot);

  double best_sol = max(sol_i);
  //cout << "sva resenja:\n" << sol_i << endl;
  //cout << "najbolje resenje:\n" << endl;
  cout << best_sol << endl;
  //printf("%.2f\n", best_sol);

  // provera da li koristiti najbolje resenje heuristike ili postoji bolje poznato iz literature
  for (int br = 0; br < 15; br++) {
    if (best_sol >= opt_sol) {
      //cout << "best:\n" << best_sol << endl;
      gap_i(br) = 100 * abs(sol_i(br) - best_sol) / abs(best_sol);
    } else if (best_sol < opt_sol) {
      //cout << "optimal:\n" << opt_sol << endl;
      gap_i(br) = 100 * abs(sol_i(br) - opt_sol) / abs(opt_sol);
    }
  }

  //cout << "gap za svako resenje:\n" << gap_i << endl; 
  double suma_gapi = sum(gap_i);

  //cout << "suma gap_i:\n" << suma_gapi << endl; 

  double agap = suma_gapi / 15;

  //cout << "agap:\n" << agap << endl;
  printf("%.4f\n", agap);

  double sigma; // standardna drvijacija

  vec razlika_gap(15);
  for (int br = 0; br < 15; br++) {
    razlika_gap(br) = (gap_i(br) - agap) * (gap_i(br) - agap);
  }

  //cout << "razlika gap:\n" << razlika_gap << endl;

  sigma = sqrt(sum(razlika_gap) / 15);
  //cout << "standardna devijacija:\n" << sigma << endl;
  printf("%.4f\n", sigma);

  cout << "" << endl;

  system("pause");
  return 0;
}
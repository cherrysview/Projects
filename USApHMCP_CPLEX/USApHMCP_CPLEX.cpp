#include<ilcplex/ilocplex.h>

ILOSTLBEGIN
typedef IloArray < IloNumVarArray > NumVarMatrix;
typedef IloArray < NumVarMatrix > NumVar3Matrix;
typedef IloArray < NumVar3Matrix > NumVar4Matrix;
typedef IloArray < IloNumArray > NumMatrix;
typedef IloArray < NumMatrix > Num3Matrix;
typedef IloArray < Num3Matrix > Num4Matrix;

int main(int, char ** ) {
  IloEnv env;
  try {
    IloModel model(env);
    IloInt n, p, check;
    IloNum alpha, beta;

    ifstream data("USApHMCP_ulaz.txt"); // napomena: promeniti ime ulaznog fajla ili kopirati instance u .txt fajl sa ovim nazivom
    data >> n;
    data >> p;

    //definisemo matrice
    IloArray < IloNumArray > C(env, n);
    for (int i = 0; i < n; i++)
      C[i] = IloNumArray(env, n);

    IloArray < IloNumArray > D(env, n);
    for (int i = 0; i < n; i++)
      D[i] = IloNumArray(env, n);

    //citamo podatke
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
        data >> D[i][j];
        D[i][j] /= 1000;
      }

    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        data >> C[i][j];

    //cout << C << endl;
    //cout << D << endl;

    data >> alpha;
    data >> beta;

    data >> check;

    if (check == 777)
      cout << "Uspesno!\n" << endl;
    else
      cout << "Neuspesno!\n" << endl;

    // promenljive:

    // promenljiva odluke X binarna 
    IloArray < IloNumVarArray > x(env, n);
    for (int i = 0; i < n; i++) {
      x[i] = IloNumVarArray(env, n, 0, 1, ILOBOOL);
    }

    // promenljiva odluke Y binarna
    NumVar4Matrix y(env, n); //4D matrica

    for (int i = 0; i < n; i++) {
      y[i] = NumVar3Matrix(env, n); //3D matrica
      for (int j = 0; j < n; j++) {
        y[i][j] = NumVarMatrix(env, n); // 2D matrica
        for (int k = 0; k < n; k++) {
          y[i][j][k] = IloNumVarArray(env, n); //niz
          for (int m = 0; m < n; m++) {
            y[i][j][k][m] = IloNumVar(env, 0, 1, ILOBOOL); // elementi
          }
        }
      }
    }

    // q je 1 ako vazi uslov vremnskog roka putovanja, a 0 inace

    Num4Matrix q(env, n); //4D matrica
    for (int i = 0; i < n; i++) {
      q[i] = Num3Matrix(env, n); //3D matrica
      for (int j = 0; j < n; j++) {
        q[i][j] = NumMatrix(env, n); // 2D matrica
        for (int k = 0; k < n; k++) {
          q[i][j][k] = IloNumArray(env, n); //niz
          for (int m = 0; m < n; m++) {
            if (C[i][k] + alpha * C[k][m] + C[j][m] <= beta) {
              q[i][j][k][m] = 1;
            } else {
              q[i][j][k][m] = 0;
            }
          }

        }
      }
    }

    // ogranicenja:

    IloExpr Ogr1(env);

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
          for (int m = 0; m < n; m++) {
            Ogr1 = 2 * y[i][j][k][m] - x[i][k] - x[j][m];
            model.add(Ogr1 <= 0);
            Ogr1.end();
            Ogr1 = IloExpr(env);
          }
        }
      }
    }

    IloExpr Ogr2(env);

    for (int i = 0; i < n; i++) {
      for (int k = 0; k < n; k++) {
        Ogr2 += x[i][k];
      }
      model.add(Ogr2 == 1);
      Ogr2.end();
      Ogr2 = IloExpr(env);
    }

    IloExpr Ogr3(env);

    for (int k = 0; k < n; k++) {
      Ogr3 += x[k][k];
    }
    model.add(Ogr3 == p);

    IloExpr Ogr4(env);

    for (int i = 0; i < n; i++) {
      for (int k = 0; k < n; k++) {
        Ogr4 = x[i][k] - x[k][k];
        model.add(Ogr4 <= 0);
        Ogr4.end();
        Ogr4 = IloExpr(env);
      }
    }

    // funkcija cilja:

    IloExpr fo(env);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
          for (int m = 0; m < n; m++) {
            fo += D[i][j] * y[i][j][k][m] * q[i][j][k][m];
          }
        }
      }
    }

    model.add(IloMaximize(env, fo));

    IloCplex cplex(env);
    //	cplex.setOut(env.getNullStream());
    cplex.extract(model);

    // VREME I POCETNI CVOR
    // parametri
    cplex.setParam(IloCplex::Param::TimeLimit, 7200);
    cplex.setParam(IloCplex::HeurFreq, -1);

    IloNum start = clock();
    if (!cplex.solve()) {
      env.error() << "Nije moguce resiti!\n" << endl;
      throw (-1);
    }
    IloNum end = clock();

    //ispisi promenljive
    ofstream izlaz("izlaz_USApHMCP.txt");
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        izlaz << cplex.getValue(x[i][j]) << " ";
      }
      izlaz << endl;
    }
    izlaz << endl;
    env.out() << "Status: " << cplex.getStatus() << endl;
    env.out() << "Vrednost funkcije cilja: " << cplex.getObjValue() << endl;
    env.out() << "Vreme izvrsavanja: " << cplex.getTime() << endl;
    env.out() << "Broj iteracija: " << cplex.getNiterations() << endl;
    env.out() << "Broj cvorova: " << cplex.getNnodes() << endl;
    env.out() << "Gap:" << cplex.getMIPRelativeGap() * 100 << endl;
    env.out() << "Vreme izvrsavanja - clock: " << (double)(end - start) / (double) CLK_TCK << endl;

  } catch (IloException & ex) {
    cerr << "Error Cplex: " << ex << endl;
  }
}
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 4

void b(int slot, int etat[N]) {
    etat[slot] = 0;
}

void bdag(int slot, int etat[N]) {
    etat[slot] = 1;
}

int puiss2(int nombre) {
    int puissance = 1;
    for (int i = 0; i < nombre; i++) {
        puissance *= 2;
    }
    return puissance;
}

void transformation_binaire(int nombre, int *etat) {
    for (int i = 0; i < N; i++) {
        etat[i] = 0;
    }

    int tour = 1;
    while (nombre > 0) {
        int reste = nombre % 2;
        etat[N - tour] = reste;
        nombre /= 2;
        tour++;
    }
}

int transformation_decimal(int *binaire) {
    int puissance = N - 1;
    int nombre = 0;

    for (int i = 0; i < N; i++) {
        nombre += binaire[i] * puiss2(puissance);
        puissance--;
    }
    return nombre;
}

double **matrice_hamiltonien(int n) {
    double **etat_h = malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        etat_h[i] = malloc(n * sizeof(double));
        for (int j = 0; j < n; j++) {
            etat_h[i][j] = 0;
        }
    }

    for (int k = 0; k < n; k++) {
        int *saver = malloc(N * sizeof(int));
        transformation_binaire(k, saver);

        for (int m = 0; m < N; m++) {
            int next = (m + N + 1) % N;

            if (saver[next] == 1 && saver[m] == 0) {
                b(next, saver);
                bdag(m, saver);
                int indice = transformation_decimal(saver);
                etat_h[indice][k] = -1;
                bdag(next, saver);
                b(m, saver);
            }

            if (saver[next] == 0 && saver[m] == 1) {
                bdag(next, saver);
                b(m, saver);
                int indice = transformation_decimal(saver);
                etat_h[indice][k] = -1;
                b(next, saver);
                bdag(m, saver);
            }
        }

        free(saver);
    }

    return etat_h;
}

void jacobi(double **matrix, double eps, int n) {
    double maximum;

    do {
        for (int p = 0; p < n; p++) {
            for (int q = 0; q < n; q++) {
                if (q != p && fabs(matrix[p][q]) >= eps) {
                    double phi = 0.5 * atan((2 * matrix[p][q]) /
                                   (matrix[q][q] - matrix[p][p]));

                    for (int r = 0; r < n; r++) {
                        if (r != p && r != q) {
                            double m_rp = matrix[r][p];
                            double m_rq = matrix[r][q];

                            matrix[r][p] = m_rp * cos(phi) - m_rq * sin(phi);
                            matrix[r][q] = m_rq * cos(phi) + m_rp * sin(phi);

                            matrix[p][r] = matrix[r][p];
                            matrix[q][r] = matrix[r][q];
                        }
                    }

                    double m_pp = matrix[p][p];
                    double m_qq = matrix[q][q];
                    double m_pq = matrix[p][q];

                    matrix[p][p] = m_pp * cos(phi)*cos(phi)
                                 + m_qq * sin(phi)*sin(phi)
                                 - 2 * m_pq * cos(phi)*sin(phi);

                    matrix[q][q] = m_pp * sin(phi)*sin(phi)
                                 + m_qq * cos(phi)*cos(phi)
                                 + 2 * m_pq * cos(phi)*sin(phi);

                    matrix[p][q] = 0.0;
                    matrix[q][p] = 0.0;
                }
            }
        }

        maximum = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j && fabs(matrix[i][j]) > maximum) {
                    maximum = fabs(matrix[i][j]);
                }
            }
        }

    } while (maximum >= eps);
}

double *valeurs_propres(int n) {
    double **h = matrice_hamiltonien(n);
    jacobi(h, 1e-12, n);

    double *vp = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        vp[i] = h[i][i];
    }

    for (int i = 0; i < n; i++) {
        free(h[i]);
    }
    free(h);

    return vp;
}

double calcul_energie(double T, int n, double vp[]) {
    double numerateur = 0.0;
    double denominateur = 0.0;

    for (int i = 0; i < n; i++) {
        double E = vp[i];
        double et = exp(-E / T);
        numerateur += E * et;
        denominateur += et;
    }

    return numerateur / denominateur;
}

int main() {
    int n = puiss2(N);
    double *vp = valeurs_propres(n);

    int nt = 100;

    for (int i = 0; i < nt; i++) {
        double tlog = log(1e-1) + i * (log(1e2) - log(1e-1)) / (nt - 1);
        double T = exp(tlog);
        double E = calcul_energie(T, n, vp);
        printf("%f\n", E);
    }

    free(vp);
    return 0;
}

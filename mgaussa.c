#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int op = 0;

void matrix1(long double ***m, long double **b, int *n) {
    *n = 8;
    *m = malloc(*n * sizeof(long double *));
    for (int i = 0; i < *n; i++) {
        (*m)[i] = malloc((*n + 1) * sizeof(long double));
    }
    *b = malloc(*n * sizeof(long double));
    
    for (int i = 0; i < *n; i++) {
        for (int j = 0; j < *n; j++) {
           (*m)[i][j] = 1.0 / (i + j + 1); 
        }
        (*b)[i] = 1.0 / (*n + i);
    }
}

void matrix2(long double ***m, long double **b, int *n) {
    *n = 50;
    *m = malloc(*n * sizeof(long double *));
    for (int i = 0; i < *n; i++) {
        (*m)[i] = malloc((*n + 1) * sizeof(long double));
    }
    *b = malloc(*n * sizeof(long double));
    
    for (int j = 0; j < *n; j++) {  
        (*m)[0][j] = 0;
    }    
    (*m)[0][0] = (*m)[0][*n - 1] = 1; 

    for (int i = 1; i < *n; i++) {
        for (int j = 0; j < i; j++) {
            (*m)[i][j] = -(i - j)*(i - j);
        }
        (*m)[i][i] = 1;
        (*m)[i][*n - 1] = 1; 
        (*b)[i] = 0.1;
    }
}

long double* m_Gaussa(long double **m, long double *b, int n) {
    int *index_x = malloc(n * sizeof(int)); 
    for (int i = 0; i < n; i++) index_x[i] = i;

    // Прямой ход метода Гаусса
    for (int i = 0; i < n; i++) {
       long double max = fabs(m[i][i]);
       int index = i;

       for (int j = i + 1; j < n; j++) {  
            if (fabs(m[j][i]) > max) {
                max = fabs(m[j][i]);
                index = j;
            }
       }
       
        if (index != i) { 
            for(int j = 0; j < n; j++) {
                long double tmp = m[j][i];
                m[j][i] = m[j][index];
                m[j][index] = tmp;
            }
            int tmp_index = index_x[i];
            index_x[i] = index_x[index];
            index_x[index] = tmp_index;
        }

        for(int j = i + 1; j < n; j++) {
            m[i][j] /= m[i][i];
            op++;
        }
        b[i] /= m[i][i];
        op++;
        m[i][i] = 1;

        for (int j = i + 1; j < n; j++) { 
            long double f = m[j][i];
            for (int k = i; k < n; k++) {
                m[j][k] -= f * m[i][k];
                op +=2;
            }
            b[j] -= f * b[i];
            op+=2;
        }
    }

    // Обратный ход
    long double *x = malloc(n * sizeof(long double));
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= m[i][j] * x[j];
            op +=2;
        }
    }

    for (int i = 0; i < n; i++) {
        long double a = x[i];
        x[i] = x[index_x[i]];
        x[index_x[i]] = a;

        int tmp_index = index_x[i];
        index_x[i] = i;
        index_x[tmp_index] = tmp_index;
    }

    free(index_x);
    return x;
}

void endm(long double **m, long double *b, int n, long double *x) {
    for (int i = 0; i < n; i++) {
        free(m[i]);
    }
    free(m);
    free(b);
    free(x);
}

void print_res(long double **m, long double *b, int n, long double *x) {
    printf("Размер матрицы: %d*%d\n", n, n);
    printf("Количество операций %d\n", op);
    printf("Решение системы х=(");
    for(int i = 0; i < n; i++) printf(" %Lf", x[i]);
    printf(")\n");
    long double *r = malloc(n * sizeof(long double));
    long double nevyazka = 0;
    for(int i = 0; i < n; i++) {
        r[i] = b[i]; 
        for (int j = 0; j < n; j++) {
            r[i] -= m[i][j]*x[j]; 
        }
    }
    
    for (int i = 0; i < n; i++) r[i] *= r[i]; 
    for (int i = 0; i < n; i++) nevyazka += r[i];
    nevyazka = sqrt(nevyazka); 
    printf("Норма невязки: %Lf\n", nevyazka);
    
    free(r); 
}

int main() {
    long double **m, *b, *x;
    int n;

    printf("Первая СЛАУ\n");
    matrix1(&m, &b, &n);
    x = m_Gaussa(m, b, n);
    print_res(m, b, n, x);
    endm(m, b, n, x);

    op = 0;
    
    printf("Вторая СЛАУ\n");
    matrix2(&m, &b, &n);
    x = m_Gaussa(m, b, n);
    print_res(m, b, n, x);
    endm(m, b, n, x);

    return 0;
}

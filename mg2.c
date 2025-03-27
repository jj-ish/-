#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int op = 0;

void matrix1(long double ***m, long double **b, int *n) {
    *n = 8;
    *m = malloc(*n * sizeof(long double *));
    for (int i = 0; i < *n; i++) {
        (*m)[i] = malloc(*n * sizeof(long double));
    }
    *b = malloc(*n * sizeof(long double));
    
    for (int i = 0; i < *n; i++) {
        for (int j = 0; j < *n; j++) {
            (*m)[i][j] = 1.0 / (i + j + 1); //+1, так как i, j перебираем с нуля
        }
        (*b)[i] = 1.0 / ((*n) + i);
    }
}

void matrix2(long double ***m, long double **b, int *n) {
    *n = 50;
    *m = malloc(*n * sizeof(long double *));
    for (int i = 0; i < *n; i++) {
        (*m)[i] = calloc(*n, sizeof(long double)); // Инициализация нулями
    }
    *b = malloc(*n * sizeof(long double));
    
    // Заполнение матрицы
    for (int i = 0; i < *n; i++) {
        for (int j = 0; j <= i; j++) {
            (*m)[i][j] = (i == j) ? 1.0 : -(i - j) * (i - j); // Диагональ равна 1
        }
        (*m)[i][*n - 1] = 1.0; // Последний столбец
        (*b)[i] = 0.1;
    }
}

long double** copy_matrix(long double **a, int n) {
    long double **m = malloc(n * sizeof(long double *));
    for (int i = 0; i < n; i++) {
        m[i] = malloc(n * sizeof(long double));
        for (int j = 0; j < n; j++) {
            m[i][j] = a[i][j];
        }
    }
    return m;
}

long double* m_Gaussa(long double **a, long double *b, int n) {
    long double **m = copy_matrix(a, n);
    long double *temp_b = malloc(n * sizeof(long double));
    for (int i = 0; i < n; i++) temp_b[i] = b[i];

    // Прямой ход метода Гаусса
    for (int i = 0; i < n; i++) {
        // Поиск ведущего элемента
        long double max = fabs(m[i][i]);
        int index = i;
        for (int j = i + 1; j < n; j++) {
            if (fabs(m[j][i]) > max) {
                max = fabs(m[j][i]);
                index = j;
            }
        }

        if (max == 0) {
            printf("Матрица вырожденная!\n");
            exit(1);
        }

        // Перестановка строк
        if (index != i) {
            long double *tmp = m[i];
            m[i] = m[index];
            m[index] = tmp;

            long double temp = temp_b[i];
            temp_b[i] = temp_b[index];
            temp_b[index] = temp;
        }

        // Прямой ход
        for (int j = i + 1; j < n; j++) {
            long double factor = m[j][i] / m[i][i];
            op++; // Операция деления
            for (int k = i; k < n; k++) {
                m[j][k] -= factor * m[i][k];
                op += 2; // 1 умножение и 1 вычитание
            }
            temp_b[j] -= factor * temp_b[i];
            op += 2; // 1 умножение и 1 вычитание
        }
    }

    // Обратный ход
    long double *x = malloc(n * sizeof(long double));
    for (int i = n - 1; i >= 0; i--) {
        x[i] = temp_b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= m[i][j] * x[j];
            op += 2; // 1 умножение и 1 вычитание
        }
        x[i] /= m[i][i];
        op++; // Операция деления
    }

    for (int i = 0; i < n; i++) free(m[i]);
    free(m);
    free(temp_b);

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
    printf("Количество операций: %d\n", op);

    long double *r = malloc(n * sizeof(long double));
    long double nevyazka = 0;

    for (int i = 0; i < n; i++) {
        r[i] = -b[i];
        for (int j = 0; j < n; j++) {
            r[i] += m[i][j] * x[j];
        }
        nevyazka += r[i] * r[i];
    }

    nevyazka = sqrt(nevyazka);
    printf("Норма невязки: %.30Lf\n", nevyazka);

    free(r);
}

int main() {
    long double **a, *b, *x;
    int n;

    printf("Первая СЛАУ\n");
    matrix1(&a, &b, &n);
    x = m_Gaussa(a, b, n);
    print_res(a, b, n, x);
    endm(a, b, n, x);

    op = 0;

    printf("\nВторая СЛАУ\n");
    matrix2(&a, &b, &n);
    x = m_Gaussa(a, b, n);
    print_res(a, b, n, x);
    endm(a, b, n, x);

    return 0;
}

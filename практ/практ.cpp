#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

// Вычисление коэффициентов аппроксимирующей прямой
void getApprox(float(&x)[2][9], double* a, double* b, int n) {
    double sumx = 0;
    double sumy = 0;
    double sumx2 = 0;
    double sumxy = 0;
    for (int i = 0; i < n; i++) {
        sumx += x[0][i];
        sumy += x[1][i];
        sumx2 += x[0][i] * x[0][i];
        sumxy += x[0][i] * x[1][i];
    }
    *a = (n * sumxy - (sumx * sumy)) / (n * sumx2 - sumx * sumx);
    *b = (sumy - *a * sumx) / n;
    return;
}
int main() {

    setlocale(LC_ALL, "RU");
    double  a, b;
    const unsigned int DIM1 = 2;
    const unsigned int DIM2 = 9;

    float arr[DIM1][DIM2] = {
        { 0.88, 0.9, 0.91, 0.93, 0.94, 0.96, 0.97, 0.99, 1 },
        { 0.029, 0.086, 0.17, 0.31, 0.43, 0.57, 0.71, 0.86, 0.97 },
    };

    cout << "x     y" << endl;
    for (int i = 0; i < DIM2; i++)
        cout << arr[0][i] << " - " << arr[1][i] << endl;
    getApprox(arr, &a, &b, DIM2);

    cout << endl << "a = " << a << " " << "b = " << b << endl;
    cout << endl << "Уравнение регрессии имеет вид: " << "y = " << b << " + " << a << "*x" << endl;

    //Формирование массивов для изображения линии регрессии на графике.
    float y2[DIM2];
    for (int i = 0; i < DIM2; i++)
        y2[i] = b + a * arr[0][i];

    // коэффцент корреляции
    float sumx = 0, sumy = 0, srx = 0, sry = 0;

    for (int i = 0; i < DIM2; i++) {
        sumx += arr[0][i];
        sumy += arr[1][i];
    }

    srx = sumx / DIM2;
    sry = sumy / DIM2;

    float SR[DIM1][DIM2] = {  //отклонения от среднего арифметического
        { .0, .0, .0, .0, .0, .0, .0, .0, .0 },
        { .0, .0, .0, .0, .0, .0, .0, .0, .0 },
    };

    float SR2[DIM1][DIM2] = {
        { .0, .0, .0, .0, .0, .0, .0, .0, .0 },
        { .0, .0, .0, .0, .0, .0, .0, .0, .0 },
    };

    float SRxy[DIM2] = { .0, .0, .0, .0, .0, .0, .0, .0, .0 };

    for (int i = 0; i < DIM2; i++)
    {
        SR[0][i] = arr[0][i] - srx;
        SR[1][i] = arr[1][i] - sry;
    }

    for (int i = 0; i < DIM1; i++)
        for (int j = 0; j < DIM2; j++)
            SR2[i][j] = SR[i][j] * SR[i][j];

    for (int i = 0; i < DIM2; i++)
        SRxy[i] = SR[0][i] * SR[1][i];

    // сумма произведений отклонений
    float sum_pr_ot_xy = 0;
    for (int i = 0; i < DIM2; i++)
        sum_pr_ot_xy += SRxy[i];

    // суммы квадратов отклонений по x и y
    float sum_kv_ot_x = 0;
    float sum_kv_ot_y = 0;
    for (int j = 0; j < DIM2; j++) {
        sum_kv_ot_x += SR2[0][j];
        sum_kv_ot_y += SR2[1][j];
    }

    //коэфф кореляции
    float kor = 0;
    kor = sum_pr_ot_xy / (sqrt(sum_kv_ot_x * sum_kv_ot_y));

    cout << endl << "Коэффицент кореляции: " << kor << endl;
    cout << endl;

    float Y_reg[DIM2], SSE = 0, MAE = 0;

    for (int i = 0; i < DIM2; i++)
        Y_reg[i] = b + a * (arr[0][i]);

    for (int i = 0; i < DIM2; i++)
        SSE += pow(arr[1][i] - Y_reg[i], 2);

    cout << "Суммарная квадратичная ошибка для линии регрессии SSE: " << SSE << endl;

    MAE = SSE / DIM2;

    cout << endl << "Средняя ошибка для линии регрессии MAE: " << MAE << endl;


    //  Подбор параметров функции y = A^4 + Dt + K
    // считаем суммы

    float  sum_x = 0, sum_x2 = 0, sum_x4 = 0, sum_x5 = 0, sum_x8 = 0, sum_y = 0, sum_yx = 0, sum_yx4 = 0;

    for (int i = 0; i < DIM2; i++)
    {
        sum_x += arr[0][i];
        sum_x2 += round(pow(arr[0][i], 2) * 1000) / 1000;
        sum_x4 += round(pow(arr[0][i], 4) * 1000) / 1000 ;
        sum_x5 += round(pow(arr[0][i], 5) * 1000) / 1000 ;
        sum_x8 += round(pow(arr[0][i], 8) * 1000) / 1000 ;
        sum_y += arr[1][i];
        sum_yx += round(arr[0][i] * arr[1][i] * 1000) / 1000 ;
        sum_yx4 = round(pow(arr[0][i], 4) * arr[1][i] * 1000) / 1000 ;

    }

    // решение системы методом гаусса
    
    int const m = 4, n = 3;
    //создаем массив

    float matrix_0[n][m] = {{9, sum_x, sum_x4, sum_y},
                       { sum_x, sum_x2, sum_x5, sum_yx},
                       { sum_x4, sum_x5, sum_x8, sum_yx4} };

    //инициализируем

    float** matrix = new float* [n];
    for (int i = 0; i < n; i++)
        matrix[i] = new float[m];

    //инициализируем

    for (int i = 0; i < n; i++)

        for (int j = 0; j < m; j++)
        {
            matrix[i][j] = matrix_0[i][j];
        }

    //выводим массив
    cout << endl << "Матрица: " << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
            cout << matrix[i][j] << " ";
        cout << endl;
    }
    cout << endl;

    //Метод Гаусса
    //Прямой ход, приведение к верхнетреугольному виду
    float  tmp, xx[m];
    int k;

    for (int i = 0; i < n; i++)
    {
        tmp = matrix[i][i];
        for (int j = n; j >= i; j--)
            matrix[i][j] /= tmp;
        for (int j = i + 1; j < n; j++)
        {
            tmp = matrix[j][i];
            

            for (k = n; k >= i; k--)
                matrix[j][k] -= tmp * matrix[i][k];
        }
    }
    /*обратный ход*/
    xx[n - 1] = matrix[n - 1][n];
    for (int i = n - 2; i >= 0; i--)
    {
        xx[i] = matrix[i][n];
        for (int j = i + 1; j < n; j++) xx[i] -= matrix[i][j] * xx[j];
    }

    //Выводим решения
    

    cout << "Корни системы: " << endl;

    double k_a = xx[2], k_d = xx[1], k_k = xx[0], SSE2 = 0, MAE2 = 0, Y_reg2[DIM2];

    cout << endl << "A: " << k_a << " D: " << k_d << " K: " << k_k << endl;
    cout << endl << "Подобранная функциональная зависимость заданного вида: y = " << k_k << "*x^4 + "<< k_d <<"*x + "<< k_a << endl;


    for (int i = 0; i < DIM2; i++)
        Y_reg2[i] = (k_a * pow(arr[0][i], 4)) + k_d*arr[0][i] + k_k;


    for (int i = 0; i < DIM2; i++)
        SSE2 += pow((arr[1][i] - Y_reg2[i]), 2);


    cout << endl << "Суммарная квадратичная ошибка для подобранной функциональной зависимости SSE: " << SSE2 << endl;

    MAE2 = SSE2 / DIM2;

    cout << endl << "Средняя ошибка для подобранной функциональной зависимости MAE: " << MAE2 << endl << endl;

    //Формирование массивов для изображения подобранной функциональной зависимости на графике.
    double y3[DIM2];
    for (int i = 0; i < DIM2; i++)
        y3[i] = (k_a * pow(arr[0][i], 4)) + k_d * arr[0][i] + k_k;

    ofstream f;
    f.open("data.txt");
    for (int i = 0; i < DIM2; i++) {
        f << arr[0][i] << " ";
        f << y2[i] << " ";
        f << y3[i] << " ";
        f << arr[1][i] << "\n";
    }
    f.close();
}
// plot "C:\\Users\\PoPo\\Desktop\\Техпрог\\Практика_Техпрог\\data.txt" u 1:2 w l, "C:\\Users\\PoPo\\Desktop\\Техпрог\\Практика_Техпрог\\data.txt" u 1:3 w l, "C:\\Users\\PoPo\\Desktop\\Техпрог\\Практика_Техпрог\\data.txt" u 1:4 w l, "C:\\Users\\PoPo\\Desktop\\Техпрог\\Практика_Техпрог\\data.txt" u 1:4 

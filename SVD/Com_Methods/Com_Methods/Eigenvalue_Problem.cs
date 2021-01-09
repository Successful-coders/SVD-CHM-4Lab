using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Com_Methods
{
    //исключение i-ой строки и i-го столбка из матрицы А
    class Eigenvalue_Problem
    {
        void Delete_Row_Column(Matrix A, int i)
        {
            //исключаем i-ые строку и столбец
            for (int I = i; I < A.M - 1; I++) 
            {
                A.Elem[I] = A.Elem[I + 1];
                //перемещаем столбцы в I-ой строке
                for (int J = i; J < A.N - 1; J++)
                {
                    A.Elem[I][J] = A.Elem[I][J + 1];
                }
            }
            A.Size_Reduction(A.M - 1, A.N - 1);
        }

        //сдвиг по Рэлею
        double Rayleigh_Shift(Matrix A)
        {
            return A.Elem[A.M - 1][A.M - 1];
        }

        //сдвиг по Уилкинсону
        double Wilcoxon_Shift(Matrix A)
        {
            int N = A.N;
            double a = A.Elem[N - 2][N - 2],
                   b = A.Elem[N - 1][N - 1],
                   c = A.Elem[N - 1][N - 2],
                   d = A.Elem[N - 2][N - 1];
            double D = (a + b) * (a + b) - 4 * (a * b - c * d);
            if (D < 0) throw new Exception("The matrix has complex eigen value...");
            return ((a + b) + Math.Sqrt(D)) * 0.5;
        }

        //реализация сдвига
        void Shift(Matrix A, double shift)
        {
            for (int i = 0; i < A.M; i++) A.Elem[i][i] += shift;
        }

        //реализация алгоритма поиска собственных знаечний: QR - итерации
        public List<double> Find_Eigenvalues_QR_Iterations (Matrix A, QR_Decomposition.QR_Algorithm Method)
        {
            //список собственных значений
            var EigenValues = new List<double>();
            int iter = 0;
            //итерационный процесс
            while(A.M != 1)
            {
                
                //обновление спсика собственных значений
                for (int i = A.M - 1; i > 0; i--)
                {
                    //если элемент A[i][i-1] == 0, A[i][i] - собственное значение
                    if (Math.Abs(A.Elem[i][i - 1]) < 1e-6)
                    {
                        
                        EigenValues.Add(A.Elem[i][i-1]);
                        //исключаем i-ые строку и столбец
                        Delete_Row_Column(A, i);
                        i = A.M;
                    }
                }

                if (A.M == 1) break;

                //прямой сдвиг
                double shift = Rayleigh_Shift(A);
                Shift(A, -shift);

                //QR-разложение
                var QR = new QR_Decomposition(A, Method);

                //новая матрица А
                A = QR.R * QR.Q;

                //обратный сдвиг
                Shift(A, shift);
            }

            //дополняем список последним оставшимся собственным значением
            //EigenValues.Add(A.Elem[0][0]);

            //сортируем по возрастанию и формируем результат
            EigenValues.Sort((double X, double Y) => { if (X > Y) return 0; else return 1; });

            return EigenValues;
        }
    }
}

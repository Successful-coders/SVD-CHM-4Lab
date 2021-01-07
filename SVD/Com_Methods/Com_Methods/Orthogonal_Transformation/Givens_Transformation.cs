using System;

namespace Com_Methods
{
    /// <summary>
    /// преобразования Гивенса
    /// </summary>
    class Givens_Transformation
    {
        /// <summary>
        /// реализация процедуры ортогонализации по методу вращений Гивенса
        /// </summary>
        /// <param name="A - исходная матрица"></param>
        /// <param name="Q - ортогональная матрица преобразований"></param>
        /// <param name="R - результат"></param>
        public static void Givens_Orthogonalization(Matrix A, Matrix Q, Matrix R)
        {
            double help1, help2;

            //косинус, синус
            double c = 0, s = 0;

            //запись матрицы А в R
            R.Copy(A);

            //алгоритм вращения Гивенса: для каждого столбца
            for (int j = 0; j < R.N - 1; j++)
            {
                //просматриваем строки в столбце
                for (int i = j + 1; i < R.M; i++)
                {
                    //если очередной элемент под диагональю не нулевой, то требуется поворот вектора
                    if (Math.Abs(R.Elem[i][j]) > CONST.EPS)
                    {
                        help1 = Math.Sqrt(Math.Pow(R.Elem[i][j], 2) + Math.Pow(R.Elem[j][j], 2));
                        c = R.Elem[j][j] / help1;
                        s = R.Elem[i][j] / help1;

                        //A_new = Gt * A
                        for (int k = j; k < R.N; k++)
                        {
                            help1 = c * R.Elem[j][k] + s * R.Elem[i][k];
                            help2 = c * R.Elem[i][k] - s * R.Elem[j][k];
                            R.Elem[j][k] = help1;
                            R.Elem[i][k] = help2;
                        }

                        //перемножаем строки матрицы Q на трансп.матрицу преобразования Q = Q * G
                        for (int k = 0; k < Q.M; k++)
                        {
                            help1 = c * Q.Elem[k][j] + s * Q.Elem[k][i];
                            help2 = c * Q.Elem[k][i] - s * Q.Elem[k][j];
                            Q.Elem[k][j] = help1;
                            Q.Elem[k][i] = help2;
                        }
                    }
                }
            }
        }

        // метод зануления j-го столбца с i-ой позиции
        public static void Column_Transformation(Matrix A, Matrix U, int j, int i)
        {
            double help1, help2;

            //косинус, синус
            double c = 0, s = 0;

            //просматриваем строки в столбце
            for (int q = i; q < A.M; q++)
            {
                //если очередной элемент под диагональю не нулевой, то требуется поворот вектора
                if (Math.Abs(A.Elem[q][j]) > CONST.EPS)
                {
                    help1 = Math.Sqrt(Math.Pow(A.Elem[q][j], 2) + Math.Pow(A.Elem[j][j], 2));
                    c = A.Elem[j][j] / help1;
                    s = A.Elem[q][j] / help1;

                    //A_new = Gt * A
                    for (int k = j; k < A.N; k++)
                    {
                        help1 = c * A.Elem[j][k] + s * A.Elem[q][k];
                        help2 = c * A.Elem[q][k] - s * A.Elem[j][k];
                        A.Elem[j][k] = help1;
                        A.Elem[q][k] = help2;
                    }

                    //перемножаем строки матрицы Q на трансп.матрицу преобразования Q = Q * G
                    for (int k = 0; k < U.M; k++)
                    {
                        help1 = c * U.Elem[k][j] + s * U.Elem[k][q];
                        help2 = c * U.Elem[k][q] - s * U.Elem[k][j];
                        U.Elem[k][j] = help1;
                        U.Elem[k][q] = help2;
                    }
                }
            }
        }

        // метод зануления i-ой строки с j-ой позиции
        public static void Row_Transformation(Matrix A, Matrix V, int i, int j1)
        {
            double help1, help2;

            //косинус, синус
            double c = 0, s = 0;

            //Ыалгоритм вращения Гивенса: для каждого столбца
            for (int j = j1 + 1; j < A.N; j++)
            {

                //если очередной элемент под диагональю не нулевой, то требуется поворот вектора
                if (Math.Abs(A.Elem[i][j]) > CONST.EPS)
                {
                    help1 = Math.Sqrt(Math.Pow(A.Elem[i][j], 2) + Math.Pow(A.Elem[j][j], 2));
                    c = A.Elem[j][j] / help1;
                    s = A.Elem[i][j] / help1;

                    //A_new = Gt * A
                    for (int k = j; k < A.N; k++)
                    {
                        help1 = c * A.Elem[j][k] + s * A.Elem[i][k];
                        help2 = c * A.Elem[i][k] - s * A.Elem[j][k];
                        A.Elem[j][k] = help1;
                        A.Elem[i][k] = help2;
                    }

                    //перемножаем строки матрицы Q на трансп.матрицу преобразования Q = Q * G
                    for (int k = 0; k < V.M; k++)
                    {
                        help1 = c * V.Elem[k][j] + s * V.Elem[k][i];
                        help2 = c * V.Elem[k][i] - s * V.Elem[k][j];
                        V.Elem[k][j] = help1;
                        V.Elem[k][i] = help2;
                    }
                }
            }
        }
    }
}


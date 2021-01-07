using System;
using System.Collections.Generic;

namespace Com_Methods
{
    //SVD на базе ортогональных преобразований
    class SVD
    {
        //левые сингулярные векторы
        public Matrix U { set; get; }
        //правые сингулярные векторы
        public Matrix V { set; get; }
        //сингулярные числа
        public Matrix Sigma { set; get; }

        //---------------------------------------------------------------------------------------------

        /// <summary>
        /// конструктор
        /// </summary>
        /// <param name="A - матрица для SVD"></param>
        public SVD(Matrix A)
        {
            //вызов метода, который выполнит построение SVD
            Start_SVD(A,1);
        }

        //---------------------------------------------------------------------------------------------

        /// <summary>
        /// проверка на положительность сингулярных чисел
        /// </summary>
        private void Check_Singular_Values()
        {
            //наименьшее измерение
            int Min_Size = Math.Min(Sigma.M, Sigma.N);

            //проверка сингулярных чисел на положительность
            for (int i = 0; i < Min_Size; i++)
            {
                if (Sigma.Elem[i][i] < 0)
                {
                    Sigma.Elem[i][i] = -Sigma.Elem[i][i];

                    for (int j = 0; j < U.M; j++)
                        U.Elem[j][i] = -U.Elem[j][i];
                }
            }
        }

        //---------------------------------------------------------------------------------------------

        /// <summary>
        /// сортировка сингулярных чисел
        /// </summary>
        private void Sort_Singular_Values()
        {
            //наименьшее измерение
            int Min_Size = Math.Min(Sigma.M, Sigma.N);

            //сортировка сингулярных чисел
            for (int I = 0; I < Min_Size; I++)
            {
                var Max_Elem = Sigma.Elem[I][I];
                int Index = I;
                for (int i = I + 1; i < Min_Size; i++)
                {
                    if (Sigma.Elem[i][i] > Max_Elem)
                    {
                        Max_Elem = Sigma.Elem[i][i];
                        Index = i;
                    }
                }
                //найден наибольший элемент
                if (I != Index)
                {
                    Sigma.Elem[Index][Index] = Sigma.Elem[I][I];
                    Sigma.Elem[I][I] = Max_Elem;
                    U.Column_Transposition(I, Index);
                    V.Column_Transposition(I, Index);
                }
            }
        }

        //---------------------------------------------------------------------------------------------

        /// <summary>
        /// редукция SVD 
        /// </summary>
        /// <param name="Reduction - порог отбрасывания сингулярных чисел"></param>
        public void Reduction_SVD(double Reduction)
        {
            //наименьшее измерение
            int Min_Size = Math.Min(Sigma.M, Sigma.N);

            //проверка на возможность редукции по сингулярным числам
            for (int i = 0; i < Min_Size; i++)
            {
                if (Math.Abs(Sigma.Elem[i][i]) < Reduction)
                {
                    Min_Size = i;
                    break;
                }
            }
            //редукция размерности матриц
            Sigma.Size_Reduction(Min_Size, Min_Size);
            U.Size_Reduction(U.M, Min_Size);
            V.Size_Reduction(V.M, Min_Size);
        }

        //---------------------------------------------------------------------------------------------

        /// <summary>
        /// двухэтапный SVD-алгоритм
        /// </summary>
        /// <param name="A - матрица для SVD"></param>
        public void Start_SVD(Matrix A, int method)
        {
            //наименьшее измерение
            int Min_Size = Math.Min(A.M, A.N);

            //размеры нижней и верхней внешних диагоналей
            int Up_Size = Min_Size - 1, Down_Size = Min_Size - 1;

            //инициализация матрицы левых сингулярных векторов 
            U = new Matrix(A.M, A.M);

            //матрица сингулярных чисел
            Sigma = new Matrix(A.M, A.N);

            //инициализация матрицы правых сингулярных векторов 
            V = new Matrix(A.N, A.N);

            //инициализация матриц для SVD
            for (int i = 0; i < A.M; i++)
            {
                U.Elem[i][i] = 1.0;
                for (int j = 0; j < A.N; j++) Sigma.Elem[i][j] = A.Elem[i][j];
            }
            for (int i = 0; i < A.N; i++) V.Elem[i][i] = 1.0;

            //**************** этап I: бидиагонализация *************************

            for (int i = 0; i < Min_Size - 1; i++)
            {
                Householder_Transformation.Column_Transformation(Sigma, U, i, i);
                Householder_Transformation.Row_Transformation(Sigma, V, i, i + 1);
            }

            //ситуация M > N - строк больше => дополнительное умножение слева 
            if (A.M > A.N)
            {
                Householder_Transformation.Column_Transformation(Sigma, U, A.N - 1, A.N - 1);
                //нижняя побочная диагональ длиннее на 1
                Down_Size += 1;
            }

            //ситуация M < N - столбцов больше => дополнительное умножение справа
            if (A.M < A.N)
            {
                Householder_Transformation.Row_Transformation(Sigma, V, A.M - 1, A.M);
                //верхняя побочная диагональ длиннее на 1
                Up_Size += 1;
            }

            //**************** этап II: преследование ************
            //********* приведение к диагональному виду **********
            switch (method)
            {
                case 1:
                    {
                        //для хранение изменяющихся элементов верхней диагонали
                        var Up = new double[Up_Size];
                        var Down = new double[Down_Size];
                        //число неизменившихся элементов над главной диагональю
                        int CountUpElements;

                        //процедура преследования
                        do
                        {


                            CountUpElements = 0;

                            //обнуление верхней диагонали
                            for (int i = 0; i < Up_Size; i++)
                            {
                                if (Math.Abs(Up[i] - Sigma.Elem[i][i + 1]) > CONST.EPS)
                                {
                                    Up[i] = Sigma.Elem[i][i + 1];
                                    Householder_Transformation.Row_Transformation(Sigma, V, i, i);
                                }
                                else
                                    CountUpElements++;
                            }

                            //обнуление нижней диагонали
                            for (int i = 0; i < Down_Size; i++)
                            {
                                if (Math.Abs(Down[i] - Sigma.Elem[i + 1][i]) > CONST.EPS)
                                {
                                    Down[i] = Sigma.Elem[i + 1][i];
                                    Householder_Transformation.Column_Transformation(Sigma, U, i, i);
                                }
                            }

                        }
                        while (CountUpElements != Up_Size);

                        //убираем отрицательные сингулярные числа
                        Check_Singular_Values();
                        //сортируем по возрастанию сингулярные числа
                        Sort_Singular_Values();
                        break;
                    }

                case 2:
                    {
                        //для хранение изменяющихся элементов верхней диагонали
                        var Up = new double[Up_Size];
                        var Down = new double[Down_Size];
                        //число неизменившихся элементов над главной диагональю
                        int CountUpElements;

                        //процедура преследования
                        do
                        {


                            CountUpElements = 0;

                            //обнуление верхней диагонали
                            for (int i = 0; i < Up_Size; i++)
                            {
                                if (Math.Abs(Up[i] - Sigma.Elem[i][i + 1]) > CONST.EPS)
                                {
                                    Up[i] = Sigma.Elem[i][i + 1];
                                    Givens_Transformation.Row_Transformation(Sigma, V, i, i);
                                }
                                else
                                    CountUpElements++;
                            }

                            //обнуление нижней диагонали
                            for (int i = 0; i < Down_Size; i++)
                            {
                                if (Math.Abs(Down[i] - Sigma.Elem[i + 1][i]) > CONST.EPS)
                                {
                                    Down[i] = Sigma.Elem[i + 1][i];
                                    Givens_Transformation.Column_Transformation(Sigma, U, i, i+1);
                                }
                            }

                        }
                        while (CountUpElements != Up_Size);

                        //убираем отрицательные сингулярные числа
                        Check_Singular_Values();
                        //сортируем по возрастанию сингулярные числа
                        Sort_Singular_Values();
                        break;
                    }

            }
        }
           
    
        //---------------------------------------------------------------------------------------------

        //ранг матрицы
        public int Rank()
        {
            return Sigma.M;
        }

        //---------------------------------------------------------------------------------------------

        //модуль определителя матрицы
        public double Abs_Det()
        {
            //ранг матрицы
            int Size = Rank();

            if (Size == 0) throw new Exception("Error in SVD.Rank: SVD is not built ...");
            double det = 1;
            for (int i = 0; i < Size; i++)
                det *= Sigma.Elem[i][i];
            return det;
        }

        //---------------------------------------------------------------------------------------------

        //число обусловленности матрицы
        public double Cond()
        {
            //ранг матрицы
            int Size = Rank();

            if (Size == 0) throw new Exception("Error in SVD.Rank: SVD is not built ...");
            return Sigma.Elem[0][0] / Sigma.Elem[Size - 1][Size - 1];
        }

        //---------------------------------------------------------------------------------------------
        

        //нормальное псевдорешение СЛАУ Ax = F => x = V * S^(-1) * Ut * F
        public Vector Start_Solver(Vector F)
        {
            //ранг матрицы
            int Size = Rank();

            if (Size == 0) throw new Exception("Error in SVD.Rank: SVD is not built ...");

            //UtF = Ut * F
            var UtF = U.Multiplication_Trans_Matrix_Vector(F);

            //UtF = S^(-1) * Ut * F
            for (int i = 0; i < UtF.N; i++) UtF.Elem[i] /= Sigma.Elem[i][i];

            //Res = V * S^(-1) * Ut * F
            var Res = V * UtF;
            return Res;
        }
    }
}

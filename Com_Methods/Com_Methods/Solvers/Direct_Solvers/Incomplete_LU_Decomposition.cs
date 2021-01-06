using System;

namespace Com_Methods
{
    /// <summary>
    /// интерфейс неполного LU-разложения
    /// </summary>
    public interface IIncomplete_LU_Decomposition
    {
        //решение СЛАУ L * x = F с нижним треугольником матрицы
        void SLAU_L(Vector X, Vector F);
        //решение СЛАУ Lt * x = F с нижним транспонированным треугольником матрицы
        void SLAU_Lt(Vector X, Vector F);
        //решение СЛАУ U * x = F с верхним треугольником матрицы
        void SLAU_U(Vector X, Vector F);
        //решение СЛАУ Ut * x = F с верхним транспонированным треугольником матрицы
        void SLAU_Ut(Vector X, Vector F);
    }
    
    /// <summary>
    /// неполное LU-разложение в формате CSlR
    /// </summary>
    class Incomplete_LU_Decomposition_CSlR : IIncomplete_LU_Decomposition
    {
        //матрица неполной LU-декомпозиции
        private CSlR_Matrix ILU;

        //конструктор
        public Incomplete_LU_Decomposition_CSlR(CSlR_Matrix A)
        {
            ILU = new CSlR_Matrix();
            ILU.N = A.N;
            ILU.iptr = A.iptr;
            ILU.jptr = A.jptr;
            int N_autr = A.autr.Length;
            ILU.autr = new double[N_autr];
            ILU.altr = new double[N_autr];
            ILU.di = new double[A.N];

            for (int i = 0; i < N_autr; i++) { ILU.altr[i] = A.altr[i]; ILU.autr[i] = A.autr[i]; }
            for (int i = 0; i < A.N; i++) { ILU.di[i] = A.di[i]; }

            //начинаем с i = 1, т.к. в первой строке нижнего треугольника только диагональный элемент
            for (int i = 1; i < A.N; i++)
            {
                for (int j = A.iptr[i] - 1; j < A.iptr[i + 1] - 1; j++)
                {
                    for (int a = A.iptr[i] - 1; a < j; a++)
                    {
                        for (int b = A.iptr[A.jptr[j] - 1] - 1; b < A.iptr[A.jptr[j]] - 1; b++)
                        {
                            if (A.jptr[a] == A.jptr[b])
                            {
                                ILU.altr[j] -= ILU.altr[a] * ILU.autr[b];
                                ILU.autr[j] -= ILU.autr[a] * ILU.altr[b];
                            }
                        }
                    }
                    ILU.autr[j] /= ILU.di[A.jptr[j] - 1];
                    ILU.di[i] -= ILU.autr[j] * ILU.altr[j];
                }
            }
        }

        //-------------------------------------------------------------------------------------------------

        //решение СЛАУ L * x = F с нижним треугольником матрицы
        public void SLAU_L(Vector X, Vector F)
        {
            for (int i = 0; i < ILU.N; i++)
            {
                X.Elem[i] = F.Elem[i];
                for (int j = ILU.iptr[i] - 1; j < ILU.iptr[i + 1] - 1; j++)
                    X.Elem[i] -= X.Elem[ILU.jptr[j] - 1] * ILU.altr[j];
                X.Elem[i] /= ILU.di[i];
            }
        }

        //-------------------------------------------------------------------------------------------------

        //решение СЛАУ Lt * x = F с нижним транспонированным треугольником матрицы
        public void SLAU_Lt(Vector X, Vector F)
        {
            double[] V = new double[ILU.N];
            for (int i = 0; i < ILU.N; i++) V[i] = F.Elem[i];
            for (int i = ILU.N - 1; i >= 0; i--)
            {
                X.Elem[i] = V[i] / ILU.di[i];
                for (int j = ILU.iptr[i] - 1; j < ILU.iptr[i + 1] - 1; j++)
                    V[ILU.jptr[j] - 1] -= X.Elem[i] * ILU.altr[j];
            }
        }

        //-------------------------------------------------------------------------------------------------

        //решение СЛАУ U * x = F с верхним треугольником матрицы
        public void SLAU_U(Vector X, Vector F)
        {
            for (int i = 0; i < ILU.N; i++) X.Elem[i] = F.Elem[i];
            for (int i = ILU.N - 1; i >= 0; i--)
            {
                for (int j = ILU.iptr[i] - 1; j < ILU.iptr[i + 1] - 1; j++)
                {
                    X.Elem[ILU.jptr[j] - 1] -= X.Elem[i] * ILU.autr[j];
                }
            }
        }

        //-------------------------------------------------------------------------------------------------

        //решение СЛАУ Ut * x = F с верхним транспонированным треугольником матрицы
        public void SLAU_Ut(Vector X, Vector F)
        {
            for (int i = 0; i < ILU.N; i++)
            {
                X.Elem[i] = F.Elem[i];
                for (int j = ILU.iptr[i] - 1; j < ILU.iptr[i + 1] - 1; j++)
                    X.Elem[i] -= X.Elem[ILU.jptr[j] - 1] * ILU.autr[j];
            }
        }
    }
}

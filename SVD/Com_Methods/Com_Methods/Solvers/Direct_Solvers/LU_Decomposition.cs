using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Com_Methods
{
    //LU-разложение
    class LU_Decomposition
    {
        public Matrix LU { set; get; }

        //получение LU матрицы
        public void CreateLU(Matrix A)
        {
            LU = new Matrix(A.M, A.N);

            //копирование матрицы A в LU, чтобы не испортить A
            for (int i = 0; i < A.M; i++)
                for (int j = 0; j < A.N; j++)
                    LU.Elem[i][j] = A.Elem[i][j];

            //прямой ход метода Гаусса для получения верхнего треугольника
            //являющегося верхней треугольной матрицей U
            Gauss_Method.Direct_Way(LU);

            //формирование нижнего треугольника
            //являющегося нижней треугольной матрицей L
            for(int i = 0; i < A.M; i++)
            {
                for(int j = 0; j < i; j++)
                {
                    double sumLikUkj = 0;
                    for (int k = 0; k < j; k++)
                        sumLikUkj += LU.Elem[i][k] * LU.Elem[k][j];

                    LU.Elem[i][j] = (A.Elem[i][j] - sumLikUkj) / LU.Elem[j][j];
                }
            }
        }

        //прямой ход, работающий с объектом классам Matrix, для обработки нижнего треугольника
        public Vector Direct_Way(Matrix A, Vector RES)
        {
            for (int i = 0; i < A.M; i++)
            {
                for (int j = 0; j < i; j++)
                    RES.Elem[i] -= A.Elem[i][j] * RES.Elem[j];
            }
            return RES;
        }

        //обратный ход, работающий с объектом классам Matrix, для обработки нижнего треугольника
        public Vector Back_Way(Matrix A, Vector RES)
        {
            for (int i = RES.N - 1; i >= 0; i--)
            {
                for (int j = i + 1; j < RES.N; j++)
                    RES.Elem[i] -= A.Elem[i][j] * RES.Elem[j];

                RES.Elem[i] /= A.Elem[i][i];
            }
            return RES;
        }

        //решение СЛАУ
        public Vector Start_Solver(Matrix A, Vector F)
        {
            var RES  = Direct_Way(LU,F);
            Back_Way(LU,RES);

            return RES;
        }
    }
}
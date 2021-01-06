using System;

namespace Com_Methods
{
    //диагональный предобусловливатель (Якоби)
    class Diagonal_Preconditioner : Preconditioner
    {
        //диагональ матрицы
        Vector Diag { get; }

        //конструктор диагонального преобусловливателя
        public Diagonal_Preconditioner(CSlR_Matrix A)
        {
            Diag = new Vector(A.N);
            for (int i = 0; i < A.N; i++)
            {
                if (Math.Abs(A.di[i]) < CONST.EPS)
                    throw new Exception("Error in Diagonal_Preconditioner: di" + (i + 1).ToString() + " = " + A.di[i].ToString());
                Diag.Elem[i] = A.di[i];
            }
        }

        //реализация преобусловливателя
        public override void Start_Preconditioner(Vector X, Vector RES)
        {
            for (int i = 0; i < Diag.N; i++)
            {
                RES.Elem[i] = X.Elem[i] / Diag.Elem[i];
            }
        }

        //реализация транспонированного преобусловливателя
        public override void Start_Tr_Preconditioner(Vector X, Vector RES)
        {
            for (int i = 0; i < Diag.N; i++)
            {
                RES.Elem[i] = X.Elem[i] / Diag.Elem[i];
            }
        }
    }
}

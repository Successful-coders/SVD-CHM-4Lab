using System;

namespace Com_Methods
{
    //LU-предобусловливатель (Incomplete LU-decomposition)
    class LU_Preconditioner : Preconditioner
    {
        //неполная LU-декомпозиция
        IIncomplete_LU_Decomposition ILU;
        
        //конструктор LU-преобусловливателя
        public LU_Preconditioner(CSlR_Matrix A)
        {
            ILU = new Incomplete_LU_Decomposition_CSlR(A);
        }

        //реализация преобусловливателя
        public override void Start_Preconditioner(Vector X, Vector RES)
        {
            //решаем СЛАУ с нижней треугольной матрицей
            ILU.SLAU_L(RES, X);
            //решаем СЛАУ с верхней треугольной матрицей
            ILU.SLAU_U(RES, RES);
        }

        //реализация транспонированного преобусловливателя
        public override void Start_Tr_Preconditioner(Vector X, Vector RES)
        {
            //решаем СЛАУ с нижней треугольной матрицей
            ILU.SLAU_Ut(RES, X);
            //решаем СЛАУ с верхней треугольной матрицей
            ILU.SLAU_Lt(RES, RES);
        }
    }
}

using System;

namespace Com_Methods
{
    //абстрактный класс предобусловливателя
    abstract class Preconditioner
    {
        //реализация предобусловливателя
        abstract public void Start_Preconditioner(Vector X, Vector RES);
        //реализация транспонированного предобусловливателя
        abstract public void Start_Tr_Preconditioner(Vector X, Vector RES);
        //тип предобусловливателей
        public enum Type_Preconditioner
        {
            Diagonal_Preconditioner = 1,
            LU_Decomposition
        }
    }
}
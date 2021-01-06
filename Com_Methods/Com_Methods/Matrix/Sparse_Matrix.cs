using System;
using System.IO;

namespace Com_Methods
{
    public interface ISparse_Matrix
    {
        //размер матрицы
        int N { set; get; }
        //умножение матрицы на вектор y = A * x
        void Mult_MV(Vector X, Vector Y);
        //умножение транспонированной матрицы на вектор y = At * x
        void Mult_MtV(Vector X, Vector Y);
    }

    //матрица в разреженном строчно-столбцовом формате CSlR
    class CSlR_Matrix : ISparse_Matrix
    {
        //размер матрицы
        public int N { set; get; }
        //диагональ матрицы
        public double[] di {set; get;}
        //нижний треугольник
        public double[] altr { set; get; }
        //верхний треугольник
        public double[] autr { set; get; }
        //номера строк (столбцов) ненулевых элементов
        public int[] jptr { set; get; }
        //номера строк (столбцов), с которых начинается jptr
        public int[] iptr { set; get; }

        //конструктор по умолчанию
        public CSlR_Matrix() { }
       
        //конструктор по файлам
        public CSlR_Matrix(string PATH)
        {
            char[] Separator = new char[] {' '};

            //размер системы
            using (var Reader = new StreamReader(File.Open(PATH + "Size.txt", FileMode.Open)))
            {
                N = Convert.ToInt32(Reader.ReadLine());
                //выделение памяти под массивы di и iptr
                iptr = new int[N + 1];
                di   = new double[N];
            }

            //диагональ матрицы
            using (var Reader = new StreamReader(File.Open(PATH + "di.txt", FileMode.Open)))
            {
                for (int i = 0; i < N; i++)
                {
                    di[i] = Convert.ToDouble(Reader.ReadLine().Split(Separator, StringSplitOptions.RemoveEmptyEntries)[0]);
                }
            }

            //массив iptr
            using (var Reader = new StreamReader(File.Open(PATH + "iptr.txt", FileMode.Open)))
            {
                for (int i = 0; i <= N; i++)
                {
                    iptr[i] = Convert.ToInt32(Reader.ReadLine().Split(Separator, StringSplitOptions.RemoveEmptyEntries)[0]);
                }
            }

            //выделение памяти под массивы jptr, altr, autr
            int Size = iptr[N] - 1;
            jptr = new int[Size];
            altr = new double[Size];
            autr = new double[Size];
            var Reader1 = new StreamReader(File.Open(PATH + "jptr.txt", FileMode.Open));
            var Reader2 = new StreamReader(File.Open(PATH + "altr.txt", FileMode.Open));
            var Reader3 = new StreamReader(File.Open(PATH + "autr.txt", FileMode.Open));
            for (int i = 0; i < Size; i++)
            {
                jptr[i] = Convert.ToInt32 (Reader1.ReadLine().Split(Separator, StringSplitOptions.RemoveEmptyEntries)[0]);
                altr[i] = Convert.ToDouble(Reader2.ReadLine().Split(Separator, StringSplitOptions.RemoveEmptyEntries)[0]);
                autr[i] = Convert.ToDouble(Reader3.ReadLine().Split(Separator, StringSplitOptions.RemoveEmptyEntries)[0]);
            }
            Reader1.Close(); Reader2.Close(); Reader3.Close();
        }

        //-------------------------------------------------------------------------------------------------

        //умножение матрицы на вектор y = A * x
        public void Mult_MV (Vector X, Vector Y)
        {
            for (int i = 0; i < N; i++) Y.Elem[i] = X.Elem[i] * di[i];
            for (int i = 0; i < N; i++)
                for (int j = iptr[i] - 1; j < iptr[i + 1] - 1; j++)
                {
                    Y.Elem[i] += X.Elem[jptr[j] - 1] * altr[j];
                    Y.Elem[jptr[j] - 1] += X.Elem[i] * autr[j];
                }
        }

        //-------------------------------------------------------------------------------------------------

        //умножение транспонированной матрицы на вектор y = At * x
        public void Mult_MtV(Vector X, Vector Y)
        {
            for (int i = 0; i < N; i++) Y.Elem[i] = X.Elem[i] * di[i];
            for (int i = 0; i < N; i++)
                for (int j = iptr[i] - 1; j < iptr[i + 1] - 1; j++)
                {
                    Y.Elem[i] += X.Elem[jptr[j] - 1] * autr[j];
                    Y.Elem[jptr[j] - 1] += X.Elem[i] * altr[j];
                }
        }
    }
}
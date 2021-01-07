using System;
using System.Diagnostics;
using System.Threading;


namespace Com_Methods
{
    class CONST
    {
        //точность решения
        public static double EPS = 1e-15;
    }

    class Tools
    {
        //замер времени
        public static string Measurement_Time(Thread thread)
        {
            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();
            thread.Start();
            while (thread.IsAlive) ;
            stopWatch.Stop();
            TimeSpan ts = stopWatch.Elapsed;
            string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds / 10);
            return ("RunTime: " + elapsedTime);
        }

        //относительная невязка
        public static double Relative_Discrepancy (Matrix A, Vector X, Vector F)
        {
            var x = A * X;
            for (int i = 0; i < x.N; i++) x.Elem[i] -= F.Elem[i];
            return x.Norm() / F.Norm();
        }
        
    }

    class Program
    {
        
        static void Main()
        {
            try
            {
                //SVD: анализ с унитарными преобразованиями
                var SVD_Thread = new Thread(() =>
                {
                    //Excel
                    var Excel_WR = new Excel_Reader_Writer();
                    
                    //плотная матрица из файла Excel
                    var A = Excel_WR.Read_Matrix("C:\\Users\\Жопчики\\Desktop\\Даша\\projects-with-git\\чм 4 лаба\\SVD\\Com_Methods\\Com_Methods\\Data\\Excel_Data\\2a.xlsx");



                    //открываем Excel
                    //Excel_WR.Write_Matrix(A, "Matrix A");

                    //Binary_Reader_Writer.Save_Object(A, "C:\\Users\\Жопчики\\Desktop\\Даша\\projects-with-git\\чм 4 лаба\\SVD\\Com_Methods\\Com_Methods\\Data\\Original.bin");

                    var SVD = new SVD(A);
                    SVD.Reduction_SVD(1e-15);

                    Matrix Res = SVD.U * SVD.Sigma * SVD.V.Transpose_Matrix();
                   // Excel_WR.Write_Matrix(SVD.U * SVD.Sigma * SVD.V.Transpose_Matrix(), "Compression Matrix");
                    //double qmax = 0;
                    //for (int i = 0; i < A.M; i++)
                    //    for (int j = 0; j < A.N; j++)
                    //    {
                    //        if (A.Elem[i][j] > CONST.EPS)
                    //        {
                    //            double q = Math.Abs(A.Elem[i][j] - Res.Elem[i][j]) / Math.Abs(A.Elem[i][j]);
                    //            if (qmax < q) qmax = q;
                    //        }
                    //    }
                    //Console.WriteLine("Quality:"+qmax);
                    
                   // Binary_Reader_Writer.Save_ObjectThree(SVD.Sigma, SVD.U, SVD.V, "C:\\Users\\Жопчики\\Desktop\\Даша\\projects-with-git\\чм 4 лаба\\SVD\\Com_Methods\\Com_Methods\\Data\\SVD.bin");

                    //Console.WriteLine("\nMatrix U:");
                    //SVD.U.Console_Write_Matrix();

                    //Console.WriteLine("\nMatrix Sigma (size {0}):", SVD.Sigma.M);
                    //for(int i = 0; i < SVD.Sigma.M; i++) Console.WriteLine(SVD.Sigma.Elem[i][i]);

                    //Console.WriteLine("\nMatrix V:");
                    //SVD.V.Console_Write_Matrix();

                    //Console.WriteLine("\nMatrix A = U * Sigma * Vt:");
                    //(SVD.U * SVD.Sigma * SVD.V.Transpose_Matrix()).Console_Write_Matrix();

                    ////ранг матрицы
                    //Console.WriteLine("\nrk(A) = " + SVD.Rank());

                    ////определитель матрицы
                    //Console.WriteLine("\n|det(A)| = " + SVD.Abs_Det());

                    //число обусловленности
                      Console.WriteLine("\nCond(A) = " + SVD.Cond());

                    //число обусловленности
                    // Console.WriteLine("\nCond2(A) = " + A.Cond_Square_Matrix_Parallel());

                    //число обусловленности
                    // Console.WriteLine("\nCondInfinity(A) = " + A.Cond_InfinityNorm());

                    //решение СЛАУ
                    //var X_True = new Vector(A.N);
                    //for (int i = 0; i < A.N; i++) X_True.Elem[i] = 1.0;
                    //var F = A * X_True;
                    //var X = SVD.Start_Solver(F);
                    //Console.WriteLine("\nResult:");
                    //foreach (var el in X.Elem) Console.WriteLine(el);
                    ////норма невязки ||Ax - f||
                    //Console.WriteLine("\nRelative_Discrepancy: {0}", Tools.Relative_Discrepancy(A, X, F));


                });

                //время работы программы
                Console.WriteLine("\n" + Tools.Measurement_Time(SVD_Thread));

            }
            catch (Exception E)
            {
                Console.WriteLine("\n*** Error! ***");
                Console.WriteLine("Method:  {0}",   E.TargetSite);
                Console.WriteLine("Message: {0}\n", E.Message);
            }
        }
    }
}
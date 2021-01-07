using System;
using System.IO;
using System.Runtime.Serialization.Formatters.Binary;

namespace Com_Methods
{
    class Binary_Reader_Writer
    {
        /// <summary>
        /// сериализация объекта
        /// Obj - сохраняемый объект
        /// PATH - путь к файлу для записи
        ///</summary>
        public static void Save_Object (Object Obj, string PATH)
        {
            BinaryFormatter binFormat = new BinaryFormatter();
            using (Stream fStream = new FileStream(PATH, FileMode.Create, FileAccess.Write, FileShare.None))
            {
                binFormat.Serialize(fStream, Obj);
            }
        }

        public static void Save_ObjectThree(Object Obj1, Object Obj2, Object Obj3, string PATH)
        {
            BinaryFormatter binFormat = new BinaryFormatter();
            using (Stream fStream = new FileStream(PATH, FileMode.Create, FileAccess.Write, FileShare.None))
            {
                binFormat.Serialize(fStream, Obj1);
                binFormat.Serialize(fStream, Obj2);
                binFormat.Serialize(fStream, Obj3);
            }
        }

        /// <summary>
        /// десериализация объекта
        /// вернёт объект из файла
        /// PATH - путь к файлу для записи
        ///</summary>
        public static Object Load_Object (string PATH)
        {
            Object Obj = null;
            BinaryFormatter binFormat = new BinaryFormatter();
            using (Stream fStream = File.OpenRead(PATH))
            {
                Obj = binFormat.Deserialize(fStream);
            }
            return Obj;
        }
    }
}

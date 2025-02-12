using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace A
{
    class Program
    {
        static void Main(string[] args)
        {
            try
            {
                string[] input = File.ReadAllLines("input.txt");
                double[] cords = input[0].Split().Select(double.Parse).ToArray;
                double w = double.Parce(input[1]);
                double finish = input[2].Split().Select(double.Parse).ToArray;
                int n = int.Parce(input[3]);
                List<Tuple<double, double, double>> turns = new List<Tuple<double, double, double>>();
                for (int i=0;i<n;i++)
                {
                    double[] turn = input[4+i].Split().Select(double.Parse).ToArray();
                    turns.Add(Tuple.Create(turn[0], turn[1], turn[2]));
                }
                double x1 = cords[0], y1 = cords[1];
                double x1 = cords[2], y1 = cords[3];
                double d = x1-x2;
                double L = y1-y2;
                if (d <= 0 || L <= 0)
                {
                    Console.WriteLine(-1); // Некорректные данные (например, колёса "перепутаны")
                    return;
                }

                List<double> alphas = new List<double>(); // Углы поворота для каждого поворота
                bool possible = true; // Флаг возможности прохождения всех поворотов

                // Обработка каждого поворота
                foreach (var turn in turns)
                {
                    double Xi = turn.Item1, Yi = turn.Item2, Ri = turn.Item3;
                    
                    // Расчёт радиуса траектории центра автомобиля:
                    // R_center = Внешний радиус поворота (Ri) - ширина трека (w) + половина ширины колеи (d/2)
                    double R_center = Ri - w + d / 2;

                    // Проверка на минимально допустимый радиус:
                    // Если R_center слишком мал, автомобиль не сможет повернуть
                    if (R_center <= d / 2)
                    {
                        possible = false;
                        break;
                    }

                    // Расчёт угла поворота левого колеса:
                    // Используется геометрия поворота: α = arctg(колесная база / (R_center - d/2))
                    double alpha = Math.Atan(L / (R_center - d / 2));
                    alphas.Add(alpha);
                }

                // Если хотя бы один поворот невозможен, вывод -1
                if (!possible)
                {
                    Console.WriteLine(-1);
                    return;
                }
                double s = 0.0;

                // Вывод общей длины траектории с точностью до 5 знаков
                Console.WriteLine(s.ToString("F5"));
                
                // Вывод углов для каждого поворота
                foreach (var alpha in alphas)
                {
                    Console.WriteLine(alpha.ToString("F5"));
                }
            }
            catch (Exception)
            {
                Console.WriteLine(-1);
            }
        }
    }
}


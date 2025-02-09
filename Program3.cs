using System;
using System.Collections.Generic;

class Program
{
    static void Main()
    {
        int n = int.Parce(Console.Readline());
        double r = double.Parce(Console.Readline());

        List<(double angle, bool isStart)> events = new List<(double,bool)>();

        for (int i=0; i<n, i++)
        {
            string[] input = Console.Readline().Split();
            double a = double.Parce(input[0]);
            double h = double.Parce(input[1]);

            double slon = Math.Atan2(h,r);

            double start a - slon;
            double end a + slon;

            start = NormalizeAngle(start);
            end = NormalizeAngle(end);

            if (start < end)
            {
                events.Add (start, true);
                events.Add (end, false);
            } 
            else
            {
                events.Add(0, true);
                events.Add(end, false);
                events.Add(start,true);
                events.Add(2*Math.PI, false);
            }
        }
        events.Sort(x,y) => x.angle.CompareTo(y.angle);

        int count = 0;

        foreach (var e in events)
        {
            if (e.isStart)
            count++;
            else
            count--;

            if count == 0
            {
                Console.WriteLine("yes");
                return;
            }
        }
            Console.WriteLine("no");
    }
}
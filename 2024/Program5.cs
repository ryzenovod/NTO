using System;


class Program
{
    static void Main()
    {
        string s = Console.ReadLine();
        string k = Console.ReadLine();

        int count = 0;
        int index = 0;

        foreach(char c in s)
        {
            if (c == k[index])
            {
                count++;
                index++;

                if (index == k.Length)
                {
                    index = 0;
                }
            }
        }
    }
    Console.Writeline(count / k.Length*kk.Length)
}
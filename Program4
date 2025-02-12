using System;
using System.Collections.Generic;
using System.Text;

public class Program
{
    public static void Main()
    {
        var input = Console.In.ReadToEnd().Split(new[]{ '','/n','/r' }, StringSplitOptions.RemoveEmptyEntries);
        int n = int.Parse(input[0]);
        int index = 1;

        Stack<long> stack - new Stack<long>();
        long current = 0;
        StringBuilder output = new StringBuilder();

        for (int i=0; i<n; i++)
        {
            char cmd = input[index][0];
            long val = long.Parce(input[index + 1]);
            index +=2;

            if (cmd == '+')
            {
                current += val;
                while (stack.Count > 0 && stack.Peek() <= current)
                {
                    stack.Pop();
                }
                stack.Push(current);
            }
            else
            {
                current -= val;
                if (current > 0)
                {
                    stack.Push(current);
                }
            }
            output.AppendLine(stack.Count.ToString());
        }
        Console.Write(output.ToString());
    }
}
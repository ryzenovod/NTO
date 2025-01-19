using System;
using System.Collections.Generic;
using System.Linq;

class Program
{
    static void Main()
    {
        int[] integers ={4,3,6,7,10,26,99,40,21,14,12,22,148,29,18,81,28,31,27,20};

        List<string> binaryString = integers.Select(n => Convert.ToString(n,2)).ToList();

        string concatenated = string.Join("", binaryString);

        string hex = ConvertBinaryToHex(concatenated);

        Console.WriteLine(hex);
    }

    static string ConvertBinaryToHex (string binary)
    {
        int bites = 4- (binary.Length % 4);
        if (bites != 4)
        {
            binary = binary.PadLeft(binary.Length + bites, '0');
        }
        List<string> hexList = new List<string>();
        for (int i=0; i< binary.Length; i+=4)
        {
            string fourBites = binary.Substring(i,4);
            int decimalValue = Convert.ToInt32(fourBites, 2);
            hexList.Add(decimalValue.ToString("x"));
        }

        return string.Join("", hexList);
    }

}

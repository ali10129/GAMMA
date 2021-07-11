using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections;

namespace GammaV1
{
    class Program
    {
        public struct Element
        {
            public float Data;
            public uint row;
            public uint column;
            public bool visited;

            public Element clear()
            {
                Data = 0;
                row = 0;
                column = 0;
                visited = false;
                return this;
            }
        }

        public class PE
        {
            public static uint R = 64, bufferSize = 64;

            Element[] Amk_ScalingFactorRegFile;
            Element[,] Bkn;
            int id;
            public static float[,] Cmn;


            public PE(int M , int N,int id)
            {
                this.id = id;
                Amk_ScalingFactorRegFile = new Element[R];
                for (uint i = 0; i < R; i++)
                {
                    Amk_ScalingFactorRegFile[i].clear();
                    Amk_ScalingFactorRegFile[i].Data = 1;
                    Amk_ScalingFactorRegFile[i].column = i;
                }

                Bkn = new Element[R, bufferSize];
                for (int i = 0; i < R; i++)
                {
                    for (int j = 0; j < bufferSize; j++)
                    {
                        Bkn[i, j].clear();
                    }
                }
            }
            public float[,] calc()
            {
                Queue<Element> OUT = merger();
                Element fromQueue, fromMul, toC;
                toC = new Element();
                toC.clear();
                toC.row = (uint)id;
                for (int i = 0; i < OUT.Count; i++)
                {
                    fromQueue = OUT.Dequeue();
                    foreach (var item in Amk_ScalingFactorRegFile)
                    {
                        if (item.column == fromQueue.row)
                        {
                            fromMul = Mul(fromQueue, item);
                            if(fromMul.column != toC.column)
                            {
                                /***********/

                                Cmn[toC.row, toC.column] += toC.Data;
                                toC = fromMul;
                            }
                            else
                            {
                                toC.Data += fromMul.Data;
                            }
                        }
                    }
                    OUT.Enqueue(fromQueue);
                }
                Cmn[toC.row, toC.column] += toC.Data;

                return Cmn;
            }

            public int loadA(Element[] Am, int start = 0)
            {
                int count = 0;
                for (int i = 0; i < R; i++)
                {
                    Amk_ScalingFactorRegFile[i] = Am[i + start];
                    count = i + start;
                }
                return count;
            }

            public int loadB(Element[,] B, int startRow = 0, int startColumn = 0)
            {
                for (int i = startRow; i < startRow + R; i++)
                {
                    for (int j = startColumn; j < startColumn + bufferSize; j++)
                    {

                        Bkn[i, j] = B[i, j];
                    }
                }
                return 0;
            }

            public Queue<Element> merger()
            {
                Queue<Element> elements = new Queue<Element>();

                uint min_index = uint.MaxValue;
                uint min_value = uint.MaxValue;
                uint j = 0;

                uint counter = R;
                

                for (int i = 0; i < bufferSize; i++)
                {
                    while (counter > 0)
                    {
                        if (j < R)
                        {
                            if(Bkn[j,i].visited == false && Bkn[j, i].column < min_value)
                            {
                                min_index = j;
                            }
                            j++;
                        }
                        else
                        {
                            Bkn[min_index, i].visited = true;
                            elements.Enqueue(Bkn[min_index, i]);

                            counter--;
                            j = 0;
                            min_index = uint.MaxValue;
                            min_value = uint.MaxValue;
                        }
                    }
                    counter = R;
                }
                return elements;
            }

            public Element Mul(Element From_Merger, Element From_Amk_ScalinFactor)
            {
                Element temp = new Element();

                temp.Data = From_Amk_ScalinFactor.Data * From_Merger.Data;
                temp.row = From_Amk_ScalinFactor.row;
                temp.column = From_Merger.column;

                return temp;
 
            }
        }
        static void Main(string[] args)
        {
            Random _random = new Random();
            int m, n, k;
            m = 64;
            n = 64;
            k = 64;

            float[,] A = new float[m, k];
            Element[][] Am = new Element[m][];

            float[,] B = new float[k, n];
            Element[,] Bn = new Element[k, n];

            PE.Cmn = new float[m, n];
            float[,] C1 = new float[m, n];
            float[,] C2 = new float[m, n];

            //step 1
            int tmp = 0;
            for (int i = 0; i < m; i++)
            {
                Am[i] = new Element[k];
                for (int j = 0; j < k; j++)
                {
                    tmp = _random.Next(5);
                    A[i, j] = tmp < 5 ? tmp : 0;
                    Am[i][j].Data = A[i, j];
                    Am[i][j].row = (uint)i;
                    Am[i][j].column = (uint)j;

                }
            }
            for (int i = 0; i < k; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    tmp = _random.Next(5);
                    B[i, j] = tmp < 5 ? tmp : 0;

                    Bn[i, j].Data = B[i, j];
                    Bn[i, j].row = (uint)i;
                    Bn[i, j].column = (uint)j;
                }
            }

            PE[] PEs = new PE[m];
            for (int i = 0; i < m; i++)
            {
                PEs[i] = new PE(m, n, i);
                PEs[i].loadA(Am[i]);
                PEs[i].loadB(Bn);
                PEs[i].calc();
            }
            C1 = PE.Cmn;

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    for (int l = 0; l < k; l++)
                    {
                        C2[i, j] += A[i, l] * B[l, j];
                    }
                }
            }

            Console.WriteLine("\n\n\n computer output:");
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Console.Write(C1[i, j] + ",");
                    if (C2[i, j] != C1[i, j])
                    {
                        Console.WriteLine("Errrrrrrrrrrrorrrrrrrr");
                    }
                }
                Console.WriteLine();
            }

            Console.ReadKey();

        }
    }
}

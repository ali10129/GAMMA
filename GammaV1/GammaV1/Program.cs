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
            uint R, bufferSize;
            Element[] Amk_ScalingFactorRegFile;
            Element[,] Bkn;
            Element[,] Cmn;


            public PE(uint M , uint N, uint R = 64, uint bufferSize = 64)
            {
                this.R = R;
                this.bufferSize = bufferSize;

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

                Cmn = new Element[M, N];
                for (int i = 0; i < M; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        Cmn[i, j].clear();
                    }
                }
            }

            public uint loadA(Element[] Am, int start = 0)
            {
                uint count = 0;
                for (int i = 0; i < R; i++)
                {
                    if(Am[i+start].column %= )
                    Amk_ScalingFactorRegFile[i]  = 
                }
                return count;
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
        }
    }
}

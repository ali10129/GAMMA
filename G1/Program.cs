﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace G1
{
    class Program
    {
         class Cache
        {
            uint data_cache_size;
            uint block_size;
            uint n_block_in_assoc;

            uint number_of_cache_access = 0;

            uint n_sets;

            ulong[,,] Tags;


            uint write_hit, write_miss, read_hit, read_miss;
            public uint Cache_hitmiss_penalty(uint read_hitCycle=1, uint write_hitCycle=1,uint read_missCycle=10, uint write_missCycle=10)
            {
                Console.WriteLine("\n\n_____________________________________________________________________________");
                Console.WriteLine("Cache with LRU replacement policy");
                Console.WriteLine("Write:\t\tHit:{0,-15}\t\t-->Write Hit Cycle:{1}", write_hit, write_hitCycle);
                Console.WriteLine("Write:\t\tMiss:{0,-15}\t\t-->Write Miss Cycle:{1}", write_miss, write_missCycle);
                Console.WriteLine("Read:\t\tHit:{0,-15}\t\t-->Read Hit Cycle:{1}", read_hit, read_hitCycle);
                Console.WriteLine("Read:\t\tMiss:{0,-15}\t\t-->Read Miss Cycle:{1}", read_miss, read_missCycle);
                uint total = read_hit* read_hitCycle +read_miss * read_missCycle + write_hit * write_hitCycle + write_miss * write_missCycle;

                Console.WriteLine("\n\n ==============\n\tTotal cycle:\t{0}\n\tNumber of Accesses to the Cache:\t{1}", total, number_of_cache_access);
                Console.WriteLine("\tAverage Access Cycle:\t(Total cycle / Number of Accesses)\t\t{0}", (total/(double)number_of_cache_access));
                Console.WriteLine("\n_____________________________________________________________________________\n\n");
                return total;

            }
            

            public Cache(uint total_cache_size_in_byte, uint block_size_in_byte, uint number_of_blocks_in_associtivity)
            {
                n_block_in_assoc = number_of_blocks_in_associtivity;
                data_cache_size = total_cache_size_in_byte;
                block_size = block_size_in_byte;
                write_hit = write_miss = read_hit = read_miss = 0;

                n_sets = total_cache_size_in_byte / (block_size_in_byte * number_of_blocks_in_associtivity);

                Tags = new ulong[n_sets,n_block_in_assoc,2];
                for (uint i = 0; i < n_sets; i++)
                {
                    for (uint j = 0; j < n_block_in_assoc; j++)
                    {
                        Tags[i, j, 0] = 0;
                        Tags[i, j, 1] = j;
                    }
                }
            }

            public bool check_and_put_Data(ulong address, bool write = false)
            {
                // return True if hit, return False if miss;
                number_of_cache_access++;

                ulong index_bits =(address/block_size) % n_sets;
                ulong tag = (address / block_size) / n_sets;

                bool hit =  putData(tag, index_bits);

                if (write && hit)
                {
                    write_hit++;
                }
                else if (write && !hit)
                {
                    write_miss++;
                }
                else if (!write && hit)
                {
                    read_hit++;
                }
                else if (!write && !hit)
                {
                    read_miss++;
                }

                return hit;
            }

            private bool putData(ulong tag, ulong index)
            {
                bool hit = false;
                ulong block_id = 0;
                for (uint i = 0; i < n_block_in_assoc; i++)
                {
                    if (Tags[index,i,0] == tag)
                    {
                        hit = true;
                        block_id = Tags[index, i, 1];
                        break;
                    }

                    block_id = block_id < Tags[index, i, 1] ? Tags[index, i, 1] : block_id;
                }

                //LRU replacement Policy
                for (int i = 0; i < n_block_in_assoc; i++)
                {
                    if (Tags[index, i, 1] < block_id)
                    {
                        Tags[index, i, 1]++;
                    }
                    else if (Tags[index, i, 1] == block_id)
                    {
                        Tags[index, i, 0] = tag;
                        Tags[index, i, 1] = 0;
                    }
                }


                return hit;
            }

        }
        public static ulong getElementAddress(float[,] t, int i , int j)
        {
            ulong k0= 0;
            unsafe
            {
                // Must pin object on heap so that it doesn't move while using interior pointers.
                fixed (float* p1 = &t[i, j])
                {
                    k0 = (ulong)p1;
                } 
            }
            return k0;
        }

        static void Main(string[] args)
        {

            // total size 8 MB , 256 B block size, 32 Associative in each set, LRU replacement Policy:
            Cache DataCache0 = new Cache(8 * 1024 * 1024, 256, 32);
            Cache DataCache1 = new Cache(8 * 1024 * 1024, 256, 32);
            Cache DataCache2 = new Cache(8 * 1024 * 1024, 256, 32);

            int M = 100;
            int K = 100;
            int N = 100;

            float[,] A = new float[M, K];
            float[,] B = new float[K, N];

            float[,] C = new float[M, N];

            ulong k0,k1,k2;

            unsafe
            {
                // Must pin object on heap so that it doesn't move while using interior pointers.
                fixed (float* p1 = &A[0, 0])
                {
                    k0 = (ulong)p1;
                }
                fixed (float* p2 = &B[0, 0])
                {
                    k1 = (ulong)p2;
                }
                fixed (float* p3 = &C[0, 0])
                {
                    k2 = (ulong)p3;
                }
            }


            int aa = A.Length * sizeof(float);
            int bb = B.Length * sizeof(float);
            int cc = C.Length * sizeof(float);

            Console.WriteLine("size of A:{0} in (0x{1:x16})", aa, k0);
            Console.WriteLine("size of B:{0} in (0x{1:x16})", bb, k1);
            Console.WriteLine("size of C:{0} in (0x{1:x16})", cc, k2);


            //inner product
            for (int m = 0; m < M; m++)
            {
                for (int n = 0; n < N; n++)
                {

                    DataCache0.check_and_put_Data(getElementAddress(C, m, n));       //read C[m,n]
                    float REGISTER = C[m, n];

                    for (int k = 0; k < K; k++)
                    {
                        DataCache0.check_and_put_Data(getElementAddress(A, m, k));       //read A[m,k]
                        DataCache0.check_and_put_Data(getElementAddress(B, k, n));       //read B[k,n]

                        REGISTER += A[m, k] * B[k, n];

                    }
                    
                    DataCache0.check_and_put_Data(getElementAddress(C, m, n),true);       //write C[m,n]
                    C[m, n] = REGISTER;

                    Console.WriteLine("Inner-product dataflow:\t m {0}\t n {1}", m, n);
                }
            }

            //outer product
            for (int k = 0; k < K; k++)
            {
                for (int m = 0; m < M; m++)
                {
                    DataCache1.check_and_put_Data(getElementAddress(A, m, k));       //read A[m,k]
                    float REGISTER = A[m, k];
                    for (int n = 0; n < N; n++)
                    {

                        DataCache1.check_and_put_Data(getElementAddress(C, m, n));       //read C[m,n]      
                        DataCache1.check_and_put_Data(getElementAddress(B, k, n));       //read B[k,n]

                        C[m, n] += REGISTER * B[k, n];

                        DataCache1.check_and_put_Data(getElementAddress(C, m, n), true);       //write C[m,n]

                    }
                    Console.WriteLine("Outer-productdataflow:\t k {0}\t m {1}", k, m);
                }
            }

            //Gustavson product
            for (int m = 0; m < M; m++)
            {
                for (int k = 0; k < K; k++)
                {
                    DataCache2.check_and_put_Data(getElementAddress(A, m, k));       //read A[m,k]
                    float REGISTER = A[m, k];
                    for (int n = 0; n < N; n++)
                    {
                        DataCache2.check_and_put_Data(getElementAddress(C, m, n));       //read C[m,n]
                        DataCache2.check_and_put_Data(getElementAddress(B, k, n));       //read B[k,n]

                        C[m, n] += REGISTER * B[k, n];

                        DataCache2.check_and_put_Data(getElementAddress(C, m, n), true);       //write C[m,n]

                    }
                    Console.WriteLine("Gustavson dataflow:\t m {0}\t k {1}", m, k);

                }
            }


            Console.WriteLine("\n\nM * K * N  = {0}\n\n", M*K*N);
            Console.WriteLine("=================inter product:=====================");
            DataCache0.Cache_hitmiss_penalty();
            Console.WriteLine("=================outer product:=====================");
            DataCache1.Cache_hitmiss_penalty();
            Console.WriteLine("=================Gustavson product:=====================");
            DataCache2.Cache_hitmiss_penalty();

            Console.ReadKey();
        }
    }
}
using System;
using System.Threading;
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
            uint n_ways;

            uint number_of_cache_access = 0;

            uint n_sets;

            ulong[,,] Tags;


            uint write_hit, write_miss, read_hit, read_miss;
            public string Cache_hitmiss_penalty(uint read_hitCycle=1, uint write_hitCycle=1,uint read_missCycle=10, uint write_missCycle=10)
            {
                string txt ="";
                //txt += string.Format("\n\n_____________________________________________________________________________\n");
                //txt += string.Format("Cache with LRU replacement policy\n");
                //txt += string.Format("Read:\t\tHit:{0,-15}\t\t-->Read Hit Cycle:{1}\n", read_hit, read_hitCycle);
                txt += string.Format("{0},", read_miss);
                //txt += string.Format("Write:\t\tHit:{0,-15}\t\t-->Write Hit Cycle:{1}\n", write_hit, write_hitCycle);
                txt += string.Format("{0},", write_miss);
                //uint total = read_hit* read_hitCycle +read_miss * read_missCycle + write_hit * write_hitCycle + write_miss * write_missCycle;

                //txt += string.Format("\n\n ==============\n\tTotal cycle:\t{0}\n\tNumber of Accesses to the Cache:\t{1}\n", total, number_of_cache_access);
                //txt += string.Format("\tAverage Access Cycle:\t(Total cycle / Number of Accesses)\t\t{0}\n", (total/(double)number_of_cache_access));
                //txt += string.Format("\n_____________________________________________________________________________\n\n");

                return txt;
            }
            

            public Cache(uint total_cache_size_in_byte, uint block_size_in_byte, uint num_ways)
            {
                n_ways = num_ways;
                data_cache_size = total_cache_size_in_byte;
                block_size = block_size_in_byte;
                write_hit = write_miss = read_hit = read_miss = 0;

                n_sets = total_cache_size_in_byte / (block_size_in_byte * num_ways);

                Tags = new ulong[n_sets,n_ways,2];
                for (uint i = 0; i < n_sets; i++)
                {
                    for (uint j = 0; j < n_ways; j++)
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
                for (uint i = 0; i < n_ways; i++)
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
                for (int i = 0; i < n_ways; i++)
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
        public static ulong getElementAddress(float[,] t, int i , int j )
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
        public static ulong getElementAddress(ulong a0, int i, int j, int column_in_row, int floatsize = sizeof(float))
        {
            ulong k0 = a0 + (ulong)((i * floatsize) + j);
            return k0;
        }
        static void delta(int dm, int dk, int dn, uint cmb, uint _bsize, uint _asso)
        {

            var Inner = new StringBuilder();
            var Outer = new StringBuilder();
            var Gustavson = new StringBuilder();
            Random rand = new Random();

            Inner.AppendFormat("{0},{1},{2},{3},{4},{5},", dm, dk, dn, cmb, _bsize, _asso);
            Outer.AppendFormat("{0},{1},{2},{3},{4},{5},", dm, dk, dn, cmb, _bsize, _asso);
            Gustavson.AppendFormat("{0},{1},{2},{3},{4},{5},", dm, dk, dn, cmb, _bsize, _asso);

            // total size 8 MB , 256 B block size, 32 Associative in each set, LRU replacement Policy:
            Cache DataCache0 = new Cache(cmb * 1024 * 1024, _bsize, _asso);
            Cache DataCache1 = new Cache(cmb * 1024 * 1024, _bsize, _asso);
            Cache DataCache2 = new Cache(cmb * 1024 * 1024, _bsize, _asso);

            int M = dm;
            int K = dk;
            int N = dn;

            ulong k0, k1, k2;
            
            
            unsafe
            {
                float[,] As = new float[4, 4];
                fixed (float* p1 = &As[0, 0])
                {
                    k0 = (ulong)p1;
                }
            }

            k1 = k0 + (ulong)(M * K * sizeof(float) + rand.Next(100, 10000) * sizeof(float));
            k2 = k1 + (ulong)(K * N * sizeof(float) + rand.Next(100, 10000) * sizeof(float));

            ulong A = k0;
            ulong B = k1;
            ulong C = k2;


            Thread t0 = new Thread(() => inner_product(DataCache0,M,K,N,A,B,C));
            Thread t1 = new Thread(() => outer_product(DataCache0, M, K, N, A, B, C));
            Thread t2 = new Thread(() => gustavson(DataCache0, M, K, N, A, B, C));

            t0.Start();
            t1.Start();
            t2.Start();

            t0.Join();
            t1.Join();
            t2.Join();

            Inner.AppendLine(DataCache0.Cache_hitmiss_penalty());
            Outer.AppendLine(DataCache1.Cache_hitmiss_penalty());
            Gustavson.AppendLine(DataCache2.Cache_hitmiss_penalty());
            System.IO.File.AppendAllText("Inner.csv", Inner.ToString());
            System.IO.File.AppendAllText("Outer.csv", Outer.ToString());
            System.IO.File.AppendAllText("Gustavson.csv", Gustavson.ToString());
        }

        static void inner_product(Cache DataCache0,int M,int K,int N,ulong A,ulong B,ulong C)
        {
            //inner product
            for (int m = 0; m < M; m++)
            {
                for (int n = 0; n < N; n++)
                {

                    DataCache0.check_and_put_Data(getElementAddress(C, m, n, N));       //read C[m,n]
                    //float REGISTER = C[m, n];

                    for (int k = 0; k < K; k++)
                    {
                        DataCache0.check_and_put_Data(getElementAddress(A, m, k, K));       //read A[m,k]
                        DataCache0.check_and_put_Data(getElementAddress(B, k, n, N));       //read B[k,n]

                        //REGISTER += A[m, k] * B[k, n];

                    }

                    DataCache0.check_and_put_Data(getElementAddress(C, m, n, N), true);       //write C[m,n]
                    //C[m, n] = REGISTER;

                    //Console.WriteLine("Inner-product dataflow:\t m {0}\t n {1}", m, n);
                }
            }

        }
        static void outer_product(Cache DataCache1, int M, int K, int N, ulong A, ulong B, ulong C)
        {
            //outer product
            for (int k = 0; k < K; k++)
            {
                for (int m = 0; m < M; m++)
                {
                    DataCache1.check_and_put_Data(getElementAddress(A, m, k, K));       //read A[m,k]
                    //float REGISTER = A[m, k];
                    for (int n = 0; n < N; n++)
                    {

                        DataCache1.check_and_put_Data(getElementAddress(C, m, n, N));       //read C[m,n]      
                        DataCache1.check_and_put_Data(getElementAddress(B, k, n, N));       //read B[k,n]

                        //C[m, n] += REGISTER * B[k, n];

                        DataCache1.check_and_put_Data(getElementAddress(C, m, n, N), true);       //write C[m,n]

                    }
                    //Console.WriteLine("Outer-productdataflow:\t k {0}\t m {1}", k, m);
                }
            }
        }
        static void gustavson(Cache DataCache2, int M, int K, int N, ulong A, ulong B, ulong C)
        {
            //Gustavson product
            for (int m = 0; m < M; m++)
            {
                for (int k = 0; k < K; k++)
                {
                    DataCache2.check_and_put_Data(getElementAddress(A, m, k, K));       //read A[m,k]
                    //float REGISTER = A[m, k];
                    for (int n = 0; n < N; n++)
                    {
                        DataCache2.check_and_put_Data(getElementAddress(C, m, n, N));       //read C[m,n]
                        DataCache2.check_and_put_Data(getElementAddress(B, k, n, N));       //read B[k,n]

                        //C[m, n] += REGISTER * B[k, n];

                        DataCache2.check_and_put_Data(getElementAddress(C, m, n, N), true);       //write C[m,n]

                    }
                    //Console.WriteLine("Gustavson dataflow:\t m {0}\t k {1}", m, k);

                }
            }
        }

        static void Main(string[] args)
        {
            System.IO.File.WriteAllText("Inner.csv", "M,K,N,cache size (MB),block size (Byte),# of ways,# of read miss,# write miss,\n");
            System.IO.File.AppendAllText("Outer.csv", "M,K,N,cache size (MB),block size (Byte),# of ways,# of read miss,# write miss,\n");
            System.IO.File.AppendAllText("Gustavson.csv", "M,K,N,cache size (MB),block size (Byte),# of ways,# of read miss,# write miss,\n");
            int[] matDims = { 256, 2048, 65536 };   
            uint[] cacheSizes = { 1 };                 //MB
            uint[] blockSizes = { 32 };          //Byte
            uint[] Ways = {8};           // sets = cachesize / (block size * Ways)
            string txt = "";
            int index = 0;
            foreach (var item1 in matDims)
                foreach (var item2 in matDims)
                    foreach (var item3 in matDims)
                    {
                        Console.Clear();

                        Console.WriteLine("{0},{1},{2}\n\t==>> {3} / {4}", item1, item2, item3,++index,matDims.Length* matDims.Length * matDims.Length);
                        foreach (var item4 in cacheSizes)
                            foreach (var item5 in blockSizes)
                                foreach (var item6 in Ways)
                                {
                                    delta(item1, item2, item3, item4, item5, item6);
                                }
                    }

            System.IO.File.AppendAllText("Inner.csv", "\n\n");
            System.IO.File.AppendAllText("Outer.csv", "\n\n");
            System.IO.File.AppendAllText("Gustavson.csv", "\n\n");
        }
    }
}

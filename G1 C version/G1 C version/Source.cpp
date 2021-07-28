#include <iostream>

#include <cstdlib>
#include <ctime>

#include <fstream>

#define Tags(i,j,k)  Tags[(i*n_ways*2)+(j*2)+k]

#include <sstream>

//#include <pthread.h>

using namespace std;

class Cache
{
public:
	int data_cache_size;
	int block_size;
	int n_ways;

	int number_of_cache_access = 0;

	int n_sets;

	long* Tags;


	int write_hit, write_miss, read_hit, read_miss;

	string Cache_hitmiss_penalty(int read_hitCycle = 1, int write_hitCycle = 1, int read_missCycle = 10, int write_missCycle = 10)
	{
		string txt = "";
		txt += to_string(read_miss) + ",";
		txt += to_string(write_miss) + ",";
		return txt;
	}


	Cache(int total_cache_size_in_byte, int block_size_in_byte, int num_ways)
	{
		n_ways = num_ways;
		data_cache_size = total_cache_size_in_byte;
		block_size = block_size_in_byte;
		write_hit = write_miss = read_hit = read_miss = 0;

		n_sets = total_cache_size_in_byte / (block_size_in_byte * num_ways);

		Tags = new long[n_sets*n_ways * 2];
		for (int i = 0; i < n_sets; i++)
		{
			for (int j = 0; j < n_ways; j++)
			{
				Tags(i, j, 0) = 0;
				Tags(i, j, 1) = j;
			}
		}
	}

	bool check_and_put_Data(long address, bool write = false)
	{
		// return True if hit, return False if miss;
		number_of_cache_access++;

		long index_bits = (address / block_size) % n_sets;
		long tag = (address / block_size) / n_sets;

		bool hit = putData(tag, index_bits);

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

private: 
	bool putData(long tag, long index)
	{
		bool hit = false;
		long block_id = 0;
		for (int i = 0; i < n_ways; i++)
		{
			if (Tags(index, i, 0) == tag)
			{
				hit = true;
				block_id = Tags(index, i, 1);
				break;
			}

			block_id = block_id < Tags(index, i, 1) ? Tags(index, i, 1) : block_id;
		}

		//LRU replacement Policy
		for (int i = 0; i < n_ways; i++)
		{
			if (Tags(index, i, 1) < block_id)
			{
				Tags(index, i, 1)++;
			}
			else if (Tags(index, i, 1) == block_id)
			{
				Tags(index, i, 0) = tag;
				Tags(index, i, 1) = 0;
			}
		}


		return hit;
	}

};


inline long getElementAddress(long a0, int i, int j, int column_in_row, int floatsize = sizeof(float))
{
	long k0 = a0 + (long)((i * floatsize*column_in_row) + j);
	return k0;
}
/*
struct ThreadArgs
{
public:
	void* cache;
	int M, K, N;
	int A, B, C;
};


//void* inner_product(Cache DataCache0, int M, int K, int N, long A, long B, long C)
void* inner_product(void* args)
{

	Cache* DataCache0;
	int M, K, N, A, B, C;

	struct ThreadArgs *Args;
	Args = (struct ThreadArgs *) args;

	DataCache0 = (Cache*)Args->cache;

	M = Args->M;
	K = Args->K;
	N = Args->N;

	A = Args->A;
	B = Args->B;
	C = Args->C;


	//inner product
	for (int m = 0; m < M; m++)
	{
		for (int n = 0; n < N; n++)
		{

			DataCache0->check_and_put_Data(getElementAddress(C, m, n, N));       //read C[m,n]
																				 //float REGISTER = C[m, n];

			for (int k = 0; k < K; k++)
			{
				DataCache0->check_and_put_Data(getElementAddress(A, m, k, K));       //read A[m,k]
				DataCache0->check_and_put_Data(getElementAddress(B, k, n, N));       //read B[k,n]

																					 //REGISTER += A[m, k] * B[k, n];

			}

			DataCache0->check_and_put_Data(getElementAddress(C, m, n, N), true);       //write C[m,n]
																					   //C[m, n] = REGISTER;

																					   //Console.WriteLine("Inner-product dataflow:\t m {0}\t n {1}", m, n);
		}
	}
	pthread_exit(NULL);
	return 0;
}
void* outer_product(void* args)
{

	Cache* DataCache1;
	int M, K, N, A, B, C;

	struct ThreadArgs *Args;
	Args = (struct ThreadArgs *) args;

	DataCache1 = (Cache*)Args->cache;

	M = Args->M;
	K = Args->K;
	N = Args->N;

	A = Args->A;
	B = Args->B;
	C = Args->C;
	//outer product
	for (int k = 0; k < K; k++)
	{
		for (int m = 0; m < M; m++)
		{
			DataCache1->check_and_put_Data(getElementAddress(A, m, k, K));       //read A[m,k]
																				 //float REGISTER = A[m, k];
			for (int n = 0; n < N; n++)
			{

				DataCache1->check_and_put_Data(getElementAddress(C, m, n, N));       //read C[m,n]      
				DataCache1->check_and_put_Data(getElementAddress(B, k, n, N));       //read B[k,n]

																					 //C[m, n] += REGISTER * B[k, n];

				DataCache1->check_and_put_Data(getElementAddress(C, m, n, N), true);       //write C[m,n]

			}
			//Console.WriteLine("Outer-productdataflow:\t k {0}\t m {1}", k, m);
		}
	}
	pthread_exit(NULL);
	return 0;
}

void *gustavson(void* args)
{

	Cache* DataCache2;
	int M, K, N, A, B, C;

	struct ThreadArgs *Args;
	Args = (struct ThreadArgs *) args;

	DataCache2 = (Cache*)Args->cache;

	M = Args->M;
	K = Args->K;
	N = Args->N;

	A = Args->A;
	B = Args->B;
	C = Args->C;
	//Gustavson product
	for (int m = 0; m < M; m++)
	{
		for (int k = 0; k < K; k++)
		{
			DataCache2->check_and_put_Data(getElementAddress(A, m, k, K));       //read A[m,k]
																				 //float REGISTER = A[m, k];
			for (int n = 0; n < N; n++)
			{
				DataCache2->check_and_put_Data(getElementAddress(C, m, n, N));       //read C[m,n]
				DataCache2->check_and_put_Data(getElementAddress(B, k, n, N));       //read B[k,n]

																					 //C[m, n] += REGISTER * B[k, n];

				DataCache2->check_and_put_Data(getElementAddress(C, m, n, N), true);       //write C[m,n]

			}
			//Console.WriteLine("Gustavson dataflow:\t m {0}\t k {1}", m, k);

		}
	}
	pthread_exit(NULL);
	return 0;
}
*/
//void* inner_product(Cache DataCache0, int M, int K, int N, long A, long B, long C)
void inner_product0(Cache* DataCache0, int M, int K, int N, long A, long B, long C)
{


	//inner product
	for (int m = 0; m < M; m++)
	{
		for (int n = 0; n < N; n++)
		{

			DataCache0->check_and_put_Data(getElementAddress(C, m, n, N));       //read C[m,n]
																				 //float REGISTER = C[m, n];

			for (int k = 0; k < K; k++)
			{
				DataCache0->check_and_put_Data(getElementAddress(A, m, k, K));       //read A[m,k]
				DataCache0->check_and_put_Data(getElementAddress(B, k, n, N));       //read B[k,n]

																					 //REGISTER += A[m, k] * B[k, n];

			}

			DataCache0->check_and_put_Data(getElementAddress(C, m, n, N), true);       //write C[m,n]
																					   //C[m, n] = REGISTER;

																					   //Console.WriteLine("Inner-product dataflow:\t m {0}\t n {1}", m, n);
		}
	}
}
void outer_product0(Cache* DataCache1, int M, int K, int N, long A, long B, long C)
{
	//outer product
	for (int k = 0; k < K; k++)
	{
		for (int m = 0; m < M; m++)
		{
			DataCache1->check_and_put_Data(getElementAddress(A, m, k, K));       //read A[m,k]
																				//float REGISTER = A[m, k];
			for (int n = 0; n < N; n++)
			{

				DataCache1->check_and_put_Data(getElementAddress(C, m, n, N));       //read C[m,n]      
				DataCache1->check_and_put_Data(getElementAddress(B, k, n, N));       //read B[k,n]

																					//C[m, n] += REGISTER * B[k, n];

				DataCache1->check_and_put_Data(getElementAddress(C, m, n, N), true);       //write C[m,n]

			}
			//Console.WriteLine("Outer-productdataflow:\t k {0}\t m {1}", k, m);
		}
	}
}
void gustavson0(Cache* DataCache2, int M, int K, int N, long A, long B, long C)
{
	//Gustavson product
	for (int m = 0; m < M; m++)
	{
		for (int k = 0; k < K; k++)
		{
			DataCache2->check_and_put_Data(getElementAddress(A, m, k, K));       //read A[m,k]
																				 //float REGISTER = A[m, k];
			for (int n = 0; n < N; n++)
			{
				DataCache2->check_and_put_Data(getElementAddress(C, m, n, N));       //read C[m,n]
				DataCache2->check_and_put_Data(getElementAddress(B, k, n, N));       //read B[k,n]

																					 //C[m, n] += REGISTER * B[k, n];

				DataCache2->check_and_put_Data(getElementAddress(C, m, n, N), true);       //write C[m,n]

			}
			//Console.WriteLine("Gustavson dataflow:\t m {0}\t k {1}", m, k);

		}
	}
}



static void delta(int dm, int dk, int dn, int cmb, int _bsize, int _asso)
{
	char Inner[100];
	snprintf(Inner, sizeof(Inner), "%d,%d,%d,%d,%d,%d,", dm, dk, dn, cmb, _bsize, _asso);
	char Outer[100];
	snprintf(Outer, sizeof(Outer), "%d,%d,%d,%d,%d,%d,", dm, dk, dn, cmb, _bsize, _asso);
	char Gustavson[100];
	snprintf(Gustavson, sizeof(Gustavson), "%d,%d,%d,%d,%d,%d,", dm, dk, dn, cmb, _bsize, _asso);


	std::srand(std::time(nullptr)); // use current time as seed for random generator

	// total size xxx Kb , 256 B block size, 32 Associative in each set, LRU replacement Policy:
	Cache DataCache0 = Cache(cmb * 1024, _bsize, _asso);
	Cache DataCache1 = Cache(cmb * 1024, _bsize, _asso);
	Cache DataCache2 = Cache(cmb * 1024, _bsize, _asso);

	int M = dm;
	int K = dk;
	int N = dn;

	long k0, k1, k2;


	float *As = new float[4, 4];

	k0 = (long)&As;


	k1 = k0 + (long)(M * K * sizeof(float) + (rand() % 16) * 256);
	k2 = k1 + (long)(K * N * sizeof(float) + (rand() % 16) * 256);

	long A = k0;
	long B = k1;
	long C = k2;

	inner_product0(&DataCache0, M, K, N, A, B, C);
	outer_product0(&DataCache1, M, K, N, A, B, C);
	gustavson0(&DataCache2, M, K, N, A, B, C);
	/*
	pthread_t t0, t1, t2;

	struct ThreadArgs Args0, Args1, Args2;

	Args0.cache = &DataCache0;

	Args0.A = A;
	Args0.B = B;
	Args0.C = C;

	Args0.M = M;
	Args0.K = K;
	Args0.N = N;

	Args1.cache = &DataCache1;

	Args1.A = A;
	Args1.B = B;
	Args1.C = C;

	Args1.M = M;
	Args1.K = K;
	Args1.N = N;

	Args2.cache = &DataCache2;

	Args2.A = A;
	Args2.B = B;
	Args2.C = C;

	Args2.M = M;
	Args2.K = K;
	Args2.N = N;
	pthread_create(&t0, NULL, inner_product, (void *)&Args0);
	pthread_create(&t1, NULL, outer_product, (void *)&Args1);
	pthread_create(&t2, NULL, gustavson, (void *)&Args2);

	void* status;
	pthread_join(t0, &status);
	pthread_join(t1, &status);
	pthread_join(t2, &status);
	*/

	string I = "", O = "", G = "";

	I = string(Inner);
	O = string(Outer);
	G = string(Gustavson);

	I += DataCache0.Cache_hitmiss_penalty() + "\n";
	O += DataCache1.Cache_hitmiss_penalty() + "\n";
	G += DataCache2.Cache_hitmiss_penalty() + "\n";


	ofstream ofs;
	ofs.open("Inner.csv", std::ofstream::out | std::ofstream::app);
	ofs << I.c_str();
	ofs.close();

	ofs.open("Outer.csv", std::ofstream::out | std::ofstream::app);
	ofs << O.c_str();
	ofs.close();

	ofs.open("Gustavson.csv", std::ofstream::out | std::ofstream::app);
	ofs << G.c_str();
	ofs.close();
}





int main(int argc, char *argv[])
{
	if (argc != 7) {
		cout << "Invalid input arguments error!\nYou should input 6 integer arguments as:\nM K N Cache_size(KB) Block_size(B) #_of_Ways\nExit -1" << endl;
		cin >> argc;
		return -1;
	}

	int Data[6];
	for (int i = 0; i < 6; i++) {
		std::stringstream ss(argv[i + 1]);
		if (ss >> Data[i])
			std::cout << "Arg" << i << " is: " << Data[i] << endl;
		else
			std::cout << "error";
	}



	ofstream ofs;
	ofs.open("Inner.csv", std::ofstream::out | std::ofstream::app);
	ofs << "\nM,K,N,cache size (KB),block size (Byte),# of ways,# of read miss,# write miss," << std::endl;
	ofs.close();

	ofs.open("Outer.csv", std::ofstream::out | std::ofstream::app);
	ofs << "\nM,K,N,cache size (KB),block size (Byte),# of ways,# of read miss,# write miss," << std::endl;
	ofs.close();

	ofs.open("Gustavson.csv", std::ofstream::out | std::ofstream::app);
	ofs << "\nM,K,N,cache size (KB),block size (Byte),# of ways,# of read miss,# write miss," << std::endl;
	ofs.close();

	delta(Data[0], Data[1], Data[2], Data[3], Data[4], Data[5]);

	ofs.open("Inner.csv", std::ofstream::out | std::ofstream::app);
	ofs << std::endl << std::endl;
	ofs.close();

	ofs.open("Outer.csv", std::ofstream::out | std::ofstream::app);
	ofs << std::endl << std::endl;
	ofs.close();

	ofs.open("Gustavson.csv", std::ofstream::out | std::ofstream::app);
	ofs << std::endl << std::endl;
	ofs.close();

	return 0;
}
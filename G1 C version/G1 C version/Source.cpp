#include <iostream>

#include <cstdlib>
#include <ctime>

#include <fstream>

#define Tags(i,j,k)  Tags[(i*n_ways*2)+(j*2)+k]

#include <sstream>

//#include <pthread.h>
#include<thread>		//To compile programs with std::thread support use:
						//g++ -std = c++11 -pthread

using namespace std;



typedef unsigned long ulong;
//typedef unsigned long long ulong;			// if using windows OS uncomment this.
typedef unsigned int uint;

class Cache
{
public:
	uint data_cache_size;
	uint block_size;
	uint n_ways;

	ulong number_of_cache_access = 0;

	uint n_sets;

	ulong* Tags;


	ulong write_hit, write_miss, read_hit, read_miss;

	string Cache_hitmiss_penalty(ulong read_hitCycle = 1, ulong write_hitCycle = 1, ulong read_missCycle = 10, ulong write_missCycle = 10)
	{
		string txt = "";
		txt += to_string(read_miss) + ",";
		txt += to_string(write_miss) + ",";
		return txt;
	}


	Cache(uint total_cache_size_in_byte, uint block_size_in_byte, uint num_ways)
	{
		n_ways = num_ways;
		data_cache_size = total_cache_size_in_byte;
		block_size = block_size_in_byte;
		write_hit = write_miss = read_hit = read_miss = 0;

		n_sets = total_cache_size_in_byte / (block_size_in_byte * num_ways);

		Tags = new ulong[n_sets*n_ways * 2];
		for (uint i = 0; i < n_sets; i++)
		{
			for (uint j = 0; j < n_ways; j++)
			{
				Tags(i, j, 0) = 0;
				Tags(i, j, 1) = j;
			}
		}
	}

	bool check_and_put_Data(ulong address, bool write = false)
	{
		// return True if hit, return False if miss;
		number_of_cache_access++;

		ulong index_bits = (address / block_size) % n_sets;
		ulong tag = (address / block_size) / n_sets;

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
	bool putData(ulong tag, ulong index)
	{
		bool hit = false;
		ulong block_id = 0;
		for (uint i = 0; i < n_ways; i++)
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
		for (uint i = 0; i < n_ways; i++)
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


inline ulong getElementAddress(ulong a0, int i, int j, int column_in_row, int floatsize = sizeof(float))
{
	ulong k0 = a0 + (ulong)((i * floatsize*column_in_row) + j);
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
void inner_product0(Cache* DataCache0, int M, int K, int N, ulong A, ulong B, ulong C)
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
void outer_product0(Cache* DataCache1, int M, int K, int N, ulong A, ulong B, ulong C)
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
void gustavson0(Cache* DataCache2, int M, int K, int N, ulong A, ulong B, ulong C)
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



static void delta(uint dm, uint dk, uint dn, uint cmb, uint _bsize, uint _asso)
{

	

	// total size xxx Kb , 256 B block size, 32 Associative in each set, LRU replacement Policy:
	Cache DataCache0 = Cache(cmb * 1024, _bsize, _asso);
	Cache DataCache1 = Cache(cmb * 1024, _bsize, _asso);
	Cache DataCache2 = Cache(cmb * 1024, _bsize, _asso);

	int M = dm;
	int K = dk;
	int N = dn;

	ulong k0, k1, k2;


	float *As = new float[4, 4];
	k0 = (ulong)&As;

	std::srand(std::time(nullptr)); // use current time as seed for random generator
	k1 = k0 + (ulong)(M * K * sizeof(float) + (rand() % 16) * 256);
	k2 = k1 + (ulong)(K * N * sizeof(float) + (rand() % 16) * 256);

	ulong A = k0;
	ulong B = k1;
	ulong C = k2;

	//inner_product0(&DataCache0, M, K, N, A, B, C);
	//outer_product0(&DataCache1, M, K, N, A, B, C);
	//gustavson0(&DataCache2, M, K, N, A, B, C);

	thread t0(inner_product0, &DataCache0, M, K, N, A, B, C);
	thread t1(outer_product0, &DataCache1, M, K, N, A, B, C);
	thread t2(gustavson0, &DataCache2, M, K, N, A, B, C);

	t0.join();
	t1.join();
	t2.join();

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



	////////////Save results:///////////////////////////////////////////////////////////////////////////////////////////

	char Inner[100];
	snprintf(Inner, sizeof(Inner), "%d,%d,%d,%d,%d,%d,", dm, dk, dn, cmb, _bsize, _asso);
	char Outer[100];
	snprintf(Outer, sizeof(Outer), "%d,%d,%d,%d,%d,%d,", dm, dk, dn, cmb, _bsize, _asso);
	char Gustavson[100];
	snprintf(Gustavson, sizeof(Gustavson), "%d,%d,%d,%d,%d,%d,", dm, dk, dn, cmb, _bsize, _asso);

	string I = "", O = "", G = "";

	I = string(Inner);
	O = string(Outer);
	G = string(Gustavson);

	I += DataCache0.Cache_hitmiss_penalty() + "\n";
	O += DataCache1.Cache_hitmiss_penalty() + "\n";
	G += DataCache2.Cache_hitmiss_penalty() + "\n";

	cout << "inner:\t" << I << endl;
	cout << "outer:\t" << O << endl;
	cout << "gustav:\t" << G << endl;

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
		//cin >> argc;
		return -1;
	}

	uint Data[6];
	for (int i = 0; i < 6; i++) {
		std::stringstream ss(argv[i + 1]);
		if (ss >> Data[i])
			std::cout << "Arg" << i << " is: " << Data[i] << endl;
		else
			std::cout << "error";
	}

	//cout << "size of int: " << sizeof(uint) << "\tsize of long: " << sizeof(ulong);
	//ofstream ofs;
	/*
	ofs.open("Inner.csv", std::ofstream::out | std::ofstream::app);
	ofs << "\nM,K,N,cache size (KB),block size (Byte),# of ways,# of read miss,# write miss," << std::endl;
	ofs.close();

	ofs.open("Outer.csv", std::ofstream::out | std::ofstream::app);
	ofs << "\nM,K,N,cache size (KB),block size (Byte),# of ways,# of read miss,# write miss," << std::endl;
	ofs.close();

	ofs.open("Gustavson.csv", std::ofstream::out | std::ofstream::app);
	ofs << "\nM,K,N,cache size (KB),block size (Byte),# of ways,# of read miss,# write miss," << std::endl;
	ofs.close();
	*/

	delta(Data[0], Data[1], Data[2], Data[3], Data[4], Data[5]);

	/*
	ofs.open("Inner.csv", std::ofstream::out | std::ofstream::app);
	ofs << std::endl << std::endl;
	ofs.close();

	ofs.open("Outer.csv", std::ofstream::out | std::ofstream::app);
	ofs << std::endl << std::endl;
	ofs.close();

	ofs.open("Gustavson.csv", std::ofstream::out | std::ofstream::app);
	ofs << std::endl << std::endl;
	ofs.close();
	*/

	cout << "finished successfully!" << endl << endl;
	return 0;
}


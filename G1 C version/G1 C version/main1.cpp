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
//typedef unsigned long long int ulong;			// if using windows OS uncomment this.
typedef unsigned int uint;

class Cache
{
private:

	uint data_cache_size;
	uint block_size;
	uint n_ways;

	ulong number_of_cache_access = 0;

	uint n_sets;

	ulong* Tags;
	ulong* n_hammer;


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
	string Cache_hitmiss_penalty()
	{
		string txt = "";
		txt += to_string(read_miss + write_miss) + ",";
		txt += to_string(number_of_cache_access) + ",";
		return txt;
	}


public:
	ulong write_hit, write_miss, read_hit, read_miss;
	ulong RWMISS = 0;



	Cache(uint setSize, uint block_size_in_byte, uint num_ways)
	{
		n_ways = num_ways;
		data_cache_size = setSize * num_ways * block_size_in_byte;
		block_size = block_size_in_byte;
		write_hit = write_miss = read_hit = read_miss = 0;

		n_sets = setSize; 

		n_hammer = new ulong[n_sets];

		Tags = new ulong[n_sets*n_ways * 2];
		for (uint i = 0; i < n_sets; i++)
		{
			for (uint j = 0; j < n_ways; j++)
			{
				Tags(i, j, 0) = 0;
				Tags(i, j, 1) = j;
			}
			n_hammer[i] = 0;
		}
	}

	uint check_and_put_Data(ulong address, bool write = false)
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

		if (!hit) {
			n_hammer[index_bits]++;
		}

		if (hit)
			return 0;
		else
			return 1;
	}


	static ulong rwmiss(Cache* a, Cache* b, Cache* c) {
		a->RWMISS = a->read_miss + a->write_miss;
		b->RWMISS = b->read_miss + b->write_miss;
		c->RWMISS = c->read_miss + c->write_miss;
		return a->RWMISS + b->RWMISS + c->RWMISS;
	}
	static void info(Cache* A0, Cache* A1, Cache* A2, Cache* B0, Cache* B1, Cache* B2, Cache* C0, Cache* C1, Cache* C2, uint dm, uint dk, uint dn, uint _bsize, uint _asso, string filename = "_Results_SpiltCache.csv")
	{
		char buffer[200];
		snprintf(buffer, sizeof(buffer), "%d,%d,%d,%d,%d,A0:%d,B0:%d,C0:%d,A1:%d,B1:%d,C1:%d,A2:%d,B2:%d,C2:%d, ,", dm, dk, dn, _bsize, _asso, A0->n_sets,B0->n_sets,C0->n_sets, A1->n_sets, B1->n_sets, C1->n_sets, A2->n_sets, B2->n_sets, C2->n_sets );

		string I = "";

		I = string(buffer);


		I += "[," + A0->Cache_hitmiss_penalty() + "," + B0->Cache_hitmiss_penalty() + "," + C0->Cache_hitmiss_penalty() + ", ](," + to_string(rwmiss(A0, B0, C0)) + ",)|[,";
		I += A1->Cache_hitmiss_penalty() + "," + B1->Cache_hitmiss_penalty() + "," + C1->Cache_hitmiss_penalty() + ", ](," + to_string(rwmiss(A1, B1, C1)) + ",)|[,";
		I += A2->Cache_hitmiss_penalty() + "," + B2->Cache_hitmiss_penalty() + "," + C2->Cache_hitmiss_penalty() + ", ](," + to_string(rwmiss(A2, B2, C2)) + ",),\n";

		cout << filename << ":\t" << I << endl;

		ofstream ofs;
		ofs.open(filename, std::ofstream::out | std::ofstream::app);
		ofs << I.c_str();
		ofs.close();
	}
};

inline ulong getElementAddress(ulong a0, int i, int j, int column_in_row, int floatsize = sizeof(float))
{
	ulong k0 = a0 + (ulong)((i*column_in_row) + j)*floatsize;
	return k0;
}

void inner_product0(Cache* A0, Cache* B0, Cache* C0, int M, int K, int N, ulong A, ulong B, ulong C)
{

	//inner product
	for (int m = 0; m < M; m++)
	{
		for (int n = 0; n < N; n++)
		{

			C0->check_and_put_Data(getElementAddress(C, m, n, N));       //read C[m,n]
																		 //float REGISTER = C[m, n];

			for (int k = 0; k < K; k++)
			{
				A0->check_and_put_Data(getElementAddress(A, m, k, K));       //read A[m,k]
				B0->check_and_put_Data(getElementAddress(B, n, k, K));       //read B[k,n] =====> B.T[n,k]

																			 //REGISTER += A[m, k] * B[k, n]; ===> REGISTER += A[m, k] * B.T[n,k];

			}

			C0->check_and_put_Data(getElementAddress(C, m, n, N), true);       //write C[m,n]
																			   //C[m, n] = REGISTER;

		}
	}
}
void outer_product0(Cache* A1, Cache* B1, Cache* C1, int M, int K, int N, ulong A, ulong B, ulong C)
{
	//outer product
	for (int k = 0; k < K; k++)
	{
		for (int m = 0; m < M; m++)
		{
			A1->check_and_put_Data(getElementAddress(A, k, m, M));       //read A[m,k] =====> A.T[k,m]
																		 //float REGISTER = A.T[k,m];
			for (int n = 0; n < N; n++)
			{

				C1->check_and_put_Data(getElementAddress(C, m, n, N));       //read C[m,n]      
				B1->check_and_put_Data(getElementAddress(B, k, n, N));       //read B[k,n]

																			 //C[m, n] += REGISTER * B[k, n];

				C1->check_and_put_Data(getElementAddress(C, m, n, N), true);       //write C[m,n]

			}
			//Console.WriteLine("Outer-productdataflow:\t k {0}\t m {1}", k, m);
		}
	}
}
void gustavson0(Cache* A2, Cache* B2, Cache* C2, int M, int K, int N, ulong A, ulong B, ulong C)
{
	//Gustavson product
	for (int m = 0; m < M; m++)
	{
		for (int k = 0; k < K; k++)
		{
			A2->check_and_put_Data(getElementAddress(A, m, k, K));       //read A[m,k]
																		 //float REGISTER = A[m, k];
			for (int n = 0; n < N; n++)
			{
				C2->check_and_put_Data(getElementAddress(C, m, n, N));       //read C[m,n]
				B2->check_and_put_Data(getElementAddress(B, k, n, N));       //read B[k,n]

																			 //C[m, n] += REGISTER * B[k, n];

				C2->check_and_put_Data(getElementAddress(C, m, n, N), true);       //write C[m,n]

			}
			//Console.WriteLine("Gustavson dataflow:\t m {0}\t k {1}", m, k);

		}
	}
}



static void delta(uint dm, uint dk, uint dn, uint _bsize, uint _asso)
{
	//std::srand(std::time(nullptr)); // use current time as seed for random generator

	// total size xxx Kb , xxx B block size, xxx Associative in each set, LRU replacement Policy:


	int M = dm;
	int K = dk;
	int N = dn;

	//float *As = new float[4, 4];
	//k0 = (ulong)&As;
	ulong k0, k1, k2;
	k0 = 0xF0000000;
	k1 = k0 + (ulong)(M * K * sizeof(float));		//k1 = k0 + (ulong)(M * K * sizeof(float) + (rand() % 16) * 256);
	k2 = k1 + (ulong)(K * N * sizeof(float));		//k2 = k1 + (ulong)(K * N * sizeof(float) + (rand() % 16) * 256);

	ulong A = k0;
	ulong B = k1;
	ulong C = k2;


	int sizeA = 50;
	int sizeB = 200;
	int sizeC = 50;

	Cache A0 = Cache(sizeA, _bsize, _asso);
	Cache B0 = Cache(sizeB, _bsize, _asso);
	Cache C0 = Cache(sizeC, _bsize, _asso);
	///////////////////////////////////////

	sizeA = 50;
	sizeB = 50;
	sizeC = 200;

	Cache A1 = Cache(sizeA, _bsize, _asso);
	Cache B1 = Cache(sizeB, _bsize, _asso);
	Cache C1 = Cache(sizeC, _bsize, _asso);
	///////////////////////////////////////

	sizeA = 50;
	sizeB = 200;
	sizeC = 50;

	Cache A2 = Cache(sizeA, _bsize, _asso);
	Cache B2 = Cache(sizeB, _bsize, _asso);
	Cache C2 = Cache(sizeC, _bsize, _asso);

	thread t0(inner_product0, &A0, &B0, &C0, M, K, N, A, B, C);
	thread t1(outer_product0, &A1, &B1, &C1, M, K, N, A, B, C);
	thread t2(gustavson0, &A2, &B2, &C2, M, K, N, A, B, C);

	t0.join();
	t1.join();
	t2.join();

	////////////Save results:///////////////////////////////////////////////////////////////////////////////////////////
	Cache::info(&A0, &A1, &A2, &B0, &B1, &B2, &C0, &C1, &C2, dm, dk, dn, _bsize, _asso);
}

int main(int argc, char *argv[])
{
	if (argc != 6) {
		cout << "Invalid input arguments error!\nYou should input 5 integer arguments as:\nM K N  Block_size(B) #_of_Ways\nExit -1" << endl;
		//cin >> argc;
		return -1;
	}

	uint Data[5];
	for (int i = 0; i < 5; i++) {
		std::stringstream ss(argv[i + 1]);
		if (ss >> Data[i])
			std::cout << "Arg" << i << " is: " << Data[i] << endl;
		else
			std::cout << "error";
	}

	delta(Data[0], Data[1], Data[2], Data[3], Data[4]);

	cout << "finished successfully!" << endl << endl;
	cin >> Data[0];
	return 0;
}


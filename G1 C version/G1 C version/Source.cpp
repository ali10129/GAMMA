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

	ulong write_hit, write_miss, read_hit, read_miss;
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
	string Cache_hitmiss_penalty(ulong read_hitCycle = 1, ulong write_hitCycle = 1, ulong read_missCycle = 10, ulong write_missCycle = 10)
	{
		string txt = "";
		txt += to_string(read_miss) + ",";
		txt += to_string(write_miss) + ",";
		return txt;
	}


public:

	Cache(uint total_cache_size_in_byte, uint block_size_in_byte, uint num_ways)
	{
		n_ways = num_ways;
		data_cache_size = total_cache_size_in_byte;
		block_size = block_size_in_byte;
		write_hit = write_miss = read_hit = read_miss = 0;

		n_sets = total_cache_size_in_byte / (block_size_in_byte * num_ways);

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

		if (!hit) {
			n_hammer[index_bits]++;
		}

		return hit;
	}

	static void info(Cache* c0, Cache* c1 , Cache* c2, uint dm, uint dk, uint dn, uint cmb, uint _bsize, uint _asso, string filename="Results.csv")
	{
		char buffer[100];
		snprintf(buffer, sizeof(buffer), "%d,%d,%d,%d,%d,%d,==>{IN|OUT|Gus}[R miss|W miss], ,", dm, dk, dn, cmb, _bsize, _asso);

		string I = "";

		I = string(buffer);


		I += c0->Cache_hitmiss_penalty() + " ," + c1->Cache_hitmiss_penalty() + " ," + c2->Cache_hitmiss_penalty() + "\n";

		cout << filename << ":\t"<< I << endl;

		ofstream ofs;
		ofs.open(filename, std::ofstream::out | std::ofstream::app);
		ofs << I.c_str();
		ofs.close();
	}

	static void Hammer_results(Cache* C0,Cache* C1,Cache* C2, uint dm, uint dk, uint dn, uint cmb, uint _bsize, uint _asso)
	{
		string txt = "\n" + to_string(dm) + "," + to_string(dk) + "," + to_string(dn) + "," + to_string(cmb) + "," + to_string(_bsize) + "," + to_string(_asso) + ",\n";
		txt += "# of set,Inner,Outer,Gustovson,\n";
		for (int i = 0; i < C0->n_sets; i++)
		{
			txt += to_string(i) + "," + to_string(C0->n_hammer[i]) + "," + to_string(C1->n_hammer[i]) + "," + to_string(C2->n_hammer[i]) + ",\n";
		}
		ofstream ofs;
		ofs.open("_Hammering Results.csv", std::ofstream::out | std::ofstream::app);
		ofs << txt.c_str();
		ofs.close();
	}

};


inline ulong getElementAddress(ulong a0, int i, int j, int column_in_row, int floatsize = sizeof(float))
{
	ulong k0 = a0 + (ulong)((i * floatsize*column_in_row) + j);
	return k0;
}

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
	//std::srand(std::time(nullptr)); // use current time as seed for random generator

	// total size xxx Kb , xxx B block size, xxx Associative in each set, LRU replacement Policy:
	Cache DataCache0 = Cache(cmb * 1024, _bsize, _asso);
	Cache DataCache1 = Cache(cmb * 1024, _bsize, _asso);
	Cache DataCache2 = Cache(cmb * 1024, _bsize, _asso);

	int M = dm;
	int K = dk;
	int N = dn;

	//float *As = new float[4, 4];
	//k0 = (ulong)&As;
	ulong k0, k1, k2;
	k0 = 0xF0000000;
	k1 = k0 + (ulong)(M * K * sizeof(float) );		//k1 = k0 + (ulong)(M * K * sizeof(float) + (rand() % 16) * 256);
	k2 = k1 + (ulong)(K * N * sizeof(float) );		//k2 = k1 + (ulong)(K * N * sizeof(float) + (rand() % 16) * 256);

	ulong A = k0;
	ulong B = k1;
	ulong C = k2;

	thread t0(inner_product0, &DataCache0, M, K, N, A, B, C);
	thread t1(outer_product0, &DataCache1, M, K, N, A, B, C);
	thread t2(gustavson0, &DataCache2, M, K, N, A, B, C);

	t0.join();
	t1.join();
	t2.join();

	////////////Save results:///////////////////////////////////////////////////////////////////////////////////////////
	Cache::info(&DataCache0, &DataCache1, &DataCache2, dm, dk, dn, cmb, _bsize, _asso);
	Cache::Hammer_results(&DataCache0, &DataCache1, &DataCache2, dm, dk, dn, cmb, _bsize, _asso);
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

	delta(Data[0], Data[1], Data[2], Data[3], Data[4], Data[5]);

	cout << "finished successfully!" << endl << endl;
	return 0;
}


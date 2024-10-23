////////////////////////////////////////////////////////
// Implementation of a Suduku solver using SPA Sep 2024
////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <time.h>
#include <limits.h>
//#include <complex.h>
#include <vector>
#include <iostream>
#include <omp.h>
#include <iomanip>
#include <fstream>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <Windows.h>
#include "../../../../Downloads/color.hpp"
#include <unordered_map>
#include <utility> // for std::pair

using namespace std;
 
#define Def_N				9	// Sudoku size of N x N
#define Def_VNs				Def_N * Def_N // Number of VNs of LDPC code equavilent
#define Def_CNs				3 * Def_N // Number of CNs of LDPC code equavilent
#define Def_dv				int(sqrt(Def_N)) // Degree of VNs
#define Def_dc				Def_N // Degree of CNs

#define Max_Iterations		40 // Max numer of Iterations for SPA decoder

vector<vector<int>> VNsNeighbors(Def_VNs, vector<int>(Def_dv, -1));
vector<vector<int>> CNsNeighbors(Def_CNs, vector<int>(Def_dc, -1));

vector<vector<double>> initial_probs(Def_VNs, vector<double>(Def_N, 0));

vector<vector<vector<double>>> mVNtoCN(Def_VNs, vector<vector<double>>(Def_dv, vector<double>(Def_N, 0)));
vector<vector<vector<double>>> mCNtoVN(Def_CNs, vector<vector<double>>(Def_dc, vector<double>(Def_N, 0)));

vector<vector<vector<int>>> Permutations(Def_N, vector<vector<int>>(40320, vector<int>(Def_N - 1, -1)));

vector<int> Sudoku(Def_VNs, 0);

static const char row_separator[] =
"+---+---+---+---+---+---+---+---+---+\n";

void print_Sudoku(int* Sud, bool type)
{
	int i, j, k;

	int N = Def_N;

	if (type == 0)
	{
		for (unsigned int board_row = 0; board_row < N; ++board_row)
		{
			std::cout << row_separator;
			cout << "| ";
			for (unsigned int board_column = 0;
				board_column < N;
				++board_column)
			{
				if (Sudoku[board_row * N + board_column] != 0)
					std::cout << Sud[board_row * N + board_column] << " | ";
				else
					std::cout << " " << " | ";
			}
			cout << "\n";
		}
	}
	else
	{
		for (unsigned int board_row = 0; board_row < N; ++board_row)
		{
			std::cout << row_separator;
			cout << "| ";
			for (unsigned int board_column = 0;
				board_column < N;
				++board_column)
			{
				if (Sudoku[board_row * N + board_column] != 0)
					std::cout << Sud[board_row * N + board_column] << " | ";
				else
					std::cout << dye::blue_on_black(Sud[board_row * N + board_column]) << " | ";
			}
			cout << "\n";
		}
	}
	std::cout << row_separator;
}

void connect_graph(int N, int VNs, int CNs, int dv, int dc)
{
	int i, j, k;

	int const1 = N;
	int const2 = 2 * N;
	int const3 = 3 * N;
	
	vector<int> VNs_count(VNs,0);

	for (i = 0; i < const1; i++)
	{
		for (j = 0; j < dc; j++)
		{
			int numb = i * N + j;

			CNsNeighbors[i][j] = numb;
			VNsNeighbors[numb][VNs_count[numb]++] = i;
		}
	}

	for (i = 0; i < const1; i++)
	{
		for (j = 0; j < dc; j++)
		{
			int numb = i + j * const1;

			CNsNeighbors[i + const1][j] = numb;
			VNsNeighbors[numb][VNs_count[numb]++] = i + const1;
		}
	}

	for (i = 0; i < const1; i++)
	{
		for (j = 0; j < dc; j++)
		{
			int numb = (i / 3 ) * const3 +  (i % 3 * 3) + (j % 3 + (const1 * (j / 3)));
			
			CNsNeighbors[i + const2][j] = numb;
			VNsNeighbors[numb][VNs_count[numb]++] = i + const2;
		}
	}
}


void permutations(vector<vector<int>>& res, vector<int> arr, int idx) {

	// Base case: if idx reaches the size of the array, 
	 // add the permutation to the result
	if (idx == arr.size()) {
		res.push_back(arr);
		return;
	}

	// Permutations made by swapping each element
	// starting from index `idx`
	for (int i = idx; i < arr.size(); i++) {

		// Swapping
		swap(arr[idx], arr[i]);

		// Recursive call to create permutations 
		// for the next element
		permutations(res, arr, idx + 1);

		// Backtracking
		swap(arr[idx], arr[i]);
	}
}

// Function to get the permutations
vector<vector<int>> permute(vector<int>& arr) {

	// Declaring result variable
	vector<vector<int>> res;

	// Calling permutations with idx starting at 0
	permutations(res, arr, 0);
	return res;
}


void compute_permutations(int N)
{
	int i, j, k, x;
	int size = N - 1;

	ifstream inFile;

	inFile.open("Permutations.txt");

	for (i = 0; i < N; i++)
	{
		inFile >> x;

		for (j = 0; j < 40320; j++)
		{
			for (k = 0; k < size; k++)
			{
				inFile >> x;
				Permutations[i][j][k] = x;
			}
		}
	}

	inFile.close();
}

void Sudoku_fill(int VNs)
{
	int i, j;

	//string initial_fill;

	int initial_fill[Def_VNs] = {
	
		0,0,2,  3,7,0,  0,5,0,
		1,0,0,  0,2,9,  0,0,0,
		0,0,4,  0,6,1,  0,0,0,

		0,0,0,  0,0,4,  0,0,6,
		6,0,0,  2,5,0,  0,1,8,
		0,2,7,  0,0,0,  0,0,5,

		4,0,0,  8,3,0,  5,0,1,
		0,0,0,  0,0,0,  9,0,0,
		7,5,0,  0,0,0,  0,4,2,
								};

	//cout << "Type down the filled numbers in your current Sudoku, going row by row (expty spaces are equal to zero):";

	//getline(cin, initial_fill);

	for (i = 0; i < VNs; i++)
		Sudoku[i] = initial_fill[i];


	printf("The Sudoku to be solved is as following:\n");

	print_Sudoku(initial_fill, 0);
}

void initialize_probabilites(int VNs, int N, int dv, int CNs, int dc)
{
	int i, j, k, CN, VN;

	for (i = 0; i < VNs; i++)
	{
		int numb = Sudoku[i];
		
		if(numb != 0)
			initial_probs[i][numb - 1] = 1.0;
		else
		{
			for (j = 0; j < N; j++)
				initial_probs[i][j] = 1.0;

			for (j = 0; j < dv; j++)
			{
				CN = VNsNeighbors[i][j];

				for (k = 0; k < dc; k++)
				{
					VN = CNsNeighbors[CN][k];

					if (VN != i)
					{
						int numb2 = Sudoku[VN];

						if (numb2 != 0)
							initial_probs[i][numb2 - 1] = 0.0;
						else
							continue;
					}
					else
						continue;
				}
			}
		}
		
	}

	double sum = 0;

	for (i = 0; i < VNs; i++)
	{
		sum = 0;

		for (j = 0; j < N; j++)
			sum += initial_probs[i][j];

		for (j = 0; j < N; j++)
			initial_probs[i][j] /= sum;
	}
}

//Mapp integer to vector over GF(q)
void int2vec_q(int* vect, int integ, int vect_length, int q)
{
	int i, j;
	int cnt = 0;

	for (i = 0; i < vect_length; i++)
	{
		vect[i] = 0;
	}

	while (integ > 0)
	{
		vect[cnt++] = integ % q;
		integ /= q;
	}

}


void dec2binarr(int * res, long n, int dim)
{
	// note: res[dim] will save the sum res[0]+...+res[dim-1]
	//int* res = (int*)calloc(dim + 1, sizeof(int));
	//vector<int> res(dim + 1, 0);
	int pos = dim - 1;

	for (int i = 0; i < dim; i++)
		res[i] = 0;

	// note: this will crash if dim < log_2(n)...
	while (n > 0)
	{
		res[pos] = n % 2;
		res[dim] += res[pos];
		n = n / 2; // integer division        
		pos--;
	}
}

struct pair_hash {
	template <class T1, class T2>
	std::size_t operator () (const std::pair<T1, T2>& pair) const {
		return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
	}
};

// Function to compute the permanent recursively with memoization (using double)
double compute_permanent(const vector<vector<double>>& matrix, int row, int bitmask, unordered_map<pair<int, int>, double, pair_hash>& memo) {
	int n = matrix.size();

	// Base case: If all rows are assigned (i.e., all columns are used)
	if (row == n) {
		return 1.0;
	}

	// Generate the key as a pair of row and bitmask
	pair<int, int> key = { row, bitmask };

	// Check if the result is already computed
	if (memo.find(key) != memo.end()) {
		return memo[key];
	}

	// Recursive case: Iterate over all columns
	double permanent = 0.0;
	for (int col = 0; col < n; ++col) {
		if (!(bitmask & (1 << col))) {  // If this column is not yet used
			permanent += matrix[row][col] * compute_permanent(matrix, row + 1, bitmask | (1 << col), memo);
		}
	}

	// Memoize and return the result
	memo[key] = permanent;
	return permanent;
}

// Wrapper function to start the computation
double permanent(const vector<vector<double>>& matrix) {
	int n = matrix.size();
	unordered_map<pair<int, int>, double, pair_hash> memo;  // Memoization table for doubles
	return compute_permanent(matrix, 0, 0, memo);
}


void SPA_sudoku(int N, int VNs, int CNs, int dv, int dc, int iter_max, int* solution)
{
	int i, j, k, ii, jj, iter, CN, VN;

	double prod[Def_VNs][Def_N];

	double sum_norm[Def_VNs];
	double prod2[Def_N];
	vector<int> Neighbor(Def_dc, 0);

	for (VN = 0; VN < VNs; VN++)
		for (i = 0; i < dv; i++)
			for(j = 0; j < N; j++)
				mVNtoCN[VN][i][j] = initial_probs[VN][j];

	for (iter = 0; iter < iter_max; iter++)
	{
		//CN to VN message update
		for (CN = 0; CN < CNs; CN++)
		{

			for (i = 0; i < dc; i++)
			{
				VN = CNsNeighbors[CN][i];

				for (ii = 0; ii < dv; ii++)
				{
					if (VNsNeighbors[VN][ii] == CN)
					{
						Neighbor[i] = ii;
						break;
					}
				}
			}

			vector<vector<double>> A(dc - 1, vector<double>(N - 1, -1));

			for (i = 0; i < dc; i++)
			{
				/*
				vector<int> CNsNeighs(Def_dc, 0);
				vector<int> Neighbor_temp(Def_dc, 0);

				CNsNeighs = CNsNeighbors[CN];
				Neighbor_temp = Neighbor;

				VN = CNsNeighs[i];

				CNsNeighs.erase(CNsNeighs.begin() + i);
				Neighbor_temp.erase(Neighbor_temp.begin() + i);

				double sum_CN = 0;

				for (j = 0; j < N; j++)
				{
					double prod_CN = 1;

					for (ii = 0; ii < 40320; ii++)
					{
						prod_CN = 1;

						for (jj = 0; jj < N - 1 && prod_CN != 0; jj++)
							prod_CN *= mVNtoCN[CNsNeighs[jj]][Neighbor_temp[jj]][Permutations[j][ii][jj]];

						mCNtoVN[CN][i][j] += prod_CN;
					}

					sum_CN += mCNtoVN[CN][i][j];
				}
				*/

				
				vector<int> CNsNeighs(Def_dc, 0);
				CNsNeighs = CNsNeighbors[CN];

				VN = CNsNeighs[i];
				
				double sum_CN = 0;

				for (j = 0; j < N; j++)
				{
					int count1, count2;

					count1 = 0;

					for (ii = 0; ii < dc; ii++)
					{
						if (ii != i)
						{
							count2 = 0;

							for (jj = 0; jj < N; jj++)
							{
								if (jj != j)
									A[ii - count1][jj - count2] = mVNtoCN[CNsNeighs[ii]][Neighbor[ii]][jj];
								else
									count2++;
							}
						}
						else
							count1++;
					}

					mCNtoVN[CN][i][j] = permanent(A);
					sum_CN += mCNtoVN[CN][i][j];

				}
			
				for (j = 0; j < N; j++)
					mCNtoVN[CN][i][j] /= sum_CN;

			}
		}//loop over all CNs

		for (VN = 0; VN < VNs; VN++)
			for (i = 0; i < N; i++)
				prod[VN][i] = 1.0;

		//VN to CN message update
		for (VN = 0; VN < VNs; VN++)
		{
			sum_norm[VN] = 0;

			for (i = 0; i < N; i++)
			{
				for (j = 0; j < dv; j++)
				{
					CN = VNsNeighbors[VN][j];

					for (k = 0; k < dc; k++)
					{
						if (CNsNeighbors[CN][k] == VN)
						{
							Neighbor[j] = k;
						}
					}
					
					double numb = mCNtoVN[CN][Neighbor[j]][i];
					
					prod[VN][i] *= numb;
				}//end for loop over VN neighbors

				prod[VN][i] *= initial_probs[VN][i];
			}//end for loop over GF symbols 

			for (i = 0; i < N; i++)
				prod2[i] = 1;

			for (i = 0; i < dv; i++)
			{
				double f_sum = 0;
				CN = VNsNeighbors[VN][i];
				
				for (j = 0; j < N; j++)
				{
					prod2[j] = prod[VN][j];

					double numb = mCNtoVN[CN][Neighbor[i]][j];

					if (numb != 0.0)
						prod2[j] /= numb;
					else
						prod2[j] = 0;

					f_sum += prod2[j];
				}//end for loop over VN neighbors

				if (f_sum != 0)
				{
					for (j = 0; j < N; j++)
					{
						mVNtoCN[VN][i][j] = prod2[j] / f_sum;
					}
				}
			}
		}//end loop over VNs

		printf("Iter = %d\n", iter);
	}

	//APP descicion
	for (VN = 0; VN < VNs; VN++)
	{
		int max = rand() % N;
		double max_prod = 0;

		for (i = 0; i < N; i++)
		{
			if (prod[VN][i] > max_prod)
			{
				max_prod = prod[VN][i];
				max = i;
			}
			else
				continue;
		}//end for loop over GF symbols

		solution[VN] = max + 1;
	}//end for loop over VNs
}


void main()
{
	int N = Def_N;
	int dv = Def_dv;
	int dc = Def_dc;
	int VNs = Def_VNs;
	int CNs = Def_CNs;

	int iter_max = Max_Iterations;

	int solution[Def_VNs] = { 0 };

	connect_graph(N, VNs, CNs, dv, dc);

	//compute_permutations(N);

	Sudoku_fill(VNs);

	initialize_probabilites(VNs, N, dv, CNs, dc);

	SPA_sudoku(N, VNs, CNs, dv, dc, iter_max, solution);

	/*
	for (int i = 0; i < VNs; i++)
		Sudoku[i] = solution[i];
		*/
	printf("The Solution to the above Sudoku is most likely the following:\n");
	print_Sudoku(solution, 1);
}
#include<iostream>
#include<windows.h>　
#include <ctime>
#include <cstdlib>
#include <fstream>
using namespace std;
void ordinary(int n, int** b, int* a, int* sum) {
	for (int i = 0; i < n; i++) {
		sum[i] = 0.0;
		for (int j = 0; j < n; j++) {
			sum[i] += b[j][i] * a[j];
		}
	}
}
void optimal(int n, int** b, int* a, int* sum) {
	for (int i = 0; i < n; i++) {
		sum[i] = 0.0;
	}
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			sum[i] += b[j][i] * a[j];
		}
	}
}
int getRand(int min, int max) {
	return (rand() % (max - min + 1)) + min;
}
int main()
{

	srand(time(0));
	int n;
	cin >> n;
	/*ofstream file("rand.txt");
	for (int i = 0; i < n * (n + 1); i++) {
		int r = getRand(0, 9);
		file << r<<endl;//生成随机数来提供实验样本
	}*/
	
	int* a = new int[n];
	for (int i = 0; i < n; i++) {
		cin>>a[i];
	}
	int** b = new int* [n];
	int* sum = new int[n];
	for (int i = 0; i < n; i++) {
		b[i] = new int[n];
		sum[i] = -1;
		for(int j = 0; j < n; j++)
		cin >> b[i][j];
	}

	LARGE_INTEGER t1, t2, tc;
	QueryPerformanceFrequency(&tc);
	QueryPerformanceCounter(&t1);
	for (int i = 0; i < 1000000; i++) { ordinary(n, b, a, sum); }
		QueryPerformanceCounter(&t2);
	double time = (double)(t2.QuadPart - t1.QuadPart) / (double)tc.QuadPart;
	cout << "time = " << time << endl;
	return 0;
}
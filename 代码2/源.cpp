#include<iostream>
#include<windows.h>¡¡
#include <cstdlib>
using namespace std;
int ordinary(int n, int* a) {
	int sum=0;
	for (int i = 0; i < n; i++) {
		sum += a[i];
	}
	return sum;
}
int optimal1(int n, int* a) {
	int sum1 = 0; int sum2 = 1;
	for (int i = 0; i < n; i += 2) {
		sum1 += a[i];
		sum2 += a[i + 1];
	}
	return sum1 + sum2;
}
int optimal2(int n, int* a) {
	if (n == 1)	return 0;
	else {
		for (int i = 0; i < n / 2; i++) {
			a[i] += a[n - i - 1];
			n = n / 2;
			optimal2(n,a);
		}
	}
}
int main() {
	srand(time(0));
	int n;
	cin >> n;
	int* a = new int[n];
	for (int i = 0; i < n; i++) {
		cin>>a[i];
	}
	LARGE_INTEGER t1, t2, tc;
	QueryPerformanceFrequency(&tc);
	QueryPerformanceCounter(&t1);
	for (int i = 0; i < 1000000; i++) { 
		ordinary(n,a);
		//optimal1(n, a);
		//optimal2(n, a);
	}
	QueryPerformanceCounter(&t2);
	double time = (double)(t2.QuadPart - t1.QuadPart) / (double)tc.QuadPart;
	cout << "time = " << time << endl;
	return 0;
}
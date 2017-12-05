#include<iostream>
#include<complex>
#include<math.h>
#include<time.h>

#define PI 3.141592653

using namespace std;

void dft(complex<double>* xn, int n, complex<double>* Xk) {
	complex<double> WN(cos(-2 * PI / n), sin(-2 * PI / n));
	complex<double> Wk = 1.0;

	for (int i = 0; i < n; i++) {
		complex<double> Wnk = 1.0;
		for (int j = 0; j < n; j++)
		{
			Xk[i] += xn[j] * Wnk;
			Wnk *= Wk;
		}
		Wk *= WN;
	}
}

//倒位序
void invert(complex<double>* xn, int n, complex<double>* inverted_xn){
	int bit_num = int(round(log2(n)));
	int half_width = int(floor(bit_num / 2.0));

	//求倒位序号
	for (int i = 0; i < n; i++) {

		int t = 0;
		int bit_pos1 = 1;
		int bit_pos2 = n >> 1;
		for (int j = 0; j < half_width; j++) {
			t = ( t | ((i & bit_pos1) << (bit_num - j*2 - 1)) );
			t = ( t | ((i & bit_pos2) >> (bit_num - j*2 - 1)) );
			bit_pos1 = bit_pos1 << 1;
			bit_pos2 = bit_pos2 >> 1;
		}

		if ((bit_num & 0x01) == 1) {
			t = t | (i & bit_pos1);
		}

		inverted_xn[i] = xn[t];
	}

}

//DIT蝶形运算单元
void dit_butterfly(complex<double> &x1, complex<double> &x2, complex<double> &rotation_factor, complex<double> &y1, complex<double> &y2) {
	complex<double> t = x2 * rotation_factor;
	y2 = x1 - t;
	//可同址运算
	y1 = x1 + t;
}

//基2时间fft
void dit_fft(complex<double>* xn, int n, complex<double>* Xk) {
	//层数
	int m = int(round(log2(n)));

	//获得倒位序
	invert(xn, n, Xk);

	//m层
	for (int i = 0; i < m; i++) {
		int N = 1 << (i + 1);
		complex<double> Wn(cos(-2 * PI / N), sin(-2 * PI / N));

		for (int j = 0; j < n / N; j++) {
			//旋转因子
			complex<double> Wr = 1.0;

			for (int k = 0; k < N / 2; k++) {
				//同址运算
				dit_butterfly(Xk[j*N + k], Xk[j*N + k + N / 2], Wr, Xk[j*N + k], Xk[j*N + k + N / 2]);
				Wr *= Wn;
			}
		}
	}
}

//DIF蝶形运算单元
void dif_butterfly(complex<double> &x1, complex<double> &x2, complex<double> &rotation_factor, complex<double> &y1, complex<double> &y2) {
	//临时变量，用于同址运算
	complex<double> t = x1;
	y1 = t + x2;
	y2 = (t - x2) * rotation_factor;
}

//基2频率fft
void dif_fft(complex<double>* xn, int n, complex<double>* Xk) {
	//层数
	int m = int(round(log2(n)));

	complex<double> *inverted_Xk = new complex<double>[n];
	for (int i = 0; i < n; i++) {
		inverted_Xk[i] = xn[i];
	}

	for (int i = 0; i < m; i++) {
		int N = n >> i;
		complex<double> Wn(cos(-2 * PI / N), sin(-2 * PI / N));

		for (int j = 0; j < n / N; j++) {
			//旋转因子
			complex<double> Wr = 1.0;

			for (int k = 0; k < N / 2; k++) {
				//同址运算
				dif_butterfly(inverted_Xk[j*N + k], inverted_Xk[j*N + k + N/2], Wr, inverted_Xk[j*N + k], inverted_Xk[j*N + k + N/2]);
				Wr *= Wn;
			}
		}
	}

	//倒位序操作
	invert(inverted_Xk, n, Xk);
}

void main()
{
	/*
	for(int i = 0; i < n; i++) {
		xn1[i] = cos(i * 2 * PI / 8);
	}
	*/
	for (int i = 10; i <= 16; i++) {
		int N = 2 << i;
		complex<double> *xn = new complex<double>[N];
		complex<double> *Xk = new complex<double>[N];

		for (int j = 0; j < N; j++) {
			xn[j] = 1; 
			Xk[j] = 0;
		}
		
		cout << "N = " << N << ":" << endl;

		//DFT
		clock_t t1 = clock();
		dft(xn, N, Xk);
		t1 = clock() - t1;
		cout << "DFT time cost: " << (double)t1 / CLOCKS_PER_SEC << ", " << "X[0], X[1], X[3]: " << Xk[0] << ", " << Xk[1] << ", " << Xk[2] << endl;

		//DIT-FFT
		clock_t t2 = clock();
		dit_fft(xn, N, Xk);
		t2 = clock() - t2;
		cout << "DIT-FFT time cost: " << (double)t2 / CLOCKS_PER_SEC << ", " << "X[0], X[1], X[3]: " << Xk[0] << ", " << Xk[1] << ", " << Xk[2] << endl;

		//DIF-FFT
		clock_t t3 = clock();
		dif_fft(xn, N, Xk);
		t3 = clock() - t3;
		cout << "DIT-FFT time cost: " << (double)t3 / CLOCKS_PER_SEC << ", " << "X[0], X[1], X[3]: " << Xk[0] << ", " << Xk[1] << ", " << Xk[2] << endl;

		cout << endl;

		delete []xn;
		delete []Xk;
	}
}